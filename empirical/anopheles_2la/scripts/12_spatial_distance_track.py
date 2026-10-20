#!/usr/bin/env python3
"""Stage 12A: Anopheles 2La local-tree patristic-distance tracks.

The input trees are NJ trees from haplotype distances.  This script treats
their branch lengths as descriptive substitution-distance units, never as
calibrated coalescence times.  The arrangement contrast is restricted to
within-species homozygous 2La/2L+a haplotypes; heterokaryotypes are retained
for arrangement-blind summaries only.
"""
from __future__ import annotations
import csv, gzip, hashlib, json, math, random, statistics, sys, io
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
DATA=ROOT/'data/anopheles_2la'; EMP=ROOT/'empirical/anopheles_2la'
OUT=EMP/'results/stage12_distance'; FIG=EMP/'figures'; OUT.mkdir(parents=True,exist_ok=True); FIG.mkdir(parents=True,exist_ok=True)
TREES=DATA/'processed/stage4r_all72_haplotype_local_trees.nwk.gz'
MAP=DATA/'processed/stage4r_haplotype_species_map.tsv'; GRID=DATA/'processed/stage2r_window_grid.tsv'
QC=EMP/'results/stage4r_window_qc.tsv'; TOPO=EMP/'results/stage2r_quartet_topology_signal.tsv'; REG=DATA/'metadata/region_manifest.tsv'
SEED=20261019

def rows(path):
    with open(path,newline='') as h: return list(csv.DictReader(h,delimiter='\t'))
def write(path,rs,fields):
    with open(path,'w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader()
        for r in rs: w.writerow({f:'' if r.get(f) is None else r.get(f,'') for f in fields})
def mean(v): return float(np.mean(v)) if v else float('nan')
def finite(x): return x is not None and math.isfinite(float(x))
def distance_matrix(tree):
    adj=defaultdict(list); names={}
    for parent in tree.find_clades(order='preorder'):
        for child in parent.clades:
            w=float(child.branch_length or 0); adj[id(parent)].append((id(child),w)); adj[id(child)].append((id(parent),w))
    for tip in tree.get_terminals(): names[tip.name]=id(tip)
    out={}
    for name,node in names.items():
        stack=[(node,-1,0.0)]
        while stack:
            u,prev,d=stack.pop()
            for other,nid in names.items():
                if nid==u: out[(name,other)]=d
            for v,w in adj[u]:
                if v!=prev: stack.append((v,u,d+w))
    return out
def all_pair_stats(D, tips):
    ds=[]
    for i,a in enumerate(tips):
        for b in tips[i+1:]: ds.append(D[(a,b)])
    return mean(ds),float(np.median(ds)) if ds else float('nan'),float(np.var(ds,ddof=1)) if len(ds)>1 else float('nan'),max(ds) if ds else float('nan')
def means_for_groups(D, groups):
    vals=[]
    for a,b in groups:
        vals += [D[(x,y)] for x in a for y in b if x!=y]
    return mean(vals),len(vals)
def balanced_contrast(D, by_species):
    opposite=[]; same=[]; audit=[]
    for sp,g in sorted(by_species.items()):
        std=g.get('2L+a/2L+a',[]); inv=g.get('2La/2La',[])
        if len(std) and len(inv):
            ov,n=means_for_groups(D,[(std,inv)]); opposite.append(ov)
            same_parts=[]
            if len(std)>1: same_parts.append(means_for_groups(D,[(std,std)])[0])
            if len(inv)>1: same_parts.append(means_for_groups(D,[(inv,inv)])[0])
            if same_parts: same.append(mean(same_parts))
            audit.append((sp,len(std),len(inv),ov,mean(same_parts) if same_parts else float('nan')))
    if not opposite or not same: return float('nan'),float('nan'),0,0,';'.join(f'{x[0]}:{x[1]}/{x[2]}' for x in audit)
    # Equal weight by species and by contrast class; this avoids sample-size dominance.
    opp=mean(opposite); sm=mean(same)
    # Pooled-pair sensitivity is calculated independently below from the audit groups.
    return opp-sm,opp,sm,len(opposite),';'.join(f'{x[0]}:{x[1]}/{x[2]}' for x in audit)
def contrast_pooled(D, by_species):
    op=[]; sm=[]
    for g in by_species.values():
        s=g.get('2L+a/2L+a',[]); i=g.get('2La/2La',[])
        op += [D[(x,y)] for x in s for y in i]
        sm += [D[(x,y)] for z in (s,i) for ix,x in enumerate(z) for y in z[ix+1:]]
    return mean(op)-mean(sm) if op and sm else float('nan')
def spatial_p(rows0, value, nperm=999):
    usable=[r for r in rows0 if r['region_class'] in ('inside','outside') and finite(r['A_distance'])]
    if len(usable)<4:return float('nan')
    obs=mean([float(r['A_distance']) for r in usable if r['region_class']=='inside'])-mean([float(r['A_distance']) for r in usable if r['region_class']=='outside'])
    vals=np.array([float(r['A_distance']) for r in usable]); k=sum(r['region_class']=='inside' for r in usable); rng=random.Random(SEED); null=[]
    for _ in range(nperm):
        rng.shuffle(vals.tolist())
        take=rng.sample(range(len(vals)),k); mask=np.zeros(len(vals),bool); mask[take]=True
        null.append(float(np.mean(vals[mask])-np.mean(vals[~mask])))
    return (1+sum(abs(x)>=abs(obs) for x in null))/(nperm+1)
def main():
    mp={r['haplotype_id']:r for r in rows(MAP)}; grid={r['window_id']:r for r in rows(GRID)}; qc={r['window_id']:r for r in rows(QC)}; topo={r['window_id']:r for r in rows(TOPO)}
    rec=[]; excluded=defaultdict(int)
    with gzip.open(TREES,'rt') as h:
        for line in h:
            wid,nwk=line.rstrip().split('\t',1)
            if wid not in grid: excluded['missing_grid']+=1; continue
            try: tree=Phylo.read(io.StringIO(nwk),'newick')
            except Exception: excluded['parse_error']+=1; continue
            by=defaultdict(lambda:defaultdict(list)); tips=[t.name for t in tree.get_terminals()]
            for tip in tips:
                m=mp.get(tip.split('|',1)[0])
                if m and m['arrangement'] in ('2La/2La','2L+a/2L+a'): by[m['species']][m['arrangement']].append(tip)
            D=distance_matrix(tree); ad,opp,sm,ns,design=balanced_contrast(D,by); pooled=contrast_pooled(D,by)
            allmean,med,var,diam=all_pair_stats(D,tips)
            internal=[float(c.branch_length or 0) for c in tree.get_nonterminals() if c!=tree.root]
            g=grid[wid]; q=topo.get(wid,{})
            rec.append({'window_id':wid,'start':g['start'],'end':g['end'],'midpoint':(int(g['start'])+int(g['end']))/2,'region_class':g['region_class'],'n_tips':len(tips),'n_eligible_species':ns,'arrangement_design':design,'opposite_mean':opp,'same_mean':sm,'A_distance':ad,'A_relative':ad/allmean if finite(ad) and allmean else float('nan'),'A_distance_pooled':pooled,'mean_pair_distance':allmean,'median_pair_distance':med,'pair_distance_variance':var,'tree_diameter':diam,'mean_internal_branch':mean(internal),'total_tree_length':sum(float(c.branch_length or 0) for c in tree.find_clades()),'q_arrangement':q.get('q_arrangement'),'q_species':q.get('q_species'),'q_third':q.get('q_third'),'D':q.get('D'),'n_post_filter_snps':qc.get(wid,{}).get('n_post_filter_snps')})
    fields=list(rec[0]); write(OUT/'stage12_window_distance.tsv',rec,fields)
    p=spatial_p(rec,0)
    summary=[{'metric':'n_windows','value':len(rec)},{'metric':'n_arrangement_eligible','value':sum(finite(r['A_distance']) for r in rec)},{'metric':'arrangement_spatial_null_p','value':p},{'metric':'branch_length_units','value':'NJ patristic substitution-distance units; not coalescence times'},{'metric':'heterokaryotype_policy','value':'excluded from primary arrangement contrast; retained in arrangement-blind summaries'}]
    write(OUT/'stage12_spatial_null.tsv',summary,['metric','value'])
    # Plot with physical gaps preserved by using points and interval bars, not interpolated lines.
    x=np.array([float(r['midpoint'])/1e6 for r in rec]); a=np.array([float(r['A_distance']) if finite(r['A_distance']) else np.nan for r in rec]); ar=np.array([float(r['A_relative']) if finite(r['A_relative']) else np.nan for r in rec]); d=np.array([float(r['D']) if finite(r['D']) else np.nan for r in rec])
    fig,ax=plt.subplots(3,1,figsize=(12,8),sharex=True,constrained_layout=True)
    ax[0].plot(x,d,'.',ms=2,color='#405a8d'); ax[0].axhline(0,color='k',lw=.5); ax[0].set_ylabel('quartet D(w)')
    ax[1].plot(x,a,'.',ms=2,color='#c44e52'); ax[1].axhline(0,color='k',lw=.5); ax[1].set_ylabel('A_distance\nNJ units')
    ax[2].plot(x,ar,'.',ms=2,color='#55a868'); ax[2].axhline(0,color='k',lw=.5); ax[2].set_ylabel('A_relative'); ax[2].set_xlabel('2L position (Mb)')
    for q in ax: q.axvspan(20.524058,42.165532,color='#8c6bb1',alpha=.14)
    fig.suptitle('Anopheles 2La local-tree topology and patristic-distance tracks')
    fig.savefig(FIG/'anopheles_2la_stage12_distance.pdf'); fig.savefig(FIG/'anopheles_2la_stage12_distance.png',dpi=220); plt.close(fig)
    report=f'''# Anopheles Stage 12 spatial distance analysis\n\nProcessed {len(rec)} Stage-4R local NJ trees. Branch lengths are patristic distances in the published NJ substitution-distance scale; no molecular clock or coalescence-time interpretation is made.\n\nThe primary arrangement contrast averages within-species 2La/2L+a pair means and equalizes species before forming `A_distance = opposite - same`. Heterokaryotypes are excluded from that contrast because haplotype phase is not established, but all mapped tips contribute to arrangement-blind summaries. The individual-pair weighted sensitivity is retained per window.\n\nThe frozen 2La interval is 20.524058--42.165532 Mb. Boundary-overlap windows remain in the track and are not silently interpolated. The observed inside/outside contrast spatial-null result is recorded in `stage12_spatial_null.tsv`; it is an alignment test against the prespecified interval, not a test that a change point exists.\n\nA limitation is that the NJ trees are local substitution-distance summaries and may be affected by informative-site density. The per-window post-filter SNP count is included for scale diagnostics.\n'''
    (OUT/'report.md').write_text(report); (OUT/'methods_text.md').write_text(report); (OUT/'results_text.md').write_text(report)
    manifest={'seed':SEED,'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [TREES,MAP,GRID,TOPO,QC]},'n_windows':len(rec),'outputs':['stage12_window_distance.tsv','stage12_spatial_null.tsv']}; (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__': main()
