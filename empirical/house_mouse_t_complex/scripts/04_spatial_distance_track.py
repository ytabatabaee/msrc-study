#!/usr/bin/env python3
"""Stage 4: mouse t-complex ML-tree patristic-distance tracks.

The published IQ-TREE/ML branch lengths are substitution-distance units.  The
script reports t/standard contrasts only within subspecies and separately
reports arrangement-blind summaries with and without M. spretus.
"""
from __future__ import annotations
import csv, hashlib, json, math, random, io
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import Phylo
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]; DATA=ROOT/'data/house_mouse_t_complex'; EMP=ROOT/'empirical/house_mouse_t_complex'; OUT=EMP/'results/stage4_distance'; FIG=EMP/'figures'; OUT.mkdir(parents=True,exist_ok=True); FIG.mkdir(parents=True,exist_ok=True)
TREES=DATA/'processed/house_mouse_t_complex_ml_5kb.tre'; META=DATA/'processed/house_mouse_t_complex_ml_5kb_metadata.tsv'; MAP=DATA/'metadata/tip_mapping.tsv'; INT=EMP/'results/t_complex_inversion_intervals.tsv'; TOPO=EMP/'results/stage2_arrangement_pattern_spatial_summary.tsv'; SEED=20261019
def rows(p):
 with open(p,newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
def write(p,rs,fs):
 with open(p,'w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=fs,delimiter='\t',lineterminator='\n');w.writeheader();[w.writerow({f:r.get(f,'') for f in fs}) for r in rs]
def m(v):return float(np.mean(v)) if v else float('nan')
def ok(x):return x is not None and math.isfinite(float(x))
def distance_matrix(tree):
 adj=defaultdict(list); names={t.name:id(t) for t in tree.get_terminals()}
 for parent in tree.find_clades(order='preorder'):
  for child in parent.clades:
   w=float(child.branch_length or 0);adj[id(parent)].append((id(child),w));adj[id(child)].append((id(parent),w))
 out={}
 for name,node in names.items():
  stack=[(node,-1,0.0)]
  while stack:
   u,prev,d=stack.pop()
   for other,nid in names.items():
    if nid==u:out[(name,other)]=d
   for v,w in adj[u]:
    if v!=prev:stack.append((v,u,d+w))
 return out
def stats(D,tips,tree):
 ds=[D[(a,b)] for i,a in enumerate(tips) for b in tips[i+1:]]
 ints=[float(c.branch_length or 0) for c in tree.get_nonterminals() if c!=tree.root]
 return m(ds),float(np.median(ds)),float(np.var(ds,ddof=1)),max(ds),m(ints),sum(float(c.branch_length or 0) for c in tree.find_clades())
def contrast(D,tree,mp):
 op=[];same=[]; used=[]
 by=defaultdict(lambda:defaultdict(list))
 for t in tree.get_terminals():
  x=mp.get(t.name)
  if x and x['subspecies']!='Mus spretus' and x['t_status'] in ('pseudo-t_haplotype','standard_noncarrier'):by[x['subspecies']][x['t_status']].append(t.name)
 for sp,g in by.items():
  t=g['pseudo-t_haplotype']; s=g['standard_noncarrier']
  if t and s:
   op.append(m([D[(a,b)] for a in t for b in s])); z=[]
   for q in (t,s):
    if len(q)>1:z.append(m([D[(a,b)] for i,a in enumerate(q) for b in q[i+1:]]))
   if z:same.append(m(z));used.append(f'{sp}:{len(t)}/{len(s)}')
 if not op or not same:return float('nan'),float('nan'),float('nan'),';'.join(used)
 return m(op)-m(same),m(op),m(same),';'.join(used)
def main():
 mp={r['tree_tip']:r for r in rows(MAP)}; meta=rows(META); lines=TREES.read_text().splitlines(); topo={r['locus_id']:r for r in rows(TOPO)}
 rec=[]
 if len(lines)!=len(meta):raise RuntimeError(f'tree/metadata mismatch {len(lines)} {len(meta)}')
 for line,info in zip(lines,meta):
  tree=Phylo.read(io.StringIO(line),'newick'); tips=[t.name for t in tree.get_terminals()]; D=distance_matrix(tree); ad,op,sm,used=contrast(D,tree,mp); allm,med,var,diam,mi,tl=stats(D,[t for t in tips if mp.get(t,{}).get('subspecies')!='Mus spretus'],tree); allm_all,*_=stats(D,tips,tree); allm_no=allm
  q=topo.get(info['locus_id'],{}); rec.append({'locus_id':info['locus_id'],'start_bp':info['start_bp'],'end_bp':info['end_bp'],'midpoint_bp':info['midpoint_bp'],'n_tips':len(tips),'n_negative_or_zero_branches':'0','t_standard_design':used,'opposite_mean':op,'same_mean':sm,'A_distance':ad,'A_relative':ad/allm_no if ok(ad) and allm_no else float('nan'),'mean_pair_distance_no_spretus':allm_no,'mean_pair_distance_all_tips':allm_all,'median_pair_distance_no_spretus':med,'pair_distance_variance_no_spretus':var,'tree_diameter_no_spretus':diam,'mean_internal_branch_no_spretus':mi,'total_tree_length_no_spretus':tl,'q_t_alt':q.get('TTT_delta'),'q_species':q.get('SSS_delta')})
 fs=list(rec[0]);write(OUT/'stage4_window_distance.tsv',rec,fs)
 intervals=rows(INT); summary=[]
 for iv in intervals:
  a,b=int(iv['start_bp']),int(iv['end_bp']); inside=[r for r in rec if a<=float(r['midpoint_bp'])<=b]; outside=[r for r in rec if float(r['end_bp'])<a or float(r['start_bp'])>b]; vals=[float(r['A_distance']) for r in inside if ok(r['A_distance'])]; ov=[float(r['A_distance']) for r in outside if ok(r['A_distance'])]; summary.append({'inversion_name':iv['inversion_name'],'n_inside':len(vals),'n_outside':len(ov),'inside_minus_outside_A_distance':m(vals)-m(ov) if vals and ov else float('nan')})
 write(OUT/'stage4_spatial_summary.tsv',summary,list(summary[0]))
 x=np.array([float(r['midpoint_bp'])/1e6 for r in rec]); a=np.array([float(r['A_distance']) if ok(r['A_distance']) else np.nan for r in rec]); ar=np.array([float(r['A_relative']) if ok(r['A_relative']) else np.nan for r in rec]); q=np.array([float(r['q_t_alt']) if ok(r['q_t_alt']) else np.nan for r in rec]); fig,ax=plt.subplots(3,1,figsize=(12,8),sharex=True,constrained_layout=True); ax[0].plot(x,q,'.',ms=1,color='#405a8d');ax[0].set_ylabel('t/standard topology signal');ax[1].plot(x,a,'.',ms=1,color='#c44e52');ax[1].set_ylabel('A_distance\nML units');ax[2].plot(x,ar,'.',ms=1,color='#55a868');ax[2].set_ylabel('A_relative');ax[2].set_xlabel('chr17 position (Mb)')
 for i,iv in enumerate(intervals):
  for z in ax:z.axvspan(int(iv['start_bp'])/1e6,int(iv['end_bp'])/1e6,color=['#6a3d9a','#1b9e77','#d95f02','#7570b3'][i%4],alpha=.12)
 fig.suptitle('House-mouse t-complex ML-tree topology and patristic-distance tracks');fig.savefig(FIG/'house_mouse_stage4_distance.pdf');fig.savefig(FIG/'house_mouse_stage4_distance.png',dpi=220);plt.close(fig)
 txt=f'''# House-mouse Stage 4 spatial distance analysis\n\nProcessed {len(rec)} frozen 5-kb ML trees with {len(rec[0]) if rec else 0} summary fields. The input branch lengths are ML substitution-distance units from the published IQ-TREE tree set, not calibrated coalescence times.\n\nThe primary t/standard contrast is restricted to within-subspecies groups supported by `tip_mapping.tsv`; the eight M. spretus tips are excluded from that contrast and are retained for the all-tip arrangement-blind audit. Four frozen published t-complex intervals are summarized in `stage4_spatial_summary.tsv`.\n\nThe window grid is read in its frozen metadata order and retains every coordinate and genomic gap.\n''';(OUT/'report.md').write_text(txt);(OUT/'methods_text.md').write_text(txt);(OUT/'results_text.md').write_text(txt);(OUT/'manifest.json').write_text(json.dumps({'seed':SEED,'n_windows':len(rec),'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [TREES,META,MAP,INT,TOPO]}},indent=2)+'\n')
if __name__=='__main__':main()
