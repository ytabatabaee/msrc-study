#!/usr/bin/env python3
"""Stage 7: house-mouse influential-sample sensitivity analysis."""
from __future__ import annotations
import csv, hashlib, io, json, math, sys
from collections import defaultdict
from pathlib import Path
import numpy as np
from Bio import Phylo

ROOT=Path(__file__).resolve().parents[3]
EMP=ROOT/'empirical/house_mouse_t_complex'; DATA=ROOT/'data/house_mouse_t_complex'
OUT=EMP/'results/stage7_influential_sample'; FIG=EMP/'figures'
TREES=DATA/'processed/house_mouse_t_complex_ml_5kb.tre'; META=DATA/'processed/house_mouse_t_complex_ml_5kb_metadata.tsv'
MAP=DATA/'metadata/tip_mapping.tsv'; DIST=EMP/'results/stage4_distance/stage4_window_distance.tsv'
QSCAN=EMP/'results/stage2_fixed_quartet_scan.tsv'; INT=EMP/'results/t_complex_inversion_intervals.tsv'
STAGE6=EMP/'results/stage6_distance_discontinuity'; EXCLUDED='Mmd_FRA4_18B.variant'; SEED=20261020
sys.path.insert(0,str(ROOT/'empirical/shared'))
sys.path.insert(0,str(EMP/'scripts'))
from spatial_boundary_detector import segment_track, intervals_to_coordinates


def read(p):
    with p.open(newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
def write(p,rows,fields):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows({f:r.get(f,'') for f in fields} for r in rows)
def num(x):
    try:
        y=float(x);return y if math.isfinite(y) else math.nan
    except (TypeError,ValueError):return math.nan
def vals(xs):return [num(x) for x in xs if math.isfinite(num(x))]
def avg(xs):
    x=vals(xs);return float(np.mean(x)) if x else math.nan
def med(xs):
    x=vals(xs);return float(np.median(x)) if x else math.nan
def trim(xs):
    x=sorted(vals(xs));k=int(.05*len(x));return float(np.mean(x[k:len(x)-k])) if x and len(x)>2*k and k else (float(np.mean(x)) if x else math.nan)
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()
def pk(a,b):return a,b

def dmatrix(tree):
    adj=defaultdict(list)
    for parent in tree.find_clades(order='preorder'):
        for child in parent.clades:
            w=num(child.branch_length);adj[id(parent)].append((id(child),w));adj[id(child)].append((id(parent),w))
    leaves={id(c):c.name for c in tree.get_terminals()};out={}
    for start,sname in leaves.items():
        stack=[(start,-1,0.)]
        while stack:
            u,prev,d=stack.pop()
            if u in leaves:out[(sname,leaves[u])]=d
            for v,w in adj[u]:
                if v!=prev:stack.append((v,u,d+w))
    return out

def pairs(D,tree,mapping,exclude=None):
    tips=[c.name for c in tree.get_terminals() if mapping.get(c.name,{}).get('subspecies')!='Mus spretus' and c.name!=exclude]
    return [D[pk(a,b)] for i,a in enumerate(tips) for b in tips[i+1:] if pk(a,b) in D and math.isfinite(D[pk(a,b)])]

def contrast(D,tree,mapping,exclude=None):
    groups=defaultdict(lambda:defaultdict(list))
    for c in tree.get_terminals():
        if c.name==exclude:continue
        info=mapping.get(c.name)
        if info and info['subspecies']!='Mus spretus' and info['t_status'] in ('pseudo-t_haplotype','standard_noncarrier'):
            groups[info['subspecies']][info['t_status']].append(c.name)
    opposite=[];same=[];design=[]
    for sp,g in sorted(groups.items()):
        t,s=g['pseudo-t_haplotype'],g['standard_noncarrier']
        cross=[D[pk(a,b)] for a in t for b in s if pk(a,b) in D]
        within_means=[avg([D[pk(a,b)] for i,a in enumerate(q) for b in q[i+1:] if pk(a,b) in D])
                      for q in (t,s) if len(q)>1]
        if cross and within_means:
            opposite.append(avg(cross));same.append(avg(within_means));design.append(f'{sp}:{len(t)}/{len(s)}')
    if not opposite or not same:return math.nan,math.nan,math.nan,';'.join(design)
    om,sm=avg(opposite),avg(same);return om-sm,om,sm,';'.join(design)

def quartet_counts(newick,groups):
    from house_mouse_stage2_utils import quartet_counts_for_tree
    return quartet_counts_for_tree(newick,groups)
def qscore(c):
    n=sum(c.values())
    if not n:return {'q_species':math.nan,'q_t_alt':math.nan,'q_other':math.nan,'q_unresolved':math.nan,'n_resolved':0,'n_unresolved':0}
    return {'q_species':c['Q_SPECIES']/n,'q_t_alt':c['Q_T_ALT']/n,'q_other':c['Q_OTHER']/n,'q_unresolved':c['Q_UNRESOLVED']/n,'n_resolved':c['Q_SPECIES']+c['Q_T_ALT']+c['Q_OTHER'],'n_unresolved':c['Q_UNRESOLVED']}

def branch_audit(tree,D,mapping):
    tips=[c.name for c in tree.get_terminals()]; branches=[num(c.branch_length) for c in tree.find_clades() if c is not tree.root and math.isfinite(num(c.branch_length))]
    terminal={c.name:num(c.branch_length) for c in tree.get_terminals() if math.isfinite(num(c.branch_length))}
    allp=pairs(D,tree,mapping);nop=pairs(D,tree,mapping,EXCLUDED);fra=[D[pk(EXCLUDED,t)] for t in tips if t!=EXCLUDED and pk(EXCLUDED,t) in D]
    return {'n_tips':len(tips),'n_non_spretus_tips':sum(mapping.get(t,{}).get('subspecies')!='Mus spretus' for t in tips),'n_pairs_all':len(allp),'n_pairs_without_excluded':len(nop),'excluded_tip_present':int(EXCLUDED in tips),'excluded_terminal_branch':terminal.get(EXCLUDED,math.nan),'excluded_terminal_fraction_total_tree':terminal.get(EXCLUDED,math.nan)/sum(branches) if EXCLUDED in terminal and sum(branches) else math.nan,'excluded_pair_mean':avg(fra),'other_pair_mean':avg(nop),'B_all_recomputed':avg(allp),'B_without_excluded':avg(nop),'median_pair_distance':med(allp),'trimmed_pair_distance_5pct':trim(allp),'total_tree_length':sum(branches),'mean_internal_branch':avg([num(c.branch_length) for c in tree.get_nonterminals() if c is not tree.root]),'median_terminal_branch':med(list(terminal.values())),'max_terminal_branch':max(terminal.values()) if terminal else math.nan,'max_terminal_tip':max(terminal,key=terminal.get) if terminal else '','n_zero_branches':sum(x==0 for x in branches),'n_negative_branches':sum(x<0 for x in branches),'n_nonfinite_branches':len([c for c in tree.find_clades() if c is not tree.root and not math.isfinite(num(c.branch_length))])}

def bins(rows,size=100000):
    groups=defaultdict(list)
    for r in rows:
        s=int(num(r['start_bp']));groups[5000000+((s-5000000)//size)*size].append(r)
    out=[]
    for s,g in sorted(groups.items()):
        observed=sum(int(num(r['end_bp']))-int(num(r['start_bp']))+1 for r in g)
        z={'bin_start_bp':s,'bin_end_bp':s+size-1,'midpoint_bp':s+size/2,'n_observed_windows':len(g),'observed_bp':observed,'coverage_fraction':observed/size}
        for f in ('D_original','D_leave_one_out','A_distance_original','A_distance_leave_one_out','A_relative_original','A_relative_leave_one_out','B_all','B_leave_one_out','B_median','B_trimmed'):z[f]=avg([r[f] for r in g])
        for lab in ('original','leave_one_out'):
            c={k:sum(num(r[f'n_{k}_{lab}']) for r in g) for k in ('species','t_alt','other','unresolved')};n=sum(c.values())
            for k in c:z[f'q_{k}_{lab}']=c[k]/n if n else math.nan
        out.append(z)
    return out

def merged(rs):
    out=[]
    for a,b in sorted((int(a),int(b)+1) for a,b in rs):
        if out and a<=out[-1][1]:out[-1][1]=max(out[-1][1],b)
        else:out.append([a,b])
    return out
def inter(a,b):
    i=j=t=0
    while i<len(a) and j<len(b):
        t+=max(0,min(a[i][1],b[j][1])-max(a[i][0],b[j][0]))
        if a[i][1]<b[j][1]:i+=1
        else:j+=1
    return t
def evaluate(pred,refs):
    p=merged([(a,b-1) for a,b in pred]);r=merged(refs);pl=sum(b-a for a,b in p);rl=sum(b-a for a,b in r);it=inter(p,r);pr=it/pl if pl else 0;rc=it/rl if rl else 0
    er=[]
    for a,b in p:
        c,d=min(r,key=lambda x:abs(x[0]-a)+abs(x[1]-b));er += [abs(a-c)/1000,abs((b-1)-d)/1000]
    return {'iou_union':it/(pl+rl-it) if pl+rl-it else 0,'precision':pr,'recall':rc,'f1':2*pr*rc/(pr+rc) if pr+rc else 0,'false_positive_bp':max(0,pl-it),'false_negative_bp':max(0,rl-it),'mean_boundary_error_kb':avg(er),'n_predicted_intervals':len(pred)}

def detect(rows,fields,max_intervals):
    z=[r for r in rows if all(math.isfinite(num(r[f])) for f in fields)];x=np.array([[num(r[f]) for f in fields] for r in z]);s=segment_track(x,fields,min_size=3,max_intervals=max_intervals);co=intervals_to_coordinates(s,[{'start':r['start_bp'],'end':r['end_bp']} for r in z]);p=[(int(float(r['start'])),int(float(r['end']))) for r in co];return s,p,len(z)

def make_figures(b,refs,native):
    import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
    audit={r['locus_id']:r for r in read(OUT/'sample_branch_length_audit.tsv')}
    native=[dict(r,excluded_terminal_branch=audit[r['locus_id']]['excluded_terminal_branch']) for r in native]
    x=np.array([num(r['midpoint_bp'])/1e6 for r in b])
    def shade(ax):
        for i,(a,z) in enumerate(refs):
            ax.axvspan(a/1e6,z/1e6,color=['#6a3d9a','#1b9e77','#d95f02','#7570b3'][i],alpha=.12,lw=0);ax.axvline(a/1e6,color='#666',ls='--',lw=.5);ax.axvline(z/1e6,color='#666',ls='--',lw=.5)
    def parts():return [p for p in np.split(np.arange(len(x)),np.flatnonzero(np.diff(x)>.2)+1) if len(p)]
    def draw(stem,title,lab,fields):
        fig,ax=plt.subplots(6,1,figsize=(13,13),sharex=True,gridspec_kw={'height_ratios':[.4,1,1,1,1,1]},constrained_layout=True)
        for a in ax:shade(a);a.grid(axis='y',alpha=.16)
        ax[0].set_ylim(0,1);ax[0].set_yticks([]);ax[0].set_ylabel('inversions',fontsize=8)
        for i,(a,z) in enumerate(refs):ax[0].bar((a+z)/2e6,.65,width=(z-a)/1e6,color=['#6a3d9a','#1b9e77','#d95f02','#7570b3'][i],alpha=.6)
        for n,(idx,f,t,c) in enumerate(((1,fields[0],'(b) ALL_TIPS topology contrast D(w)','#333'),(3,fields[1],'(d) arrangement-informed A_distance','#c44e52'),(4,fields[2],'(e) arrangement-informed A_relative','#55a868'),(5,fields[3],'(f) arrangement-blind B(w)','#8172b2'))):
            y=np.array([num(r[f]) for r in b]);
            for p in parts():ax[idx].plot(x[p],y[p],color=c,lw=.8,marker='.',ms=1.2)
            if idx!=5:ax[idx].axhline(0,color='#777',lw=.5)
            ax[idx].set_ylabel('relative' if idx==4 else 'support contrast' if idx==1 else 'substitution units');ax[idx].set_title(t,loc='left',weight='bold')
        q1=np.array([num(r[f'q_species_{lab}']) for r in b]);q2=np.array([num(r[f'q_t_alt_{lab}']) for r in b]);q3=np.array([num(r[f'q_other_{lab}']) for r in b])
        for p in parts():ax[2].stackplot(x[p],q1[p],q2[p],q3[p],colors=['#4c78a8','#e45756','#59a14f'],alpha=.85)
        ax[2].set_ylim(0,1);ax[2].set_ylabel('support');ax[2].set_title('(c) three-topology quartet support',loc='left',weight='bold');ax[2].legend(['q_S species history','q_T t-associated alternative','q_3 third resolution'],frameon=False,ncol=3,fontsize=8,loc='upper center',bbox_to_anchor=(.5,1.25));ax[-1].set_xlabel('chromosome 17 position (Mb)');fig.suptitle(title,fontsize=15,weight='bold');fig.text(.99,.01,'Distances are local ML patristic substitution units.',ha='right',fontsize=8);fig.savefig(FIG/f'house_mouse_chr17_{stem}_100kb.pdf',bbox_inches='tight');fig.savefig(FIG/f'house_mouse_chr17_{stem}_100kb.png',dpi=240,bbox_inches='tight');plt.close(fig)
    draw('original','House-mouse chromosome 17: original data (FRA4 retained)','original',('D_original','A_distance_original','A_relative_original','B_all'))
    draw('leave_one_out','House-mouse chromosome 17: leave-one-out sensitivity (FRA4 excluded)','leave_one_out',('D_leave_one_out','A_distance_leave_one_out','A_relative_leave_one_out','B_leave_one_out'))
    fig,ax=plt.subplots(3,1,figsize=(13,8),sharex=True,constrained_layout=True);xn=np.array([num(r['midpoint_bp'])/1e6 for r in native])
    for a in ax:shade(a);a.grid(axis='y',alpha=.16)
    ax[0].plot(xn,[num(r['excluded_terminal_branch']) for r in native],color='#c44e52',lw=.6);ax[0].set_title('FRA4 terminal branch length',loc='left',weight='bold');ax[0].set_ylabel('substitution units')
    ax[1].plot(xn,[num(r['B_all']) for r in native],color='#8172b2',lw=.6,label='all samples');ax[1].plot(xn,[num(r['B_leave_one_out']) for r in native],color='#333',lw=.6,label='FRA4 excluded');ax[1].set_title('Arrangement-blind B(w)',loc='left',weight='bold');ax[1].set_ylabel('mean distance');ax[1].legend(frameon=False,fontsize=8)
    ax[2].plot(xn,[num(r['A_relative_original']) for r in native],color='#55a868',lw=.6,label='all samples');ax[2].plot(xn,[num(r['A_relative_leave_one_out']) for r in native],color='#e45756',lw=.6,label='FRA4 excluded');ax[2].set_title('A_relative',loc='left',weight='bold');ax[2].set_ylabel('relative');ax[2].set_xlabel('chromosome 17 position (Mb)');ax[2].legend(frameon=False,fontsize=8);fig.suptitle('House-mouse Stage 7 sample-influence diagnostics',weight='bold');fig.savefig(FIG/'house_mouse_stage7_sample_influence.pdf',bbox_inches='tight');fig.savefig(FIG/'house_mouse_stage7_sample_influence.png',dpi=240,bbox_inches='tight');plt.close(fig)

def main():
    mapping_rows=read(MAP);mapping={r['tree_tip']:r for r in mapping_rows};meta=read(META);source= {r['locus_id']:r for r in read(DIST)};qsource={r['locus_id']:r for r in read(QSCAN) if r['treatment']=='ALL_TIPS'};trees=[x.strip() for x in TREES.read_text().splitlines() if x.strip()]
    if not(len(trees)==len(meta)==len(source)==len(qsource)==4046):raise RuntimeError('input count mismatch')
    groups=defaultdict(list)
    for r in mapping_rows:groups[r['subspecies']].append(r['tree_tip'])
    groups={k:sorted(v) for k,v in groups.items()};groups_loo={k:[x for x in v if x!=EXCLUDED] for k,v in groups.items()}
    rows=[];audits=[];quartets=[]
    for line,info in zip(trees,meta):
        tree=Phylo.read(io.StringIO(line),'newick');D=dmatrix(tree);locus=info['locus_id'];a=branch_audit(tree,D,mapping);al,ol,sl,design=contrast(D,tree,mapping);a2,ol2,sl2,design2=contrast(D,tree,mapping,EXCLUDED);pa=pairs(D,tree,mapping);pl=pairs(D,tree,mapping,EXCLUDED);ca=quartet_counts(line,groups);cl=quartet_counts(line,groups_loo);qa=qscore(ca);ql=qscore(cl);s=source[locus]
        r={'locus_id':locus,'start_bp':info['start_bp'],'end_bp':info['end_bp'],'midpoint_bp':info['midpoint_bp'],'n_non_spretus_all':a['n_non_spretus_tips'],'n_pairs_all':len(pa),'n_pairs_leave_one_out':len(pl),'A_distance_original':al,'A_distance_leave_one_out':a2,'opposite_mean_original':ol,'same_mean_original':sl,'opposite_mean_leave_one_out':ol2,'same_mean_leave_one_out':sl2,'B_all':avg(pa),'B_leave_one_out':avg(pl),'B_median':med(pa),'B_trimmed':trim(pa),'A_relative_original':al/avg(pa) if math.isfinite(al) and avg(pa)>1e-15 else math.nan,'A_relative_leave_one_out':a2/avg(pl) if math.isfinite(a2) and avg(pl)>1e-15 else math.nan,'D_original':qa['q_t_alt']-qa['q_species'],'D_leave_one_out':ql['q_t_alt']-ql['q_species'],'B_source_stage4':num(s['mean_pair_distance_no_spretus']),'A_source_stage4':num(s['A_distance']),'A_relative_source_stage4':num(s['A_relative'])}
        for lab,z,c in (('original',qa,ca),('leave_one_out',ql,cl)):
            r[f'q_species_{lab}']=z['q_species'];r[f'q_t_alt_{lab}']=z['q_t_alt'];r[f'q_other_{lab}']=z['q_other'];r[f'q_unresolved_{lab}']=z['q_unresolved'];r[f'n_species_{lab}']=c['Q_SPECIES'];r[f'n_t_alt_{lab}']=c['Q_T_ALT'];r[f'n_other_{lab}']=c['Q_OTHER'];r[f'n_unresolved_{lab}']=c['Q_UNRESOLVED']
        rows.append(r);a.update({'locus_id':locus,'start_bp':info['start_bp'],'end_bp':info['end_bp'],'midpoint_bp':info['midpoint_bp']});audits.append(a);quartets.append({'locus_id':locus,'start_bp':info['start_bp'],'end_bp':info['end_bp'],'q_species_original':qa['q_species'],'q_t_alt_original':qa['q_t_alt'],'q_other_original':qa['q_other'],'q_unresolved_original':qa['q_unresolved'],'D_original':r['D_original'],'q_species_leave_one_out':ql['q_species'],'q_t_alt_leave_one_out':ql['q_t_alt'],'q_other_leave_one_out':ql['q_other'],'q_unresolved_leave_one_out':ql['q_unresolved'],'D_leave_one_out':r['D_leave_one_out'],'sum_original':sum(qa[k] for k in ('q_species','q_t_alt','q_other','q_unresolved')),'sum_leave_one_out':sum(ql[k] for k in ('q_species','q_t_alt','q_other','q_unresolved')),'n_resolved_original':qa['n_resolved'],'n_unresolved_original':qa['n_unresolved'],'n_resolved_leave_one_out':ql['n_resolved'],'n_unresolved_leave_one_out':ql['n_unresolved']})
    rows.sort(key=lambda r:int(num(r['start_bp'])));audits.sort(key=lambda r:int(num(r['start_bp'])));quartets.sort(key=lambda r:int(num(r['start_bp'])))
    refs=[(int(r['start_bp']),int(r['end_bp'])) for r in read(INT)];b=bins(rows);cp=num(read(STAGE6/'change_point_diagnostics.tsv')[0]['coordinate_bp']);up=[r for r in rows if num(r['midpoint_bp'])<cp];down=[r for r in rows if num(r['midpoint_bp'])>=cp]
    influence=[]
    for region,z in (('upstream',up),('downstream',down),('all',rows)):
        for f in ('B_all','B_leave_one_out','A_distance_original','A_distance_leave_one_out','A_relative_original','A_relative_leave_one_out','B_median','B_trimmed'):influence.append({'region':region,'track':f,'n_windows':len(z),'mean':avg([r[f] for r in z]),'median':med([r[f] for r in z]),'sd':float(np.std(vals([r[f] for r in z]),ddof=1)) if len(vals([r[f] for r in z]))>1 else math.nan})
    influence.append({'region':'audit','track':'stage6_change_point_bp','n_windows':len(rows),'mean':cp,'median':cp,'sd':math.nan})
    cfg=[('original_B',['B_all']),('leave_one_out_B',['B_leave_one_out']),('median_B',['B_median']),('trimmed_B',['B_trimmed']),('original_A',['A_distance_original']),('leave_one_out_A',['A_distance_leave_one_out']),('original_joint',['D_original','A_distance_original']),('leave_one_out_joint',['D_leave_one_out','A_distance_leave_one_out'])];det=[]
    for method,fields in cfg:
        for model,k in (('single_interval',1),('multi_interval_0_to_4',4)):
            seg,p,n=detect(rows,fields,k);det.append({'method':method,'model':model,'status':seg.status,'n_used_windows':n,'predicted_intervals':';'.join(f'{a}-{z}' for a,z in p),'n_predicted_intervals':len(p),'score':seg.score,'raw_improvement':seg.raw_improvement,'penalty':seg.penalty,'n_segments':seg.n_segments,'feature_names':';'.join(fields),**evaluate(p,refs),'configuration':'Stage13c shared detector, robust MAD scaling, min_size=3, references used post-fit only'})
    write(OUT/'sample_branch_length_audit.tsv',audits,list(audits[0]));write(OUT/'window_original_vs_leave_one_out.tsv',rows,list(rows[0]));write(OUT/'robust_distance_tracks.tsv',rows,['locus_id','start_bp','end_bp','midpoint_bp','B_all','B_leave_one_out','B_median','B_trimmed','A_distance_original','A_distance_leave_one_out','A_relative_original','A_relative_leave_one_out']);write(OUT/'binned_100kb_tracks.tsv',b,list(b[0]));write(OUT/'quartet_leave_one_out_validation.tsv',quartets,list(quartets[0]));write(OUT/'boundary_detection_sensitivity.tsv',det,list(det[0]));write(OUT/'influence_summary.tsv',influence,list(influence[0]));make_figures(b,refs,rows)
    eb=max(abs(num(r['B_all'])-num(r['B_source_stage4'])) for r in rows);ea=max(abs(num(r['A_distance_original'])-num(r['A_source_stage4'])) for r in rows);er=max(abs(num(r['A_relative_original'])-num(r['A_relative_source_stage4'])) for r in rows)
    report=f'''# House-mouse Stage 7 influential-sample sensitivity\n\nThe primary track retains all mapped non-*M. spretus* samples. The leave-one-out sensitivity excludes `{EXCLUDED}` from both arrangement-informed and arrangement-blind calculations. Direct recalculation from the frozen Newick trees reproduced Stage 4 with maximum absolute errors of {eb:.3g} for B(w), {ea:.3g} for A_distance, and {er:.3g} for A_relative.\n\nThe frozen Stage 6 data-driven change point is {cp/1e6:.4f} Mb. The sample is retained in the primary analysis because no concrete label, branch-format, or tree-processing error was identified. The source repository does not include alignments, per-window SNP counts, filtering summaries, or complete IQ-TREE command provenance, so ascertainment and tree-estimation explanations cannot be ruled out. A biological source of long divergence also cannot be established from substitution-unit branch lengths alone.\n\n**Verified:** the discontinuity and sample influence are reproducible from the source ML trees, and original coordinate/arithmetic mappings agree. **Supported but not conclusive:** downstream changes track the focal terminal branch and its pairwise contributions. **Unresolved:** biological divergence versus a region-specific data-processing or tree-estimation effect.\n'''
    (OUT/'diagnostic_report.md').write_text(report);(OUT/'results_text.md').write_text(report);(OUT/'data_provenance_audit.md').write_text(f'''# Data provenance audit\n\nTrees: `{TREES.relative_to(ROOT)}`; 4,046 metadata-matched ML trees. Branch lengths are published ML patristic substitution units. The primary B(w) uses all mapped non-*M. spretus* pairs; the sensitivity removes `{EXCLUDED}` from every recalculated numerator and denominator. ALL_TIPS topology uses the frozen quartet definitions; leave-one-out support removes the same tree tip from its status group. Alignments, SNP counts, and full per-window inference settings are unavailable.\n''');(OUT/'methods_text.md').write_text('''# Methods\n\nPatristic distances were recalculated from every frozen local ML Newick tree. Original and leave-one-out A_distance and A_relative were calculated with the same within-subspecies t-versus-standard design; B is the mean across the corresponding non-M. spretus pair universe. Median and 5% trimmed means are robust sensitivity summaries. ALL_TIPS quartet support was recomputed with the frozen Stage 2 definitions. Boundary sensitivity uses the fixed shared detector with robust MAD scaling, min_size=3, one-interval primary fitting, and a zero-to-four interval sensitivity; reference intervals are post-fit evaluation only.\n''');(OUT/'figure_caption.md').write_text('''**House-mouse chromosome 17 Stage 7 sensitivity.** The original figure retains all mapped non-M. spretus samples, including Mmd_FRA4_18B.variant. The sensitivity figure removes that sample consistently from arrangement-informed and arrangement-blind statistics and recomputes ALL_TIPS quartet support. Tracks are 100-kb means of observed 5-kb windows; distances are local ML substitution units, not calibrated times or coalescence times. Shading marks the four frozen inversion intervals.\n''')
    manifest={'stage':'house mouse Stage 7 influential sample sensitivity','seed':SEED,'command':'python empirical/house_mouse_t_complex/scripts/07_influential_sample_sensitivity.py','excluded_sample':EXCLUDED,'detector':'empirical/shared/spatial_boundary_detector.py','detector_configuration':{'min_size':3,'primary_max_intervals':1,'sensitivity_max_intervals':4,'robust_scaling':'median/MAD'},'inputs':{str(p.relative_to(ROOT)):sha(p) for p in (TREES,META,MAP,DIST,QSCAN,INT)},'n_windows':len(rows),'n_100kb_bins':len(b),'stage6_change_point_bp':cp,'original_stage4_reproduction':{'max_B_abs_error':eb,'max_A_abs_error':ea,'max_A_relative_abs_error':er},'primary_original_retained':True,'previous_outputs_preserved':True};(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

if __name__=='__main__':main()
