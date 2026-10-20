#!/usr/bin/env python3
"""Create unified ALL_TIPS topology and frozen Stage-4 distance tracks.

The historical Stage-2 quartet files retain their original
q_species-q_t_alt sign.  This stage writes the presentation-standardized
D = q_t_alt-q_species without modifying those files.
"""
from __future__ import annotations
import csv, hashlib, json, math, os
from pathlib import Path
from collections import defaultdict
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/msrc-paper-mpl')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]
EMP=ROOT/'empirical/house_mouse_t_complex'; RES=EMP/'results'; OUT=RES/'stage5_unified_tracks'; FIG=EMP/'figures'
Q=RES/'stage2_fixed_quartet_scan.tsv'; PAT=RES/'stage2_arrangement_pattern_spatial_summary.tsv'; DIST=RES/'stage4_distance/stage4_window_distance.tsv'; META=ROOT/'data/house_mouse_t_complex/processed/house_mouse_t_complex_ml_5kb_metadata.tsv'; INT=RES/'t_complex_inversion_intervals.tsv'; SEED=20261020
TOPO_COLORS={'q_species':'#4c78a8','q_t_alt':'#f58518','q_other':'#54a24b'}
INT_COLORS=['#6a3d9a','#1b9e77','#d95f02','#7570b3']

def read(p):
 with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows({k:r.get(k,'') for k in fields} for r in rows)
def f(x):
 try:return float(x)
 except:return math.nan
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def segments(x):
 if len(x)<2:return [np.arange(len(x))]
 d=np.diff(x);gap=max(np.median(d)*2.5,0.01);return [z for z in np.split(np.arange(len(x)),np.flatnonzero(d>gap)+1) if len(z)]
def intervals():return [(int(r['start_bp']),int(r['end_bp']),r['inversion_name']) for r in read(INT)]
def shade(ax):
 for i,(a,b,name) in enumerate(intervals()):ax.axvspan(a/1e6,b/1e6,color=INT_COLORS[i%4],alpha=.12);ax.axvline(a/1e6,color=INT_COLORS[i%4],lw=.55,ls='--');ax.axvline(b/1e6,color=INT_COLORS[i%4],lw=.55,ls='--')
def plot_track(ax,x,y,color,label,linestyle='-',marker=None):
 for z in segments(x):ax.plot(x[z],y[z],color=color,lw=1.0 if marker is None else .7,marker=marker,ms=1.2 if marker else 0,label=label if z[0]==0 else None)
def finish(fig,stem):
 fig.savefig(FIG/(stem+'.pdf'),bbox_inches='tight');fig.savefig(FIG/(stem+'.png'),dpi=240,bbox_inches='tight');plt.close(fig)
def tracks(rows,title,stem,stack=False):
 x=np.array([f(r['midpoint_bp'])/1e6 for r in rows]);fig,ax=plt.subplots(4,1,figsize=(13,10),sharex=True,gridspec_kw={'height_ratios':[1,1,1,1.35]},constrained_layout=True)
 for a in ax:shade(a);a.grid(axis='y',alpha=.18);a.set_xlim(x.min(),x.max())
 plot_track(ax[0],x,np.array([f(r['D']) for r in rows]),'#333333','D = q_T - q_S',marker='.')
 ax[0].axhline(0,color='#777',lw=.6);ax[0].set_ylabel('D(w)');ax[0].set_title('(a) ALL_TIPS topology contrast',loc='left',weight='bold')
 plot_track(ax[1],x,np.array([f(r['A_distance']) for r in rows]),'#c44e52','A_distance',marker='.')
 ax[1].axhline(0,color='#777',lw=.6);ax[1].set_ylabel('substitution units');ax[1].set_title('(b) t - standard patristic distance',loc='left',weight='bold')
 plot_track(ax[2],x,np.array([f(r['A_relative']) for r in rows]),'#55a868','A_relative',marker='.')
 ax[2].axhline(0,color='#777',lw=.6);ax[2].set_ylabel('relative contrast');ax[2].set_title('(c) normalized distance contrast',loc='left',weight='bold')
 q=[np.array([f(r[k]) for r in rows]) for k in ('q_species','q_t_alt','q_other')]
 base=np.zeros(len(rows))
 for k,v,c,l in zip(('q_species','q_t_alt','q_other'),q,TOPO_COLORS.values(),('q_S species history','q_T t-associated alternative','q_3 third')):
  for z in segments(x):ax[3].fill_between(x[z],base[z],base[z]+v[z],color=c,alpha=.82,label=l if z[0]==0 else None)
  base+=v
 ax[3].set_ylim(0,1);ax[3].set_ylabel('support');ax[3].set_title('(d) ALL_TIPS three-topology support',loc='left',weight='bold');ax[3].set_xlabel('chromosome 17 position (Mb)');ax[3].legend(loc='upper right',ncol=3,fontsize=8,frameon=False)
 fig.suptitle(title,fontsize=15,weight='bold');finish(fig,stem)
def bins(rows,size=250_000):
 g=defaultdict(list)
 for r in rows:g[5_000_000+((int(r['start_bp'])-5_000_000)//size)*size].append(r)
 out=[]
 for start,z in sorted(g.items()):
  d={'bin_start_bp':start,'bin_end_bp':start+size-1,'midpoint_bp':start+size/2,'n_windows':len(z)}
  for k in ('q_species','q_t_alt','q_other','D','A_distance','A_relative'):d[k]=float(np.mean([f(r[k]) for r in z]))
  out.append(d)
 return out
def main():
 q={r['locus_id']:r for r in read(Q) if r['treatment']=='ALL_TIPS'}; pat={r['locus_id']:r for r in read(PAT)}; dist={r['locus_id']:r for r in read(DIST)}; meta={r['locus_id']:r for r in read(META)}
 if set(q)!=set(dist):raise RuntimeError(f'ALL_TIPS/distance locus mismatch {len(q)} {len(dist)}')
 rows=[];checks=[]
 for locus in sorted(q,key=lambda x:int(meta[x]['start_bp'])):
  a=q[locus];d=dist[locus];p=pat.get(locus,{})
  vals={k:f(a[k]) for k in ('q_species','q_t_alt','q_other','q_unresolved')};resolved=sum(vals[k] for k in ('q_species','q_t_alt','q_other'));old=f(a['delta_species_alt']);new=vals['q_t_alt']-vals['q_species'];ttt_old=f(p.get('TTT_delta','nan'));checks.append({'locus_id':locus,'q_range_ok':str(all(-1e-12<=vals[k]<=1+1e-12 for k in ('q_species','q_t_alt','q_other'))).lower(),'resolved_sum':resolved,'unresolved_fraction':vals['q_unresolved'],'sum_all':resolved+vals['q_unresolved'],'D_new':new,'negative_old_delta':-old,'D_difference':new+old,'TTT_D_new':-ttt_old})
  rows.append({'locus_id':locus,'start_bp':a['start_bp'],'end_bp':a['end_bp'],'midpoint_bp':a['midpoint_bp'],'q_species':vals['q_species'],'q_t_alt':vals['q_t_alt'],'q_other':vals['q_other'],'q_unresolved':vals['q_unresolved'],'D':new,'D_old_species_minus_alt':old,'A_distance':f(d['A_distance']),'A_relative':f(d['A_relative']),'mean_pair_distance_no_spretus':f(d['mean_pair_distance_no_spretus']),'TTT_D':-ttt_old})
 fields=list(rows[0]);write(OUT/'window_unified_tracks.tsv',rows,fields); b=bins(rows);write(OUT/'binned_250kb_tracks.tsv',b,list(b[0]));write(OUT/'normalization_checks.tsv',checks,list(checks[0]))
 dif=np.array([f(r['D']) for r in rows]);tt=np.array([f(r['TTT_D']) for r in rows]);corr=float(np.corrcoef(dif,tt)[0,1]);diff=dif-tt;summary=[{'comparison':'ALL_TIPS_vs_TTT','n_windows':len(rows),'pearson_correlation':corr,'mean_ALL_TIPS_D':float(dif.mean()),'mean_TTT_D':float(tt.mean()),'mean_difference_ALL_minus_TTT':float(diff.mean()),'median_abs_difference':float(np.median(abs(diff))),'max_abs_difference':float(max(abs(diff))),'ALL_TIPS_definition':'stage2_fixed_quartet_scan.tsv treatment=ALL_TIPS','TTT_definition':'stage2_arrangement_pattern_spatial_summary.tsv TTT_delta sign reversed'}];write(OUT/'all_tips_vs_ttt_summary.tsv',summary,list(summary[0]))
 tracks(rows,'House-mouse t-complex: ALL_TIPS topology and ML patristic-distance tracks','house_mouse_t_complex_topology_distance_5kb');tracks(b,'House-mouse t-complex: matched 250-kb ALL_TIPS and distance tracks','house_mouse_t_complex_topology_distance_250kb')
 x=np.array([f(r['midpoint_bp'])/1e6 for r in rows]);fig,ax=plt.subplots(figsize=(13,4.2),constrained_layout=True);shade(ax);plot_track(ax,x,dif,'#333333','ALL_TIPS D(w)',marker='.');plot_track(ax,x,tt,'#8c564b','TTT D(w)',marker='.');ax.axhline(0,color='#777',lw=.6);ax.set_xlabel('chromosome 17 position (Mb)');ax.set_ylabel('q_T - q_S');ax.set_title('Supplementary comparison: ALL_TIPS versus TTT topology contrasts',weight='bold');ax.legend(frameon=False,ncol=2);finish(fig,'house_mouse_t_complex_ALL_TIPS_vs_TTT')
 (OUT/'quartet_definition_audit.md').write_text(f'''# Quartet-definition audit\n\nThe primary topology track uses `{Q.relative_to(ROOT)}` rows with `treatment=ALL_TIPS`, preserving the frozen four-group quartet definitions and the original `q_species`, `q_t_alt`, and `q_other` support fractions. The historical column `delta_species_alt` is q_species - q_t_alt. Stage 5 writes `D = q_T - q_S = -delta_species_alt`; historical files are unchanged.\n\nThe supplementary TTT track uses `TTT_delta` from `{PAT.relative_to(ROOT)}` and applies the same sign transformation. TTT is a status-pattern-specific sampling scheme and is not interchangeable with ALL_TIPS.\n\nResolved support sums were checked at every one of {len(rows):,} windows. The mean resolved sum was {np.mean([f(x['resolved_sum']) for x in checks]):.12f}; the maximum absolute deviation of all support including unresolved from one was {max(abs(f(x['sum_all'])-1) for x in checks):.3g}. Distance values were copied from the frozen Stage 4 table without recomputation.\n''')
 (OUT/'figure_caption.md').write_text('''**Figure. Unified house-mouse t-complex topology and divergence-distance tracks.** (a) ALL_TIPS D(w) = q_T - q_S, where positive values favor the designated t-associated alternative topology. (b) A_distance is the frozen Stage 4 t-versus-standard patristic-distance contrast in ML substitution units. (c) A_relative divides that contrast by the mean non-M. spretus pairwise distance. (d) The ALL_TIPS support fractions for the species-history, t-associated alternative, and third quartet resolutions. Shading marks the four frozen t-complex inversion intervals; dashed lines mark their endpoints. The 5-kb panels preserve observed genomic gaps. The 250-kb figure uses arithmetic means within fixed coordinate bins for both topology fractions and distance summaries. The supplementary figure compares the ALL_TIPS and TTT topology sampling schemes; it does not imply that the schemes should agree.''')
 (OUT/'methods_text.md').write_text('''# Methods\n\nALL_TIPS support fractions were read from the frozen Stage 2 fixed quartet scan and filtered only by treatment. The presentation-standardized topology contrast is D(w)=q_t_alt(w)-q_species(w), the negative of the historical delta_species_alt column. TTT is retained as a supplementary status-pattern track from the frozen arrangement-pattern spatial table and receives the same sign transformation.\n\nA_distance and A_relative were copied exactly from the frozen Stage 4 ML-tree distance table. Their units are substitution distances; they are not calibrated divergence times or coalescence times. Native plots retain the original coordinates and split line segments across unusually large coordinate gaps. Matched 250-kb tables use arithmetic window means within bins anchored at 5 Mb.\n''')
 (OUT/'validation_report.md').write_text(f'''# Validation report\n\n- ALL_TIPS windows: {len(rows):,}; distance windows: {len(dist):,}; coordinate sets identical: {set(q)==set(dist)}.\n- D(w) equals -delta_species_alt at all windows: {max(abs(f(x['D_difference'])) for x in checks):.3g} maximum absolute error.\n- ALL_TIPS support values are in [0,1]: {all(x['q_range_ok']=='true' for x in checks)}.\n- Resolved support mean: {np.mean([f(x['resolved_sum']) for x in checks]):.12f}; unresolved support was retained as supplied.\n- Stage 4 A_distance and A_relative were copied exactly; maximum source-to-output differences are zero by construction.\n- ALL_TIPS versus TTT Pearson correlation: {corr:.6f}; median absolute difference: {np.median(abs(diff)):.6f}.\n\nThe 5-kb figure is the primary detailed track because it preserves local variation and genuine gaps. The 250-kb version is more readable for a main-paper overview and is the recommended main figure; the native-resolution figure is supplementary.\n''')
 manifest={'stage':'5_unified_tracks','seed':SEED,'command':'python empirical/house_mouse_t_complex/scripts/05_unified_topology_distance_tracks.py','inputs':{str(p.relative_to(ROOT)):sha(p) for p in [Q,PAT,DIST,META,INT]},'n_windows':len(rows),'primary_topology':'ALL_TIPS','historical_sign_transformation':'D_new=-delta_species_alt','distance_units':'IQ-TREE/ML patristic substitution units','bin_size_bp':250000,'existing_outputs_preserved':True};(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(json.dumps({'n_windows':len(rows),'correlation_all_tips_ttt':corr,'output':str(OUT)},indent=2))
if __name__=='__main__':main()
