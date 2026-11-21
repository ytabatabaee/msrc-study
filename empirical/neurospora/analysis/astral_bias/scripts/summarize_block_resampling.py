#!/usr/bin/env python3
"""Summarize spatial bootstrap and matched-removal ASTRAL output."""
import csv, json, math
from collections import defaultdict
from pathlib import Path
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parents[1]; RES=HERE/'results'
rows=list(csv.DictReader((RES/'block_resampling.tsv').open(),delimiter='\t'))
def num(x): return None if x in ('',None,'None') else float(x)
def q(x,p):
 x=sorted(x)
 if not x:return None
 i=(len(x)-1)*p; lo=int(i); hi=min(lo+1,len(x)-1); return x[lo]*(hi-i)+x[hi]*(i-lo)
branches=defaultdict(dict); controls=defaultdict(dict); splitsets=defaultdict(set); counts=defaultdict(dict)
for r in rows:
 d=r['design']; rep=int(r['replicate']); typ=r['analysis']; c=r['condition']
 if typ=='bootstrap_branch':
  branches[(d,rep,c)][r['split']]={'cu':num(r['CU']),'su':num(r['SU_per_site'])}
  splitsets[(d,rep,c)].add(r['split'])
 elif typ=='matched_removal_branch':
  controls[(d,rep)][r['split']]={'cu':num(r['CU']),'su':num(r['SU_per_site']),'dcu':num(r['delta_CU_vs_all']),'dsu':num(r['delta_SU_vs_all'])}
  splitsets[(d,rep,'control')].add(r['split'])
 elif typ=='matched_removal': counts[(d,rep)]={'n_retained':int(r['n_trees']),'n_removed':int(r['n_removed']),'target':int(r['target_region_n'])}
# load point estimates
point=list(csv.DictReader((RES/'branch_comparison.tsv').open(),delimiter='\t'))
summary=[]; fullsummary={}
for d in ('species_level','background_lineage'):
 shared=[r for r in point if r['design']==d and r['status']=='shared']
 for r in shared:
  split=r['split_canonical']; obs_cu=num(r['delta_CU_excluded_minus_all']); obs_su=num(r['delta_SU_excluded_minus_all'])
  boot=[]
  reps=sorted({rep for (dd,rep,cond) in branches if dd==d})
  for rep in reps:
   a=branches.get((d,rep,'boot_all'),{}).get(split); b=branches.get((d,rep,'boot_excluded'),{}).get(split)
   if a and b: boot.append((b['cu']-a['cu'],b['su']-a['su']))
  null=[v.get(split) for (dd,rep),vals in controls.items() if dd==d for v in [vals] if split in vals and vals[split]['dcu'] is not None]
  nullcu=[v['dcu'] for v in null]; nullsu=[v['dsu'] for v in null if v['dsu'] is not None]
  pcu=(1+sum(abs(x)>=abs(obs_cu) for x in nullcu))/(1+len(nullcu)) if nullcu else None
  psu=(1+sum(abs(x)>=abs(obs_su) for x in nullsu))/(1+len(nullsu)) if nullsu else None
  summary.append({'design':d,'split':split,'observed_delta_CU_excluded_minus_all':obs_cu,'bootstrap_n':len(boot),'bootstrap_CU_median':q([x[0] for x in boot],.5),'bootstrap_CU_95_low':q([x[0] for x in boot],.025),'bootstrap_CU_95_high':q([x[0] for x in boot],.975),'observed_delta_SU_excluded_minus_all':obs_su,'bootstrap_SU_median':q([x[1] for x in boot],.5),'bootstrap_SU_95_low':q([x[1] for x in boot],.025),'bootstrap_SU_95_high':q([x[1] for x in boot],.975),'matched_control_n_CU':len(nullcu),'matched_control_CU_95_low':q(nullcu,.025),'matched_control_CU_95_high':q(nullcu,.975),'matched_control_CU_abs_empirical_p':pcu,'matched_control_n_SU':len(nullsu),'matched_control_SU_95_low':q(nullsu,.025),'matched_control_SU_95_high':q(nullsu,.975),'matched_control_SU_abs_empirical_p':psu,'matched_control_removals_mean':sum(counts[(d,rep)]['n_removed'] for (dd,rep) in counts if dd==d)/max(1,sum(1 for dd,rep in counts if dd==d))})
 # topology bootstrap stability and RF, reference splits
 for cond in ('boot_all','boot_excluded'):
  vals=[splitsets.get((d,rep,cond),set()) for rep in reps]
  ref={r['split_canonical'] for r in shared}
  fullsummary[f'{d}_{cond}_exact_topology_fraction']=sum(v==ref for v in vals)/len(vals) if vals else None
  fullsummary[f'{d}_{cond}_mean_RF']=sum(len(v^ref) for v in vals)/len(vals) if vals else None
 # control topology deviation from all point tree
 ref={r['split_canonical'] for r in shared}; vals=[splitsets.get((d,rep,'control'),set()) for (dd,rep) in controls if dd==d]
 fullsummary[f'{d}_matched_control_exact_topology_fraction']=sum(v==ref for v in vals)/len(vals) if vals else None
 fullsummary[f'{d}_matched_control_mean_RF']=sum(len(v^ref) for v in vals)/len(vals) if vals else None
RES.joinpath('block_resampling_summary.json').write_text(json.dumps({'branch_results':summary,'topology_stability':fullsummary,'interpretation':'Empirical two-sided tail-area estimate uses absolute matched-removal branch shifts, with plus-one correction. Twenty replicates provide coarse uncertainty; block-width sensitivity is not included.'},indent=2)+'\n')
with RES.joinpath('block_resampling_summary.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summary)
# plot matched null distributions and observed shifts
fig,axs=plt.subplots(2,1,figsize=(12,9))
for ax,metric,obskey,title in [(axs[0],'dcu','observed_delta_CU_excluded_minus_all','CU'),(axs[1],'dsu','observed_delta_SU_excluded_minus_all','SU/site')]:
 xs=[]; labels=[]; obs=[]
 for d in ('species_level','background_lineage'):
  for r in [x for x in point if x['design']==d and x['status']=='shared']:
   split=r['split_canonical']; vals=[v[split][metric] for (dd,rep),v in controls.items() if dd==d and split in v and v[split][metric] is not None]
   if not vals:continue
   x=len(xs); xs.append(vals); labels.append(('sp: ' if d=='species_level' else 'lin: ')+split.replace('|',' / ')); obs.append(float(r['delta_'+('CU' if metric=='dcu' else 'SU')+'_excluded_minus_all']))
 if xs:
  ax.boxplot(xs,positions=range(len(xs)),vert=False,showfliers=False,patch_artist=True,boxprops={'facecolor':'#c7d8e8'},medianprops={'color':'black'})
  for i,v in enumerate(obs):ax.plot(v,i,'D',color='#d62728',label='Observed focal-region shift' if i==0 else None)
  ax.axvline(0,color='black',lw=.7);ax.set_yticks(range(len(xs)),labels,fontsize=8);ax.grid(axis='x',alpha=.2);ax.set_xlabel('Matched outside-region block removal shift');ax.set_title(title+' branch changes: observed focal exclusion vs removal controls')
fig.text(.5,.008,'Red diamonds = observed focal-region shifts; boxes = matched outside-region removal distributions (20 replicates).',ha='center',fontsize=9);fig.tight_layout(rect=[0,.035,1,1])
F=HERE/'figures';
for ext in ('png','pdf','svg'):fig.savefig(F/f'matched_removal_controls.{ext}',dpi=300,bbox_inches='tight')
plt.close(fig)
print(json.dumps({'topology_stability':fullsummary,'branch_results':summary},indent=2))
