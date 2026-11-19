#!/usr/bin/env python3
"""Summarize Anopheles tracks, block uncertainty, and legacy-D regression."""
import csv,math
from pathlib import Path
from collections import defaultdict
import numpy as np
ROOT=Path(__file__).resolve().parents[5];AN=ROOT/'empirical/anopheles_2la/analysis/joint_genomic_tracks';R=AN/'results';OLD=ROOT/'empirical/anopheles_2la/results'
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def fnum(x):
 try:
  z=float(x);return z if math.isfinite(z) else 0.
 except: return 0.
def ratio(x):
 if isinstance(x,dict):a=sum(v[0] for v in x.values());b=sum(v[1] for v in x.values())
 else:a,b=x
 return a/b if b else None
def jk_ratio(blocks):
 vals=[(float(a),float(b)) for a,b in blocks.values() if b>0]
 n=len(vals);num=sum(a for a,b in vals);den=sum(b for a,b in vals);est=num/den if den else None
 if est is None or n<2:return est,None,None,None,n
 leave=[(num-a)/(den-b) for a,b in vals if den!=b];m=float(np.mean(leave));se=math.sqrt((len(leave)-1)/len(leave)*sum((v-m)**2 for v in leave)) if len(leave)>1 else None
 return est,se,est-1.96*se if se is not None else None,est+1.96*se if se is not None else None,n
def jk_difference(regblocks):
 keys=set().union(*(set(x) for x in regblocks.values()));full={k:ratio(v) for k,v in regblocks.items()};delta=full['inside']-full['outside'] if full['inside'] is not None and full['outside'] is not None else None
 if delta is None or len(keys)<2:return delta,None,None,None,len(keys)
 leave=[]
 for k in keys:
  vals={reg:ratio((sum(v[0] for b,v in blocks.items() if b!=k),sum(v[1] for b,v in blocks.items() if b!=k))) for reg,blocks in regblocks.items()}
  if vals['inside'] is not None and vals['outside'] is not None:leave.append(vals['inside']-vals['outside'])
 if len(leave)<2:return delta,None,None,None,len(keys)
 m=float(np.mean(leave));se=math.sqrt((len(leave)-1)/len(leave)*sum((v-m)**2 for v in leave))
 return delta,se,delta-1.96*se,delta+1.96*se,len(keys)
def main():
 w=read(R/'anopheles_2la_fdm_D_corrected_windows.tsv');old=read(OLD/'patterson_d/window_results.tsv')
 oldix={(x['window_id'],x['test_id'],x['filter_label']):x for x in old}
 comp=[];diffs=[]
 for r in w:
  o=oldix.get((r['window_id'],r['test_id'],r['filter_label']))
  if not o:continue
  a,b=fnum(r['abba']),fnum(r['baba']);aa,bb=fnum(o['abba']),fnum(o['baba']);d=fnum(r['D_patterson_corrected']);do=fnum(o['D_patterson'])
  comp.append(dict(window_id=r['window_id'],test_id=r['test_id'],filter_label=r['filter_label'],D_corrected=r['D_patterson_corrected'],D_frozen_legacy=o['D_patterson'],absolute_D_difference=abs(d-do),ABBA_corrected=a,ABBA_legacy=aa,BABA_corrected=b,BABA_legacy=bb,n_sites_pass_corrected=r['n_sites_pass'],n_sites_pass_legacy=o['n_sites_pass'],n_biallelic_or_monomorphic=r['n_biallelic_or_monomorphic'],n_legacy_code_le1=r['n_legacy_code_le1']))
  if r['D_patterson_corrected'] not in ('',None) and o['D_patterson'] not in ('',None):diffs.append(abs(d-do))
 write(R/'anopheles_2la_corrected_vs_frozen_D.tsv',comp)
 regression=[dict(n_matched=len(diffs),n_exact_matches=sum(x<1e-12 for x in diffs),mean_absolute_difference=float(np.mean(diffs)),median_absolute_difference=float(np.median(diffs)),max_absolute_difference=max(diffs),interpretation='Non-regression is expected where legacy code treated any allele index >1 as multiallelic; corrected analysis counts observed nucleotide states at each site. Frozen result file is preserved.')]
 write(R/'anopheles_2la_D_regression_summary.tsv',regression)
 # Ratio components for f_dM and corrected D, aggregated within/outside with physical block jackknife.
 summ=[];diffrows=[]
 for test in sorted(set(r['test_id'] for r in w)):
  for filt in ('primary','MAC2'):
   subset=[r for r in w if r['test_id']==test and r['filter_label']==filt]
   regblocks={x:defaultdict(lambda:[0.,0.]) for x in ('inside','outside')}
   for x in ('inside','outside'):
    rr=[r for r in subset if r['region_class']==x]
    blocks=regblocks[x]
    for r in rr:
     bid=((int(r['start'])-1)//500000)*500000+1;b=blocks[bid]
     b[0]+=fnum(r['f_dM_numerator']);b[1]+=fnum(r['f_dM_denominator'])
    n=sum(fnum(r['f_dM_numerator']) for r in rr);d=sum(fnum(r['f_dM_denominator']) for r in rr);est,se,lo,hi,nb=jk_ratio(blocks)
    summ.append(dict(statistic='f_dM',test_id=test,filter_label=filt,region=x,n_windows=len(rr),numerator=n,denominator=d,estimate=n/d if d else '',jackknife_se=se,ci_low=lo,ci_high=hi,n_physical_blocks=nb))
   dlt,se,lo,hi,nb=jk_difference(regblocks)
   diffrows.append(dict(statistic='f_dM',test_id=test,filter_label=filt,inside_minus_outside=dlt,jackknife_se=se,ci_low=lo,ci_high=hi,n_physical_blocks=nb))
 write(R/'anopheles_2la_fdm_inside_outside.tsv',summ);write(R/'anopheles_2la_fdm_inside_outside_contrast.tsv',diffrows)
 dsum=[];ddiff=[]
 for test in sorted(set(r['test_id'] for r in w)):
  for filt in ('primary','MAC2'):
   subset=[r for r in w if r['test_id']==test and r['filter_label']==filt];regs={x:defaultdict(lambda:[0.,0.]) for x in ('inside','outside')}
   for x in ('inside','outside'):
    rr=[r for r in subset if r['region_class']==x];
    for r in rr:
     bid=((int(r['start'])-1)//500000)*500000+1;b=regs[x][bid];a=fnum(r['abba']);bb=fnum(r['baba']);b[0]+=a-bb;b[1]+=a+bb
    num=sum(fnum(r['abba'])-fnum(r['baba']) for r in rr);den=sum(fnum(r['abba'])+fnum(r['baba']) for r in rr);est,se,lo,hi,nb=jk_ratio(regs[x])
    dsum.append(dict(statistic='Patterson_D_corrected',test_id=test,filter_label=filt,region=x,n_windows=len(rr),numerator=num,denominator=den,estimate=num/den if den else '',jackknife_se=se,ci_low=lo,ci_high=hi,n_physical_blocks=nb))
   delta,se,lo,hi,nb=jk_difference(regs);ddiff.append(dict(statistic='Patterson_D_corrected',test_id=test,filter_label=filt,inside_minus_outside=delta,jackknife_se=se,ci_low=lo,ci_high=hi,n_physical_blocks=nb))
 write(R/'anopheles_2la_D_corrected_inside_outside.tsv',dsum);write(R/'anopheles_2la_D_corrected_inside_outside_contrast.tsv',ddiff)
 # Descriptive pi/dXY/FST summaries and paired physical-block delete-one uncertainty.
 pi=read(R/'anopheles_2la_pi_windows.tsv');div=read(R/'anopheles_2la_dxy_fst_windows.tsv');ds=[]
 for group in sorted(set(r['group'] for r in pi)):
  sub=[r for r in pi if r['group']==group];regs={x:defaultdict(lambda:[0.,0.]) for x in ('inside','outside')}
  for r in sub:
   x=r['region_class']
   if x not in regs or not r['callable_sites']:continue
   b=regs[x][((int(r['start'])-1)//500000)*500000+1];b[0]+=fnum(r['pi_sum']);b[1]+=int(r['callable_sites'])
  for x in ('inside','outside'):
   z=[r for r in sub if r['region_class']==x];val=sum(fnum(r['pi_sum']) for r in z)/sum(int(r['callable_sites']) for r in z) if sum(int(r['callable_sites']) for r in z) else None
   est,se,lo,hi,nb=jk_ratio(regs[x]);ds.append(dict(statistic='pi',group=group,comparison='',region=x,estimate=val,jackknife_se=se,ci_low=lo,ci_high=hi,n_windows=len(z),callable_sites=sum(int(r['callable_sites']) for r in z),effective_pairs='',n_physical_blocks=nb))
  delta,se,lo,hi,nb=jk_difference(regs);ds.append(dict(statistic='pi',group=group,comparison='',region='inside_minus_outside',estimate=delta,jackknife_se=se,ci_low=lo,ci_high=hi,n_windows='',callable_sites='',effective_pairs='',n_physical_blocks=nb))
 for pair in sorted(set(r['comparison'] for r in div)):
  sub=[r for r in div if r['comparison']==pair];regs={x:defaultdict(lambda:[0.,0.]) for x in ('inside','outside')}
  for r in sub:
   x=r['region_class']
   if x not in regs or not r['callable_sites']:continue
   bid=((int(r['start'])-1)//500000)*500000+1;b=regs[x][bid]
   # dXY is averaged over callable sites; Hudson F_ST is ratio of summed numerators/denominators.
   dsite=int(r['hudson_fst_sites']);dxy_sum=fnum(r['dxy_sum']);fstden=dxy_sum if dsite==int(r['callable_sites']) else 0.
   if dsite and fstden:
    b[0]+=fnum(r['hudson_fst'])*fstden;b[1]+=fstden
  for stat in ('dXY','Hudson_FST'):
   for x in ('inside','outside'):
    z=[r for r in sub if r['region_class']==x];sites=sum(int(r['callable_sites']) for r in z);npairs=sum(int(r['effective_allele_pairs']) for r in z)
    if stat=='dXY':
     blocks=defaultdict(lambda:[0.,0.])
     for r in z:
      bid=((int(r['start'])-1)//500000)*500000+1;b=blocks[bid];b[0]+=fnum(r['dxy_sum']);b[1]+=int(r['callable_sites'])
     den=sum(v[1] for v in blocks.values());num=sum(v[0] for v in blocks.values());value=num/den if den else None
     est,se,lo,hi,nb=jk_ratio(blocks)
    else:
     eligible=[r for r in z if int(r['hudson_fst_sites'])==int(r['callable_sites'])]
     den=sum(fnum(r['dxy_sum']) for r in eligible);num=sum(fnum(r['hudson_fst'])*fnum(r['dxy_sum']) for r in eligible);value=num/den if den else None
    if stat=='Hudson_FST':est,se,lo,hi,nb=jk_ratio(regs[x])
    ds.append(dict(statistic=stat,group='',comparison=pair,region=x,estimate=value,jackknife_se=se,ci_low=lo,ci_high=hi,n_windows=len(z),callable_sites=sites,effective_pairs=npairs,n_physical_blocks=nb))
   if stat=='Hudson_FST':
    delta,se,lo,hi,nb=jk_difference(regs);ds.append(dict(statistic=stat,group='',comparison=pair,region='inside_minus_outside',estimate=delta,jackknife_se=se,ci_low=lo,ci_high=hi,n_windows='',callable_sites='',effective_pairs='',n_physical_blocks=nb))
   elif stat=='dXY':
    db={x:defaultdict(lambda:[0.,0.]) for x in ('inside','outside')}
    for r in sub:
     if r['region_class'] in db:
      bid=((int(r['start'])-1)//500000)*500000+1;b=db[r['region_class']][bid];b[0]+=fnum(r['dxy_sum']);b[1]+=int(r['callable_sites'])
    delta,se,lo,hi,nb=jk_difference(db);ds.append(dict(statistic=stat,group='',comparison=pair,region='inside_minus_outside',estimate=delta,jackknife_se=se,ci_low=lo,ci_high=hi,n_windows='',callable_sites='',effective_pairs='',n_physical_blocks=nb))
 write(R/'anopheles_2la_diversity_inside_outside.tsv',ds)
 # Merge genetic and genealogical tracks and calculate Spearman correlations.
 topo=read(OLD/'stage2r_quartet_topology_signal.tsv');dist=read(OLD/'stage12_distance/stage12_window_distance.tsv');top={x['window_id']:x for x in topo};di={x['window_id']:x for x in dist}
 try:from scipy.stats import spearmanr
 except Exception:spearmanr=None
 correlations=[]
 metrics=[('f_dM_A_merus','quartet:f_dM:A_merus'),('f_dM_B_merus','quartet:f_dM:B_merus'),('D_Patterson_corrected_A_merus','quartet:D_patterson_corrected:A_merus'),('D_Patterson_corrected_B_merus','quartet:D_patterson_corrected:B_merus'),('pi_G_s','pi:G_s'),('pi_G_i','pi:G_i'),('pi_C_s','pi:C_s'),('pi_C_i','pi:C_i'),('dXY_Gs_Gi','dxy:primary_within_gambiae'),('dXY_Cs_Ci','dxy:primary_within_coluzzii'),('FST_Gs_Gi','fst:primary_within_gambiae'),('FST_Cs_Ci','fst:primary_within_coluzzii')]
 for label,spec in metrics:
  vals=[]
  if spec.startswith('quartet:'):
   _,kind,design=spec.split(':')
   field='f_dM' if kind=='f_dM' else 'D_patterson_corrected'
   for r in w:
    if r['test_id']!=design or r['filter_label']!='primary':continue
    tr=top.get(r['window_id']);dt=di.get(r['window_id'])
    if tr and dt:vals.append((r['window_id'],fnum(r[field]),tr,dt))
  elif spec.startswith('pi:'):
   group=spec.split(':')[1];g={x['window_id']:x for x in pi if x['group']==group}
   for wid,x in g.items():
    if int(x['callable_sites'])>=25000 and wid in top and wid in di:vals.append((wid,fnum(x['pi_per_callable_site']),top[wid],di[wid]))
  elif spec.startswith(('dxy:','fst:')):
   kind,pair=spec.split(':');key='dxy_per_callable_site' if kind=='dxy' else 'hudson_fst';g={x['window_id']:x for x in div if x['comparison']==pair}
   for wid,x in g.items():
    calls=int(x['hudson_fst_sites'] if kind=='fst' else x['callable_sites'])
    if calls>=25000 and wid in top and wid in di:vals.append((wid,fnum(x[key]),top[wid],di[wid]))
  if not vals:continue
  for outcome in ('q_species','q_arrangement','q_third','D','mean_internal_branch','tree_diameter','mean_pair_distance'):
   pairs=[(v[1],fnum(v[2].get(outcome,'')) if outcome in v[2] else fnum(v[3].get(outcome,''))) for v in vals]
   pairs=[x for x in pairs if math.isfinite(x[0]) and math.isfinite(x[1])]
   if len(pairs)<4:continue
   rho,pv=(spearmanr([x[0] for x in pairs],[x[1] for x in pairs]) if spearmanr else (float('nan'),float('nan')))
   correlations.append(dict(popgen_track=label,tree_or_branch_measure=outcome,n_windows=len(pairs),spearman_rho=float(rho),p_value=float(pv),interpretation='descriptive spatial association; not causal; p-values are not adjusted for spatial autocorrelation'))
 if correlations:write(R/'anopheles_2la_track_correlations.tsv',correlations)
 low=[]
 for r in pi:
  n=int(r['callable_sites']);low.append(dict(window_id=r['window_id'],metric='pi',group=r['group'],comparison='',callable_sites=n,window_bp=int(r['end'])-int(r['start'])+1,callability_fraction=n/(int(r['end'])-int(r['start'])+1),plot_eligible=n>=25000,threshold_callable_sites=25000))
 for r in div:
  n=int(r['hudson_fst_sites']);low.append(dict(window_id=r['window_id'],metric='dXY_FST',group='',comparison=r['comparison'],callable_sites=n,window_bp=int(r['end'])-int(r['start'])+1,callability_fraction=n/(int(r['end'])-int(r['start'])+1),plot_eligible=n>=25000,threshold_callable_sites=25000))
 write(R/'anopheles_2la_low_callability_windows.tsv',low)
 print('matched D rows',len(comp),'mean abs delta',float(np.mean(diffs)),'max',max(diffs),'summaries',len(summ),'diversity rows',len(ds),'correlations',len(correlations))
if __name__=='__main__':main()
