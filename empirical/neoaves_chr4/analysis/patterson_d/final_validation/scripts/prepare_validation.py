#!/usr/bin/env python3
"""Build additive Neoaves D validation tables from frozen locus-level outputs."""
import csv, math, random
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[6]
BASE=ROOT/'empirical/neoaves_chr4/analysis/patterson_d'
OUT=BASE/'final_validation/results'
INTERVALS=[('PNAS_1',25030000,32670000),('PNAS_2',33510000,34470000),('PNAS_3',44130000,56810000)]
def rows(path):
 with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(path,rs):
 path.parent.mkdir(parents=True,exist_ok=True)
 if not rs:return
 with open(path,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rs[0]),delimiter='\t');w.writeheader();w.writerows(rs)
def D(a,b):return (a-b)/(a+b) if a+b else None
def jk_d(vals,inside_value,bs=1000000):
 g=defaultdict(lambda:[0,0])
 for st,ins,a,b in vals:
  if ins==inside_value:g[st//bs][0]+=a;g[st//bs][1]+=b
 aa=sum(v[0] for v in g.values());bb=sum(v[1] for v in g.values());theta=D(aa,bb)
 reps=[D(aa-v[0],bb-v[1]) for v in g.values()]
 reps=[x for x in reps if x is not None]
 if theta is None or len(reps)<2:return theta,None
 m=sum(reps)/len(reps);return theta,math.sqrt((len(reps)-1)/len(reps)*sum((x-m)**2 for x in reps))
def jk_delta(vals,bs=1000000):
 # vals are (start,inside,a,b); each genomic block is retained, never site-resampled.
 g=defaultdict(lambda:[0,0,0,0])
 for s,inside,a,b in vals:
  x=g[s//bs];k=0 if inside else 2;x[k]+=a;x[k+1]+=b
 ia=sum(v[0] for v in g.values());ib=sum(v[1] for v in g.values());oa=sum(v[2] for v in g.values());ob=sum(v[3] for v in g.values())
 di,do=D(ia,ib),D(oa,ob);diff=di-do if di is not None and do is not None else None
 reps=[]
 for v in g.values():
  x,y=D(ia-v[0],ib-v[1]),D(oa-v[2],ob-v[3])
  if x is not None and y is not None:reps.append(x-y)
 if len(reps)<2:return diff,None
 m=sum(reps)/len(reps);se=math.sqrt((len(reps)-1)/len(reps)*sum((x-m)**2 for x in reps))
 return diff,se

def main():
 summary=rows(BASE/'results/patterson_d_summary.tsv')
 locus=rows(BASE/'results/patterson_d_locus_results.tsv')
 manuscript=rows(BASE/'results/manuscript_results_table.tsv')
 designs=rows(BASE/'config/frozen_quartets.tsv')
 # final pooled values, 1-Mb genomic delete-block SE; separate H0s are explicit.
 final=[]
 for d in [x['design'] for x in designs]:
  for filt in ('all','transversions_only'):
   sub=[r for r in summary if r['design']==d and r['filter']==filt and r['region']=='all_PNAS' and r['block_size']=='1000000']
   ins=next((r for r in sub if r['region_class']=='inside'),None);out=next((r for r in sub if r['region_class']=='outside'),None)
   if not ins or not out:continue
   final.append(dict(design=d,filter=filt,interval_scope='all_PNAS_union',inside_D=ins['D'],inside_se=ins['jackknife_se'],inside_CI95_low=float(ins['D'])-1.96*float(ins['jackknife_se']),inside_CI95_high=float(ins['D'])+1.96*float(ins['jackknife_se']),outside_D=out['D'],outside_se=out['jackknife_se'],outside_CI95_low=float(out['D'])-1.96*float(out['jackknife_se']),outside_CI95_high=float(out['D'])+1.96*float(out['jackknife_se']),inside_minus_outside_D=ins['inside_minus_outside_D'],inside_minus_outside_se=ins['inside_minus_outside_se'],delta_CI95_low=float(ins['inside_minus_outside_D'])-1.96*float(ins['inside_minus_outside_se']),delta_CI95_high=float(ins['inside_minus_outside_D'])+1.96*float(ins['inside_minus_outside_se']),test_Doutside_equals_0='separate null; use outside_CI95',test_Dinside_equals_Doutside='separate null; use delta_CI95',n_inside_loci=ins['n_loci'],n_outside_loci=out['n_loci'],inside_informative_sites=ins['informative_sites'],outside_informative_sites=out['informative_sites'],ABBA_inside=ins['ABBA'],BABA_inside=ins['BABA'],ABBA_outside=out['ABBA'],BABA_outside=out['BABA']))
 # Derive transition-only counts exactly as all-sites minus transversions for matching loci.
 bykey={(r['design'],r['locus_id'],r['filter']):r for r in locus}
 for d in [x['design'] for x in designs]:
  trans=[]
  for key,r in bykey.items():
   dd,lid,filt=key
   if dd!=d or filt!='all':continue
   tx=bykey.get((dd,lid,'transversions_only'))
   if not tx:continue
   a=int(r['ABBA'])-int(tx['ABBA']);b=int(r['BABA'])-int(tx['BABA'])
   if a<0 or b<0:raise ValueError('transversion counts exceed all-site counts')
   trans.append((int(r['start']),r['class']=='inside',a,b))
  if not trans:continue
  inside,ise=jk_d(trans,True);outside,ose=jk_d(trans,False);delta,dse=jk_delta(trans)
  ia=sum(a for s,i,a,b in trans if i);ib=sum(b for s,i,a,b in trans if i);oa=sum(a for s,i,a,b in trans if not i);ob=sum(b for s,i,a,b in trans if not i)
  final.append(dict(design=d,filter='transitions_only',interval_scope='all_PNAS_union',inside_D=inside,inside_se=ise,inside_CI95_low=inside-1.96*ise if ise is not None else '',inside_CI95_high=inside+1.96*ise if ise is not None else '',outside_D=outside,outside_se=ose,outside_CI95_low=outside-1.96*ose if ose is not None else '',outside_CI95_high=outside+1.96*ose if ose is not None else '',inside_minus_outside_D=delta,inside_minus_outside_se=dse,delta_CI95_low=delta-1.96*dse if dse is not None else '',delta_CI95_high=delta+1.96*dse if dse is not None else '',test_Doutside_equals_0='separate null; use outside_CI95',test_Dinside_equals_Doutside='separate null; use delta_CI95',n_inside_loci=sum(i for s,i,a,b in trans),n_outside_loci=sum(not i for s,i,a,b in trans),inside_informative_sites=ia+ib,outside_informative_sites=oa+ob,ABBA_inside=ia,BABA_inside=ib,ABBA_outside=oa,BABA_outside=ob))
 write(OUT/'final_results.tsv',final)
 # taxa/outgroup occupancy and sensitivity; groups not extrapolated from focal outcomes.
 tax=[]
 loci_by=defaultdict(set)
 for r in locus:
  if r['filter']=='all':loci_by[r['design']].add(r['locus_id'])
 for d in designs:
  did=d['design'];x=next((r for r in final if r['design']==did and r['filter']=='all'),None)
  tax.append(dict(design=did,relationship=d['relationship'],P1=d['P1'],P2=d['P2'],P3=d['P3'],outgroup=d['O'],outgroup_role=d['outgroup_role'],loci_with_all_four_taxa=len(loci_by[did]),chr4_manifest_loci=6247,quartet_occupancy_fraction=round(len(loci_by[did])/6247,6),inside_loci=x['n_inside_loci'] if x else '',outside_loci=x['n_outside_loci'] if x else '',inside_D=x['inside_D'] if x else '',outside_D=x['outside_D'] if x else '',delta_D=x['inside_minus_outside_D'] if x else '',sensitivity_type='Galliform outgroup alternative only' if 'alt_outgroup' in did else 'frozen representative quartet',independence='Correlated designs share taxa and loci; not independent replicates',ingroup_representative_sensitivity='Not evaluated: only frozen ingroup quartets are represented in saved nucleotide counts; additional ingroup taxa require a new alignment scan.',biological_rationale=d['rationale']))
 write(OUT/'taxon_sensitivity.tsv',tax)
 # interval-specific and exclusion estimates from saved per-locus counts; block SE by source coordinate.
 reg=[]
 for d in ['Columbea_primary','N61_primary','N62_primary']:
  rr=[r for r in locus if r['filter']=='all' and r['design']==d]
  scenarios=[]
  for name,lo,hi in INTERVALS:scenarios.append((name,lambda s,lo=lo,hi=hi:lo<=s<=hi))
  scenarios += [('all_PNAS',lambda s:any(a<=s<=b for _,a,b in INTERVALS)),('exclude_PNAS_2_inside_PNAS_1+3',lambda s:any(a<=s<=b for n,a,b in INTERVALS if n!='PNAS_2')),('exclude_PNAS_3_inside_PNAS_1+2',lambda s:any(a<=s<=b for n,a,b in INTERVALS if n!='PNAS_3'))]
  for name,pred in scenarios:
   vals=[(int(r['start']),pred(int(r['start'])),int(r['ABBA']),int(r['BABA'])) for r in rr]
   ia=sum(a for s,i,a,b in vals if i);ib=sum(b for s,i,a,b in vals if i);oa=sum(a for s,i,a,b in vals if not i);ob=sum(b for s,i,a,b in vals if not i)
   delta,se=jk_delta(vals)
   reg.append(dict(design=d,scenario=name,inside_D=D(ia,ib),outside_D=D(oa,ob),inside_minus_outside_D=delta,delta_se_1Mb=se,inside_ABBA=ia,inside_BABA=ib,outside_ABBA=oa,outside_BABA=ob,inside_informative_sites=ia+ib,outside_informative_sites=oa+ob,interpretation='Block-jackknife descriptive interval; PNAS_2 has low counts.' if name=='PNAS_2' else 'Pooled region estimate; outside means all complementary chr4 loci.',coordinate_axis='Stiller source window-start coordinate used by frozen GalGal4 predicate'))
  # Boundary movement on the same source coordinate used by the frozen predicate.
  for mode,amount in [('expand',250000),('shrink',250000),('expand',1000000),('shrink',1000000)]:
   intervals=[]
   for _,lo,hi in INTERVALS:
    intervals.append((max(0,lo-amount),hi+amount) if mode=='expand' else (lo+amount,hi-amount))
   pred=lambda s,iv=intervals:any(a<=s<=b for a,b in iv)
   vals=[(int(r['start']),pred(int(r['start'])),int(r['ABBA']),int(r['BABA'])) for r in rr]
   ia=sum(a for s,i,a,b in vals if i);ib=sum(b for s,i,a,b in vals if i);oa=sum(a for s,i,a,b in vals if not i);ob=sum(b for s,i,a,b in vals if not i)
   delta,se=jk_delta(vals)
   reg.append(dict(design=d,scenario=f'{mode}_{amount//1000}kb',inside_D=D(ia,ib),outside_D=D(oa,ob),inside_minus_outside_D=delta,delta_se_1Mb=se,inside_ABBA=ia,inside_BABA=ib,outside_ABBA=oa,outside_BABA=ob,inside_informative_sites=ia+ib,outside_informative_sites=oa+ob,interpretation='Boundary perturbation; descriptive because GalGal4-to-GalGal6 liftover is not made here.',coordinate_axis='Stiller source window-start coordinate used by frozen GalGal4 predicate'))
 write(OUT/'regional_robustness.tsv',reg)
 # Recompute pooled sums from locus counts and compare to stored pooled rows; also check coordinate mask.
 checks=[]
 for d in [r['design'] for r in designs]:
  for filt in ('all','transversions_only'):
   sub=[r for r in locus if r['design']==d and r['filter']==filt]
   for cls in ('inside','outside'):
    chosen=[r for r in sub if r['class']==cls]
    aa=sum(int(r['ABBA']) for r in chosen);bb=sum(int(r['BABA']) for r in chosen)
    old=next((r for r in summary if r['design']==d and r['filter']==filt and r['region']=='all_PNAS' and r['region_class']==cls and r['block_size']=='1000000'),None)
    checks.append(dict(check='pooled_ABBA_BABA_identity',design=d,filter=filt,scope=cls,raw_locus_ABBA=aa,raw_locus_BABA=bb,stored_ABBA=old['ABBA'] if old else '',stored_BABA=old['BABA'] if old else '',status='PASS' if old and aa==int(old['ABBA']) and bb==int(old['BABA']) else 'FAIL'))
   bad=0
   for r in sub:
    s=int(r['start']); expected=any(a<=s<=b for _,a,b in INTERVALS)
    if (r['class']=='inside')!=expected:bad+=1
   checks.append(dict(check='coordinate_membership_vs_frozen_predicate',design=d,filter=filt,scope='all_chr4',raw_locus_ABBA='',raw_locus_BABA='',stored_ABBA='',stored_BABA='',status='PASS' if not bad else f'FAIL:{bad}'))
 for d in [x['design'] for x in designs]:
  for cls in ('inside','outside'):
   allr=next(x for x in final if x['design']==d and x['filter']=='all');txr=next(x for x in final if x['design']==d and x['filter']=='transversions_only');trr=next(x for x in final if x['design']==d and x['filter']=='transitions_only')
   aa=sum(int(x[f'ABBA_{cls}']) for x in (txr,trr));bb=sum(int(x[f'BABA_{cls}']) for x in (txr,trr))
   checks.append(dict(check='transition_transversion_partition_pooled',design=d,filter='transitions_only+transversions_only',scope=cls,raw_locus_ABBA=aa,raw_locus_BABA=bb,stored_ABBA=allr[f'ABBA_{cls}'],stored_BABA=allr[f'BABA_{cls}'],status='PASS' if aa==int(allr[f'ABBA_{cls}']) and bb==int(allr[f'BABA_{cls}']) else 'FAIL'))
 write(OUT/'estimator_validation.tsv',checks)
 # Spatial summaries at the existing physical-block scales and exact matched gene-tree bins.
 spatial=[]
 for r in summary:
  if r['region']=='all_PNAS' and r['filter'] in ('all','transversions_only'):
   spatial.append(dict(kind='inside_outside_block_jackknife',design=r['design'],filter=r['filter'],block_size_bp=r['block_size'],region_class=r['region_class'],D=r['D'],SE=r['jackknife_se'],inside_minus_outside_D=r['inside_minus_outside_D'],inside_minus_outside_SE=r['inside_minus_outside_se'],n_loci=r['n_loci'],informative_sites=r['informative_sites'],note='Existing 500-kb/1-Mb/2-Mb genomic delete-block results; blocks are physical, sites are not independent replicates.'))
 # 1-Mb exact matched site/tree Pearson correlations and 2-Mb block bootstrap CI.
 comp=rows(OUT/'nucleotide_topology_comparison.tsv')
 for d in ('Columbea_primary','N61_primary','N62_primary'):
  rr=[r for r in comp if r['design']==d and int(r['informative_sites'])>0 and int(r['n_exact_quartet_trees'])>0]
  pairs=[]
  for r in rr:
   s=int(r['window_start']); inside=any(a<=s<=b for _,a,b in INTERVALS)
   pairs.append((s,float(r['D_site_pooled']),float(r['D_topology']),inside))
  def pearson(p):
   if len(p)<3:return None
   x=[v[1] for v in p];y=[v[2] for v in p];mx=sum(x)/len(x);my=sum(y)/len(y);den=math.sqrt(sum((z-mx)**2 for z in x)*sum((z-my)**2 for z in y))
   return sum((a-mx)*(b-my) for a,b in zip(x,y))/den if den else None
  for label,ps in [('whole_chr4',pairs),('inside',[p for p in pairs if p[3]]),('outside',[p for p in pairs if not p[3]])]:
   spatial.append(dict(kind='Dsite_Dtopology_correlation_1Mb',design=d,filter='all_sites',block_size_bp='1000000',region_class=label,D='',SE='',inside_minus_outside_D='',inside_minus_outside_SE='',n_loci=len(ps),informative_sites='',note=f'Pearson r={pearson(ps)}; exact same-locus quartet trees and nucleotide site counts; descriptive.'))
  # partial correlation after demeaning inside/outside categories
  means={flag:(sum(p[1] for p in pairs if p[3]==flag)/sum(p[3]==flag for p in pairs),sum(p[2] for p in pairs if p[3]==flag)/sum(p[3]==flag for p in pairs)) for flag in (False,True)}
  residual=[(s,x-means[i][0],y-means[i][1],i) for s,x,y,i in pairs]
  rr2=pearson([(s,x,y,i) for s,x,y,i in residual])
  spatial.append(dict(kind='Dsite_Dtopology_partial_correlation_1Mb',design=d,filter='all_sites',block_size_bp='1000000',region_class='demeaned_inside_outside',D='',SE='',inside_minus_outside_D='',inside_minus_outside_SE='',n_loci=len(pairs),informative_sites='',note=f'Pearson residual r={rr2}; descriptive, not a significance test.'))
 write(OUT/'spatial_validation.tsv',spatial)
 print('wrote final_results, taxon_sensitivity, regional_robustness, estimator_validation, spatial_validation')
if __name__=='__main__':main()
