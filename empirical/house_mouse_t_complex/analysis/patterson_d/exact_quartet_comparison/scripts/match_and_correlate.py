#!/usr/bin/env python3
"""Join finished polarization and exact-tree tables; aggregate and test spatial correspondence."""
import csv,math,statistics
from pathlib import Path
from collections import defaultdict
from run_analysis import GROUPS,INV,classify,association,write
BASE=Path(__file__).resolve().parents[1];P=BASE.parent;POL=P/'polarization_validation'/'results';OUT=BASE/'results'
def read(p):
 with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def main():
 corrected=read(POL/'corrected_window_D.tsv');trees=read(OUT/'matched_gene_tree_quartets.tsv');oldtop=read(P/'results'/'patterson_d_topology_comparison.tsv');top={(x['design'],int(x['start'])):x for x in oldtop};nt={(x['design'],int(x['start'])):x for x in corrected};matched=[]
 for t in trees:
  x=nt.get((t['design'],int(t['start'])));z=top.get((t['design'],int(t['start'])),{});matched.append({'chromosome':'chr17','start':int(t['start']),'end':int(t['end']),'design':t['design'],'ABBA':x['ABBA'] if x else '', 'BABA':x['BABA'] if x else '', 'D_site':x['D'] if x else '', 'informative_sites':x['informative_sites'] if x else '', 'q1':t['q1'],'q2':t['q2'],'q3':t['q3'],'D_topology':t['D_topology'],'eligible_quartets':t['eligible_resolved_quartets'],'possible_quartets':t['possible_induced_quartets'],'unresolved_quartets':t['unresolved_quartets'],'missing_samples':t['missing_samples'],'mean_internal_branch_length':t['mean_internal_branch_length'],'n_P1':t['n_P1'],'n_P2':t['n_P2'],'n_P3':t['n_P3'],'n_O':t['n_O'],'requested_P1':t['requested_P1'],'requested_P2':t['requested_P2'],'requested_P3':t['requested_P3'],'requested_O':t['requested_O'],'coverage_P1':int(t['n_P1'])/int(t['requested_P1']),'coverage_P2':int(t['n_P2'])/int(t['requested_P2']),'coverage_P3':int(t['n_P3'])/int(t['requested_P3']),'coverage_O':int(t['n_O'])/int(t['requested_O']),'recombination_inversion_class':classify(int(t['start']),int(t['end'])),'source_recombination_class':z.get('recombination_class','NA')})
 write(OUT/'matched_window_dataset.tsv',matched)
 # Equal-locus genealogy weights: pool nucleotide counts; average each local-tree q once per locus.
 multi=[]
 for res in (50000,100000,250000,500000):
  bins=defaultdict(list)
  for r in matched:
   if r['D_site'] not in ('',None) and r['D_topology'] not in ('',None):bins[(r['design'],r['start']//res*res)].append(r)
  for (d,st),xs in sorted(bins.items()):
   a=sum(float(x['ABBA']) for x in xs);b=sum(float(x['BABA']) for x in xs);cov=statistics.mean(int(x['eligible_quartets'])/int(x['possible_quartets']) for x in xs if int(x['possible_quartets']))
   mean=lambda field:statistics.mean(float(x[field]) for x in xs if x[field] not in ('',None))
   multi.append({'design':d,'resolution_bp':res,'start':st,'end':min(st+res,37055000),'ABBA':a,'BABA':b,'D_site':(a-b)/(a+b) if a+b else '', 'informative_sites':sum(int(x['informative_sites']) for x in xs),'q1':mean('q1'),'q2':mean('q2'),'q3':mean('q3'),'D_topology':mean('D_topology'),'n_loci':len(xs),'eligible_quartets_sum_diagnostic_only':sum(int(x['eligible_quartets']) for x in xs),'mean_locus_quartet_resolution_fraction':cov,'mean_internal_branch_length':mean('mean_internal_branch_length'),'region':classify(st,min(st+res,37055000))})
 write(OUT/'matched_nucleotide_genealogy_tracks.tsv',multi)
 stats=[]
 for res in (50000,100000,250000,500000):
  for d in GROUPS:
   base=[x for x in multi if x['resolution_bp']==res and x['design']==d]
   for min_sites in (1,5,10,25):
    xs=[x for x in base if int(x['informative_sites'])>=min_sites]
    for region in ('matched_all','inversion','outside_inversion'):
     if region=='matched_all':sub=xs
     elif region=='inversion':sub=[x for x in xs if any(a<=x['start'] and x['end']<=b for _,a,b in INV)]
     else:sub=[x for x in xs if not any(x['start']<b and x['end']>a for _,a,b in INV)]
     v=association(sub);stats.append({'design':d,'resolution_bp':res,'minimum_informative_sites':min_sites,'minimum_quartet_resolution_fraction':'','region':region,'n_windows':len(sub),'pearson_r':v[0],'spearman_rho':v[1],'sign_agreement':v[2],'block_bootstrap_low':v[3],'block_bootstrap_high':v[4],'circular_shift_null_abs_r_p':v[5],'block_level_pearson_r':v[6]})
   for min_cov in (.8,.95,.999):
    xs=[x for x in base if float(x['mean_locus_quartet_resolution_fraction'])>=min_cov and int(x['informative_sites'])>=1];v=association(xs);stats.append({'design':d,'resolution_bp':res,'minimum_informative_sites':1,'minimum_quartet_resolution_fraction':min_cov,'region':'matched_all','n_windows':len(xs),'pearson_r':v[0],'spearman_rho':v[1],'sign_agreement':v[2],'block_bootstrap_low':v[3],'block_bootstrap_high':v[4],'circular_shift_null_abs_r_p':v[5],'block_level_pearson_r':v[6]})
 write(OUT/'spatial_correlations.tsv',stats)
 print(f'Wrote {len(matched)} matched windows; {len(multi)} resolution tracks; {len(stats)} correlation/sensitivity rows')
if __name__=='__main__':main()
