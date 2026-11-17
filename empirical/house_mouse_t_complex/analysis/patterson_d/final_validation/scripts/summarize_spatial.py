#!/usr/bin/env python3
"""Assemble spatial correlations and block uncertainty with explicit units."""
from pathlib import Path
import csv, math, random
import numpy as np
P=Path(__file__).resolve().parents[2]; OUT=P/'final_validation/results'
def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rows):
 fields=['design','analysis','metric','region','resolution_bp','estimate','ci_low','ci_high','n_windows','n_blocks','circular_shift_p_descriptive','note']
 with open(p,'w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
rows=[]
# Whole-span source correlations, primary support threshold only.
for x in read(P/'exact_quartet_comparison/results/spatial_correlations.tsv'):
 if x['region']=='matched_all' and x['minimum_informative_sites']=='1' and x['minimum_quartet_resolution_fraction']=='' and int(x['resolution_bp']) in (50000,100000,250000,500000):
  rows.append({'design':x['design'],'analysis':'whole_span','metric':'Pearson_r_Dsite_vs_Dtopology','region':'matched_all','resolution_bp':x['resolution_bp'],'estimate':x['pearson_r'],'ci_low':x['block_bootstrap_low'],'ci_high':x['block_bootstrap_high'],'n_windows':x['n_windows'],'n_blocks':math.ceil(32_055_000/500_000),'circular_shift_p_descriptive':x['circular_shift_null_abs_r_p'],'note':'Whole analyzed tree/SNP matched span; original 500-kb block bootstrap. Circular shifts are descriptive only under irregular coverage.'})
# Raw correlations stratified into inversion and collinear coverage; outside has very few blocks.
for x in read(P/'exact_quartet_comparison/results/spatial_correlations.tsv'):
 if x['region'] in ('inversion','outside_inversion') and x['minimum_informative_sites']=='1' and x['minimum_quartet_resolution_fraction']=='' and int(x['resolution_bp']) in (50000,100000,250000,500000):
  rows.append({'design':x['design'],'analysis':'raw_region_subset','metric':'Pearson_r_Dsite_vs_Dtopology','region':x['region'],'resolution_bp':x['resolution_bp'],'estimate':x['pearson_r'],'ci_low':x['block_bootstrap_low'],'ci_high':x['block_bootstrap_high'],'n_windows':x['n_windows'],'n_blocks':'','circular_shift_p_descriptive':x['circular_shift_null_abs_r_p'],'note':'Region-specific raw correlation. Outside-inversion coverage is sparse; do not interpret its circular-shift p-value as general spatial inference.'})
tracks=read(P/'exact_quartet_comparison/results/matched_nucleotide_genealogy_tracks.tsv')
def corr(x,y):
 x=np.asarray(x,float);y=np.asarray(y,float)
 if len(x)<3 or np.std(x)==0 or np.std(y)==0:return float('nan')
 return float(np.corrcoef(x,y)[0,1])
def residualize(rows,group_field):
 xs=np.array([float(r['D_site']) for r in rows]);ys=np.array([float(r['D_topology']) for r in rows]);gs=[r[group_field] for r in rows]
 xr=xs.copy();yr=ys.copy()
 for g in set(gs):
  ix=[i for i,a in enumerate(gs) if a==g]
  xr[ix]-=xs[ix].mean();yr[ix]-=ys[ix].mean()
 return corr(xr,yr)
def bootstrap_rows(data,center_group,seed,reps=1000):
 blocks=sorted(set(int(x['start'])//500000 for x in data))
 if len(blocks)<4:return float('nan'),float('nan')
 by={b:[x for x in data if int(x['start'])//500000==b] for b in blocks};rng=random.Random(seed);vals=[]
 for _ in range(reps):
  chosen=[rng.choice(blocks) for i in blocks];sample=[]
  for b in chosen:sample.extend(by[b])
  v=residualize(sample,center_group)
  if math.isfinite(v):vals.append(v)
 if len(vals)<reps*.5:return float('nan'),float('nan')
 return float(np.quantile(vals,.025)),float(np.quantile(vals,.975))
# Region-centered point correlations and physical 500-kb cluster bootstrap.
for d in sorted(set(x['design'] for x in tracks)):
 for res in (50000,100000,250000,500000):
  base=[x for x in tracks if x['design']==d and int(x['resolution_bp'])==res]
  groups=[('all_region_centered',base,'region'),('inversion_only_centered',[x for x in base if x['region'].startswith('In(17)')],'region'),('collinear_flanks_gaps_centered',[x for x in base if not x['region'].startswith('In(17)')],'region')]
  for name,data,gfield in groups:
   if len(data)<4:continue
   est=residualize(data,gfield);lo,hi=bootstrap_rows(data,gfield,seed=20261117+res+len(d)+len(name))
   rows.append({'design':d,'analysis':'region_centered','metric':'Pearson_r_Dsite_vs_Dtopology_residual','region':name,'resolution_bp':res,'estimate':est,'ci_low':lo,'ci_high':hi,'n_windows':len(data),'n_blocks':len(set(int(x['start'])//500000 for x in data)),'circular_shift_p_descriptive':'','note':'Residualized D_site and D_topology within source interval class; cluster bootstrap samples 500-kb blocks. Outside-only estimates have few blocks and remain descriptive.'})
# Independent delete-block jackknife from the new corrected contribution table.
for x in read(OUT/'block_uncertainty_validation.tsv'):
 rows.append({'design':x['design'],'analysis':'delete_block_jackknife','metric':'D_inside_minus_D_outside','region':x['scope'],'resolution_bp':x['block_bp'],'estimate':x['estimate_inside_minus_outside_D'],'ci_low':x['normal_95_low'],'ci_high':x['normal_95_high'],'n_windows':'','n_blocks':x['delete_blocks'],'circular_shift_p_descriptive':'','note':'Normal 95% interval from delete-block jackknife SE; D recomputed after each physical block deletion.'})
# Published internal bootstrap results for corrected nucleotide and exact topology.
for x in read(P/'polarization_validation/results/corrected_inside_outside_block_uncertainty.tsv'):
 rows.append({'design':x['design'],'analysis':'existing_block_bootstrap','metric':x['region'],'region':'full_contig_inside_minus_complement','resolution_bp':x['block_bp'],'estimate':x['point_difference'],'ci_low':x['block_bootstrap_low'],'ci_high':x['block_bootstrap_high'],'n_windows':x['inside_windows']+x['complement_windows'],'n_blocks':x['inside_blocks']+x['complement_blocks'],'circular_shift_p_descriptive':'','note':'Corrected fixed-REF+ALT nucleotide block bootstrap; source output retained.'})
for x in read(P/'exact_quartet_comparison/results/genomic_block_uncertainty.tsv'):
 rows.append({'design':x['design'],'analysis':'existing_block_bootstrap','metric':x['statistic'],'region':x['region'],'resolution_bp':x['block_bp'],'estimate':x['point_difference'],'ci_low':x['block_bootstrap_low'],'ci_high':x['block_bootstrap_high'],'n_windows':'','n_blocks':str(x['blocks_region'])+'/'+str(x['blocks_reference']),'circular_shift_p_descriptive':'','note':'Exact matched coverage; source block bootstrap; induced quartets are not independent units.'})
write(OUT/'spatial_validation.tsv',rows)
print('wrote',len(rows),'spatial rows')
