#!/usr/bin/env python3
"""2-Mb physical-block bootstrap intervals for matched Dsite/Dtopology correlations."""
import csv, math, random
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[6]
BASE=ROOT/'empirical/neoaves_chr4/analysis/patterson_d'
P=BASE/'final_validation/results/nucleotide_topology_comparison.tsv'
OUT=BASE/'final_validation/results/spatial_validation.tsv'
INTERVALS=((25030000,32670000),(33510000,34470000),(44130000,56810000))
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def corr(p):
 if len(p)<3:return None
 x=[v[1] for v in p];y=[v[2] for v in p];mx=sum(x)/len(x);my=sum(y)/len(y);den=math.sqrt(sum((z-mx)**2 for z in x)*sum((z-my)**2 for z in y))
 return sum((a-mx)*(b-my) for a,b in zip(x,y))/den if den else None
def main():
 result=read(OUT); rows=read(P); rng=random.Random(20261117)
 for d in ('Columbea_primary','N61_primary','N62_primary'):
  vals=[]
  for r in rows:
   if r['design']!=d or int(r['informative_sites'])==0 or int(r['n_exact_quartet_trees'])==0:continue
   s=int(r['window_start']); inside=any(a<=s<=b for a,b in INTERVALS)
   vals.append((s,float(r['D_site_pooled']),float(r['D_topology']),inside))
  grouped=defaultdict(list)
  for v in vals:grouped[v[0]//2_000_000].append(v)
  keys=list(grouped); names=('whole_chr4','inside','outside','demeaned_inside_outside'); boot={k:[] for k in names}
  means={flag:(sum(v[1] for v in vals if v[3]==flag)/sum(v[3]==flag for v in vals),sum(v[2] for v in vals if v[3]==flag)/sum(v[3]==flag for v in vals)) for flag in (False,True)}
  for _ in range(1000):
   sample=[]
   for k in [rng.choice(keys) for __ in keys]:sample.extend(grouped[k])
   subsets={'whole_chr4':sample,'inside':[v for v in sample if v[3]],'outside':[v for v in sample if not v[3]],'demeaned_inside_outside':[(v[0],v[1]-means[v[3]][0],v[2]-means[v[3]][1],v[3]) for v in sample]}
   for name,subset in subsets.items():
    c=corr(subset)
    if c is not None and math.isfinite(c):boot[name].append(c)
  for name,bs in boot.items():
   obs={'whole_chr4':vals,'inside':[v for v in vals if v[3]],'outside':[v for v in vals if not v[3]],'demeaned_inside_outside':[(v[0],v[1]-means[v[3]][0],v[2]-means[v[3]][1],v[3]) for v in vals]}[name]
   c=corr(obs); bs.sort(); lo=bs[int(.025*(len(bs)-1))] if bs else '';hi=bs[int(.975*(len(bs)-1))] if bs else ''
   result.append(dict(kind='Dsite_Dtopology_2Mb_block_bootstrap',design=d,filter='all_sites',block_size_bp='2000000',region_class=name,D='',SE='',inside_minus_outside_D='',inside_minus_outside_SE='',n_loci=len(obs),informative_sites='',note=f'Observed r={c}; 2-Mb physical-block bootstrap percentile 95% interval [{lo}, {hi}], {len(bs)} valid replicates; 1000 resamples, seed 20261117. Descriptive only.'))
 with OUT.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter='\t');w.writeheader();w.writerows(result)
 print('appended 12 physical-block bootstrap rows')
if __name__=='__main__':main()
