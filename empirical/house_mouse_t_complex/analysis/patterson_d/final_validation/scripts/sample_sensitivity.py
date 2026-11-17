#!/usr/bin/env python3
"""Indexed leave-one-out and balanced P1/P2 subsampling on exact tree windows."""
import os, random, csv, math
from pathlib import Path
import certifi, numpy as np, pysam
os.environ['SSL_CERT_FILE']=certifi.where();os.environ['REQUESTS_CA_BUNDLE']=certifi.where()
P=Path(__file__).resolve().parents[2]
ROOT=P.parents[2]
URL='https://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz'
import sys
sys.path.insert(0,str(P/'scripts'))
from run_patterson_d import GROUPS
OUT=P/'final_validation/results'
MATCH=P/'exact_quartet_comparison/results/matched_window_dataset.tsv'
INV=[(5747331,6957919),(7435308,13016510),(14068229,18051507),(19401000,37055000)]
def complete(gt):return gt is not None and len(gt)==2 and None not in gt
def region(st):
 return 'inside' if any(a<=st and st+5000<=b for a,b in INV) else 'outside'
def dval(a,b):return (a-b)/(a+b) if a+b else float('nan')
def main():
 trees={}
 for x in csv.DictReader(open(MATCH),delimiter='\t'):
  trees.setdefault(x['design'],{})[int(x['start'])]=x['recombination_inversion_class']
 v=pysam.VariantFile(URL); idx={s:i for i,s in enumerate(v.header.samples)}
 designs={d:tuple(tuple(sorted(g)) for g in gs) for d,gs in GROUPS.items()}
 needed=set(s for gs in designs.values() for g in gs for s in g)
 if not needed<=set(v.header.samples):raise RuntimeError('frozen samples missing from VCF')
 # Generate fixed deterministic balanced P1/P2 subsets. Keep the smaller group in full.
 subs={}; reps=100
 for d,(g1,g2,g3,og) in designs.items():
  rng=random.Random(20261117+list(designs).index(d)); k=min(len(g1),len(g2)); large=g1 if len(g1)>len(g2) else g2
  small=g2 if len(g1)>len(g2) else g1
  picks=[rng.sample(large,k) for _ in range(reps)]
  # Fixed seed subset matrix: rows are members of the larger P1/P2 group.
  selection=np.zeros((len(large),reps),dtype=bool); pos={s:i for i,s in enumerate(large)}
  for j,ss in enumerate(picks):
   for s in ss:selection[pos[s],j]=True
  subs[d]=(small,large,selection,k)
 # Accumulators for leave-one-out values, and balanced sample D vectors.
 loo={}; bal={}
 for d,(g1,g2,g3,og) in designs.items():
  loo[d]={}
  for g in (g1,g2,g3):
   for s in g:loo[d][s]={z:[0.,0.] for z in ('all','inside','outside')}
  bal[d]={z:np.zeros((reps,2),dtype=float) for z in ('all','inside','outside')}
 # Fetch only the bounded 5.0–37.055 Mb range using indexed 500-kb requests.
 for lo in range(5_000_000,37_055_000,500_000):
  hi=min(lo+500_000,37_055_000)
  for rec in v.fetch('chr17',lo,hi):
   if len(rec.alleles)!=2 or any(a is None or len(a)!=1 or a not in 'ACGT' for a in rec.alleles):continue
   st=(rec.pos-1)//5000*5000
   active=[d for d in designs if st in trees.get(d,{})]
   if not active:continue
   for d in active:
    g1,g2,g3,og=designs[d]
    gt={s:rec.samples[s].get('GT') for s in (*g1,*g2,*g3,*og)}
    ogts=[gt[s] for s in og]
    if not all(complete(x) for x in ogts):continue
    alleles={a for x in ogts for a in x}
    if alleles=={0}:flip=False
    elif alleles=={1}:flip=True
    else:continue
    arrays=[]; vals=[]; valid=[]
    for g in (g1,g2,g3):
     vv=np.array([sum(gt[s]) if complete(gt[s]) else 0 for s in g],dtype=float)
     mm=np.array([complete(gt[s]) for s in g],dtype=bool)
     arrays.append(vv);valid.append(mm);vals.append(vv.sum()/(2*mm.sum()) if mm.any() else np.nan)
    if any(not math.isfinite(x) for x in vals):continue
    p1,p2,p3=([1-x for x in vals] if flip else vals)
    rg='inside' if trees[d][st].startswith('In(17)') else 'outside'
    # Per-individual leave-one-out contributions in the population role of that sample.
    for role,g in enumerate((g1,g2,g3)):
     vv,mm=arrays[role],valid[role]
     for j,s in enumerate(g):
      n=int(mm.sum())-int(mm[j]); alt=vv.sum()-vv[j]
      if n<1 or not mm[j]:continue
      ps=list((1-p1,1-p2,1-p3) if flip else (p1,p2,p3))
      ps[role]=alt/(2*n)
      if flip:ps=[1-x for x in ps]
      a=(1-ps[0])*ps[1]*ps[2];b=ps[0]*(1-ps[1])*ps[2]
      for z in ('all',rg):loo[d][s][z][0]+=a;loo[d][s][z][1]+=b
    # Balanced P1/P2 subsampling, retaining all P3 members.
    small,large,selection,k=subs[d]; small_role=0 if len(g1)<=len(g2) else 1; large_role=1-small_role
    def freq_subset(g,role):
     if role==small_role:
      vv,mm=arrays[role],valid[role]
      return vv.sum()/(2*mm.sum()) if mm.any() else np.full(reps,np.nan)
     vv,mm=arrays[role],valid[role];mat=selection
     counts=mm.astype(float)@mat
     sums=vv@mat
     return np.divide(sums,2*counts,out=np.full(reps,np.nan),where=counts>0)
    pp=list(vals)
    pp[small_role]=freq_subset((g1,g2)[small_role],small_role)
    pp[large_role]=freq_subset((g1,g2)[large_role],large_role)
    if flip:pp=[1-x for x in pp]
    aa=(1-pp[0])*pp[1]*pp[2];bb=pp[0]*(1-pp[1])*pp[2]
    for z in ('all',rg):
     bal[d][z][:,0]+=np.nan_to_num(aa,nan=0.);bal[d][z][:,1]+=np.nan_to_num(bb,nan=0.)
 # Emit a tidy per-individual table and replicate distributions, not inferential p-values.
 loo_rows=[]
 for d,ss in loo.items():
  for sample,rr in ss.items():
   for z,(a,b) in rr.items():
    loo_rows.append({'design':d,'sample_excluded':sample,'region':z,'ABBA':a,'BABA':b,'D_leave_one_out':dval(a,b),'unit':'sample-level descriptive; not block uncertainty'})
 with open(OUT/'leave_one_individual_out.tsv','w') as f:
  w=csv.DictWriter(f,fieldnames=list(loo_rows[0]),delimiter='\t');w.writeheader();w.writerows(loo_rows)
 bal_rows=[]
 for d,zs in bal.items():
  for z,ab in zs.items():
   ds=np.array([dval(a,b) for a,b in ab]);finite=ds[np.isfinite(ds)]
   bal_rows.append({'design':d,'region':z,'replicates':len(finite),'seed':20261117+list(designs).index(d),'balanced_n_P1_P2':min(len(designs[d][0]),len(designs[d][1])),'D_median':np.median(finite),'D_2.5pct':np.quantile(finite,.025),'D_97.5pct':np.quantile(finite,.975),'D_min':finite.min(),'D_max':finite.max(),'unit':'deterministic sample subsampling; not a population-structure correction'})
 with open(OUT/'balanced_population_subsampling.tsv','w') as f:
  w=csv.DictWriter(f,fieldnames=list(bal_rows[0]),delimiter='\t');w.writeheader();w.writerows(bal_rows)
 print(f'wrote {len(loo_rows)} leave-one-out rows and {len(bal_rows)} balanced subsampling summaries')
if __name__=='__main__':main()
