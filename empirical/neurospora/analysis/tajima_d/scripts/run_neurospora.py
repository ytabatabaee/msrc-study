#!/usr/bin/env python3
"""Run group-specific Tajima D on the audited haploid Neurospora L10 panel."""
import csv, gzip, json, math, sys
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[5]
SHARED=ROOT/'empirical/analysis/tajimas_d';sys.path.insert(0,str(SHARED))
from tajimas import tajima_from_counts, weighted_delta, constants
GENO=ROOT/'data/neurospora/raw/neurospora/geno/cras_sit_hist_tet_L10.min20MAC1.geno.gz'
INV=ROOT/'empirical/neurospora/results/gene_tree_inventory.tsv'
ASSIGN=ROOT/'empirical/neurospora/analysis/patterson_d/results/cj_mating_type_resolution.tsv'
OUT=ROOT/'empirical/neurospora/analysis/tajima_d'
POOLED='L10_resolved_pooled'; A='L10_mat_A'; a='L10_mat_a'; GROUPS=(POOLED,A,a)
class Acc:
 def __init__(self): self.panel=0;self.complete=0;self.S=0;self.pi=0.;self.multi=0;self.missing=0

def main():
 with ASSIGN.open() as f: ar=list(csv.DictReader(f,delimiter='\t'))
 amap={r['tree_tip']:r['mating_type'] for r in ar if r['L10_membership']=='yes' and r['mating_type'] in ('A','a')}
 agroups={A:sorted([s for s,t in amap.items() if t=='A']),a:sorted([s for s,t in amap.items() if t=='a'])}
 agroups[POOLED]=sorted(agroups[A]+agroups[a])
 with gzip.open(GENO,'rt') as f: header=f.readline().split();samples=header[2:]
 if len(samples)!=26 or len(set(samples))!=26: raise ValueError('Neurospora genotype header mismatch/duplicates')
 for k,ss in agroups.items():
  if not set(ss)<=set(samples): raise ValueError(f'{k} samples missing from geno matrix')
 ix={s:i for i,s in enumerate(samples)};gix={k:[ix[s] for s in ss] for k,ss in agroups.items()}
 wins=[]
 with INV.open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   r['window']=int(r['window']);r['start']=int(r['start']);r['end']=int(r['end']);r['mid']=int(r['mid']);wins.append(r)
 per_chr=defaultdict(list)
 for i,w in enumerate(wins):per_chr[w['linkage_group']].append((w['start'],w['end'],i))
 for c in per_chr:per_chr[c].sort()
 starts={c:[x[0] for x in z] for c,z in per_chr.items()};max_end={c:max(x[1] for x in z) for c,z in per_chr.items()}
 specs={}
 for size in (50000,100000):
  res=f'{size//1000}kb';specs[res]={}
  for c,end in max_end.items():
   for lo in range(1,end+1,size):specs[res][(c,lo,lo+size-1)]=None
 specs['gene_tree']={(w['linkage_group'],w['start'],w['end']):i for i,w in enumerate(wins)}
 acc={res:defaultdict(Acc) for res in specs}
 n_records=bad_rows=n_duplicate=multiallelic_global=missing_total=0;seen=set();missing_hist=defaultdict(int)
 with gzip.open(GENO,'rt') as f:
  next(f)
  for line in f:
   x=line.split()
   if len(x)!=len(header):bad_rows+=1;continue
   chrom=x[0];pos=int(x[1]);calls=x[2:];n_records+=1
   kpos=(chrom,pos)
   if kpos in seen:n_duplicate+=1
   seen.add(kpos);valid=[v.upper() for v in calls]
   miss=sum(v not in ('A','C','G','T') for v in valid);missing_total+=miss;missing_hist[miss]+=1
   if len({v for v in valid if v in ('A','C','G','T')})>2:multiallelic_global+=1
   keys=[]
   if chrom in per_chr:
    j=bisect_right(starts[chrom],pos)-1
    if j>=0 and pos<=per_chr[chrom][j][1]:
     w=wins[per_chr[chrom][j][2]];keys.append(('gene_tree',(chrom,w['start'],w['end'])))
   for size in (50000,100000):
    res=f'{size//1000}kb';lo=((pos-1)//size)*size+1;k=(chrom,lo,lo+size-1)
    if k in specs[res]:keys.append((res,k))
   for res,wkey in keys:
    for g in GROUPS:
     ac=acc[res][(wkey,g)];ac.panel+=1
     counts=defaultdict(int);ok=True
     for j in gix[g]:
      v=valid[j]
      if v not in ('A','C','G','T'):ok=False;break
      counts[v]+=1
     if not ok:ac.missing+=1;continue
     ac.complete+=1
     if len(counts)>1:ac.S+=1
     if len(counts)>2:ac.multi+=1
     n=len(gix[g]);ac.pi+=(n*n-sum(v*v for v in counts.values()))/(n*(n-1))
 if bad_rows or n_duplicate:raise ValueError(f'bad rows {bad_rows}, duplicate coordinates {n_duplicate}')
 def finalize(wkey,res,g):
  ac=acc[res].get((wkey,g),Acc());n=len(gix[g]);S=ac.S;pi=ac.pi
  a1,_,e1,e2=constants(n);theta=S/a1;var=e1*S+e2*S*(S-1);reason='';D=None
  if n<=3:reason='n_haplotypes<=3'
  elif S==0:reason='no_segregating_sites'
  elif var<=0 or not math.isfinite(var):reason='nonpositive_or_nonfinite_classical_variance'
  else:D=(pi-theta)/math.sqrt(var)
  if ac.panel==0 and S==0:reason='no_variant_records_in_window'
  return dict(n_individuals=len(agroups[g]),ploidy='haploid',n_haplotypes=n,n_panel_variant_records=ac.panel,n_complete_variant_records=ac.complete,complete_call_fraction=ac.complete/ac.panel if ac.panel else None,missing_call_fraction=ac.missing/ac.panel if ac.panel else None,n_multiallelic_within_group=ac.multi,S=S,pi_total=pi,theta_W=theta,tajimas_D=D,classical_variance=var,D_undefined_reason=reason,callable_bases=None,pi_per_callable_base=None,thetaW_per_callable_base=None,per_base_undefined_reason='variant-only genotype matrix; no invariant callable-base mask')
 rows={res:[] for res in specs}
 for res,keys in specs.items():
  for wkey in keys:
   if res=='gene_tree':
    i=keys[wkey];w=wins[i];chrom=w['linkage_group'];lo=w['start'];hi=w['end'];mid=w['mid'];window_id=w['window']
   else:chrom,lo,hi=wkey;mid=(lo+hi)//2;window_id=f'{chrom}_{lo}_{hi}';w=None
   stats={g:finalize(wkey,res,g) for g in GROUPS}
   delta,weighted=weighted_delta(stats[POOLED]['tajimas_D'],stats[A]['tajimas_D'],stats[a]['tajimas_D'],len(agroups[A]),len(agroups[a]))
   inside=chrom=='LGI' and 1500000<=mid<=8000000
   row=dict(resolution=res,window_id=window_id,linkage_group=chrom,start=lo,end=hi,mid=mid,interval_class='LGI_suppressed_region' if inside else ('LGI_flank' if chrom=='LGI' else 'other_linkage_group'),D_Tajima_pooled=stats[POOLED]['tajimas_D'],D_Tajima_X=stats[A]['tajimas_D'],D_Tajima_Y=stats[a]['tajimas_D'],Delta_D_Tajima=delta,Delta_D_Tajima_weighted=weighted,n_pooled=stats[POOLED]['n_haplotypes'],n_X=stats[A]['n_haplotypes'],n_Y=stats[a]['n_haplotypes'],S_pooled=stats[POOLED]['S'],S_X=stats[A]['S'],S_Y=stats[a]['S'],pi_total_pooled=stats[POOLED]['pi_total'],pi_total_X=stats[A]['pi_total'],pi_total_Y=stats[a]['pi_total'],thetaW_pooled=stats[POOLED]['theta_W'],thetaW_X=stats[A]['theta_W'],thetaW_Y=stats[a]['theta_W'],panel_records_pooled=stats[POOLED]['n_panel_variant_records'],complete_records_pooled=stats[POOLED]['n_complete_variant_records'],missing_fraction_pooled=stats[POOLED]['missing_call_fraction'],multiallelic_pooled=stats[POOLED]['n_multiallelic_within_group'],callable_bases=None,pi_per_callable_base=None,thetaW_per_callable_base=None,per_base_undefined_reason='No invariant callable-site mask; variant-only matrix',undefined_pooled=stats[POOLED]['D_undefined_reason'],undefined_X=stats[A]['D_undefined_reason'],undefined_Y=stats[a]['D_undefined_reason'],gene_tree_q1=float(w['q1_species']) if w else None,gene_tree_q2=float(w['q2_other']) if w else None,gene_tree_q3=float(w['q3_hisp_A']) if w else None,D_topology=float(w['D_q3_minus_q1']) if w else None,ascertainment='published globally variable-site matrix; min20 MAC1; fixed haplotype set, complete calls within group/window')
   rows[res].append(row)
 OUT.mkdir(parents=True,exist_ok=True);(OUT/'results').mkdir(exist_ok=True);(OUT/'figures').mkdir(exist_ok=True)
 for res,rr in rows.items():
  with (OUT/'results'/f'neurospora_tajima_d_{res}.tsv').open('w') as f:
   wr=csv.DictWriter(f,fieldnames=list(rr[0]),delimiter='\t',lineterminator='\n');wr.writeheader();wr.writerows(rr)
 import numpy as np
 rng=np.random.default_rng(20521);sums=[]
 for res,rr in rows.items():
  if res not in ('50kb','100kb'):continue
  for stat in ('D_Tajima_pooled','D_Tajima_X','D_Tajima_Y','Delta_D_Tajima','Delta_D_Tajima_weighted'):
   inside=[r for r in rr if r['linkage_group']=='LGI' and r['interval_class']=='LGI_suppressed_region' and r[stat] is not None]
   outside=[r for r in rr if r['linkage_group']=='LGI' and r['interval_class']=='LGI_flank' and r[stat] is not None]
   def blockmeans(z):
    b=defaultdict(list)
    for r in z:b[(r['linkage_group'],(r['start']-1)//100000)].append(float(r[stat]))
    return [float(np.mean(v)) for v in b.values()]
   ib,ob=blockmeans(inside),blockmeans(outside);obs=float(np.mean([r[stat] for r in inside])-np.mean([r[stat] for r in outside])) if inside and outside else None;boot=[]
   if len(ib)>1 and len(ob)>1:
    for _ in range(10000):boot.append(float(rng.choice(ib,len(ib),replace=True).mean()-rng.choice(ob,len(ob),replace=True).mean()))
   sums.append(dict(resolution=res,statistic=stat,n_inside_windows=len(inside),n_outside_windows=len(outside),inside_mean=float(np.mean([r[stat] for r in inside])) if inside else None,outside_mean=float(np.mean([r[stat] for r in outside])) if outside else None,difference_inside_minus_flank=obs,block_width_bp=100000,n_inside_blocks=len(ib),n_outside_blocks=len(ob),block_bootstrap_95_low=float(np.quantile(boot,.025)) if boot else None,block_bootstrap_95_high=float(np.quantile(boot,.975)) if boot else None,interpretation='descriptive block interval; not a neutrality-test p-value'))
 with (OUT/'results/neurospora_tajima_d_region_summary.tsv').open('w') as f:
  wr=csv.DictWriter(f,fieldnames=list(sums[0]),delimiter='\t',lineterminator='\n');wr.writeheader();wr.writerows(sums)
 fig,axs=plt.subplots(4,1,figsize=(13,10),sharex=True,constrained_layout=True);labels=[('D_Tajima_pooled','Pooled resolved L10'),('D_Tajima_X','mat-A only'),('D_Tajima_Y','mat-a only'),('Delta_D_Tajima','Pooled minus mean backgrounds')]
 lgi=[r for r in rows['gene_tree'] if r['linkage_group']=='LGI']
 for ax,(col,label) in zip(axs,labels):
  ax.axvspan(1.5,8,color='#e9a63b',alpha=.20,lw=0);x=[r['mid']/1e6 for r in lgi];y=[float(r[col]) if r[col] is not None else float('nan') for r in lgi];ax.plot(x,y,marker='.',markersize=3,lw=.6,color='#245b78');ax.axhline(0,color='black',lw=.6);ax.set_ylabel(label);ax.grid(axis='y',alpha=.2)
 axs[-1].set_xlabel('LGI midpoint (Mb; 2509 study-reference coordinate axis)');fig.suptitle('Neurospora L10 Tajima’s D on the filtered SNP matrix\nHaploid isolates; group-specific complete calls; no per-base denominators',fontsize=13)
 fig.savefig(OUT/'figures/neurospora_tajima_d_LGI.pdf');fig.savefig(OUT/'figures/neurospora_tajima_d_LGI.png',dpi=300);fig.savefig(OUT/'figures/neurospora_tajima_d_LGI.svg');plt.close(fig)
 summary={'geno_records':n_records,'duplicate_coordinates':n_duplicate,'bad_rows':bad_rows,'global_multiallelic_records':multiallelic_global,'global_missing_rate':missing_total/(n_records*len(samples)),'missing_by_n':dict(missing_hist),'samples':{k:agroups[k] for k in (A,a)},'n_samples':{k:len(agroups[k]) for k in GROUPS},'maximum_window_end':max_end,'n_windows':{k:len(v) for k,v in rows.items()},'callability':'No invariant records or callable-site mask; per-base pi/theta omitted.'}
 (OUT/'results/neurospora_tajima_d_audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2));print('regional rows',len(sums))
if __name__=='__main__':main()
