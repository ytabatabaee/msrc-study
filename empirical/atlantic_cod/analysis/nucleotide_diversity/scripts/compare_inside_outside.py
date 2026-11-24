#!/usr/bin/env python3
"""Compare inversion-associated cohorts inside and outside suppression.

Inside the inversion, groups are interpreted as homozygous arrangement
backgrounds. Outside, these same sample cohorts are population-grouped only;
there is no claim that alleles remain linked to the inversion arrangement.
"""
import csv,gzip
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[5]
BASE=ROOT/'empirical/atlantic_cod/analysis/nucleotide_diversity'
R=BASE/'results';FIG=BASE/'figures'
VCF=ROOT/'data/atlantic_cod/raw/zenodo/local_downloads/gadus_morhua_outside_supergenes_snapp.vcf.gz'
WINDOWS=ROOT/'data/atlantic_cod/processed/cod_window_trees.tsv'
REGIONS=ROOT/'data/atlantic_cod/metadata/region_manifest.tsv'
from estimate_tracks import read_tsv,read_assignment,call_ok,site_stats,hudson_fst

def main():
 windows=[r for r in read_tsv(WINDOWS) if r['lg'] in {'LG01','LG02','LG07','LG12'}]
 regions={r['lg']:r for r in read_tsv(REGIONS) if r['lg'] in {'LG01','LG02','LG07','LG12'}}
 arrang=read_assignment()
 statuses={}
 for w in windows:
  lg=w['lg'];a=int(w['start']);b=int(w['end']);reg=regions[lg]
  inside=a>=int(reg['start']) and b<=int(reg['end'])
  overlap=a<=int(reg['end']) and b>=int(reg['start'])
  statuses[(lg,w['window_id'])]='inside' if inside else 'boundary' if overlap else 'outside'
 groups={}
 for s in []:pass
 for lg in ('LG01','LG02','LG07','LG12'):
  groups[lg]={'X':[],'Y':[]}
  for pop,(state,conf,source) in arrang.items():
   prefix=pop
   bg='X' if state=='ancestral' else 'Y'
   groups[lg][bg].append(prefix)
 # VCF sample ID is prefix plus an individual suffix; map to population state.
 accum={}
 for w in windows:
  lg=w['lg'];wid=w['window_id'];scope=statuses[(lg,wid)]
  accum[(lg,wid)]={'scope':scope,'total':0,'eligible':0,'sumX':0.,'sumY':0.,'sumD':0.,
                   'callX':0,'callY':0,'expectX':0,'expectY':0,'copiesX':0,'copiesY':0}
 with gzip.open(VCF,'rt') as f:
  for line in f:
   if line.startswith('#CHROM'):
    header=line.rstrip().split('\t')[9:]
    idx={s:i for i,s in enumerate(header)}
    sample_groups={}
    for lg in groups:
     sample_groups[lg]={bg:[idx[s] for s in header if s.startswith(pop) for pop in []] for bg in []}
     sample_groups[lg]={'X':[],'Y':[]}
     for s,j in idx.items():
      if not s.startswith('Gadmor_'):continue
      pop=s[:10]
      match=arrang.get(pop)
      if match:
       bg='X' if match[0]=='ancestral' else 'Y'
       sample_groups[lg][bg].append(j)
    continue
   if line.startswith('#'):continue
   a=line.rstrip().split('\t');lg=a[0]
   if lg not in groups or len(a[3])!=1 or len(a[4].split(','))!=1 or len(a[4])!=1:continue
   if a[6] not in ('PASS','.') or a[5] in ('','.') or float(a[5])<30:continue
   pos=int(a[1]);hit=None
   for w in windows:
    if w['lg']==lg and int(w['start'])<=pos<=int(w['end']) and statuses[(lg,w['window_id'])]=='outside':
     hit=w['window_id'];break
   if hit is None:continue
   z=accum[(lg,hit)];z['total']+=1;fmt=a[8].split(':')
   vals=[call_ok(v,fmt) for v in a[9:]]
   gx=[vals[j] for j in sample_groups[lg]['X']];gy=[vals[j] for j in sample_groups[lg]['Y']]
   cx=[x for x in gx if x is not None];cy=[x for x in gy if x is not None]
   z['expectX']+=2*len(gx);z['expectY']+=2*len(gy);z['callX']+=2*len(cx);z['callY']+=2*len(cy)
   nx=2*len(cx);ny=2*len(cy);z['copiesX']+=nx;z['copiesY']+=ny
   stats=site_stats(sum(cx),nx,sum(cy),ny)
   if stats is None:continue
   px,py,d=stats;z['sumX']+=px;z['sumY']+=py;z['sumD']+=d;z['eligible']+=1
 # Inside values use the arrangement-specific scanner (same genotypes / filters / thresholds).
 with (R/'atlantic_cod_pi_dxy_fst_windows.tsv').open() as f: existing={(r['chromosome'],r['window_id']):r for r in csv.DictReader(f,delimiter='\t')}
 out=[]
 for w in windows:
  lg=w['lg'];wid=w['window_id'];key=(lg,wid);scope=statuses[key];base=accum[key]
  row={'chromosome':lg,'window_id':wid,'start':w['start'],'end':w['end'],'comparison_scope':scope}
  if scope=='inside':
   e=existing[key]
   for dst,src in [('pi_X_group_SNP','pi_X_SNP'),('pi_Y_group_SNP','pi_Y_SNP'),('dXY_group_SNP','dXY_SNP'),('FST_group_Hudson','FST_Hudson')]:row[dst]=e[src]
   for dst,src in [('joint_eligible_SNPs','joint_eligible_SNPs'),('mean_haploid_copies_X','mean_haploid_copies_X'),('mean_haploid_copies_Y','mean_haploid_copies_Y'),('missingness_X','missingness_X'),('missingness_Y','missingness_Y')]:row[dst]=e[src]
   row['interpretation']='homozygous arrangement-background comparison'
  elif scope=='outside':
   n=base['eligible'];row['joint_eligible_SNPs']=n
   row['mean_haploid_copies_X']=base['copiesX']/base['total'] if base['total'] else 'NA'
   row['mean_haploid_copies_Y']=base['copiesY']/base['total'] if base['total'] else 'NA'
   row['missingness_X']=1-base['callX']/base['expectX'] if base['expectX'] else 'NA'
   row['missingness_Y']=1-base['callY']/base['expectY'] if base['expectY'] else 'NA'
   if n:
    px=base['sumX']/n;py=base['sumY']/n;d=base['sumD']/n
    row.update({'pi_X_group_SNP':px,'pi_Y_group_SNP':py,'dXY_group_SNP':d,'FST_group_Hudson':hudson_fst(px,py,d)})
   else:row.update({'pi_X_group_SNP':'NA','pi_Y_group_SNP':'NA','dXY_group_SNP':'NA','FST_group_Hudson':'NA'})
   row['interpretation']='population-group comparison only; no linked arrangement assignment outside inversion'
  else:
   row.update({'pi_X_group_SNP':'NA','pi_Y_group_SNP':'NA','dXY_group_SNP':'NA','FST_group_Hudson':'NA','joint_eligible_SNPs':'NA','mean_haploid_copies_X':'NA','mean_haploid_copies_Y':'NA','missingness_X':'NA','missingness_Y':'NA','interpretation':'boundary-mixed window excluded'})
  out.append(row)
 fields=list(out[0])
 with (R/'inside_outside_cohort_group_windows.tsv').open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=fields,delimiter='\t');wr.writeheader();wr.writerows(out)
 # Equal-window means for interpretable comparison; report variable site counts separately.
 summary=[]
 for lg in ('LG01','LG02','LG07','LG12'):
  g=[r for r in out if r['chromosome']==lg]
  for scope in ('inside','outside'):
   q=[r for r in g if r['comparison_scope']==scope and r['pi_X_group_SNP']!='NA']
   summary.append({'lg':lg,'scope':scope,'valid_windows':len(q),'mean_pi_X_group_SNP':_mean(q,'pi_X_group_SNP'),'mean_pi_Y_group_SNP':_mean(q,'pi_Y_group_SNP'),'mean_dXY_group_SNP':_mean(q,'dXY_group_SNP'),'mean_FST_group_Hudson':_mean(q,'FST_group_Hudson'),'min_eligible_SNPs':min([int(x['joint_eligible_SNPs']) for x in q]) if q else 'NA','max_eligible_SNPs':max([int(x['joint_eligible_SNPs']) for x in q]) if q else 'NA','interpretation':'arrangement backgrounds' if scope=='inside' else 'population groups; arrangement linkage not assumed'})
 with (R/'inside_outside_cohort_group_summary.tsv').open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t');wr.writeheader();wr.writerows(summary)
 contrast=[]
 for lg in ('LG01','LG02','LG07','LG12'):
  a=next(x for x in summary if x['lg']==lg and x['scope']=='inside');b=next(x for x in summary if x['lg']==lg and x['scope']=='outside')
  contrast.append({'lg':lg,'inside_windows':a['valid_windows'],'outside_windows':b['valid_windows'],**{('outside_minus_inside_'+c):float(b['mean_'+c])-float(a['mean_'+c]) for c in ('pi_X_group_SNP','pi_Y_group_SNP','dXY_group_SNP','FST_group_Hudson')}})
 with (R/'inside_outside_cohort_contrast.tsv').open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=list(contrast[0]),delimiter='\t');wr.writeheader();wr.writerows(contrast)
 _figure(out,regions)
 with (R/'inside_outside_comparison.md').open('w') as f:
  f.write('# Inside versus outside comparison\n\n')
  f.write('The primary arrangement-specific table continues to report pi_X, pi_Y, dXY, and FST only within fully contained inversion windows. This additional comparison uses the same per-LG ancestral-associated (17 individuals) and derived-associated (7 individuals) cohorts throughout the chromosome. Inside suppressed regions, those cohorts correspond to homozygous arrangement backgrounds based on the published study selection. Outside, values are population-grouped summaries only; recombination means their alleles cannot be assigned to the inversion arrangement. The outside results therefore describe group differentiation, not continued X/Y chromosomal backgrounds.\n\n')
  f.write('| LG | inside windows | outside windows | inside πX | outside πX group | inside πY | outside πY group | inside dXY | outside dXY group | inside FST | outside FST group |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n')
  for lg in ('LG01','LG02','LG07','LG12'):
   a=next(x for x in summary if x['lg']==lg and x['scope']=='inside');b=next(x for x in summary if x['lg']==lg and x['scope']=='outside')
   f.write(f"| {lg} | {a['valid_windows']} | {b['valid_windows']} | {float(a['mean_pi_X_group_SNP']):.5f} | {float(b['mean_pi_X_group_SNP']):.5f} | {float(a['mean_pi_Y_group_SNP']):.5f} | {float(b['mean_pi_Y_group_SNP']):.5f} | {float(a['mean_dXY_group_SNP']):.5f} | {float(b['mean_dXY_group_SNP']):.5f} | {float(a['mean_FST_group_Hudson']):.5f} | {float(b['mean_FST_group_Hudson']):.5f} |\n")
  f.write('\nAll estimates use the same biallelic PASS/dot, QUAL≥30, GQ≥20, DP≥8 filters and at least four observed haploid copies in each group at a site. The same eligible SNP set defines within-group pi and between-group dXY at each window. Values are SNP-conditional because the source VCFs are variant-only and no callable invariant mask is available. Boundary-crossing windows remain NA. Summary values are equal-window averages; do not interpret the contrast causally because population structure and ecotype are confounded with the background cohorts outside inversions. Numeric outside-minus-inside differences are in inside_outside_cohort_contrast.tsv. Full window values are in inside_outside_cohort_group_summary.tsv and inside_outside_cohort_group_windows.tsv.\n')

def _mean(rows,col):
 a=[float(x[col]) for x in rows if x[col] not in ('NA','nan','')]
 return sum(a)/len(a) if a else 'NA'
def _figure(rows,regions):
 import matplotlib.pyplot as plt
 from matplotlib.lines import Line2D
 FIG.mkdir(parents=True,exist_ok=True)
 cols=['pi_X_group_SNP','pi_Y_group_SNP','dXY_group_SNP','FST_group_Hudson']
 titles=['πX-associated group','πY-associated group','dXY between groups','Hudson FST']
 fig,axs=plt.subplots(4,4,figsize=(16,10),sharex='col')
 for i,lg in enumerate(('LG01','LG02','LG07','LG12')):
  g=[r for r in rows if r['chromosome']==lg];reg=regions[lg]
  for j,(col,title) in enumerate(zip(cols,titles)):
   ax=axs[i,j];inside=[r for r in g if r['comparison_scope']=='inside' and r[col]!='NA'];outside=[r for r in g if r['comparison_scope']=='outside' and r[col]!='NA']
   xi=[(int(r['start'])+int(r['end']))/2/1e6 for r in inside];yi=[float(r[col]) for r in inside]
   xo=[(int(r['start'])+int(r['end']))/2/1e6 for r in outside];yo=[float(r[col]) for r in outside]
   ax.plot(xi,yi,'.-',color='#215a7a',ms=2.5,lw=.7);ax.plot(xo,yo,'.',color='#cf6d3d',ms=2.5,alpha=.65)
   ax.axvspan(int(reg['start'])/1e6,int(reg['end'])/1e6,color='#e6ad4b',alpha=.2);ax.grid(alpha=.15)
   if i==0:ax.set_title(title)
   if j==0:ax.set_ylabel(lg)
   if i==3:ax.set_xlabel('Position (Mb)')
 fig.legend([Line2D([0],[0],color='#215a7a',marker='o',lw=.8),Line2D([0],[0],color='#cf6d3d',marker='o',linestyle='')],['Inside: arrangement backgrounds','Outside: population cohorts'],loc='upper center',bbox_to_anchor=(.5,.965),ncol=2,frameon=False)
 fig.suptitle('Inside/outside cohort comparison; X=ancestral-associated, Y=derived-associated\nOutside values are population-grouped, not arrangement-linked',y=1.02);fig.tight_layout(rect=[0,0,1,.94]);fig.savefig(FIG/'atlantic_cod_inside_outside_cohort_comparison.png',dpi=180,bbox_inches='tight');fig.savefig(FIG/'atlantic_cod_inside_outside_cohort_comparison.pdf',bbox_inches='tight');plt.close(fig)

if __name__=='__main__':main()
