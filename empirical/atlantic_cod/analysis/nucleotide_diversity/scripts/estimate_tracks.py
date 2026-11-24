#!/usr/bin/env python3
"""Estimate within-arrangement diversity in cod's four inversion regions.

X denotes the ancestral arrangement; Y denotes derived. Sample-level labels
follow the published homokaryotype-selected cohort and published population /
ecotype mapping, not a forced PCA cluster assignment. Windows crossing frozen
recombination boundaries or outside inversions remain NA.
"""
import csv,gzip,math
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[5]
BASE=ROOT/'empirical/atlantic_cod/analysis/nucleotide_diversity'
VCFDIR=ROOT/'data/atlantic_cod/raw/zenodo/local_downloads'
WINDOWS=ROOT/'data/atlantic_cod/processed/cod_window_trees.tsv'
REGIONS=ROOT/'data/atlantic_cod/metadata/region_manifest.tsv'
POPS=ROOT/'data/atlantic_cod/processed/population_arrangements.tsv'
OUT=BASE/'results'
LGFILES={lg:VCFDIR/f'gadus_morhua_supergenes_{lg.lower()}_snapp.vcf.gz' for lg in ['LG01','LG02','LG07','LG12']}
MIN_GQ=20;MIN_DP=8;MIN_QUAL=30


def read_tsv(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))

def gt_dosage(v):
 g=v.split(':',1)[0]
 if '.' in g:return None
 a=g.replace('|','/').split('/')
 if len(a)!=2 or any(x not in ('0','1') for x in a):return None
 return int(a[0])+int(a[1])

def call_ok(v,fmt):
 d=dict(zip(fmt,v.split(':')));dos=gt_dosage(v)
 try:gq=float(d.get('GQ','.'));dp=float(d.get('DP','.'))
 except ValueError:return None
 return dos if dos is not None and gq>=MIN_GQ and dp>=MIN_DP else None

def corrected_pi(k,n):
 if n<2:return None
 return 2*k*(n-k)/(n*(n-1))
def site_stats(kx,nx,ky,ny):
 if nx<4 or ny<4:return None
 px=kx/nx;py=ky/ny
 return corrected_pi(kx,nx),corrected_pi(ky,ny),px*(1-py)+py*(1-px)
def hudson_fst(pi_x,pi_y,d_xy):
 return 'NA' if d_xy==0 else 1-(pi_x+pi_y)/(2*d_xy)

def read_assignment():
 states={}
 for r in read_tsv(POPS):
  name=r['tree_label'].replace('_spc','')
  for lg in LGFILES:
   state=r[f'{lg}_state']; conf=r[f'{lg}_confidence']
   if state in ('ancestral','derived'):
    states[name]=(state,conf,r[f'{lg}_source'])
 return states

def main():
 windows=[r for r in read_tsv(WINDOWS) if r['lg'] in LGFILES]
 regions={r['lg']:r for r in read_tsv(REGIONS) if r['lg'] in LGFILES}
 assign=read_assignment();sample_rows=[];sites_by_lg={};qc=[]
 for lg,path in LGFILES.items():
  sites=[];sample_header=[];counts={};total_records=passed=0;gq_sum=dp_sum=gq_n=dp_n=0;raw_miss=0
  with gzip.open(path,'rt') as fh:
   for line in fh:
    if line.startswith('#CHROM'):
     sample_header=line.rstrip().split('\t')[9:]
     for s in sample_header:
      pop=s[:10] if s.startswith('Gadmor_') else ''
      if not pop:continue
      st=assign.get(pop)
      if not st:continue
      state,conf,source=st;bg='X' if state=='ancestral' else 'Y'
      counts.setdefault(bg,[]).append(s)
      sample_rows.append({'lg':lg,'sample_id':s,'population':pop,'arrangement_state':state,'arrangement_genotype':'XX' if bg=='X' else 'YY','assignment_confidence':'high' if conf=='direct_outgroup_alignment' else 'moderate','evidence':f'published homokaryotype-selected cohort; population/ecotype-to-arrangement mapping; orientation: {conf}; {source}'})
     continue
    if line.startswith('#'):continue
    a=line.rstrip().split('\t')
    if a[0]!=lg or len(a[3])!=1 or len(a[4].split(','))!=1 or len(a[4])!=1:continue
    total_records+=1
    if a[6] not in ('PASS','.') or a[5] in ('.','') or float(a[5])<MIN_QUAL:continue
    passed+=1;fmt=a[8].split(':');gts=[]
    for v in a[9:]:
     d=dict(zip(fmt,v.split(':')))
     try:
      q=float(d.get('GQ','.'));dp=float(d.get('DP','.'))
      gq_sum+=q;gq_n+=1;dp_sum+=dp;dp_n+=1
     except ValueError:pass
     gt=call_ok(v,fmt)
     if gt is None:raw_miss+=1
     gts.append(gt)
    sites.append((int(a[1]),gts))
  sites_by_lg[lg]=(sample_header,counts,sites)
  qc.append({'lg':lg,'n_vcf_samples':len(sample_header),'n_codon_samples_assigned':sum(map(len,counts.values())),'n_X_individuals':len(counts.get('X',[])),'n_Y_individuals':len(counts.get('Y',[])),'variant_records_in_region':total_records,'PASS_biallelic_SNPs_QUAL30':passed,'GQ_DP_filtered_missing_fraction':raw_miss/(passed*len(sample_header)) if passed and sample_header else 'NA','mean_GQ':gq_sum/gq_n if gq_n else 'NA','mean_DP':dp_sum/dp_n if dp_n else 'NA'})
 # individual-class table, de-duplicate any repeated records by (LG,sample)
 unique={(r['lg'],r['sample_id']):r for r in sample_rows};sample_rows=list(unique.values())
 OUT.mkdir(parents=True,exist_ok=True)
 with (OUT/'individual_arrangement_assignments.tsv').open('w',newline='') as f:
  cols=['lg','sample_id','population','arrangement_state','arrangement_genotype','assignment_confidence','evidence'];w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(sample_rows)
 with (OUT/'genotype_qc.tsv').open('w',newline='') as f:
  cols=list(qc[0]);w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(qc)
 out=[]
 for w in windows:
  lg=w['lg'];a=int(w['start']);b=int(w['end']);r=regions[lg];fully=a>=int(r['start']) and b<=int(r['end'])
  row={'chromosome':lg,'window_id':w['window_id'],'start':a,'end':b,'pi_X_SNP': 'NA','pi_Y_SNP':'NA','dXY_SNP':'NA','FST_Hudson':'NA','joint_eligible_SNPs':'NA','callable_sites_X':'NA','callable_sites_Y':'NA','n_individuals_X':'NA','n_individuals_Y':'NA','mean_haploid_copies_X':'NA','mean_haploid_copies_Y':'NA','missingness_X':'NA','missingness_Y':'NA','normalization':'SNP-conditional','assignment_status':'outside_or_boundary_overlapping_suppressed_region' if not fully else 'no_jointly_eligible_sites','inside_suppressed_recombination':str(fully).lower(),'overlaps_suppressed_recombination':str(a<=int(r['end']) and b>=int(r['start'])).lower()}
  if fully:
   header,groups,sites=sites_by_lg[lg];ix={s:i for i,s in enumerate(header)};xs=groups.get('X',[]);ys=groups.get('Y',[]);xi=[ix[s] for s in xs];yi=[ix[s] for s in ys]
   sums=[0.,0.,0.];n_sites=0;calledX=calledY=expectedX=expectedY=0;copiesX=copiesY=0;window_site_count=0
   for pos,gts in sites:
    if not a<=pos<=b:continue
    window_site_count+=1
    gx=[gts[j] for j in xi];gy=[gts[j] for j in yi]
    expectedX+=2*len(gx);expectedY+=2*len(gy)
    cx=[z for z in gx if z is not None];cy=[z for z in gy if z is not None]
    calledX+=2*len(cx);calledY+=2*len(cy)
    nx=2*len(cx);ny=2*len(cy);copiesX+=nx;copiesY+=ny
    kx=sum(cx);ky=sum(cy);stats=site_stats(kx,nx,ky,ny)
    if stats is None:continue
    hx,hy,d=stats
    sums[0]+=hx;sums[1]+=hy;sums[2]+=d;n_sites+=1
   row.update({'n_individuals_X':len(xs),'n_individuals_Y':len(ys),'joint_eligible_SNPs':n_sites,'missingness_X':1-calledX/expectedX if expectedX else 'NA','missingness_Y':1-calledY/expectedY if expectedY else 'NA','mean_haploid_copies_X':copiesX/window_site_count if window_site_count else 'NA','mean_haploid_copies_Y':copiesY/window_site_count if window_site_count else 'NA','callable_sites_X':'NA','callable_sites_Y':'NA'})
   if n_sites:
    pxw=sums[0]/n_sites;pyw=sums[1]/n_sites;dw=sums[2]/n_sites
    fst=hudson_fst(pxw,pyw,dw)
    row.update({'pi_X_SNP':pxw,'pi_Y_SNP':pyw,'dXY_SNP':dw,'FST_Hudson':fst,'assignment_status':'estimated'})
  out.append(row)
 cols=list(out[0]);
 with (OUT/'atlantic_cod_pi_dxy_fst_windows.tsv').open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=cols,delimiter='\t');wr.writeheader();wr.writerows(out)
 # validation arithmetic, inside/outside availability, and concise report
 ins=[r for r in out if r['inside_suppressed_recombination']=='true'];valid=[r for r in ins if r['pi_X_SNP']!='NA']
 with (OUT/'track_validation_report.md').open('w') as f:
  f.write('# Arrangement-specific diversity tracks: validation and interpretation\n\n')
  f.write('The SNAPP VCFs contain diploid GT/DP/GQ calls for 29 samples, of which 24 are Atlantic cod. Matschiner et al. report that the 22 additional cod individuals were selected because each was homokaryotypic across all four supergenes; the two reference/assembly individuals were complementary homozygotes. Population/ecotype labels are mapped to ancestral/derived arrangements using the paper’s classification convention and the frozen orientation table. X=ancestral and Y=derived. This supports two homozygous classes in this deliberately ascertained cohort; it cannot test a third heterokaryotype class because such individuals were selected out.\n\n')
  f.write('The assignment file records source/confidence per individual and LG. PCA remains diagnostic and is not used to force classes. LG12 orientation is less certain because ancestral/derived orientation was inferred demographically, and both Bornholm samples have a published ~275-kb heterozygous introgressed tract inside LG12; local ancestry there is mixed even though the chromosome arrangement assignment is ancestral.\n\n')
  f.write(f'For windows fully inside frozen suppressed-recombination regions, {len(valid)}/{len(ins)} windows have at least one jointly eligible SNP. Windows outside or overlapping region boundaries are NA because X/Y ancestry is not extended across recombining sequence. SNP panel is all biallelic PASS (or missing FILTER) records with QUAL≥30 from the shared all-sample VCF; genotype copies require GQ≥20 and DP≥8. A site is jointly eligible when each background contributes at least four called haploid copies (at least two diploid individuals). Within-background π uses the unbiased finite-sample gene-diversity correction. dXY uses the mean cross-background pairwise difference at the same sites. Each quantity is the mean across that same eligible SNP set; it is SNP-conditional, not per-base. No callable invariant mask was located. Hudson FST is calculated from those matched window estimates and is not truncated when negative; it is NA when dXY is zero.\n\n')
  f.write('No sample-level arrangement heterozygotes are used. The 250-kb output preserves all windows and NA values where arrangement labels are biologically undefined. See genotype_qc.tsv and individual_genotype_classification.tsv for call QC and assignments.\n')

def _selfcheck():
 assert abs(corrected_pi(3,4)-.5)<1e-12
 assert corrected_pi(0,1) is None
 assert site_stats(1,2,1,2) is None
 assert hudson_fst(.1,.1,0)=='NA'

if __name__=='__main__':_selfcheck();main()
