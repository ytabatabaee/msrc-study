#!/usr/bin/env python3
"""Additive f_dM scan for frozen Atlantic-cod trios and Dsuite VCF inputs."""
from __future__ import annotations
import csv,gzip,json,math,sys
from collections import defaultdict
from pathlib import Path
ROOT=Path(__file__).resolve().parents[5]
sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT/'empirical/atlantic_cod/analysis/patterson_d/scripts'))
from patterson_d_core import parse_gt,read_tree_windows,window_id
from empirical.scripts.joint_genomic_tracks.fdm import fdm_site,jackknife_ratio
DATA=ROOT/'data/atlantic_cod';AN=ROOT/'empirical/atlantic_cod/analysis/joint_genomic_tracks';RAW=DATA/'raw/zenodo/local_downloads'
WINDOW=250_000
TREE_WINDOWS={}

def tsv(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def freq(gt_strings, require_all=True):
 alleles=[];n=0
 for s in gt_strings:
  g=parse_gt(s or '.')
  if g is None:
   if require_all:return None,n
   continue
  n+=1;alleles.extend(g)
 return (sum(alleles)/len(alleles),n) if alleles else (None,n)
def scan(path,designs,region_label):
 by_lg={d['lg']:d for d in designs}
 acc=defaultdict(lambda:{'num':0.,'den':0.,'abba':0.,'baba':0.,'n_sites':0,'informative':0,'records':0,'missing':0,'out_missing':0,'out_het':0,'out_fixed_ref':0,'out_fixed_alt':0})
 samples=[];ix={}
 with gzip.open(path,'rt') as fh:
  for line in fh:
   if line.startswith('##'):continue
   if line.startswith('#CHROM'):
    samples=line.rstrip('\n').split('\t')[9:];ix={x:i for i,x in enumerate(samples)}
    all_ids={sid for d in designs for key in ('P1_sample_ids','P2_sample_ids','P3_sample_ids','O_sample_ids') for sid in d[key].split(',')}
    missing=sorted(all_ids-set(ix))
    if missing:raise ValueError(f'{path.name} missing frozen sample IDs {missing}')
    continue
   fields=line.rstrip('\n').split('\t'); chrom=fields[0]
   if chrom not in by_lg:continue
   design=by_lg[chrom];p1=design['P1_sample_ids'].split(',');p2=design['P2_sample_ids'].split(',');p3=design['P3_sample_ids'].split(',');outid=design['O_sample_ids'].split(',')
   try:pos=int(fields[1])
   except ValueError:continue
   wid,start,end=window_id(chrom,pos,WINDOW)
   if wid not in TREE_WINDOWS.get(chrom,{}):continue
   winmeta=TREE_WINDOWS[chrom][wid]
   lo,hi=BOUNDS[chrom]
   if region_label=='inside' and not (winmeta['start']>=lo and winmeta['end']<=hi):continue
   if region_label=='outside' and not (winmeta['end']<lo or winmeta['start']>hi):continue
   start=winmeta['start'];w=acc[(design['design_id'],start)];w['records']+=1
   ref,alt=fields[3],fields[4]
   if fields[6] not in ('PASS','.'):continue
   if ',' in alt or len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT':continue
   fm=fields[8].split(':')
   try:gi=fm.index('GT')
   except ValueError:w['missing']+=1;continue
   gt=[s.split(':')[gi] if gi<len(s.split(':')) else '.' for s in fields[9:]]
   def get(ids):return [gt[ix[s]] for s in ids]
   a,n1=freq(get(p1));b,n2=freq(get(p2));c,n3=freq(get(p3))
   if None in (a,b,c):w['missing']+=1;continue
   outs=get(outid)
   og=[]
   for s in outs:
    g=parse_gt(s or '.')
    if g is not None:og.extend(g)
   if not og:w['out_missing']+=1;continue
   if len(set(og))!=1:w['out_het']+=1;continue
   state=og[0]
   if state==0:w['out_fixed_ref']+=1;p1d,p2d,p3d=a,b,c
   else:w['out_fixed_alt']+=1;p1d,p2d,p3d=1-a,1-b,1-c
   comp=fdm_site(p1d,p2d,p3d,0.)
   if not comp.included:continue
   w['n_sites']+=1;w['num']+=comp.numerator;w['den']+=comp.denominator;w['abba']+=comp.abba;w['baba']+=comp.baba;w['informative']+=int(comp.abba+comp.baba>0)
 rows=[]
 for (design_id,start),v in sorted(acc.items(),key=lambda kv:(kv[0][0],kv[0][1])):
  design=next(d for d in designs if d['design_id']==design_id)
  end=start+WINDOW-1
  if region_label=='inside' and (start<summary_start(design['lg']) or end>summary_end(design['lg'])):continue
  if region_label=='outside' and not (end<summary_start(design['lg']) or start>summary_end(design['lg'])):continue
  den=v['den'];d_den=v['abba']+v['baba']
  rows.append(dict(lg=design['lg'],design_id=design['design_id'],region=region_label,window_start=start,window_end=end,
   f_dM_numerator=v['num'],f_dM_denominator=den,f_dM=v['num']/den if den else 'NA',
   ABBA=v['abba'],BABA=v['baba'],D_Patterson=(v['abba']-v['baba'])/d_den if d_den else 'NA',
   n_fdm_sites=v['n_sites'],n_informative_sites=v['informative'],n_variant_records=v['records'],
   missing_group_sites=v['missing'],outgroup_missing_sites=v['out_missing'],outgroup_heterozygous_sites=v['out_het'],
   outgroup_fixed_ref_sites=v['out_fixed_ref'],outgroup_fixed_alt_sites=v['out_fixed_alt']))
 return rows

def summary_start(lg):return BOUNDS[lg][0]
def summary_end(lg):return BOUNDS[lg][1]
BOUNDS={}
def main():
 global BOUNDS,TREE_WINDOWS
 TREE_WINDOWS=read_tree_windows(DATA/'processed/cod_window_trees.tsv')
 BOUNDS={r['lg']:(int(r['start']),int(r['end'])) for r in tsv(DATA/'metadata/region_manifest.tsv')}
 designs=tsv(ROOT/'empirical/atlantic_cod/analysis/patterson_d/config/frozen_quartets.tsv')
 files=[];allrows=[]
 # Shared outside VCF stream is scanned once per frozen design; inside VCFs are the deposited SNAPP inputs (exact frozen sample membership).
 outpath=RAW/'gadus_morhua_outside_supergenes_snapp.vcf.gz'
 for d in designs:
  lg=d['lg']; inside=RAW/f'gadus_morhua_supergenes_{lg.lower()}_snapp.vcf.gz'
  if not inside.exists():raise FileNotFoundError(inside)
  rr=scan(inside,[d],'inside');allrows.extend(rr);files.append((d['design_id'],'inside',inside.name,len(rr)))
 if not outpath.exists():raise FileNotFoundError(outpath)
 rr=scan(outpath,designs,'outside');allrows.extend(rr);files.append(('all_designs','outside',outpath.name,len(rr)))
 out=AN/'results/atlantic_cod_fdm_windows.tsv';out.parent.mkdir(parents=True,exist_ok=True)
 with out.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(allrows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(allrows)
 # Pool windows by 500-kb physical block and report jackknife; estimate regions from the sums, never averaging window ratios.
 summary=[]
 for d in designs:
  for region in ('inside','outside'):
   subset=[r for r in allrows if r['design_id']==d['design_id'] and r['region']==region]
   num=sum(float(r['f_dM_numerator']) for r in subset);den=sum(float(r['f_dM_denominator']) for r in subset)
   blocks=defaultdict(lambda:[0.,0.])
   for r in subset:
    bid=((int(r['window_start'])-1)//500000)*500000+1;q=blocks[bid];q[0]+=float(r['f_dM_numerator']);q[1]+=float(r['f_dM_denominator'])
   jk=jackknife_ratio(blocks)
   summary.append(dict(design_id=d['design_id'],lg=d['lg'],region=region,n_windows=len(subset),numerator=num,denominator=den,fdm=num/den if den else 'NA',
      n_sites=sum(int(r['n_fdm_sites']) for r in subset),n_informative_sites=sum(int(r['n_informative_sites']) for r in subset),
      jackknife_fdm=jk['fdm'],jackknife_se=jk['se'],ci_low=jk['ci_low'],ci_high=jk['ci_high'],n_blocks=jk['n_blocks']))
 p=AN/'results/atlantic_cod_fdm_regions.tsv'
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(summary)
 print('wrote',out,'and',p,'rows=',len(allrows))
if __name__=='__main__':main()
