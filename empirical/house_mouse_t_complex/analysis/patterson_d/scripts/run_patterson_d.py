#!/usr/bin/env python3
"""Indexed chr17 Patterson D on Harr et al. AllMouse PASS VCF.
Uses genotype allele frequencies and fixed ancestral polarization by M. spretus.
"""
import os, ssl, certifi
os.environ['SSL_CERT_FILE']=certifi.where(); os.environ['REQUESTS_CA_BUNDLE']=certifi.where()
import csv, json, math, hashlib, sys
from pathlib import Path
import pysam
ROOT=Path(__file__).resolve().parents[5]
BASE=Path(__file__).resolve().parents[1]
RESULTS=BASE/'results'; FIG=BASE/'figures'
URL='https://wwwuser.gwdguser.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz'
TIPS=['Mmd_FRA5_B2C.variant','Mmd_FRA6_C1.variant','Mmd_FRA8_F1B.variant','Mmd_GER6_TP51D.variant2','Mmm_AFG1_396.variant7','Mmm_AFG3_416.variant7','Mmm_AFG6_444.variant7','Mmm_CZE1_CR12.variant5','Mmm_CZE2_CR13.variant5','Mmm_CZE4_CR23.variant5','Mmm_KAZ3_AL19.variant6','Mmm_KAZ7_AL41.variant6','Mmc_CAST1_H12.variant8','Mmc_CAST2_H14.variant8','Mmc_CAST6_H27.variant8']
# Variant IDs are exact VCF sample names, confirmed against the header.
SPRETUS={f'Ms_SPRE{i}_SP{n}.variant9' for i,n in enumerate([36,39,41,51,62,68,69,70],1)}
MMD_STD={'Mmd_FRA1_14.variant','Mmd_FRA2_15B.variant','Mmd_FRA3_16B.variant','Mmd_FRA4_18B.variant','Mmd_FRA7_E1.variant','Mmd_GER1_TP1.variant2','Mmd_GER2_TP121B.variant2','Mmd_GER3_TP17-2.variant2','Mmd_GER4_TP3-92.variant2','Mmd_GER5_TP4a.variant2','Mmd_GER7_TP7-10F1A2.variant2','Mmd_GER8_TP81B.variant2'}
MMM_STD={'Mmm_AFG2_413.variant7','Mmm_AFG4_424.variant7','Mmm_AFG5_435.variant7','Mmm_CZE3_CR16.variant5','Mmm_CZE5_CR25.variant5','Mmm_CZE7_CR46.variant5','Mmm_CZE8_CR59.variant5','Mmm_KAZ1_AL1.variant6','Mmm_KAZ2_AL16.variant6','Mmm_KAZ4_AL33.variant6','Mmm_KAZ5_AL38.variant6','Mmm_KAZ6_AL40.variant6','Mmm_KAZ8_AL42.variant6'}
MMD_T={'Mmd_FRA5_B2C.variant','Mmd_FRA6_C1.variant','Mmd_FRA8_F1B.variant','Mmd_GER6_TP51D.variant2'}
MMM_T={'Mmm_AFG1_396.variant7','Mmm_AFG3_416.variant7','Mmm_AFG6_444.variant7','Mmm_CZE1_CR12.variant5','Mmm_CZE2_CR13.variant5','Mmm_CZE4_CR23.variant5','Mmm_KAZ3_AL19.variant6','Mmm_KAZ7_AL41.variant6'}
MMC_T={'Mmc_CAST1_H12.variant8','Mmc_CAST2_H14.variant8','Mmc_CAST6_H27.variant8'}
MMC_STD={'Mmc_CAST3_H15.variant8','Mmc_CAST4_H24.variant8','Mmc_CAST5_H26.variant8','Mmc_CAST7_H28.variant8','Mmc_CAST8_H30.variant8','Mmc_CAST9_H34.variant8','Mmc_CAST10_H36.variant8'}
GROUPS={
'D1_Mmd_t_vs_standard':(MMD_T,MMD_STD,MMM_STD,SPRETUS),
'D2_Mmm_t_vs_standard':(MMM_T,MMM_STD,MMD_STD,SPRETUS),
'D3_Mmc_t_vs_standard':(MMC_T,MMC_STD,MMM_STD,SPRETUS),
}
INV=[('In(17)1',5747331,6957919),('In(17)2',7435308,13016510),('In(17)3',14068229,18051507),('In(17)4',19401000,40049964)]
def D(a,b): return (a-b)/(a+b) if a+b else ''
def write(path,rows):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def freq(call):
 if call is None or call[0] is None or len(call)<2 or call[1] is None:return None
 return (int(call[0]==1)+int(call[1]==1))/2
def outgroup_fixed_ref(calls,names):
 return all(calls.get(name) is not None and calls[name]==0 for name in names)
def site_contrib(p1,p2,p3):
 return (1-p1)*p2*p3, p1*(1-p2)*p3

def main():
 v=pysam.VariantFile(URL); samples=set(v.header.samples); designs={}
 for name,groups in GROUPS.items():
  for g in groups:
   if not g<=samples: raise RuntimeError(f'{name} missing samples: {g-samples}')
  designs[name]=groups
 sample_idx={s:i for i,s in enumerate(v.header.samples)}
 # Cache and checkpoint in 500-kb indexed chunks to minimize remote HTTP range requests.
 bins={}; start,end=0,v.header.contigs['chr17'].length; ck=RESULTS/'checkpoint_counts.json'
 if ck.exists(): bins=json.loads(ck.read_text())
 for chunk in range(start,end,500_000):
  pending=[st for st in range(chunk,min(chunk+500_000,end),5000) if str(st) not in bins]
  if not pending: continue
  upper=min(chunk+500_000,end)
  for st in pending:
   ct={'variants':0,'callable_sites':0,'missing_calls':0}
   for des in designs: ct[des]=[0.,0.,0,0,0,0]
   bins[str(st)]=ct
  try: it=v.fetch('chr17',chunk,upper)
  except Exception:
   ck.write_text(json.dumps(bins)); continue
  for rec in it:
   if len(rec.alleles)>2 or rec.alleles[1] is None: continue
   if not rec.ref or len(rec.ref)!=1 or rec.ref not in {'A','C','G','T'} or not rec.alts or any(len(alt)!=1 or alt not in {'A','C','G','T'} for alt in rec.alts): continue
   st=(rec.pos-1)//5000*5000; ct=bins.get(str(st))
   if ct is None: continue
   ct['variants']+=1
   calls={sample:freq(rec.samples[sample].get('GT')) for sample in samples}
   for des,(p1s,p2s,p3s,oset) in designs.items():
    if not outgroup_fixed_ref(calls,oset): ct[des][4]+=1; continue
    vals=[]
    for ss in (p1s,p2s,p3s,oset):
     x=[calls[name] for name in ss if calls[name] is not None]
     vals.append(sum(x)/len(x) if x else None)
    p1,p2,p3,po=vals
    if po is None or p1 is None or p2 is None or p3 is None: ct[des][4]+=1; continue
    if po!=0: continue
    ab,ba=site_contrib(p1,p2,p3)
    ct[des][0]+=ab; ct[des][1]+=ba; ct[des][2]+=1; ct[des][3]+=int(ab>0 or ba>0)
  ck.write_text(json.dumps(bins))
 raw=[]
 for st in range(start,end,5000):
  b=bins[str(st)]
  for des in designs:
   a,ba,polarized_n,informative,missing,zero=b[des]
   raw.append({'design':des,'chromosome':'chr17','start':st,'end':min(st+5000,end),'variants':b['variants'],'callable_site_proxy':'NA_not_in_SNP_only_VCF','sites_excluded_for_missing_population':missing,'polarized_callable_sites':polarized_n,'ABBA':a,'BABA':ba,'informative_sites':informative,'D':D(a,ba)})
 write(RESULTS/'patterson_d_window_results.tsv',raw)
 # Pooled contiguous physical bins at requested resolutions; aggregate in one pass.
 multi=[]; resolutions=(5000,50000,100000,250000,500000); pooled={}
 for r in raw:
  for res in resolutions:
   st=int(r['start'])//res*res; key=(r['design'],res,st)
   g=pooled.setdefault(key,[0.,0.,0])
   g[0]+=float(r['ABBA']);g[1]+=float(r['BABA']);g[2]+=int(r['informative_sites'])
 for res in resolutions:
  for des in designs:
   for st in range(start,end,res):
    a,ba,informative=pooled.get((des,res,st),(0.,0.,0))
    multi.append({'design':des,'resolution_bp':res,'start':st,'end':min(st+res,end),'ABBA':a,'BABA':ba,'informative_sites':informative,'D':D(a,ba)})
 write(RESULTS/'patterson_d_multiresolution.tsv',multi)
 # Paired inside/outside using physical blocks for delete-one-block jackknife.
 inside=[]
 for des in designs:
  for label,lo,hi in INV+[('combined_inversions',INV[0][1],INV[-1][2])]:
   for mode in ('overlap_union','individual_interval'):
    ir=[r for r in raw if r['design']==des and ((lo<=r['start'] and r['end']<=hi) if mode=='individual_interval' else any(r['start']<b and r['end']>a for _,a,b in INV) if label=='combined_inversions' else r['start']<hi and r['end']>lo)]
    inside_starts={r['start'] for r in ir}
    out=[r for r in raw if r['design']==des and r['start'] not in inside_starts]
    ai=sum(float(r['ABBA']) for r in ir);bi=sum(float(r['BABA']) for r in ir);ao=sum(float(r['ABBA']) for r in out);bo=sum(float(r['BABA']) for r in out)
    inside.append({'design':des,'region':label,'selection':mode,'inside_ABBA':ai,'inside_BABA':bi,'inside_informative':sum(int(r['informative_sites']) for r in ir),'inside_D':D(ai,bi),'outside_ABBA':ao,'outside_BABA':bo,'outside_informative':sum(int(r['informative_sites']) for r in out),'outside_D':D(ao,bo),'inside_minus_outside_D':(D(ai,bi)-D(ao,bo) if D(ai,bi)!='' and D(ao,bo)!='' else '')})
 write(RESULTS/'patterson_d_inside_outside.tsv',inside)
 audit={'source_url':URL,'chr17':v.header.contigs['chr17'].length,'samples_total':len(samples),'sample_ids_used':len(set.union(*(set(g) for gs in designs.values() for g in gs))),'designs':{k:[len(g) for g in gs] for k,gs in designs.items()},'bins':len(bins),'region':[start,end],'coordinate':'mm10/GRCm38 chr17 one-based VCF POS','site_estimator':'unphased diploid ALT allele frequencies; only sites with all outgroup ALT frequencies exactly zero; SNP-only; ABBA=(1-p1)*p2*p3; BABA=p1*(1-p2)*p3'}
 (RESULTS/'analysis_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps(audit,indent=2))
if __name__=='__main__':main()
