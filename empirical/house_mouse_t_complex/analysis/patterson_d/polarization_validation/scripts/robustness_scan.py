#!/usr/bin/env python3
"""Indexed matched-region transversion and ingroup-call-rate sensitivity scan."""
import os,certifi,csv,json,sys,math
from pathlib import Path
os.environ['SSL_CERT_FILE']=certifi.where();os.environ['REQUESTS_CA_BUNDLE']=certifi.where()
import pysam
BASE=Path(__file__).resolve().parents[1];P=BASE.parent;OUT=BASE/'results';EX=P/'exact_quartet_comparison'/'results';URL='https://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz'
sys.path.insert(0,str(P/'scripts'));from run_patterson_d import GROUPS
THRESHOLDS=(.8,.9,1.0);FILTERS=('all_biallelic_SNPs','transversions_only');LO,HI=5_000_000,37_055_000
CK=OUT/'robustness_checkpoint.json'
def complete(gt):return gt is not None and len(gt)==2 and all(a is not None for a in gt)
def fixed(gts):
 if any(not complete(x) for x in gts):return 'missing_outgroup'
 a={i for gt in gts for i in gt}
 if a=={0}:return 'fixed_ref'
 if a=={1}:return 'fixed_alt'
 return 'polymorphic_outgroup'
def transversion(ref,alt):return frozenset((ref,alt)) in {frozenset(('A','C')),frozenset(('A','T')),frozenset(('G','C')),frozenset(('G','T'))}
def main():
 v=pysam.VariantFile(URL); samples=set(v.header.samples); ck=json.loads(CK.read_text()) if CK.exists() else {};used=sorted(set.union(*(set(g) for gs in GROUPS.values() for g in gs)))
 if not set(used)<=samples:raise RuntimeError('sample ID mismatch')
 for chunk in range(LO,HI,500_000):
  k=str(chunk)
  if k in ck and ck[k].get('complete'):continue
  local={};excluded={d:{'missing_outgroup':0,'polymorphic_outgroup':0,'accepted_fixed_outgroup':0,'ingroup_callrate_failures':0} for d in GROUPS}
  for rec in v.fetch('chr17',chunk,min(chunk+500_000,HI)):
   if len(rec.alleles)!=2 or any(a is None or len(a)!=1 or a not in 'ACGT' for a in rec.alleles):continue
   gts={s:rec.samples[s].get('GT') for s in used}; st=(rec.pos-1)//5000*5000;tv=transversion(rec.ref,rec.alts[0]); filt_ids=['all_biallelic_SNPs']+(['transversions_only'] if tv else [])
   for d,gs in GROUPS.items():
    cl=fixed([gts[s] for s in gs[3]])
    if cl!='fixed_ref' and cl!='fixed_alt':excluded[d][cl]+=1;continue
    excluded[d]['accepted_fixed_outgroup']+=1
    vals=[];rates=[]
    for group in gs[:3]:
     calls=[gts[s] for s in group];valid=[gt for gt in calls if complete(gt)];rates.append(len(valid)/len(group));vals.append(sum((gt[0]==1)+(gt[1]==1) for gt in valid)/(2*len(valid)) if valid else None)
    for threshold in THRESHOLDS:
     if min(rates)<threshold or any(p is None for p in vals):excluded[d]['ingroup_callrate_failures']+=1;continue
     p=vals if cl=='fixed_ref' else [1-x for x in vals];p1,p2,p3=p;ab=(1-p1)*p2*p3;ba=p1*(1-p2)*p3
     for filt in filt_ids:
      key=(d,threshold,filt,st);x=local.setdefault(key,{'ABBA':0.,'BABA':0.,'informative_sites':0,'fixed_ref_sites':0,'fixed_alt_sites':0,'included_sites':0});x['ABBA']+=ab;x['BABA']+=ba;x['informative_sites']+=int(ab>0 or ba>0);x['fixed_ref_sites']+=int(cl=='fixed_ref');x['fixed_alt_sites']+=int(cl=='fixed_alt');x['included_sites']+=1
  ck[k]={'complete':True,'data':[{'design':z[0],'callrate_threshold':z[1],'filter':z[2],'start':z[3],**x} for z,x in local.items()],'exclusions':excluded};CK.write_text(json.dumps(ck))
 data=[x for c in ck.values() if c.get('complete') for x in c['data']];
 def write(name,rows):
  with (OUT/name).open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 pooled={};
 for x in data:
  for res in (5000,50000,100000,250000,500000):
   key=(x['design'],x['callrate_threshold'],x['filter'],res,int(x['start'])//res*res);g=pooled.setdefault(key,{'ABBA':0.,'BABA':0.,'informative_sites':0,'fixed_ref_sites':0,'fixed_alt_sites':0,'included_sites':0})
   for col in g:g[col]+=x[col]
 tracks=[]
 for key,g in sorted(pooled.items(),key=lambda z:(z[0][0],z[0][1],z[0][2],z[0][3],z[0][4])):
  d,th,filt,res,st=key;den=g['ABBA']+g['BABA'];tracks.append({'design':d,'callrate_threshold':th,'filter':filt,'resolution_bp':res,'start':st,'end':min(st+res,HI),**g,'D':(g['ABBA']-g['BABA'])/den if den else ''})
 write('robustness_window_counts.tsv',data);write('robustness_multiresolution.tsv',tracks)
 write('robustness_exclusion_counts.tsv',[{'chunk_start':int(k),'design':d,**v} for k,z in ck.items() if z.get('complete') for d,v in z['exclusions'].items()])
 # Region contrasts (matched tree range, inversion union vs complement) and correlations at 100kb.
 trees=list(csv.DictReader(open(EX/'matched_gene_tree_quartets.tsv'),delimiter='\t'));matched={}
 for t in trees:
  if t['D_topology'] in ('',None):continue
  k=(t['design'],int(t['start'])//100000*100000);matched.setdefault(k,[]).append(float(t['D_topology']))
 corrs=[];contrasts=[]
 for d in GROUPS:
  for th in THRESHOLDS:
   for filt in FILTERS:
    tt=[x for x in tracks if x['design']==d and x['callrate_threshold']==th and x['filter']==filt and x['resolution_bp']==100000 and x['start']>=LO and x['end']<=HI and x['D']!=''];xs=[]
    for x in tt:
     q=matched.get((d,int(x['start']))); 
     if q:xs.append({'design':d,'resolution_bp':100000,'start':x['start'],'ABBA':x['ABBA'],'BABA':x['BABA'],'informative_sites':x['informative_sites'],'D_site':x['D'],'D_topology':sum(q)/len(q)})
    sys.path.insert(0,str(P/'exact_quartet_comparison'/'scripts'));from run_analysis import association
    v=association(xs);corrs.append({'design':d,'callrate_threshold':th,'filter':filt,'resolution_bp':100000,'matched_windows':len(xs),'pearson_r':v[0],'spearman_rho':v[1],'sign_agreement':v[2],'block_bootstrap_low':v[3],'block_bootstrap_high':v[4],'circular_shift_null_abs_r_p':v[5],'block_level_pearson_r':v[6]})
    inv=[(5747331,6957919),(7435308,13016510),(14068229,18051507),(19401000,40049964)]
    ins=[x for x in tt if any(a<=x['start'] and x['end']<=b for a,b in inv)];outs=[x for x in tt if not any(x['start']<b and x['end']>a for a,b in inv)]
    pool=lambda rr:(sum(float(x['ABBA']) for x in rr),sum(float(x['BABA']) for x in rr));ai,bi=pool(ins);ao,bo=pool(outs);D=lambda a,b:(a-b)/(a+b) if a+b else float('nan')
    contrasts.append({'design':d,'callrate_threshold':th,'filter':filt,'inside_windows':len(ins),'outside_windows':len(outs),'inside_D':D(ai,bi),'outside_D':D(ao,bo),'inside_minus_outside_D':D(ai,bi)-D(ao,bo),'inside_informative':sum(int(x['informative_sites']) for x in ins),'outside_informative':sum(int(x['informative_sites']) for x in outs)})
 write('robustness_correlations.tsv',corrs);write('robustness_inside_outside.tsv',contrasts)
 print(json.dumps({'chunks':sum(x.get('complete',False) for x in ck.values()),'window_rows':len(data),'correlation_rows':len(corrs)},indent=2))
if __name__=="__main__":main()
