#!/usr/bin/env python3
"""Correctly polarize biallelic chr17 SNPs with fixed-ref or fixed-alt outgroups.
Checkpointed remote indexed scan. Original analysis artifacts are read-only inputs.
"""
import os, ssl, certifi, csv, json, math, sys, hashlib
from pathlib import Path
os.environ['SSL_CERT_FILE']=certifi.where(); os.environ['REQUESTS_CA_BUNDLE']=certifi.where()
import pysam
BASE=Path(__file__).resolve().parents[1]
ROOT=BASE.parents[3]
OLD=BASE.parent/'results'
URL='https://wwwuser.gwdguser.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz'
sys.path.insert(0,str(BASE.parent/'scripts'))
from run_patterson_d import GROUPS, D
OUT=BASE/'results'; OUT.mkdir(parents=True,exist_ok=True)
CK=OUT/'checkpoint.json'
def gt_freq(gt):
    if gt is None or len(gt)!=2 or any(a is None for a in gt): return None
    return (int(gt[0]==1)+int(gt[1]==1))/2
def polarization_state(gts):
    if any(x is None or len(x)!=2 or any(a is None for a in x) for x in gts): return 'missing'
    alleles=[a for gt in gts for a in gt]
    if all(a==0 for a in alleles): return 'fixed_ref'
    if all(a==1 for a in alleles): return 'fixed_alt'
    return 'polymorphic'
def contributions(p1,p2,p3): return (1-p1)*p2*p3, p1*(1-p2)*p3
def main():
    v=pysam.VariantFile(URL); sample_set=set(v.header.samples); names=sorted(set.union(*(set(s) for gs in GROUPS.values() for s in gs)))
    if not set(names)<=sample_set: raise RuntimeError(f'missing samples: {set(names)-sample_set}')
    end=v.header.contigs['chr17'].length
    resolutions=(5000,50000,100000,250000,500000)
    # Scan once, checkpointing compact per-5-kb class contributions after each 500-kb remote fetch.
    siteck=OUT/'site_checkpoint.json'; sites=json.loads(siteck.read_text()) if siteck.exists() else {}
    for chunk in range(0,end,500_000):
        key=str(chunk)
        if key in sites and sites[key].get('complete'): continue
        local={}
        for rec in v.fetch('chr17',chunk,min(chunk+500_000,end)):
            if len(rec.alleles)!=2 or any(a is None or len(a)!=1 or a not in 'ACGT' for a in rec.alleles): continue
            gts={s:rec.samples[s].get('GT') for s in names}; st=(rec.pos-1)//5000*5000
            for des,(p1s,p2s,p3s,oset) in GROUPS.items():
                cl=polarization_state([gts[s] for s in oset]); k=(des,st,cl); x=local.setdefault(k,{'ABBA':0.,'BABA':0.,'informative_sites':0,'variants':0});x['variants']+=1
                if cl not in ('fixed_ref','fixed_alt'): continue
                ps=[]
                for group in (p1s,p2s,p3s):
                    vals=[gt_freq(gts[s]) for s in group];vals=[z for z in vals if z is not None];ps.append(sum(vals)/len(vals) if vals else None)
                if any(z is None for z in ps):
                    x['variants']-=1;k=(des,st,'missing');x=local.setdefault(k,{'ABBA':0.,'BABA':0.,'informative_sites':0,'variants':0});x['variants']+=1;continue
                if cl=='fixed_alt': ps=[1-z for z in ps]
                ab,ba=contributions(*ps);x['ABBA']+=ab;x['BABA']+=ba;x['informative_sites']+=int(ab>0 or ba>0)
        rows=[{'design':k[0],'start':k[1],'class':k[2],**val} for k,val in local.items()]
        sites[key]={'complete':True,'rows':rows};siteck.write_text(json.dumps(sites))
    # flatten per 5kb window and class
    grouped={}
    for block in sites.values():
        for x in block['rows']:
            k=(x['design'],x['start'],x['class']);g=grouped.setdefault(k,{'ABBA':0.,'BABA':0.,'informative_sites':0,'variants':0})
            for f in g:g[f]+=x[f]
    original={}
    with (OLD/'patterson_d_window_results.tsv').open() as f:
        for x in csv.DictReader(f,delimiter='\t'):original[(x['design'],int(x['start']))]=x
    corrected=[];classrows=[];resagg={}
    for des in GROUPS:
      for st in range(0,end,5000):
        allc=[]
        for cl in ('fixed_ref','fixed_alt','polymorphic','missing'):
            c=grouped.get((des,st,cl),{'ABBA':0.,'BABA':0.,'informative_sites':0,'variants':0});allc.append(c)
            classrows.append({'design':des,'start':st,'end':min(st+5000,end),'polarization_class':cl,**c})
        a=sum(c['ABBA'] for c in allc);b=sum(c['BABA'] for c in allc);info=sum(c['informative_sites'] for c in allc)
        o=original[(des,st)]; od=float(o['D']) if o['D'] else float('nan'); nd=D(a,b)
        corrected.append({'design':des,'chromosome':'chr17','start':st,'end':min(st+5000,end),'fixed_ref_sites':allc[0]['variants'],'fixed_alt_sites':allc[1]['variants'],'polymorphic_outgroup_sites':allc[2]['variants'],'missing_outgroup_or_ingroup_sites':allc[3]['variants'],'fixed_ref_informative':allc[0]['informative_sites'],'fixed_alt_informative':allc[1]['informative_sites'],'ABBA_fixed_ref':allc[0]['ABBA'],'BABA_fixed_ref':allc[0]['BABA'],'ABBA_fixed_alt':allc[1]['ABBA'],'BABA_fixed_alt':allc[1]['BABA'],'ABBA':a,'BABA':b,'informative_sites':info,'D':nd,'original_D':od,'delta_D':(nd-od if math.isfinite(od) and nd!='' else '')})
        for res in resolutions:
            key=(des,res,st//res*res);g=resagg.setdefault(key,{'ABBA':0.,'BABA':0.,'informative_sites':0,'fixed_ref_informative':0,'fixed_alt_informative':0})
            g['ABBA']+=a;g['BABA']+=b;g['informative_sites']+=info;g['fixed_ref_informative']+=allc[0]['informative_sites'];g['fixed_alt_informative']+=allc[1]['informative_sites']
    def write(name,rows):
        with (OUT/name).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    write('corrected_window_D.tsv',corrected);write('polarization_class_contributions.tsv',classrows)
    write('original_vs_corrected_D.tsv',[{'design':x['design'],'start':x['start'],'end':x['end'],'original_D':x['original_D'],'corrected_D':x['D'],'delta_D':x['delta_D'],'original_informative_sites':original[(x['design'],int(x['start']))]['informative_sites'],'corrected_informative_sites':x['informative_sites'],'fixed_alt_informative':x['fixed_alt_informative']} for x in corrected])
    multi=[]
    for res in resolutions:
      for des in GROUPS:
       for st in range(0,end,res):
        g=resagg.get((des,res,st),{'ABBA':0.,'BABA':0.,'informative_sites':0,'fixed_ref_informative':0,'fixed_alt_informative':0})
        multi.append({'design':des,'resolution_bp':res,'start':st,'end':min(st+res,end),**g,'D':D(g['ABBA'],g['BABA'])})
    write('corrected_multiresolution_D.tsv',multi)
    print(json.dumps({'length':end,'five_kb_rows':len(corrected),'fixed_alt_informative_total':{d:sum(int(x['fixed_alt_informative']) for x in corrected if x['design']==d) for d in GROUPS},'bins_complete':sum(1 for x in sites.values() if x.get('complete'))},indent=2))
if __name__=='__main__':main()
