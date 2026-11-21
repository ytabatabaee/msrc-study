#!/usr/bin/env python3
"""Mouse t-vs-standard SNP-panel stats with sitewise pseudo-t allele inference.
Carrier diploid genotypes are never directly counted as pure t chromosomes.
"""
import os, certifi, json
os.environ['SSL_CERT_FILE']=certifi.where(); os.environ['REQUESTS_CA_BUNDLE']=certifi.where()
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[4]
OUT=ROOT/'empirical/house_mouse_t_complex/window_diversity'
URL='https://wwwuser.gwdguser.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz'
TBI=ROOT/'AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz.tbi'
MEM=ROOT/'empirical/house_mouse_t_complex/analysis/patterson_d/config/frozen_quartet_members.tsv'
TREE=ROOT/'empirical/house_mouse_t_complex/results/stage5_unified_tracks/window_unified_tracks.tsv'
INV=ROOT/'empirical/house_mouse_t_complex/results/t_complex_inversion_intervals.tsv'
CHECK=OUT/'results/mouse_scan_checkpoint.json'

def allele_counts(gt):
    if gt is None or len(gt)!=2 or any(a is None for a in gt): return None
    return (int(gt[0]),int(gt[1]))

def pi_unbiased(alleles):
    n=len(alleles)
    if n<2:return np.nan
    p=np.bincount(np.asarray(alleles,dtype=int),minlength=2)/n
    return n/(n-1)*(1-float((p*p).sum()))

def main():
    mem=pd.read_csv(MEM,sep='\t'); tree=pd.read_csv(TREE,sep='\t'); inv=pd.read_csv(INV,sep='\t')
    designs={}
    for des,g in mem[mem.design.str.startswith(('D1_Mmd','D2_Mmm','D3_Mmc'))].groupby('design'):
        designs[des]=(g[g.role=='P1'].sample_id.tolist(),g[g.role=='P2'].sample_id.tolist())
    v=pysam.VariantFile(URL,index_filename=str(TBI)); samples=set(v.header.samples)
    for des,(t,std) in designs.items():
        if not set(t+std)<=samples: raise ValueError(f'{des}: sample mismatch')
    sample_union=sorted(set(s for groups in designs.values() for g in groups for s in g))
    windows={int(r.start_bp):r for r in tree.itertuples()}
    start=min(windows); end=max(int(r.end_bp)+1 for r in tree.itertuples())
    state=json.loads(CHECK.read_text()) if CHECK.exists() else {'processed':[],'acc':{},'records':0,'ambiguous':{d:0 for d in designs},'carrier_opportunities':{d:0 for d in designs}}
    done=set(state['processed']); chunk_size=500_000
    chunk_starts=list(range(start,end,chunk_size))
    for chunk in chunk_starts:
        if chunk in done: continue
        print(f'Fetching chr17:{chunk:,}-{min(chunk+chunk_size,end):,}',flush=True)
        for rec in v.fetch('chr17',chunk,min(chunk+chunk_size,end)):
            if len(rec.alleles)!=2 or not rec.ref or len(rec.ref)!=1 or rec.ref not in 'ACGT' or not rec.alts or len(rec.alts[0])!=1 or rec.alts[0] not in 'ACGT': continue
            st=((rec.pos-1)//5000)*5000
            if st not in windows: continue
            state['records']+=1
            gtmap={s:allele_counts(rec.samples[s].get('GT')) for s in sample_union}
            for des,(ts,ss) in designs.items():
                std_calls=[gtmap[x] for x in ss]
                if any(x is None for x in std_calls): continue
                std_alleles=[a for gt in std_calls for a in gt]
                present=set(std_alleles); t_alleles=[]; amb=0
                for x in ts:
                    gt=gtmap[x]
                    if gt is None: amb+=1; continue
                    if gt[0]==gt[1]: t_alleles.append(gt[0]); continue
                    unique=[a for a in set(gt) if a not in present]
                    if len(unique)==1: t_alleles.append(unique[0])
                    else: amb+=1
                state['ambiguous'][des]+=amb; state['carrier_opportunities'][des]+=len(ts)
                if len(t_alleles)<2: continue
                n=len(t_alleles); p=np.bincount(np.asarray(t_alleles,dtype=int),minlength=2)/n
                ns=len(std_alleles); q=np.bincount(np.asarray(std_alleles,dtype=int),minlength=2)/ns
                pi_t=n/(n-1)*(1-float((p*p).sum())); pi_s=ns/(ns-1)*(1-float((q*q).sum()))
                dxy=1-float((p*q).sum()); num=dxy-(pi_t+pi_s)/2
                key=f'{des}|{st}'; a=state['acc'].setdefault(key,[0.0,0.0,0.0,0.0,0,0,0])
                a[0]+=pi_t; a[1]+=pi_s; a[2]+=dxy; a[3]+=num; a[4]+=1; a[5]+=n; a[6]+=amb
        state['processed'].append(chunk); done.add(chunk)
        CHECK.write_text(json.dumps(state,separators=(',',':')))
    v.close()
    outputs=[]
    for des,(ts,ss) in designs.items():
        d=tree.copy(); rows=[]
        for st in d.start_bp.astype(int):
            a=state['acc'].get(f'{des}|{st}')
            if a and a[4]:
                fst=a[3]/a[2] if a[2]>0 else np.nan
                rows.append((a[4],fst,a[2]/a[4],a[0]/a[4],a[1]/a[4],a[5]/a[4],a[6]/a[4]))
            else: rows.append((0,np.nan,np.nan,np.nan,np.nan,0,np.nan))
        d['design']=des; d['n_t_carriers']=len(ts); d['n_standard_diploids']=len(ss)
        d[['n_variant_sites_reconstructed','Hudson_FST_SNP_panel','dXY_per_observed_SNP','pi_t_pseudohap_per_observed_SNP','pi_standard_per_observed_SNP','mean_reconstructed_t_carriers_per_site','mean_ambiguous_carrier_calls_per_site']]=pd.DataFrame(rows,index=d.index)
        d['call_policy']='PASS biallelic SNPs; standard GT complete; carrier hom GT direct; het GT used only if one allele absent from all same-subspecies standards; otherwise masked'
        outputs.append(d)
    out=pd.concat(outputs,ignore_index=True); out.to_csv(OUT/'results/mouse_t_pseudohap_snp_panel_window_stats.tsv',sep='\t',index=False,na_rep='NA')
    invs=list(zip(inv.start_bp.astype(int),inv.end_bp.astype(int))); summ=[]; rng=np.random.default_rng(20261120)
    metrics=['Hudson_FST_SNP_panel','dXY_per_observed_SNP','pi_t_pseudohap_per_observed_SNP','pi_standard_per_observed_SNP']
    for des,g in out.groupby('design'):
        g=g.copy(); g['inside']=g.apply(lambda r:any(r.start_bp<=hi and r.end_bp>=lo for lo,hi in invs),axis=1); g['block']=g.midpoint_bp//500_000
        for col in metrics:
            x=g[g.inside][col].dropna(); y=g[~g.inside][col].dropna(); delta=x.mean()-y.mean()
            xb=g[g.inside].groupby('block')[col].mean().dropna().to_numpy(); yb=g[~g.inside].groupby('block')[col].mean().dropna().to_numpy(); boot=[]
            if len(xb)>1 and len(yb)>1:
                for _ in range(10000): boot.append(rng.choice(xb,len(xb),replace=True).mean()-rng.choice(yb,len(yb),replace=True).mean())
            ci=np.quantile(boot,[.025,.975]) if boot else [np.nan,np.nan]
            summ.append(dict(design=des,statistic=col,n_inside_windows=len(x),n_outside_windows=len(y),inside_mean=x.mean(),outside_mean=y.mean(),difference=delta,block_width_bp=500000,n_inside_blocks=len(xb),n_outside_blocks=len(yb),block_bootstrap_95_low=ci[0],block_bootstrap_95_high=ci[1]))
    pd.DataFrame(summ).to_csv(OUT/'results/mouse_t_pseudohap_inside_outside.tsv',sep='\t',index=False,na_rep='NA')
    corr=[]
    for des,g in out.groupby('design'):
        for stat in metrics:
            for track in ['A_distance','A_relative','mean_pair_distance_no_spretus','q_t_alt','D']:
                z=g[[stat,track]].dropna()
                if len(z)>2: corr.append(dict(design=des,statistic=stat,track=track,n=len(z),pearson=z[stat].corr(z[track]),spearman=z[stat].rank().corr(z[track].rank())))
    pd.DataFrame(corr).to_csv(OUT/'results/mouse_t_pseudohap_tree_track_correlations.tsv',sep='\t',index=False,na_rep='NA')
    for des,g in out.groupby('design'):
        fig,axs=plt.subplots(4,1,figsize=(12,9),sharex=True,constrained_layout=True)
        for ax,metric,label in zip(axs,metrics,['Hudson FST (SNP panel)','dXY per observed SNP','π inferred t allele per observed SNP','π standard per observed SNP']):
            for lo,hi in invs: ax.axvspan(lo/1e6,hi/1e6,color='#e7a34a',alpha=.14,lw=0)
            ax.scatter(g.midpoint_bp/1e6,g[metric],s=5,color='#245b78',alpha=.7,rasterized=True); ax.set_ylabel(label); ax.grid(axis='y',alpha=.2)
        axs[-1].set_xlabel('mm10/GRCm38 chromosome 17 coordinate (Mb)')
        fig.suptitle(f'{des}: t pseudohaplotype vs standard\nIndexed PASS SNPs; ambiguous carrier phase masked; not per callable base',fontsize=13)
        fig.savefig(OUT/f'figures/{des}_snp_panel_tracks.pdf'); fig.savefig(OUT/f'figures/{des}_snp_panel_tracks.png',dpi=300)
    audit={'source_url':URL,'index':str(TBI.name),'samples_total':len(samples),'chr17_length':94987271,'windows':len(tree),'window_coordinate_axis':'mm10/GRCm38 chr17, 5-kb windows; same axis as existing stage5 tracks','variant_records_in_tree_span':state['records'],'processed_500kb_chunks':len(state['processed']),'designs':{d:{'n_t':len(t),'n_standard':len(s),'ambiguous_or_missing_carrier_calls':state['ambiguous'][d],'total_carrier_site_opportunities':state['carrier_opportunities'][d]} for d,(t,s) in designs.items()},'genotype_call_policy':'PASS biallelic SNPs only; standard diploid calls complete; carrier homozygotes contribute their allele; carrier heterozygotes contribute the unique allele only if absent from all same-subspecies standard chromosomes, otherwise masked','limitations':['variant-only VCF; no invariant callable denominator','sitewise pseudo-t alleles do not resolve phase across loci','PASS VCF analysis uses GT only; no read-level DP/GQ filters','standard groups include geographic structure within subspecies']}
    (OUT/'results/mouse_source_and_reconstruction_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2)); print(pd.DataFrame(summ).to_string(index=False))
if __name__=='__main__': main()
