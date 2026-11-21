#!/usr/bin/env python3
"""SB vs Sb SNP-panel diversity on mapped fire-ant chr16 BUSCO intervals.
All-male S. invicta Argentine cohort; outputs are per observed SNP, not per base."""
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[4]
OUT=ROOT/'empirical/fire_ants_chr16/window_diversity'
VCF=ROOT/'data/fire_ants_chr16/raw/genotypes/gt.vcf.gz'
MAP=ROOT/'empirical/fire_ants_chr16/analysis/patterson_d/results/window_gene_coordinate_mapping.tsv'
MEM=ROOT/'empirical/fire_ants_chr16/analysis/patterson_d/config/frozen_sample_membership.tsv'
INV=ROOT/'empirical/fire_ants_chr16/analysis/patterson_d/results/vcf_sample_inventory.tsv'
DP_MIN=3
MIN_CALL_FRACTION=0.80

def main():
    mp=pd.read_csv(MAP,sep='\t'); mp=mp[mp.tree_chromosome=='chr16'].copy()
    win=mp.drop_duplicates('window_index')[['window_index','tree_chromosome','tree_window_start','tree_window_end','tree_window_mid','supplement_genomic_region']]
    memb=pd.read_csv(MEM,sep='\t'); inv=pd.read_csv(INV,sep='\t')
    eligible=memb[memb.quartet_id=='Q1_primary_Argentina'].merge(inv,left_on='sample_id',right_on='vcf_sample_id')
    elig=eligible[(eligible.species=='S.invicta')&(eligible.caste=='Male')]
    sb=elig[elig.supergene_variant=='SB'].vcf_sample_id.tolist(); sb2=elig[elig.supergene_variant=='Sb'].vcf_sample_id.tolist()
    vf=pysam.VariantFile(str(VCF)); samples=list(vf.header.samples); sidx=[samples.index(x) for x in sb]; tidx=[samples.index(x) for x in sb2]
    if len(sb)!=78 or len(sb2)!=30: raise ValueError(f'Unexpected sample counts: SB={len(sb)}, Sb={len(sb2)}')
    if len(set(sb+sb2))!=len(sb)+len(sb2): raise ValueError('Duplicate sample membership')
    # Exact BUSCO gene intervals assign observed VCF SNPs to local tree windows.
    targets=defaultdict(list)
    for _,r in mp.iterrows():
        if pd.notna(r.vcf_contig) and int(r.vcf_contig_present)==1:
            targets[str(r.vcf_contig)].append((int(r.vcf_gene_start_1based),int(r.vcf_gene_end_1based_inclusive),int(r.window_index)))
    acc=defaultdict(list); n_het=0; n_dp=0; n_snp=0; n_matched=0
    # The VCF is variant-only; stream once, decoding sample fields only in mapped intervals.
    for rec in vf:
        iv=targets.get(rec.contig)
        if not iv: continue
        wins={w for lo,hi,w in iv if lo<=rec.pos<=hi}
        if not wins: continue
        if len(rec.alleles)!=2 or any(a not in ('A','C','G','T') for a in rec.alleles): continue
        group_calls=[]
        for group in (sidx,tidx):
            calls=[]
            for j in group:
                sm=rec.samples[samples[j]]; gt=sm.get('GT'); dp=sm.get('DP')
                if not gt or any(a is None for a in gt): continue
                if len(set(gt))!=1: n_het+=1; continue
                if dp is None or dp<DP_MIN: n_dp+=1; continue
                calls.append(gt[0])
            group_calls.append(calls)
        if len(group_calls[0])<int(np.ceil(MIN_CALL_FRACTION*len(sb))) or len(group_calls[1])<int(np.ceil(MIN_CALL_FRACTION*len(sb2))): continue
        n_snp+=1
        a=np.asarray(group_calls[0],dtype=int); b=np.asarray(group_calls[1],dtype=int)
        # Include only biallelic SNPs observed among sampled arrangement chromosomes.
        p=np.array([(a==k).mean() for k in (0,1)]); q=np.array([(b==k).mean() for k in (0,1)])
        pi1=len(a)/(len(a)-1)*(1-(p*p).sum()); pi2=len(b)/(len(b)-1)*(1-(q*q).sum())
        dxy=1-float((p*q).sum()); num=dxy-(pi1+pi2)/2
        val=(pi1,pi2,dxy,num,len(a),len(b))
        for w in wins: acc[w].append(val)
        n_matched+=1
    vf.close()
    out=win.copy(); rows=[]
    for _,r in out.iterrows():
        arr=np.asarray(acc.get(int(r.window_index),[]),float)
        if arr.size:
            den=arr[:,2].sum(); fst=arr[:,3].sum()/den if den>0 else np.nan
            rows.append((len(arr),fst,arr[:,2].mean(),arr[:,0].mean(),arr[:,1].mean(),arr[:,4].mean(),arr[:,5].mean()))
        else: rows.append((0,np.nan,np.nan,np.nan,np.nan,np.nan,np.nan))
    out[['n_variable_SNPs_within_80pct_call_rule','Hudson_FST_per_observed_SNP_panel','dXY_per_observed_SNP','pi_SB_per_observed_SNP','pi_Sb_per_observed_SNP','mean_n_SB_called','mean_n_Sb_called']]=pd.DataFrame(rows,index=out.index)
    out['n_SB']=len(sb); out['n_Sb']=len(sb2); out['call_policy']=f'homozygous GT, DP>={DP_MIN}, at least {MIN_CALL_FRACTION:.0%} samples called per group; males treated as haploid; SNPs only'
    out.to_csv(OUT/'results/sinvicta_chr16_snp_panel_window_stats.tsv',sep='\t',index=False,na_rep='NA')
    # Description-only contrast by mapped chr16r region labels; resample 1-Mb blocks.
    reg=out[out.n_variable_SNPs_within_80pct_call_rule>0].copy(); reg['inside']=reg.supplement_genomic_region.eq('chr16nr'); reg['block']=reg.tree_window_mid//1_000_000
    summ=[]; rng=np.random.default_rng(20261120)
    for col in ['Hudson_FST_per_observed_SNP_panel','dXY_per_observed_SNP','pi_SB_per_observed_SNP','pi_Sb_per_observed_SNP']:
        x=reg[reg.inside][col].dropna(); y=reg[~reg.inside][col].dropna(); d=x.mean()-y.mean()
        xb=reg[reg.inside].groupby('block')[col].mean().dropna().to_numpy(); yb=reg[~reg.inside].groupby('block')[col].mean().dropna().to_numpy(); boot=[]
        if len(xb)>1 and len(yb)>1:
            for _ in range(10000): boot.append(rng.choice(xb,len(xb),replace=True).mean()-rng.choice(yb,len(yb),replace=True).mean())
        lo,hi=np.quantile(boot,[.025,.975]) if boot else (np.nan,np.nan)
        summ.append(dict(statistic=col,n_inside_windows=len(x),n_outside_windows=len(y),inside_mean=x.mean(),outside_mean=y.mean(),difference=d,block_width_bp=1_000_000,n_inside_blocks=len(xb),n_outside_blocks=len(yb),block_bootstrap_95_low=lo,block_bootstrap_95_high=hi))
    pd.DataFrame(summ).to_csv(OUT/'results/sinvicta_chr16_inside_outside.tsv',sep='\t',index=False,na_rep='NA')
    # Track plots ordered on local chromosome coordinate; shaded chr16r spans are annotated from mapped window labels.
    fig,axs=plt.subplots(4,1,figsize=(12,9),sharex=True,constrained_layout=True)
    metrics=['Hudson_FST_per_observed_SNP_panel','dXY_per_observed_SNP','pi_SB_per_observed_SNP','pi_Sb_per_observed_SNP']
    labels=['Hudson FST (SNP panel)','dXY per observed SNP','π SB per observed SNP','π Sb per observed SNP']
    for ax,col,label in zip(axs,metrics,labels):
        ax.axvspan(11.680438,27.917498,color='#e7a34a',alpha=.14,lw=0)
        ax.axvline(11.680438,color='#a75c0b',lw=.9); ax.axvline(27.917498,color='#a75c0b',lw=.9)
        ax.scatter(out.tree_window_mid/1e6,out[col],s=12,color='#245b78',alpha=.8,rasterized=True); ax.set_ylabel(label); ax.grid(axis='y',alpha=.2)
    axs[-1].set_xlabel('Chromosome 16 local tree coordinate (Mb; Stolle / gng20170922wFex axis)')
    fig.suptitle('S. invicta Argentine males: SB vs Sb\nMapped BUSCO intervals; SNP-ascertained, not per callable base',fontsize=13)
    fig.savefig(OUT/'figures/sinvicta_chr16_snp_panel_tracks.pdf'); fig.savefig(OUT/'figures/sinvicta_chr16_snp_panel_tracks.png',dpi=300)
    print(f'Samples SB={len(sb)} Sb={len(sb2)}; windows={len(out)} with_sites={(out.n_variable_SNPs_within_80pct_call_rule>0).sum()}; complete_sites={n_matched}; het_site_calls_rejected={n_het}; low_DP_rejections={n_dp}')
    print(pd.DataFrame(summ).to_string(index=False))
if __name__=='__main__': main()
