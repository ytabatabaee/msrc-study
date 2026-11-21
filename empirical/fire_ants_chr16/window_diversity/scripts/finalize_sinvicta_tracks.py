#!/usr/bin/env python3
"""Rebuild fire-ant spatial summaries/figure from saved window estimates.
chr16nr is the nonrecombining supergene; chr16r is the recombining background.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[4]
OUT=ROOT/'empirical/fire_ants_chr16/window_diversity'

def main():
    p=OUT/'results/sinvicta_chr16_snp_panel_window_stats.tsv'
    d=pd.read_csv(p,sep='\t'); d['inside']=d.supplement_genomic_region.eq('chr16nr'); d['block']=d.tree_window_mid//1_000_000
    metrics=['Hudson_FST_per_observed_SNP_panel','dXY_per_observed_SNP','pi_SB_per_observed_SNP','pi_Sb_per_observed_SNP']
    labels=['Hudson FST (SNP panel)','dXY per observed SNP','π SB per observed SNP','π Sb per observed SNP']
    rng=np.random.default_rng(20261120); rows=[]
    for col in metrics:
        x=d[d.inside][col].dropna(); y=d[~d.inside][col].dropna(); delta=x.mean()-y.mean()
        xb=d[d.inside].groupby('block')[col].mean().dropna().to_numpy(); yb=d[~d.inside].groupby('block')[col].mean().dropna().to_numpy(); boot=[]
        if len(xb)>1 and len(yb)>1:
            for _ in range(10000): boot.append(rng.choice(xb,len(xb),replace=True).mean()-rng.choice(yb,len(yb),replace=True).mean())
        ci=np.quantile(boot,[.025,.975]) if boot else [np.nan,np.nan]
        rows.append(dict(statistic=col,n_inside_windows=len(x),n_outside_windows=len(y),inside_mean=x.mean(),outside_mean=y.mean(),difference_inside_minus_outside=delta,block_width_bp=1_000_000,n_inside_blocks=len(xb),n_outside_blocks=len(yb),block_bootstrap_95_low=ci[0],block_bootstrap_95_high=ci[1]))
    pd.DataFrame(rows).to_csv(OUT/'results/sinvicta_chr16_inside_outside.tsv',sep='\t',index=False,na_rep='NA')
    tree=pd.read_csv(ROOT/'empirical/fire_ants_chr16/results/stage9_distance/stage9_window_patristic_distance.tsv',sep='\t')
    quartet=pd.read_csv(ROOT/'data/fire_ants_chr16/processed/stage4a_window_quartet_support.tsv',sep='\t')
    tracks=tree[tree.chrom=='chr16'][['window_index','A_distance','A_relative','mean_all_focal_pair_distances']].merge(quartet[['window_index','q_S','q_H','q_3']],on='window_index',how='left')
    joined=d.merge(tracks,on='window_index',how='left')
    joined.to_csv(OUT/'results/sinvicta_chr16_tree_diversity_join.tsv',sep='\t',index=False,na_rep='NA')
    corr=[]
    for region,g in joined.groupby('supplement_genomic_region'):
        for stat in metrics:
            for track in ['A_distance','A_relative','mean_all_focal_pair_distances','q_S','q_H','q_3']:
                z=g[[stat,track]].dropna()
                if len(z)>2: corr.append(dict(region=region,statistic=stat,track=track,n=len(z),pearson=z[stat].corr(z[track]),spearman=z[stat].rank().corr(z[track].rank())))
    pd.DataFrame(corr).to_csv(OUT/'results/sinvicta_chr16_tree_track_correlations.tsv',sep='\t',index=False,na_rep='NA')
    fig,axs=plt.subplots(4,1,figsize=(12,9),sharex=True,constrained_layout=True)
    for ax,col,label in zip(axs,metrics,labels):
        ax.axvspan(11.680438,27.917498,color='#e7a34a',alpha=.14,lw=0)
        ax.axvline(11.680438,color='#a75c0b',lw=.9); ax.axvline(27.917498,color='#a75c0b',lw=.9)
        ax.scatter(d.tree_window_mid/1e6,d[col],s=12,color='#245b78',alpha=.8,rasterized=True); ax.set_ylabel(label); ax.grid(axis='y',alpha=.2)
    axs[-1].set_xlabel('Chromosome 16 local tree coordinate (Mb; Stolle / gng20170922wFex axis)')
    fig.suptitle('S. invicta Argentine males: SB vs Sb\nchr16nr supergene shaded; mapped BUSCO intervals; SNP panel, not per callable base',fontsize=13)
    fig.savefig(OUT/'figures/sinvicta_chr16_snp_panel_tracks.pdf'); fig.savefig(OUT/'figures/sinvicta_chr16_snp_panel_tracks.png',dpi=300)
    print(pd.DataFrame(rows).to_string(index=False))
if __name__=='__main__': main()
