#!/usr/bin/env python3
"""Plot full-chr17 window means inside versus outside inversion intervals."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[4]
BASE=ROOT/'empirical/house_mouse_t_complex/window_diversity'
INPUT=BASE/'results/mouse_t_pseudohap_chr17_full_window_stats.tsv'
OUT_TSV=BASE/'results/mouse_t_pseudohap_chr17_inversion_barplot_summary.tsv'
OUT_PNG=BASE/'figures/mouse_t_pseudohap_chr17_inversion_barplots.png'
OUT_PDF=BASE/'figures/mouse_t_pseudohap_chr17_inversion_barplots.pdf'
BLOCK=500_000
NBOOT=10_000
SEED=20261120
METRICS=[
    ('Hudson_FST_SNP_panel','Hudson FST'),
    ('dXY_per_observed_SNP','dXY per observed SNP'),
    ('pi_t_pseudohap_per_observed_SNP','π t pseudo-haplotype per observed SNP'),
    ('pi_standard_per_observed_SNP','π standard per observed SNP'),
]
REGIONS=['inversion_union','outside_inversions']
COLORS={'inversion_union':'#D98123','outside_inversions':'#4C78A8'}
LABELS={'inversion_union':'Inside inversions','outside_inversions':'Outside inversions'}

def main():
    d=pd.read_csv(INPUT,sep='\t')
    d['block']=np.floor(d.midpoint_bp/BLOCK).astype(int)
    rng=np.random.default_rng(SEED)
    summary=[]
    boot_ci={}
    for design,g in d.groupby('design',sort=True):
        for metric,label in METRICS:
            for region in REGIONS:
                sub=g[g.region==region]
                values=sub[metric].dropna()
                block_means=sub.groupby('block')[metric].mean().dropna().to_numpy()
                point=float(values.mean()) if len(values) else np.nan
                if len(block_means)>1:
                    sims=rng.choice(block_means,size=(NBOOT,len(block_means)),replace=True).mean(axis=1)
                    low,high=np.quantile(sims,[.025,.975])
                else:
                    low=high=np.nan
                summary.append(dict(design=design,metric=metric,label=label,region=region,
                    n_windows=len(values),n_spatial_blocks=len(block_means),
                    mean_window_estimate=point,block_bootstrap_95_low=low,
                    block_bootstrap_95_high=high,block_width_bp=BLOCK,n_bootstrap=NBOOT))
                boot_ci[(design,metric,region)]=(point,float(low),float(high))
    s=pd.DataFrame(summary)
    s.to_csv(OUT_TSV,sep='\t',index=False,na_rep='NA',float_format='%.8g')

    designs=list(d.design.drop_duplicates())
    ticklabels=[x.split('_')[1] for x in designs]
    x=np.arange(len(designs)); width=.34
    fig,axs=plt.subplots(2,2,figsize=(12,8),constrained_layout=True)
    for ax,(metric,label) in zip(axs.flat,METRICS):
        for offset,region in [(-width/2,REGIONS[0]),(width/2,REGIONS[1])]:
            means=[]; lows=[]; highs=[]
            for design in designs:
                point,lo,hi=boot_ci[(design,metric,region)]
                means.append(point); lows.append(max(point-lo,0) if np.isfinite(lo) else 0); highs.append(max(hi-point,0) if np.isfinite(hi) else 0)
            ax.bar(x+offset,means,width=width,color=COLORS[region],edgecolor='white',linewidth=.5,
                   label=LABELS[region],zorder=2)
            ax.errorbar(x+offset,means,yerr=np.vstack([lows,highs]),fmt='none',ecolor='#222222',
                        elinewidth=1.05,capsize=3,zorder=3)
        ax.set_title(label,loc='left',fontsize=11,fontweight='bold')
        ax.set_xticks(x,ticklabels)
        ax.set_ylabel('Mean window estimate')
        ax.grid(axis='y',alpha=.22,zorder=0)
        ax.spines[['top','right']].set_visible(False)
    axs[0,0].legend(frameon=False,fontsize=9,loc='upper right')
    fig.suptitle('House mouse chr17: mean diversity and differentiation inside vs outside inversions',fontsize=14)
    fig.text(.5,-.015,'Bars: arithmetic mean of valid 5-kb window estimates. Error bars: 95% percentile intervals from 10,000 resamples of 500-kb spatial block means. SNP-panel values are per observed SNP, not per callable base.',ha='center',fontsize=8.5)
    fig.savefig(OUT_PNG,dpi=300,bbox_inches='tight')
    fig.savefig(OUT_PDF,bbox_inches='tight')
    print(s.to_string(index=False))

if __name__=='__main__': main()
