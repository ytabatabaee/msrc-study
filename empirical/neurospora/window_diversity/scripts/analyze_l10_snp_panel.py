#!/usr/bin/env python3
"""SNP-panel-conditioned diversity summaries for existing L10 tree windows.

These are explicitly per observed filtered variable SNP, never per callable base.
The input matrix is a globally filtered variant panel and lacks invariant sites.
"""
from pathlib import Path
import gzip
from bisect import bisect_left, bisect_right
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / 'empirical/neurospora/window_diversity'
GENO = ROOT / 'data/neurospora/raw/neurospora/geno/cras_sit_hist_tet_L10.min20MAC1.geno.gz'
WINDOWS = ROOT / 'empirical/neurospora/results/gene_tree_inventory.tsv'
A = ['CJ01','CJ03','CJ05','CJ09','CJ13','CJ15','CJ19','CJ23']
A = ['CJ01','CJ03','CJ05','CJ09','CJ13','CJ15','CJ19','CJ23']
a = ['CJ02','CJ04','CJ06','CJ08','CJ10','CJ14','CJ16','CJ20','CJ22','CJ24']
# These are the directly labeled/pair-resolved matrix samples: 8 A, 10 a. Four
# additional CJ isolates have unresolved mating type and are excluded.


def load_matrix():
    with gzip.open(GENO, 'rt') as fh:
        header = fh.readline().strip().split()
        cols = header[2:]
        rows=[]
        for line in fh:
            x=line.split()
            if len(x)==len(header): rows.append((x[0],int(x[1]),x[2:]))
    return cols, rows


def stats(x, y):
    """Unbiased within diversity and Hudson numerator/denominator for SNPs."""
    n1,n2=len(x),len(y)
    alleles=sorted(set(x)|set(y))
    p=np.array([sum(v==z for v in x)/n1 for z in alleles],float)
    q=np.array([sum(v==z for v in y)/n2 for z in alleles],float)
    # unbiased gene diversity for haploid chromosomes, conditional on observed SNP
    pi1=n1/(n1-1)*(1-(p*p).sum()) if n1>1 else np.nan
    pi2=n2/(n2-1)*(1-(q*q).sum()) if n2>1 else np.nan
    dxy=1-float((p*q).sum())
    # Hudson: dXY - (piX+piY)/2; window ratio of summed site components
    num=dxy-(pi1+pi2)/2
    return pi1,pi2,dxy,num


def main():
    cols, rows=load_matrix(); ix={c:i for i,c in enumerate(cols)}
    w=pd.read_csv(WINDOWS,sep='\t')
    records=[]
    # gene_tree_inventory coordinates define min/max positions of 50 SNPs per window;
    # join directly by linkage group and site coordinate.
    by={}
    for chrom,pos,calls in rows:
        by.setdefault(chrom,[]).append((pos,calls))
    positions={chrom:[p for p,_ in sites] for chrom,sites in by.items()}
    for _,win in w.iterrows():
        chrom_sites=by.get(win.linkage_group,[]); chrom_pos=positions.get(win.linkage_group,[])
        lo=bisect_left(chrom_pos, int(win.start)); hi=bisect_right(chrom_pos, int(win.end))
        sites=chrom_sites[lo:hi]
        if not sites: continue
        vals=np.array([[c[ix[s]] for s in A] + [c[ix[s]] for s in a] for _,c in sites])
        mask=np.all(np.isin(vals, ['A','C','G','T']),axis=1)
        vals=vals[mask]
        if len(vals)==0: continue
        site=[]
        for row in vals:
            pi1,pi2,dxy,num=stats(list(row[:len(A)]),list(row[len(A):]))
            site.append((pi1,pi2,dxy,num))
        z=np.asarray(site)
        den=z[:,2].sum()
        fst=z[:,3].sum()/den if den>0 else np.nan
        records.append(dict(window=int(win.window),linkage_group=win.linkage_group,
            start=int(win.start),end=int(win.end),mid=int(win.mid),
            n_variable_sites=len(z),n_mat_A=len(A),n_mat_a=len(a),
            pi_mat_A_per_observed_SNP=z[:,0].mean(),pi_mat_a_per_observed_SNP=z[:,1].mean(),
            dXY_per_observed_SNP=z[:,2].mean(),Hudson_FST_SNP_panel=fst,
            inside_LGI_1p5_8Mb=int(win.linkage_group=='LGI' and 1_500_000<=win.mid<=8_000_000),
            q1_species=win.q1_species,q2_other=win.q2_other,q3_hisp_A=win.q3_hisp_A))
    out=pd.DataFrame(records)
    out.to_csv(OUT/'results/l10_snp_panel_window_stats.tsv',sep='\t',index=False,na_rep='NA')
    # Spatial block bootstrap of inside-minus-flank means; resample 100-kb blocks.
    lgi=out[out.linkage_group=='LGI'].copy(); lgi['block']=lgi.mid//100_000
    metrics=['Hudson_FST_SNP_panel','dXY_per_observed_SNP','pi_mat_A_per_observed_SNP','pi_mat_a_per_observed_SNP']
    rng=np.random.default_rng(20261120); summary=[]
    for col in metrics:
        inside=lgi[lgi.inside_LGI_1p5_8Mb==1]; outside=lgi[lgi.inside_LGI_1p5_8Mb==0]
        obs=inside[col].mean()-outside[col].mean()
        ib=inside.groupby('block')[col].mean().dropna().to_numpy(); ob=outside.groupby('block')[col].mean().dropna().to_numpy()
        boot=[]
        if len(ib)>1 and len(ob)>1:
            for _ in range(10000): boot.append(rng.choice(ib,len(ib),replace=True).mean()-rng.choice(ob,len(ob),replace=True).mean())
        lo,hi=np.quantile(boot,[.025,.975]) if boot else (np.nan,np.nan)
        summary.append(dict(statistic=col,inside_windows=len(inside),outside_windows=len(outside),
            inside_mean=inside[col].mean(),outside_mean=outside[col].mean(),difference_inside_minus_outside=obs,
            spatial_100kb_block_bootstrap_95_low=lo,spatial_100kb_block_bootstrap_95_high=hi,
            n_inside_blocks=len(ib),n_outside_blocks=len(ob)))
    pd.DataFrame(summary).to_csv(OUT/'results/l10_inside_outside_spatial.tsv',sep='\t',index=False,na_rep='NA')
    # Publication tracks, separate panels, missing windows remain gaps.
    fig,axs=plt.subplots(4,1,figsize=(12,9),sharex=True,constrained_layout=True)
    for ax,col,label in zip(axs,metrics,['Hudson FST (SNP panel)','dXY per observed SNP','π mat-A per observed SNP','π mat-a per observed SNP']):
        x=lgi.mid/1e6; y=lgi[col]
        ax.scatter(x,y,s=5,color='#245b78',alpha=.7,rasterized=True)
        ax.axvspan(1.5,8.0,color='#e7a34a',alpha=.20,lw=0)
        ax.axvline(1.5,color='#a75c0b',lw=.8); ax.axvline(8.0,color='#a75c0b',lw=.8)
        ax.set_ylabel(label); ax.grid(axis='y',alpha=.2)
    axs[-1].set_xlabel('LGI coordinate (Mb; local tree-window axis)')
    fig.suptitle('Neurospora L10 diversity on filtered variable-SNP panel\nNot callable-base-normalized',fontsize=13)
    fig.savefig(OUT/'figures/l10_snp_panel_tracks.pdf'); fig.savefig(OUT/'figures/l10_snp_panel_tracks.png',dpi=300)
    # Machine-readable per-chromosome tracks, retaining q-frequency covariates.
    print(f'windows={len(out)}; mat-A={len(A)} mat-a={len(a)}; LGI={len(lgi)}')
    print(pd.DataFrame(summary).to_string(index=False))

if __name__=='__main__': main()
