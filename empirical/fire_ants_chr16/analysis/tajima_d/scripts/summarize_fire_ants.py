"""Regenerate fire-ant regional block summaries and matched-panel plot."""
from pathlib import Path
import numpy as np,pandas as pd,matplotlib.pyplot as plt
B=Path(__file__).resolve().parents[5]/'empirical/fire_ants_chr16/analysis/tajima_d'; d=pd.read_csv(B/'results/fire_ant_chr16_tajima_d_windows.tsv',sep='\t');d['inside']=d.interval_class.eq('chr16nr');d['block']=(d.mid//1_000_000).astype(int);rng=np.random.default_rng(20261121)
metrics=['pooled_tajimas_D','SB_tajimas_D','Sb_tajimas_D','matched_pooled_tajimas_D','matched_SB_tajimas_D','matched_Sb_tajimas_D','Delta_D_pooled_minus_background_mean','Delta_D_sample_weighted','matched_Delta_D_unweighted','matched_Delta_D_weighted'];rows=[]
for m in metrics:
 a=d[d.inside&d[m].notna()];b=d[(~d.inside)&d[m].notna()];av=a.groupby('block')[m].mean().values;bv=b.groupby('block')[m].mean().values;boot=[]
 if len(av)>1 and len(bv)>1:
  for _ in range(10000):boot.append(rng.choice(av,len(av),replace=True).mean()-rng.choice(bv,len(bv),replace=True).mean())
 rows.append({'statistic':m,'n_inside_windows':len(a),'n_outside_windows':len(b),'inside_mean':a[m].mean(),'outside_mean':b[m].mean(),'difference_inside_minus_outside':a[m].mean()-b[m].mean(),'block_bp':1_000_000,'n_inside_blocks':len(av),'n_outside_blocks':len(bv),'block_bootstrap_95_low':np.quantile(boot,.025) if boot else np.nan,'block_bootstrap_95_high':np.quantile(boot,.975) if boot else np.nan,'interpretation':'exploratory observed-panel contrast; no neutrality-test p-value'})
pd.DataFrame(rows).to_csv(B/'results/fire_ant_chr16_tajima_d_inside_outside.tsv',sep='\t',index=False,na_rep='NA')
fig,axs=plt.subplots(4,1,figsize=(12,9),sharex=True,constrained_layout=True)
metrics=['matched_pooled_tajimas_D','matched_SB_tajimas_D','matched_Sb_tajimas_D','matched_Delta_D_unweighted'];labels=['Pooled D, matched sites (n=108)','SB D, matched sites (n=78)','Sb D, matched sites (n=30)','Pooled − mean SB/Sb D, matched sites']
for ax,col,label in zip(axs,metrics,labels):
 ax.axvspan(11.680438,27.917498,color='#e7a34a',alpha=.16,lw=0);ax.scatter(d.mid/1e6,d[col],s=14,color='#28627b',alpha=.82);ax.axhline(0,color='.5',lw=.6);ax.set_ylabel(label);ax.grid(axis='y',alpha=.2)
axs[-1].set_xlabel('Local chr16 tree-window coordinate (Mb; BUSCO mapped VCF intervals)');fig.suptitle("S. invicta Argentine males: common-site Tajima's D contrasts\nExploratory; variant-only selected BUSCO gene panel; no callable-base denominator",fontsize=12)
for ext in ('pdf','svg','png'):fig.savefig(B/'figures'/f'fire_ant_chr16_tajima_d.{ext}',dpi=300)
print(pd.DataFrame(rows).to_string(index=False))
