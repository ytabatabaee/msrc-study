"""Descriptive pooled-minus-background contrasts from completed Anopheles tracks."""
from pathlib import Path
import pandas as pd,numpy as np,matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[5]; B=ROOT/'empirical/anopheles_2la/analysis/tajima_d'
rows=[];summary=[];rng=np.random.default_rng(20261121)
for res in (50000,100000):
 d=pd.read_csv(B/f'results/anopheles_tajima_d_{res//1000}kb.tsv',sep='\t')
 for sp in ('gambiae','coluzzii'):
  g=d[d.group.isin(['species_pooled_'+sp,'species_standard_'+sp,'species_inverted_'+sp])]
  p=g[g.group=='species_pooled_'+sp].set_index('window_id');x=g[g.group=='species_standard_'+sp].set_index('window_id');y=g[g.group=='species_inverted_'+sp].set_index('window_id')
  z=p[['start','end','mid','interval_class','tajimas_D','n_haplotypes']].join(x[['tajimas_D','n_haplotypes']],rsuffix='_standard',how='inner').join(y[['tajimas_D','n_haplotypes']],rsuffix='_inverted',how='inner')
  z=z.rename(columns={'tajimas_D':'D_pooled','n_haplotypes':'n_pooled','tajimas_D_standard':'D_standard','n_haplotypes_standard':'n_standard','tajimas_D_inverted':'D_inverted','n_haplotypes_inverted':'n_inverted'})
  z['resolution_bp']=res;z['species']=sp;z['group_definition']='pooled species across geography (includes heterokaryotypes) contrasted descriptively with homozygous arrangement groups; not the same sample panel'
  z['Delta_D_unweighted']=z.D_pooled-(z.D_standard+z.D_inverted)/2
  z['Delta_D_sample_weighted']=z.D_pooled-(z.n_standard*z.D_standard+z.n_inverted*z.D_inverted)/(z.n_standard+z.n_inverted)
  z['status']='VALIDATED within-group SFS; arrangement contrast is geographically confounded and descriptive'
  rows.extend(z.reset_index().to_dict('records'))
  for m in ('Delta_D_unweighted','Delta_D_sample_weighted'):
    a=z[z.interval_class=='inside_2La'].dropna(subset=[m]).copy();b=z[z.interval_class=='outside_2La'].dropna(subset=[m]).copy();a['block']=a.mid//500000;b['block']=b.mid//500000
    av=a.groupby('block')[m].mean().values;bv=b.groupby('block')[m].mean().values;boots=[]
    if len(av)>1 and len(bv)>1:
     for _ in range(5000):boots.append(rng.choice(av,len(av),replace=True).mean()-rng.choice(bv,len(bv),replace=True).mean())
    summary.append({'resolution_bp':res,'species':sp,'statistic':m,'inside_mean':a[m].mean(),'outside_mean':b[m].mean(),'difference_inside_minus_outside':a[m].mean()-b[m].mean(),'n_inside_windows':len(a),'n_outside_windows':len(b),'n_inside_500kb_blocks':len(av),'n_outside_500kb_blocks':len(bv),'block_bootstrap_95_low':np.quantile(boots,.025) if boots else np.nan,'block_bootstrap_95_high':np.quantile(boots,.975) if boots else np.nan,'interpretation':'descriptive spatial uncertainty; groups span geography'})
out=pd.DataFrame(rows); out.to_csv(B/'results/anopheles_tajima_d_background_contrasts.tsv',sep='\t',index=False,na_rep='NA')
pd.DataFrame(summary).to_csv(B/'results/anopheles_tajima_d_background_contrast_inside_outside.tsv',sep='\t',index=False,na_rep='NA')
fig,axs=plt.subplots(2,1,figsize=(12,6.5),sharex=True,constrained_layout=True)
for ax,sp in zip(axs,('gambiae','coluzzii')):
 z=out[(out.resolution_bp==50000)&(out.species==sp)].sort_values('mid');ax.axvspan(20.524058,42.165532,color='#e8aa46',alpha=.18);ax.scatter(z.mid/1e6,z.Delta_D_unweighted,s=10,color='#39758f',label='unweighted contrast');ax.plot(z.mid/1e6,z.Delta_D_sample_weighted,color='#ad612a',lw=.7,label='sample-size weighted');ax.axhline(0,color='.5',lw=.6);ax.set_ylabel(sp+"\nΔ Tajima's D");ax.grid(axis='y',alpha=.2);ax.legend(frameon=False,fontsize=8)
axs[-1].set_xlabel('AgamP4 2L coordinate (Mb)');fig.suptitle("Descriptive pooled-minus-background Tajima's D\nPooled species includes geographic and heterokaryotype mixture",fontsize=11)
for ext in ('pdf','svg','png'):fig.savefig(B/'figures'/f'anopheles_tajima_d_background_contrast.{ext}',dpi=300)
print(pd.DataFrame(summary).to_string(index=False))
