"""Build cross-dataset status/effect table and descriptive block-interval plot."""
from pathlib import Path
import pandas as pd,numpy as np,matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[4]; RES=ROOT/'empirical/analysis/tajimas_d/results'; rows=[]
def add(ds,group,status,n,inside,outside,diff,lo,hi,notes):rows.append({'dataset':ds,'group':group,'status':status,'n_individuals_or_haplotypes':n,'inside_mean_D':inside,'outside_mean_D':outside,'inside_minus_outside_D':diff,'block_95_low':lo,'block_95_high':hi,'notes':notes})
p=ROOT/'empirical/neurospora/analysis/tajima_d/results/neurospora_tajima_d_region_summary.tsv'
if p.exists():
 d=pd.read_csv(p,sep='\t');d=d[(d.resolution=='50kb') & d.statistic.isin(['D_Tajima_pooled','D_Tajima_X','D_Tajima_Y','Delta_D_Tajima'])];names={'D_Tajima_pooled':'pooled','D_Tajima_X':'mat-A','D_Tajima_Y':'mat-a','Delta_D_Tajima':'pooled − background mean'}
 for r in d.itertuples():add('Neurospora L10',names[r.statistic],'EXPLORATORY','18 / 8 / 10',r.inside_mean,r.outside_mean,r.difference_inside_minus_flank,r.block_bootstrap_95_low,r.block_bootstrap_95_high,'variable-site matrix; 100-kb spatial blocks')
p=ROOT/'empirical/fire_ants_chr16/analysis/tajima_d/results/fire_ant_chr16_tajima_d_inside_outside.tsv'
if p.exists():
 d=pd.read_csv(p,sep='\t');names={'pooled_tajimas_D':'pooled (108)','SB_tajimas_D':'SB (78)','Sb_tajimas_D':'Sb (30)','matched_Delta_D_unweighted':'matched pooled − background mean'}
 for r in d[d.statistic.isin(names)].itertuples():add('Fire ants chr16',names[r.statistic],'EXPLORATORY','78 SB / 30 Sb',r.inside_mean,r.outside_mean,r.difference_inside_minus_outside,r.block_bootstrap_95_low,r.block_bootstrap_95_high,'variant-only gene-ascertained panel; 1-Mb blocks; matched Δ row uses common complete sites')
p=ROOT/'empirical/anopheles_2la/analysis/tajima_d/results/anopheles_tajima_d_inside_outside.tsv'
if p.exists():
 d=pd.read_csv(p,sep='\t');d=d[(d.resolution_bp==50000)&(d.statistic=='tajimas_D')&d.group.isin(['species_pooled_gambiae','species_standard_gambiae','species_inverted_gambiae','species_pooled_coluzzii','species_standard_coluzzii','species_inverted_coluzzii'])]
 for r in d.itertuples():add('Anopheles 2La',r.group.replace('species_',''),'VALIDATED',r.group,r.inside_mean,r.outside_mean,r.inside_minus_outside,r.block_bootstrap_95_low,r.block_bootstrap_95_high,'Ag3.10 all-site; 500-kb spatial blocks; arrangement groups span geography')
p=ROOT/'empirical/anopheles_2la/analysis/tajima_d/results/anopheles_tajima_d_background_contrast_inside_outside.tsv'
if p.exists():
 d=pd.read_csv(p,sep='\t');d=d[(d.resolution_bp==50000)&(d.statistic=='Delta_D_unweighted')]
 for r in d.itertuples():add('Anopheles 2La','pooled − background mean '+r.species,'VALIDATED',r.species,r.inside_mean,r.outside_mean,r.difference_inside_minus_outside,r.block_bootstrap_95_low,r.block_bootstrap_95_high,'descriptive contrast; pooled species includes structure/heterokaryotypes; separate group masks')
p=ROOT/'empirical/atlantic_cod/analysis/tajima_d/results/atlantic_cod_tajima_d_inside_outside.tsv'
if p.exists():
 d=pd.read_csv(p,sep='\t');d=d[d.population_code=='avc']
 for r in d.itertuples():add('Atlantic cod',f'{r.LG} More_stationary','EXPLORATORY','2 diploid fish',r.inside_mean,np.nan,np.nan,np.nan,np.nan,'250-kb inside windows only; outside Relate VCF has zero records on target LGs')
feas=pd.read_csv(RES/'dataset_feasibility.tsv',sep='\t')
for r in feas.itertuples():
 if r.dataset not in ['Anopheles 2La','Atlantic cod','Fire ants chr16','Neurospora L10']:
  add(r.dataset,'not estimated',r.status,'',np.nan,np.nan,np.nan,np.nan,np.nan,'See repository-wide feasibility inventory')
out=pd.DataFrame(rows);out.to_csv(RES/'cross_dataset_tajima_d_summary.tsv',sep='\t',index=False,na_rep='NA')
plot=out[out.inside_minus_outside_D.notna()].copy();fig,ax=plt.subplots(figsize=(10,max(6,.35*len(plot))));y=np.arange(len(plot));x=plot.inside_minus_outside_D.to_numpy();lo=x-plot.block_95_low.to_numpy();hi=plot.block_95_high.to_numpy()-x;colors=['#287a51' if s=='VALIDATED' else '#d58b27' for s in plot.status]
for yi,xi,li,hii,c in zip(y,x,lo,hi,colors):ax.errorbar(xi,yi,xerr=[[li],[hii]],fmt='o',color=c,capsize=2,lw=1.2)
ax.axvline(0,color='.4',lw=.8);ax.set_yticks(y);ax.set_yticklabels(plot.dataset+' — '+plot.group,fontsize=8);ax.invert_yaxis();ax.set_xlabel('Inside − outside Tajima’s D (block-bootstrap 95% interval)');ax.set_title('Descriptive Tajima’s D regional contrasts\nDifferent sampling designs and ascertainment; not a pooled meta-analysis',loc='left',weight='bold');ax.grid(axis='x',alpha=.2);fig.tight_layout()
outdir=ROOT/'empirical/analysis/tajimas_d/figures'
for ext in ('pdf','svg','png'):fig.savefig(outdir/f'tajimas_d_cross_dataset_effects.{ext}',dpi=300,bbox_inches='tight')
print('rows',len(out),'contrasts',len(plot))
