"""Spatial block-bootstrap correlations between Tajima D and frozen tracks."""
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[5]; B=ROOT/'empirical/anopheles_2la/analysis/tajima_d'; J=ROOT/'empirical/anopheles_2la/analysis/joint_genomic_tracks/results'
d=pd.read_csv(B/'results/anopheles_tajima_d_50kb.tsv',sep='\t'); t=pd.read_csv(J/'anopheles_2la_integrated_genomic_tracks.tsv',sep='\t')
d=d.merge(t,on='window_id',how='left',suffixes=('','_existing'));d['block']=(d.mid//500000).astype(int)
cols=['pi_G_s','pi_G_i','pi_C_s','pi_C_i','dXY_Gs_Gi','Hudson_FST_Gs_Gi','dXY_Cs_Ci','Hudson_FST_Cs_Ci','q_species','q_arrangement','q_third','mean_internal_branch','tree_diameter','mean_pair_distance','D_Patterson_corrected_A_merus','D_Patterson_corrected_B_merus']
rng=np.random.default_rng(20261121); out=[]
for group,g in d.groupby('group'):
  for col in cols:
    if col not in g: continue
    z=g[['tajimas_D',col,'block']].copy();z['x']=pd.to_numeric(z.tajimas_D,errors='coerce');z['y']=pd.to_numeric(z[col],errors='coerce');z=z.dropna(subset=['x','y'])
    if len(z)<10: continue
    rho=spearmanr(z.x,z.y).statistic; blks=[x for _,x in z.groupby('block')];boots=[]
    if len(blks)>2:
      for _ in range(500):
        pick=rng.integers(0,len(blks),len(blks)); xx=np.concatenate([blks[i].x.to_numpy() for i in pick]); yy=np.concatenate([blks[i].y.to_numpy() for i in pick])
        if len(np.unique(xx))>1 and len(np.unique(yy))>1: boots.append(spearmanr(xx,yy).statistic)
    out.append({'group':group,'track':col,'n_windows':len(z),'n_500kb_blocks':len(blks),'spearman_rho':rho,'block_bootstrap_95_low':np.quantile(boots,.025) if boots else np.nan,'block_bootstrap_95_high':np.quantile(boots,.975) if boots else np.nan,'interpretation':'spatial block bootstrap; descriptive association, no neutrality p-value; linked overlapping estimands'})
pd.DataFrame(out).to_csv(B/'results/anopheles_tajima_d_track_relationships.tsv',sep='\t',index=False,na_rep='NA')
d.to_csv(B/'results/anopheles_tajima_d_joined_tracks.tsv',sep='\t',index=False,na_rep='NA')
print('relationships',len(out),'groups',d.group.nunique(),'windows per group',d.groupby('group').size().describe().to_dict())
