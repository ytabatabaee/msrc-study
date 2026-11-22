"""Join cod Tajima D on existing source window ID to linkage and tree inventory."""
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[5]; B=ROOT/'empirical/atlantic_cod/analysis/tajima_d'; D=ROOT/'data/atlantic_cod/processed'
x=pd.read_csv(B/'results/atlantic_cod_tajima_d_250kb.tsv',sep='\t')
l=pd.read_csv(D/'cod_linkage_250kb.tsv',sep='\t');t=pd.read_csv(D/'cod_window_trees.tsv',sep='\t')[['window_id','n_taxa','tree_type','posterior_support_available']]
x=x.merge(l,left_on=['window_id','LG'],right_on=['window_id','lg'],how='left',suffixes=('','_linkage')).merge(t,on='window_id',how='left',suffixes=('','_tree'))
x.to_csv(B/'results/atlantic_cod_tajima_d_joined_tracks.tsv',sep='\t',index=False,na_rep='NA')
print('rows',len(x),'linkage matched',x.median_linkage_score.notna().sum(),'tree metadata matched',x.n_taxa.notna().sum())
