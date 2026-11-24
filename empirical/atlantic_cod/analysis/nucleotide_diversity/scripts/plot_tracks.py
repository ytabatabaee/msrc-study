#!/usr/bin/env python3
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import pearsonr,spearmanr
ROOT=Path(__file__).resolve().parents[5];BASE=ROOT/'empirical/atlantic_cod/analysis/nucleotide_diversity';RES=BASE/'results';FIG=BASE/'figures'
def tsv(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def num(x):
 try:return float(x)
 except:return np.nan
def main():
 rows=tsv(RES/'atlantic_cod_pi_dxy_fst_windows.tsv');regions={r['lg']:r for r in tsv(ROOT/'data/atlantic_cod/metadata/region_manifest.tsv') if r['lg'] in ['LG01','LG02','LG07','LG12']}
 fq={r['lg']+'|'+r['window_id']:r for r in tsv(ROOT/'data/atlantic_cod/processed/stage4a_window_topology_scores.tsv')}
 bt=tsv(ROOT/'data/atlantic_cod/processed/stage7_pairwise_time_shifts.tsv');bagg={}
 for r in bt:
  if r['pair_class']=='opposite_arrangement':bagg.setdefault(r['lg']+'|'+r['window_id'],[]).append(num(r['delta_time']))
 corr=[];FIG.mkdir(parents=True,exist_ok=True)
 for lg in ['LG01','LG02','LG07','LG12']:
  g=[r for r in rows if r['chromosome']==lg];fig,axs=plt.subplots(4,1,figsize=(12,10),sharex=True)
  for ax,col,title in zip(axs,['pi_X_SNP','pi_Y_SNP','dXY_SNP','FST_Hudson'],['πX SNP-conditional','πY SNP-conditional','dXY SNP-conditional','Hudson FST']):
   x=np.array([num(r[col]) for r in g]);pos=np.array([(int(r['start'])+int(r['end']))/2/1e6 for r in g]);ok=np.isfinite(x)
   ax.plot(pos[ok],x[ok],'.-',lw=.8,ms=3,color='#215a7a');ax.set_ylabel(title);ax.grid(alpha=.2)
  r=regions[lg];axs[0].axvspan(int(r['start'])/1e6,int(r['end'])/1e6,color='#e6ad4b',alpha=.23)
  axs[-1].set_xlabel('gadMor2 position (Mb)');fig.suptitle(f'{lg}: arrangement-specific SNP-conditional diversity\nX = ancestral; Y = derived; shading = suppressed-recombination region');fig.tight_layout();fig.savefig(FIG/f'{lg}_pi_dxy_fst_tracks.png',dpi=180);fig.savefig(FIG/f'{lg}_pi_dxy_fst_tracks.pdf');plt.close(fig)
  g=[r for r in g if r['inside_suppressed_recombination']=='true' and r['pi_X_SNP']!='NA']
  for track in ['pi_X_SNP','pi_Y_SNP','dXY_SNP','FST_Hudson']:
   for name,get in [('quartet_arrangement_fraction',lambda z:num(fq[z['chromosome']+'|'+z['window_id']]['fraction_arrangement'])),('opposite_arrangement_mean_delta_time',lambda z:np.nanmean(bagg.get(z['chromosome']+'|'+z['window_id'],[np.nan])))]:
    pairs=[(num(z[track]),get(z)) for z in g if z['chromosome']+'|'+z['window_id'] in fq and np.isfinite(get(z)) and np.isfinite(num(z[track]))]
    if len(pairs)>=3:
     x,y=map(np.array,zip(*pairs))
     if np.std(x)==0 or np.std(y)==0:
      pr=(np.nan,np.nan);sr=(np.nan,np.nan)
     else:
      pr=pearsonr(x,y);sr=spearmanr(x,y)
     corr.append({'lg':lg,'diversity_track':track,'existing_track':name,'n_windows':len(x),'pearson_r':pr[0],'pearson_p':pr[1],'spearman_rho':sr[0],'spearman_p':sr[1]})
 # Overview figure, one row per LG and one column per statistic.
 fig,axs=plt.subplots(4,4,figsize=(16,10),sharex='col')
 for i,lg in enumerate(['LG01','LG02','LG07','LG12']):
  sub=[r for r in rows if r['chromosome']==lg];reg=regions[lg]
  for j,col in enumerate(['pi_X_SNP','pi_Y_SNP','dXY_SNP','FST_Hudson']):
   ax=axs[i,j];pos=np.array([(int(r['start'])+int(r['end']))/2/1e6 for r in sub]);val=np.array([num(r[col]) for r in sub]);ok=np.isfinite(val)
   ax.axvspan(int(reg['start'])/1e6,int(reg['end'])/1e6,color='#e6ad4b',alpha=.22);ax.plot(pos[ok],val[ok],'.-',lw=.7,ms=2.5,color='#215a7a');ax.grid(alpha=.15)
   if i==0:ax.set_title(['πX','πY','dXY','Hudson FST'][j])
   if j==0:ax.set_ylabel(lg)
   if i==3:ax.set_xlabel('Position (Mb)')
 fig.suptitle('Cod arrangement-specific SNP-conditional diversity; X=ancestral, Y=derived\nValues only within fully contained suppressed-recombination windows',y=1.01);fig.tight_layout();fig.savefig(FIG/'atlantic_cod_pi_tracks.png',dpi=180,bbox_inches='tight');fig.savefig(FIG/'atlantic_cod_pi_tracks.pdf',bbox_inches='tight');plt.close(fig)
 with (RES/'comparison_with_quartet_branch_tracks.tsv').open('w',newline='') as f:
  cols=['lg','diversity_track','existing_track','n_windows','pearson_r','pearson_p','spearman_rho','spearman_p'];w=csv.DictWriter(f,fieldnames=cols,delimiter='\t');w.writeheader();w.writerows(corr)
if __name__=='__main__':main()
