#!/usr/bin/env python3
"""Plot corrected Anopheles D/f_dM and all-sites population-diversity tracks."""
from pathlib import Path
import csv,math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[5]
AN=ROOT/'empirical/anopheles_2la/analysis/joint_genomic_tracks'
RES=AN/'results';FIG=AN/'figures';INV=(20_524_058,42_165_532)
def tsv(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def num(x):
 try:
  y=float(x);return y if math.isfinite(y) else np.nan
 except (TypeError,ValueError):return np.nan
def save(fig,name):
 FIG.mkdir(parents=True,exist_ok=True)
 for ext in ('pdf','svg','png'):fig.savefig(FIG/f'{name}.{ext}',dpi=300 if ext=='png' else None)
 plt.close(fig)
def setup(ax,lo=0,hi=49.4):
 ax.axvspan(INV[0]/1e6,INV[1]/1e6,color='#8b78b6',alpha=.16,lw=0,zorder=0)
 ax.set_xlim(lo,hi);ax.grid(axis='y',color='.88',lw=.6);ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False)
def line(ax,rows,xkey,ykey,label,color=None,min_calls=None,call_key='callable_sites'):
 rows=sorted(rows,key=lambda r:int(r[xkey]));x=[(int(r['start'])+int(r['end']))/2/1e6 for r in rows];y=[num(r.get(ykey,'')) if min_calls is None or int(r.get(call_key,'0') or 0)>=min_calls else np.nan for r in rows]
 ax.plot(x,y,'.-',ms=2.3,lw=.7,label=label,color=color)
def main():
 FIG.mkdir(parents=True,exist_ok=True)
 win=tsv(RES/'anopheles_2la_fdm_D_corrected_windows.tsv')
 primary=[r for r in win if r['filter_label']=='primary']
 # Corrected D and f_dM are separated from the frozen legacy D result.
 for stat,key,name,designs in [('f_dM','f_dM','anopheles_2la_fdm_tracks',['A_merus','B_merus']),('Patterson D (allele-state corrected)','D_patterson_corrected','anopheles_2la_patterson_D_corrected_tracks',['A_merus','B_merus'])]:
  fig,axs=plt.subplots(len(designs),1,figsize=(11,5.3),sharex=True,constrained_layout=True)
  for ax,did in zip(np.atleast_1d(axs),designs):
   sub=[r for r in primary if r['test_id']==did];line(ax,sub,'start',key,did);setup(ax);ax.set_ylabel(stat);ax.set_title(did,loc='left',fontsize=9,fontweight='bold');ax.axhline(0,color='.4',lw=.6)
  axs[-1].set_xlabel('AgamP4 2L coordinate (Mb)');fig.suptitle(f'Anopheles 2La: {stat}',fontsize=13,fontweight='bold');save(fig,name)
 # Sensitivity estimates remain distinct from primary designs.
 sens=[r for r in primary if r['test_id'] not in ('A_merus','B_merus')]
 for key,name,label in [('f_dM','anopheles_2la_fdm_sensitivity_tracks','f_dM'),('D_patterson_corrected','anopheles_2la_patterson_D_sensitivity_corrected_tracks','Patterson D')]:
  fig,ax=plt.subplots(figsize=(11,4.8),constrained_layout=True)
  for d in sorted(set(r['test_id'] for r in sens)):
   line(ax,[r for r in sens if r['test_id']==d],'start',key,d)
  setup(ax);ax.axhline(0,color='.4',lw=.6);ax.set_ylabel(label);ax.set_xlabel('AgamP4 2L coordinate (Mb)');ax.set_title(f'Anopheles 2La: {label} sensitivity designs',fontweight='bold');ax.legend(ncol=3,fontsize=7,frameon=False);save(fig,name)
 # Four group-specific pi tracks.
 pi=tsv(RES/'anopheles_2la_pi_windows.tsv');groups=('G_s','G_i','C_s','C_i')
 fig,axs=plt.subplots(4,1,figsize=(11,8),sharex=True,sharey=True,constrained_layout=True)
 for ax,g in zip(axs,groups):
  line(ax,[r for r in pi if r['group']==g],'start','pi_per_callable_site',g,min_calls=25000);setup(ax);ax.set_ylabel('π');ax.set_title({'G_s':'gambiae standard','G_i':'gambiae inverted','C_s':'coluzzii standard','C_i':'coluzzii inverted'}[g],loc='left',fontsize=9,fontweight='bold')
 axs[-1].set_xlabel('AgamP4 2L coordinate (Mb)');fig.suptitle('Anopheles 2La: within-group nucleotide diversity',fontsize=13,fontweight='bold');save(fig,'anopheles_2la_pi_tracks')
 dxy=tsv(RES/'anopheles_2la_dxy_fst_windows.tsv')
 pairs=('primary_within_gambiae','primary_within_coluzzii')
 for key,name,ylab,title in [('dxy_per_callable_site','anopheles_2la_dxy_tracks','dXY per callable site','Pairwise nucleotide divergence'),('hudson_fst','anopheles_2la_hudson_fst_tracks',"Hudson's F_ST",'Hudson F_ST')]:
  fig,axs=plt.subplots(2,1,figsize=(11,6),sharex=True,constrained_layout=True)
  for ax,pair in zip(axs,pairs):
   line(ax,[r for r in dxy if r['comparison']==pair],'start',key,pair,min_calls=25000,call_key='hudson_fst_sites' if key=='hudson_fst' else 'callable_sites');setup(ax);ax.set_ylabel(ylab);ax.set_title(pair.replace('_',' '),loc='left',fontsize=9,fontweight='bold')
  axs[-1].set_xlabel('AgamP4 2L coordinate (Mb)');fig.suptitle(f'Anopheles 2La: {title}',fontsize=13,fontweight='bold');save(fig,name)
 # Secondary comparisons are a separate figure for transparency.
 for key,name,ylab,title in [('dxy_per_callable_site','anopheles_2la_dxy_sensitivity_tracks','dXY per callable site','Secondary cross-species/background dXY'),('hudson_fst','anopheles_2la_hudson_fst_sensitivity_tracks',"Hudson's F_ST",'Secondary cross-species/background Hudson F_ST')]:
  fig,axs=plt.subplots(6,1,figsize=(11,10),sharex=True,constrained_layout=True)
  for ax,pair in zip(axs,('cross_species_standard','cross_species_inverted','cross_species_mixed_1','cross_species_mixed_2','primary_within_gambiae','primary_within_coluzzii')):
   line(ax,[r for r in dxy if r['comparison']==pair],'start',key,pair,min_calls=25000,call_key='hudson_fst_sites' if key=='hudson_fst' else 'callable_sites');setup(ax);ax.set_ylabel(ylab);ax.set_title(pair.replace('_',' '),loc='left',fontsize=8,fontweight='bold')
  axs[-1].set_xlabel('AgamP4 2L coordinate (Mb)');fig.suptitle(f'Anopheles 2La: {title}',fontsize=13,fontweight='bold');save(fig,name)
 print('figures generated in',FIG)
if __name__=='__main__':main()

