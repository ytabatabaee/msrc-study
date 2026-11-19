#!/usr/bin/env python3
"""Join genealogy and population-genetic tracks by frozen 50-kb window ID."""
import csv,math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[5];AN=ROOT/'empirical/anopheles_2la/analysis/joint_genomic_tracks';R=AN/'results';FIG=AN/'figures';OLD=ROOT/'empirical/anopheles_2la/results';INV=(20_524_058,42_165_532)
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def val(r,k):
 try:
  x=float(r.get(k,''));return x if math.isfinite(x) else ''
 except:return ''
def main():
 topo={x['window_id']:x for x in read(OLD/'stage2r_quartet_topology_signal.tsv')};dist={x['window_id']:x for x in read(OLD/'stage12_distance/stage12_window_distance.tsv')}
 fdm=read(R/'anopheles_2la_fdm_D_corrected_windows.tsv');pi=read(R/'anopheles_2la_pi_windows.tsv');div=read(R/'anopheles_2la_dxy_fst_windows.tsv')
 base=read(ROOT/'data/anopheles_2la/processed/stage2r_window_grid.tsv');out=[]
 for w in base:
  wid=w['window_id'];t=topo.get(wid,{});d=dist.get(wid,{});r=dict(window_id=wid,chrom=w['chrom'],start=w['start'],end=w['end'],midpoint=(int(w['start'])+int(w['end']))/2,region_class=w['region_class'],q_species=t.get('q_species',''),q_arrangement=t.get('q_arrangement',''),q_third=t.get('q_third',''),D_topology=t.get('D',''),mean_internal_branch=d.get('mean_internal_branch',''),tree_diameter=d.get('tree_diameter',''),mean_pair_distance=d.get('mean_pair_distance',''))
  for design in ('A_merus','B_merus'):
   x=next((z for z in fdm if z['window_id']==wid and z['test_id']==design and z['filter_label']=='primary'),{})
   r['D_Patterson_corrected_'+design]=x.get('D_patterson_corrected','');r['f_dM_'+design]=x.get('f_dM','')
  for g in ('G_s','G_i','C_s','C_i'):
   x=next((z for z in pi if z['window_id']==wid and z['group']==g),{});r['pi_'+g]=x.get('pi_per_callable_site','');r['pi_sites_'+g]=x.get('callable_sites','')
  for pair in ('primary_within_gambiae','primary_within_coluzzii','cross_species_standard','cross_species_inverted','cross_species_mixed_1','cross_species_mixed_2'):
   x=next((z for z in div if z['window_id']==wid and z['comparison']==pair),{});short={'primary_within_gambiae':'Gs_Gi','primary_within_coluzzii':'Cs_Ci','cross_species_standard':'Gs_Cs','cross_species_inverted':'Gi_Ci','cross_species_mixed_1':'Gs_Ci','cross_species_mixed_2':'Gi_Cs'}[pair]
   r['dXY_'+short]=x.get('dxy_per_callable_site','');r['Hudson_FST_'+short]=x.get('hudson_fst','');r['dXY_sites_'+short]=x.get('callable_sites','');r['effective_pairs_'+short]=x.get('effective_allele_pairs','')
  out.append(r)
 write(R/'anopheles_2la_integrated_genomic_tracks.tsv',out)
 fig,axs=plt.subplots(9,1,figsize=(12,14),sharex=True,constrained_layout=True)
 xs=np.array([float(r['midpoint'])/1e6 for r in out]);
 def setup(ax):
  ax.axvspan(INV[0]/1e6,INV[1]/1e6,color='#8b78b6',alpha=.16,lw=0,zorder=0);ax.grid(axis='y',color='.88',lw=.55);ax.spines['top'].set_visible(False);ax.spines['right'].set_visible(False);ax.set_xlim(xs.min(),xs.max())
 def draw(ax,key,label,color=None):
  def gety(r):
   ck=('pi_sites_'+key[3:]) if key.startswith('pi_') else ('dXY_sites_'+key[4:]) if key.startswith('dXY_') else ('dXY_sites_'+key[11:]) if key.startswith('Hudson_FST_') else None
   if ck and int(r.get(ck,'0') or 0)<25000:return np.nan
   try:return float(r[key]) if r.get(key,'') not in ('',None) and math.isfinite(float(r[key])) else np.nan
   except:return np.nan
  y=np.array([gety(r) for r in out]);ax.plot(xs,y,'.-',ms=2,lw=.65,label=label,color=color)
 ax=axs[0]
 for k,l,c in [('q_species','q species','#286090'),('q_arrangement','q arrangement','#d07a25'),('q_third','q third','#58946c')]:draw(ax,k,l,c)
 setup(ax);ax.set_ylabel('Topology\nsupport');ax.legend(ncol=3,fontsize=7,frameon=False,loc='upper right')
 specs=[(1,[('D_topology','D_top','#333333')],'D_top'),(2,[('mean_internal_branch','Mean internal branch','#3a7d44')],'Mean internal\nbranch'),(3,[('tree_diameter','Tree diameter','#8064a2')],'Tree diameter'),(4,[(f'D_Patterson_corrected_{d}',d,'#34699a' if d=='A_merus' else '#c0504d') for d in ('A_merus','B_merus')],"Patterson's D"),(5,[(f'f_dM_{d}',d,'#34699a' if d=='A_merus' else '#c0504d') for d in ('A_merus','B_merus')],'f_dM'),(6,[(f'pi_{g}',g,{'G_s':'#34699a','G_i':'#c0504d','C_s':'#4f9563','C_i':'#bf8c2d'}[g]) for g in ('G_s','G_i','C_s','C_i')],'π per callable\nsite'),(7,[('dXY_Gs_Gi','G_s vs G_i','#34699a'),('dXY_Cs_Ci','C_s vs C_i','#c0504d')],'dXY per callable\nsite'),(8,[("Hudson_FST_Gs_Gi",'G_s vs G_i','#34699a'),('Hudson_FST_Cs_Ci','C_s vs C_i','#c0504d')],"Hudson F_ST")]
 for i,tracks,label in specs:
  ax=axs[i]
  for k,l,c in tracks:draw(ax,k,l,c)
  setup(ax);ax.set_ylabel(label);ax.axhline(0,color='.5',lw=.5);ax.legend(ncol=4,fontsize=7,frameon=False,loc='upper right')
 axs[-1].set_xlabel('AgamP4 2L coordinate (Mb)');fig.suptitle('Anopheles 2La: genealogical and population-genetic tracks',fontsize=13,fontweight='bold')
 FIG.mkdir(parents=True,exist_ok=True)
 for ext in ('pdf','svg','png'):fig.savefig(FIG/f'anopheles_2la_integrated_genomic_tracks.{ext}',dpi=300 if ext=='png' else None)
 plt.close(fig);print('joined',len(out),'windows; integrated figure saved')
if __name__=='__main__':main()

