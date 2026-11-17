#!/usr/bin/env python3
"""Manuscript figures A-D from the frozen nucleotide and exact-quartet outputs."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[6]
BASE=ROOT/'empirical/neoaves_chr4/analysis/patterson_d'
OUT=BASE/'final_validation/figures'
INTERVALS=((25.03,32.67),(33.51,34.47),(44.13,56.81))
COL={'Columbea_primary':'#0072B2','N61_primary':'#D55E00','N62_primary':'#009E73'}
LABEL={'Columbea_primary':'Columbea','N61_primary':'N61','N62_primary':'N62'}
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def style(ax):
 for a,b in INTERVALS:ax.axvspan(a,b,color='#999999',alpha=.13,lw=0)
 ax.axhline(0,color='#555555',lw=.7)
def save(fig,name):
 OUT.mkdir(parents=True,exist_ok=True)
 for ext in ('png','pdf','svg'):fig.savefig(OUT/f'{name}.{ext}',dpi=300,bbox_inches='tight')
 plt.close(fig)
def main():
 # A: nucleotide tracks plus per-bin informative-site density on the same source axis.
 win=read(BASE/'results/patterson_d_windows_250kb.tsv')
 fig,axs=plt.subplots(3,1,figsize=(11,7.2),sharex=True)
 for ax,d in zip(axs,COL):
  rr=sorted([r for r in win if r['design']==d and r['filter']=='all'],key=lambda r:int(r['window_start']))
  x=np.array([float(r['midpoint'])/1e6 for r in rr]); y=np.array([float(r['Patterson_D']) if r['Patterson_D'] and r['low_information']!='True' else np.nan for r in rr]); n=np.array([int(r['informative_sites']) for r in rr])
  style(ax);ax.plot(x,y,color=COL[d],lw=1.05);ax.set_ylabel('D');ax.set_title(LABEL[d],loc='left',fontsize=10,fontweight='bold');ax.set_ylim(-.65,.65)
  ax2=ax.twinx();ax2.fill_between(x,0,n,color='#777777',alpha=.15,step='mid');ax2.set_ylim(0,max(500,n.max()*1.15));ax2.set_yticks([0,int(max(500,n.max()))]);ax2.set_ylabel('Informative sites / 250 kb',fontsize=7);ax2.tick_params(axis='y',labelsize=7)
 axs[-1].set_xlabel('Stiller chr4 source window coordinate (Mb; frozen GalGal4 window-start predicate)')
 fig.suptitle('Nucleotide Patterson’s D across chr4',y=.995,fontsize=13);fig.text(.01,.01,'Gray bands: published PNAS interval membership. Source coordinates are retained; GalGal6 structural tracks are not overlaid.',fontsize=8)
 fig.tight_layout(rect=(0,.035,1,.97));save(fig,'figure_A_nucleotide_D_tracks')
 # B: exact matched Dsite and Dtopology at 1-Mb bins.
 comp=read(BASE/'final_validation/results/nucleotide_topology_comparison.tsv')
 fig,axs=plt.subplots(3,1,figsize=(11,7.2),sharex=True)
 for ax,d in zip(axs,COL):
  rr=sorted([r for r in comp if r['design']==d],key=lambda r:int(r['window_start']))
  x=np.array([int(r['window_start'])/1e6+.5 for r in rr]); ds=np.array([float(r['D_site_pooled']) if r['D_site_pooled'] else np.nan for r in rr]);dt=np.array([float(r['D_topology']) if r['D_topology'] else np.nan for r in rr])
  style(ax);ax.plot(x,ds,color=COL[d],lw=1.2,label='D_site (matched loci)');ax.plot(x,dt,color='#333333',lw=1,ls='--',label='D_topology = q3−q2')
  ax.set_ylabel('contrast');ax.set_title(LABEL[d],loc='left',fontsize=10,fontweight='bold');ax.legend(frameon=False,ncol=2,loc='upper right',fontsize=8)
 axs[-1].set_xlabel('Stiller chr4 source window coordinate (Mb)');fig.suptitle('Nucleotide and exact quartet genealogy statistics',y=.995,fontsize=13);fig.text(.01,.01,'Each 1-Mb bin uses only loci with all four exact species in the archived tree and saved nucleotide counts. Linked bins; no site-independence p-values.',fontsize=8)
 fig.tight_layout(rect=(0,.035,1,.97));save(fig,'figure_B_matched_site_and_genealogy')
 # C: inside/outside D and their difference, with 1-Mb delete-block normal intervals.
 fr=read(BASE/'final_validation/results/final_results.tsv');fr=[r for r in fr if r['filter']=='all' and r['design'] in COL]
 fig,(ax,dx)=plt.subplots(1,2,figsize=(10.5,4.6),gridspec_kw={'width_ratios':[1.25,1]})
 yloc=np.arange(3);offset={'inside':.13,'outside':-.13}
 for i,d in enumerate(COL):
  r=next(x for x in fr if x['design']==d);c=COL[d]
  for group,key,sekey,marker in [('inside','inside_D','inside_se','o'),('outside','outside_D','outside_se','s')]:
   y=i+offset[group];v=float(r[key]);se=float(r[sekey]);ax.errorbar(v,y,xerr=1.96*se,fmt=marker,color=c,ms=5,capsize=2,lw=1.1,label=group if i==0 else None)
  v=float(r['inside_minus_outside_D']);se=float(r['inside_minus_outside_se']);dx.errorbar(v,i,xerr=1.96*se,fmt='o',color=c,ms=5,capsize=2,lw=1.1)
 ax.axvline(0,color='#555',lw=.7);dx.axvline(0,color='#555',lw=.7)
 for a in (ax,dx):a.set_yticks(yloc);a.set_yticklabels([LABEL[d] for d in COL]);a.invert_yaxis();a.grid(axis='x',alpha=.2)
 ax.set_xlabel('Patterson’s D (estimate ± 1.96 × 1-Mb block-jackknife SE)');ax.legend(frameon=False,fontsize=8);dx.set_xlabel('D_inside − D_outside')
 fig.suptitle('Inside/outside estimates and regional contrasts',fontsize=13);fig.tight_layout();save(fig,'figure_C_inside_outside')
 # D: outgroup and transversion sensitivity; independent design labels, no pooling of replicate quartets.
 rows=read(BASE/'final_validation/results/final_results.tsv');rows=[r for r in rows if r['filter'] in ('all','transitions_only','transversions_only')]
 order=[]
 for d in ('Columbea_primary','Columbea_alt_outgroup','N61_primary','N61_alt_outgroup','N62_primary','N62_alt_outgroup'):
  for filt in ('all','transitions_only','transversions_only'):
   r=next((x for x in rows if x['design']==d and x['filter']==filt),None)
   if r:order.append((d,filt,r))
 fig,ax=plt.subplots(figsize=(10,7.2));ys=np.arange(len(order))
 for y,(d,filt,r) in zip(ys,order):
  val=float(r['inside_minus_outside_D']);se=float(r['inside_minus_outside_se']);col=COL[d.replace('_alt_outgroup','_primary')];marker={'all':'o','transitions_only':'^','transversions_only':'D'}[filt]
  ax.errorbar(val,y,xerr=1.96*se,fmt=marker,color=col,ms=5,capsize=2,lw=1.1)
 labels=[f"{LABEL[d.replace('_alt_outgroup','_primary')]} — {'Meleagris' if 'alt_outgroup' in d else 'Gallus'} — {'all sites' if f=='all' else ('transitions' if f=='transitions_only' else 'transversions')}" for d,f,r in order]
 ax.set_yticks(ys);ax.set_yticklabels(labels,fontsize=8);ax.invert_yaxis();ax.axvline(0,color='#555',lw=.8);ax.set_xlabel('D_inside − D_outside ± 1.96 × 1-Mb block-jackknife SE');ax.grid(axis='x',alpha=.2)
 fig.suptitle('Outgroup and transversion sensitivity',fontsize=13);fig.text(.01,.01,'Quartets share taxa/loci and are correlated sensitivity designs; transition/transversion counts are derived by paired subtraction.',fontsize=8);fig.tight_layout(rect=(0,.035,1,.96));save(fig,'figure_D_sensitivity')
if __name__=='__main__':main()
