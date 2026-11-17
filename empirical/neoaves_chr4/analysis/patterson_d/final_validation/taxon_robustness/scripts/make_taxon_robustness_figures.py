#!/usr/bin/env python3
"""Build Figures E-G for the N61 taxon sensitivity validation."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
ROOT=Path(__file__).resolve().parents[7]
PKG=ROOT/'empirical/neoaves_chr4/analysis/patterson_d/final_validation/taxon_robustness'
OUT=PKG/'figures'
ORDER=['N61_primary','N61_P1_Cuculus','N61_P1_Ardeotis','N61_P2_Pterocles','N61_P2_Mesitornis','N61_P3_Balaeniceps','N61_P3_Nipponia','N61_P3_Aquila','N61_combo_Cuculus_Pterocles','N61_combo_Ardeotis_Mesitornis','N61_combo_Cuculus_Pterocles_Balaeniceps','N61_primary_Meleagris']
LABELS={'N61_primary':'Original: Tauraco / Columba / Ciconia / Gallus','N61_P1_Cuculus':'P1=Cuculus; rest original','N61_P1_Ardeotis':'P1=Ardeotis; rest original','N61_P2_Pterocles':'P2=Pterocles; rest original','N61_P2_Mesitornis':'P2=Mesitornis; rest original','N61_P3_Balaeniceps':'P3=Balaeniceps; rest original','N61_P3_Nipponia':'P3=Nipponia; rest original','N61_P3_Aquila':'P3=Aquila (Telluraves control)','N61_combo_Cuculus_Pterocles':'P1=Cuculus + P2=Pterocles','N61_combo_Ardeotis_Mesitornis':'P1=Ardeotis + P2=Mesitornis','N61_combo_Cuculus_Pterocles_Balaeniceps':'P1=Cuculus + P2=Pterocles + P3=Balaeniceps','N61_primary_Meleagris':'Original ingroup; Meleagris outgroup'}

def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def save(fig,name):
 OUT.mkdir(parents=True,exist_ok=True)
 for ext in ('png','pdf','svg'):fig.savefig(OUT/f'{name}.{ext}',dpi=300,bbox_inches='tight')
 plt.close(fig)
def main():
 sums=read(PKG/'results/taxon_robustness_block_summary.tsv')
 # E: all designs, all-sites and transversion-only Delta D with matched 1-Mb block CIs.
 fig,ax=plt.subplots(figsize=(10.4,7.6));ys=np.arange(len(ORDER));colors={'all':'#D55E00','transversions_only':'#0072B2'}
 for j,d in enumerate(ORDER):
  for filt,dy,marker in [('all',-.12,'o'),('transversions_only',.12,'D')]:
   r=next(x for x in sums if x['design']==d and x['filter']==filt and x['scope']=='all_PNAS' and x['block_size_bp']=='1000000')
   v=float(r['delta_D']);lo=float(r['delta_CI_low']);hi=float(r['delta_CI_high'])
   ax.errorbar(v,j+dy,xerr=[[v-lo],[hi-v]],fmt=marker,color=colors[filt],ms=4.5,capsize=2,lw=1.1)
 ax.set_yticks(ys);ax.set_yticklabels([LABELS[x] for x in ORDER],fontsize=8);ax.invert_yaxis();ax.axvline(0,color='#555',lw=.8);ax.grid(axis='x',alpha=.2)
 ax.set_xlabel('D inside − D outside (estimate and 95% 1-Mb genomic block-jackknife interval)')
 ax.set_title('N61 regional contrast across the full prespecified quartet panel',loc='left',fontweight='bold')
 ax.scatter([],[],color=colors['all'],marker='o',label='All informative sites');ax.scatter([],[],color=colors['transversions_only'],marker='D',label='Transversions only');ax.legend(frameon=False,loc='lower right',fontsize=8)
 fig.text(.01,.01,'The original Gallus quartet is retained at top; Meleagris changes only the outgroup. Designs share taxa and loci and are correlated sensitivity analyses.',fontsize=7.5)
 fig.tight_layout(rect=(0,.035,1,1));save(fig,'figure_E_alternative_quartet_delta_D')
 # F: spatial bin scatter for representative replacements. Mixed-membership boundary bins are shown separately.
 bins=read(PKG/'results/nucleotide_topology_matched_bins.tsv')
 rs=read(PKG/'results/nucleotide_topology_robustness.tsv')
 reps=['N61_primary','N61_P1_Cuculus','N61_P2_Pterocles','N61_P3_Balaeniceps','N61_P3_Aquila']
 fig,axs=plt.subplots(2,3,figsize=(11,6.8),sharex=True,sharey=True);axs=axs.ravel()
 rcolor={'inside':'#D55E00','outside':'#0072B2','mixed':'#777777'}
 for ax,d in zip(axs,reps):
  dat=[r for r in bins if r['design']==d and r['filter']=='all' and r['resolution_bp']=='1000000']
  for region in ('inside','outside','mixed'):
   sub=[r for r in dat if r['region']==region]
   if sub:ax.scatter([float(r['D_site']) for r in sub],[float(r['D_topology']) for r in sub],s=16,alpha=.65,color=rcolor[region],label=region)
  rr=next((r for r in rs if r['design']==d and r['filter']=='all' and r['resolution_bp']=='1000000' and r['scope']=='whole'),None)
  if rr:ax.set_title(f"{LABELS[d]}\nr={float(rr['pearson_Dsite_Dtopology']):.2f}; centered={float(rr['within_region_centered_pearson']):.2f}",fontsize=8)
  ax.axhline(0,color='#888',lw=.6);ax.axvline(0,color='#888',lw=.6);ax.grid(alpha=.15)
 for ax in axs[len(reps):]:ax.axis('off')
 axs[3].set_xlabel('Binned nucleotide D_site');axs[4].set_xlabel('Binned nucleotide D_site');axs[0].set_ylabel('Exact local-tree D_topology = q3 − q2');axs[3].set_ylabel('Exact local-tree D_topology = q3 − q2')
 h,l=axs[0].get_legend_handles_labels();axs[5].axis('off');axs[5].legend(h,l,frameon=False,loc='center',fontsize=9,title='PNAS status at bin start')
 fig.suptitle('Matched nucleotide and exact genealogy contrasts, 1-Mb bins',fontweight='bold');fig.text(.01,.01,'Exact four-taxon trees are matched by locus ID. Bins are classified by the frozen interval predicate applied to each bin start; correlations use 2-Mb physical-block bootstrap intervals.',fontsize=7.5)
 fig.tight_layout(rect=(0,.05,1,.94));save(fig,'figure_F_site_vs_topology')
 # G: compact sensitivity matrix (all 12 designs, including outgroup sensitivity).
 corrrows={(r['design'],r['filter'],r['resolution_bp'],r['scope']):r for r in rs}
 def summary(d,filt):return next(x for x in sums if x['design']==d and x['filter']==filt and x['scope']=='all_PNAS' and x['block_size_bp']=='1000000')
 fig,ax=plt.subplots(figsize=(12,7.5));ax.set_xlim(0,6.9);ax.set_ylim(0,len(ORDER)+2);ax.invert_yaxis();ax.axis('off')
 headers=['Quartet design','ΔD all sites','ΔD transversions','r D_site vs D_topology','r after region-centering']
 xpos=[0,3.15,4.15,5.15,6.05]
 for x,h in zip(xpos,headers):ax.text(x,.55,h,ha='left',va='center',fontsize=8,fontweight='bold')
 for i,d in enumerate(ORDER,1):
  ax.text(xpos[0],i+.5,LABELS[d],ha='left',va='center',fontsize=7.5,fontweight='bold' if d=='N61_primary' else 'normal')
  for j,filt,x in [(1,'all',xpos[1]),(2,'transversions_only',xpos[2])]:
   r=summary(d,filt);v=float(r['delta_D']);lo=-.45;hi=.05;frac=max(0,min(1,(v-lo)/(hi-lo)));cmap=plt.get_cmap('RdBu_r');col=cmap(frac)
   ax.add_patch(Rectangle((x,i+.15),.82,.65,facecolor=col,edgecolor='white'))
   ax.text(x+.41,i+.48,f'{v:+.3f}',ha='center',va='center',fontsize=7,color='black')
  r=corrrows.get((d,'all','1000000','whole'));rc=corrrows.get((d,'all','1000000','whole'))
  val=float(r['pearson_Dsite_Dtopology']) if r and r['pearson_Dsite_Dtopology']!='' else float('nan')
  cent=float(r['within_region_centered_pearson']) if r and r['within_region_centered_pearson']!='' else float('nan')
  ax.text(xpos[3]+.1,i+.5,f'{val:+.2f}',ha='left',va='center',fontsize=8)
  ax.text(xpos[4]+.1,i+.5,f'{cent:+.2f}',ha='left',va='center',fontsize=8)
 ax.text(.0,len(ORDER)+1.35,'Both ΔD columns share one diverging scale (−0.45 to +0.05); spatial-correlation columns are descriptive 1-Mb estimates.',fontsize=7.5)
 ax.set_title('N61 representative and outgroup sensitivity matrix',loc='left',fontweight='bold')
 fig.tight_layout();save(fig,'figure_G_sensitivity_matrix')
 print('figures written to',OUT)
if __name__=='__main__':main()
