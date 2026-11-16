#!/usr/bin/env python3
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[5]
RES=ROOT/'empirical/atlantic_cod/analysis/patterson_d/validation/results'
OUT=ROOT/'empirical/atlantic_cod/analysis/patterson_d/validation/figures'
OUT.mkdir(parents=True,exist_ok=True)
def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=300,bbox_inches='tight')
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    fig.savefig(OUT/(name+'.svg'),bbox_inches='tight')
    plt.close(fig)
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
# Primary source robustness: SNAPP complete-call and Dsuite missing-call-allowed VCFs.
rows=[r for r in read(RES/'inside_outside_validation.tsv') if r['design_id'].endswith('_PRIMARY') and r.get('inside_D')]
fig,axs=plt.subplots(1,2,figsize=(9,3.6),sharey=True)
for ax,region in zip(axs,['inside','outside']):
    for j,lg in enumerate(['LG01','LG02','LG07','LG12']):
        for k,src in enumerate(['SNAPP','DSUITE']):
            r=next(x for x in rows if x['lg']==lg and x['source']==src)
            ax.scatter(j+(k-.5)*.18,float(r[region+'_D']),s=45,marker='o' if src=='SNAPP' else 's',label=src if j==0 else None)
    ax.axhline(0,color='0.5',lw=.8);ax.set_xticks(range(4),['LG01','LG02','LG07','LG12']);ax.set_title(region.title());ax.set_ylabel("Patterson's D")
axs[0].legend(frameon=False,ncol=2);fig.suptitle('Primary quartet: source/filter sensitivity\nSNAPP complete-call vs Dsuite missing-call-allowed input')
fig.tight_layout();save(fig,'validation_A_source_filter_sensitivity')
# Outgroup sensitivity: primary SNAPP estimates.
rows=read(RES/'outgroup_sensitivity.tsv')
fig,axs=plt.subplots(1,2,figsize=(9,3.6),sharey=True)
for ax,region in zip(axs,['inside','outside']):
    for j,lg in enumerate(['LG01','LG02','LG07','LG12']):
        for k,og in enumerate(['Borsai','Gadoga','Gadcha','Gadmac']):
            r=next(x for x in rows if x['lg']==lg and x['region']==region and x['outgroup']==og)
            ax.scatter(j+(k-1.5)*.14,float(r['D']),s=35,label=og if j==0 else None)
    ax.axhline(0,color='0.5',lw=.8);ax.set_xticks(range(4),['LG01','LG02','LG07','LG12']);ax.set_title(region.title());ax.set_ylabel("Patterson's D")
axs[0].legend(frameon=False,ncol=2);fig.suptitle('Outgroup sensitivity (primary SNAPP quartet)')
fig.tight_layout();save(fig,'validation_B_outgroup_sensitivity')
# Quartet robustness as descriptive panel.
rows=[r for r in read(RES/'additional_quartet_designs.tsv') if r['design_class']!='primary_arrangement_concordant']
fig,axs=plt.subplots(1,4,figsize=(11,3.6),sharey=True)
for ax,lg in zip(axs,['LG01','LG02','LG07','LG12']):
    rr=[r for r in rows if r['lg']==lg]
    labs=[r['design_class'].replace('same_arrangement_q3','q3 arrangement').replace('alternative_arrangement_q2','q2 arrangement').replace('baseline_pair_control','baseline pair') for r in rr]
    vals=[float(r['inside_D']) for r in rr]
    ax.bar(range(len(rr)),vals,color=['#4C78A8','#F58518','#72B7B2']);ax.axhline(0,color='0.5',lw=.8)
    ax.set_xticks(range(len(rr)),labs,rotation=45,ha='right');ax.set_title(lg);ax.set_ylim(-1.05,1.05)
axs[0].set_ylabel("Inside-inversion D")
fig.suptitle('Frozen alternative population contrasts\nDescriptive; each bar is one predeclared quartet')
fig.tight_layout();save(fig,'validation_C_frozen_quartet_robustness')
