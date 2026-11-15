#!/usr/bin/env python3
"""Create manuscript-style spatial Patterson D figures (PNG/PDF/SVG)."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent; REPO=HERE.parents[4]
AN=REPO/'empirical/atlantic_cod/analysis/patterson_d'; RES=AN/'results'; FIG=AN/'figures'
REGIONS={}
with (REPO/'data/atlantic_cod/metadata/region_manifest.tsv').open() as f:
    for r in csv.DictReader(f,delimiter='\t'): REGIONS[r['lg']]=(int(r['start']),int(r['end']))
def rows(name):
    with (RES/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def num(v):
    try:return float(v) if v not in ('','NA','None') else np.nan
    except:return np.nan
def save(fig,name):
    FIG.mkdir(parents=True,exist_ok=True)
    for ext in ('png','pdf','svg'):fig.savefig(FIG/(name+'.'+ext),dpi=300,bbox_inches='tight')
    plt.close(fig)

def main():
    win=rows('window_patterson_d.tsv');summ=rows('inside_outside_summary.tsv');top=rows('topology_D_comparison.tsv')
    lgs=('LG01','LG02','LG07','LG12');colors={'LG01':'#1769aa','LG02':'#d17b0f','LG07':'#31936a','LG12':'#9a4c9b'}
    fig,axs=plt.subplots(4,1,figsize=(12,10),sharex=False)
    for ax,lg in zip(axs,lgs):
        rr=sorted((r for r in win if r['lg']==lg),key=lambda r:int(r['midpoint']))
        x=np.array([int(r['midpoint'])/1e6 for r in rr]);y=np.array([num(r['D_site']) for r in rr])
        a,b=REGIONS[lg];ax.axvspan(a/1e6,b/1e6,color='#e8a6a6',alpha=.35,label='Frozen inversion')
        ax.scatter(x,y,s=17,c=colors[lg],alpha=.8);ax.axhline(0,color='black',lw=.7)
        ax.set_ylabel(lg+'\nD');ax.set_ylim(-1.05,1.05);ax.grid(axis='y',alpha=.2)
    axs[-1].set_xlabel('Position on linkage group (Mb; 1-based reference coordinates)')
    fig.suptitle('Population-frequency Patterson D in published 250-kb tree windows\nQuartets frozen from outside-inversion baseline; Borsai outgroup')
    save(fig,'figure_A_D_spatial_tracks')

    fig,ax=plt.subplots(figsize=(9,5));xx=np.arange(4);width=.34
    ins=[num(next(r['inside_D'] for r in summ if r['lg']==lg)) for lg in lgs]
    out=[num(next(r['outside_D'] for r in summ if r['lg']==lg)) for lg in lgs]
    ax.bar(xx-width/2,ins,width,label='Inside inversion');ax.bar(xx+width/2,out,width,label='Outside inversion')
    ax.axhline(0,color='black',lw=.8);ax.set_xticks(xx,lgs);ax.set_ylabel('Pooled Patterson D');ax.legend(frameon=False)
    ax.set_title('Pooled inside/outside D (boundary-overlapping windows excluded)');save(fig,'figure_B_inside_outside')

    fig,axs=plt.subplots(4,1,figsize=(12,10),sharex=False)
    for ax,lg in zip(axs,lgs):
        rr=sorted((r for r in top if r['lg']==lg),key=lambda r:int(r['start']))
        x=np.array([int(r['start'])/1e6 for r in rr]);y=np.array([num(r['D_site']) for r in rr]);q=np.array([num(r['rooted_triple_contrast_q3_minus_q2']) for r in rr])
        a,b=REGIONS[lg];ax.axvspan(a/1e6,b/1e6,color='#e8a6a6',alpha=.35)
        ax.plot(x,y,color=colors[lg],lw=1,label='Nucleotide D')
        ax.step(x,q,where='mid',color='#444444',lw=.8,alpha=.65,label='Rooted 3-population SNAPP contrast')
        ax.axhline(0,color='black',lw=.6);ax.set_ylabel(lg);ax.set_ylim(-1.1,1.1);ax.legend(loc='upper right',fontsize=7,frameon=False)
    axs[-1].set_xlabel('Position (Mb)');fig.suptitle('Matched nucleotide D and rooted three-tip genealogy contrast\n(Not the exact four-taxon D_topology: published trees omit Borsai)')
    save(fig,'figure_C_D_and_topology')

    fig,axs=plt.subplots(2,2,figsize=(10,8))
    for ax,lg in zip(axs.flat,lgs):
        rr=[r for r in top if r['lg']==lg and num(r['D_site'])==num(r['D_site']) and num(r['A_opposite_minus_same'])==num(r['A_opposite_minus_same'])]
        ax.scatter([num(r['A_opposite_minus_same']) for r in rr],[num(r['D_site']) for r in rr],s=13,alpha=.65,c=colors[lg])
        ax.axhline(0,color='grey',lw=.5);ax.axvline(0,color='grey',lw=.5);ax.set_title(lg);ax.set_xlabel('Published SNAPP divergence-time contrast');ax.set_ylabel('Nucleotide D')
    fig.suptitle('Nucleotide asymmetry versus available divergence-time signal');save(fig,'figure_D_D_vs_time')

    fig,axs=plt.subplots(4,1,figsize=(12,10),sharex=False)
    for ax,lg in zip(axs,lgs):
        rr=sorted((z for z in win if z['lg']==lg),key=lambda z:int(z['midpoint']))
        x=np.array([int(z['midpoint'])/1e6 for z in rr]);density=np.array([int(z['n_accessible_biallelic_snps']) for z in rr])
        records=np.array([int(z['n_vcf_records']) for z in rr]);fail=np.array([int(z['n_outgroup_polarization_failures']) for z in rr])
        failure_rate=np.divide(fail,records,out=np.zeros_like(fail,dtype=float),where=records>0)
        a,b=REGIONS[lg];ax.axvspan(a/1e6,b/1e6,color='#e8a6a6',alpha=.35)
        ax.bar(x,density,width=.22,color=colors[lg],alpha=.65,label='Biallelic SNPs per 250-kb window')
        axr=ax.twinx();axr.plot(x,failure_rate,color='#333333',lw=.8,alpha=.8,label='Outgroup polarization failure rate')
        ax.set_ylabel(lg+'\nSNPs/window');axr.set_ylim(0,1);axr.set_ylabel('O missing/heterozygous')
        ax.grid(axis='y',alpha=.15)
    axs[-1].set_xlabel('Position (Mb; 1-based reference coordinates)')
    fig.suptitle('Callable variant density and outgroup polarization failures per matched tree window')
    save(fig,'figure_F_coverage')

    try:
        rob=rows('robustness_analysis.tsv');fig,ax=plt.subplots(figsize=(10,5))
        scenarios=list(dict.fromkeys(r['scenario'] for r in rob));base={r['lg']:num(r['inside_minus_outside_D']) for r in rob if r['scenario']=='primary_fixed_samples_Borsai'}
        marks=['o','s','^','D','v','P','X','*'];
        for j,sc in enumerate(scenarios):
            vals=[next((num(r['inside_minus_outside_D']) for r in rob if r['scenario']==sc and r['lg']==lg),np.nan) for lg in lgs]
            ax.scatter(np.arange(4)+(j-(len(scenarios)-1)/2)*.045,vals,label=sc,marker=marks[j%len(marks)],s=34)
        ax.axhline(0,color='black',lw=.7);ax.set_xticks(np.arange(4),lgs);ax.set_ylabel('Inside minus outside D');ax.set_title('Frozen-quartet outgroup, sample-size, and SNP-filter sensitivity');ax.legend(fontsize=7,ncol=2,frameon=False);save(fig,'figure_E_sensitivity')
    except FileNotFoundError: pass

if __name__=='__main__':main()
