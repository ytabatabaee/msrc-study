#!/usr/bin/env python3
"""Render six reproducible figures from the frozen D analysis tables."""
import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[5]
BASE=ROOT/'empirical/neurospora/analysis/patterson_d'
FIG=BASE/'figures'; FIG.mkdir(exist_ok=True)
def rows(name):
    with (BASE/'results'/name).open() as f: return list(csv.DictReader(f,delimiter='\t'))
def save(fig, stem):
    for ext in ('png','pdf','svg'): fig.savefig(FIG/f'{stem}.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})

# A: all linkage groups at original tree-window scale
r=[x for x in rows('window_patterson_d.tsv') if int(x['informative_sites'])>0]
fig,axs=plt.subplots(7,1,figsize=(10,13),sharey=True)
for ax,lg in zip(axs,['LGI','LGII','LGIII','LGIV','LGV','LGVI','LGVII']):
 z=[x for x in r if x['linkage_group']==lg]; x=np.array([float(t['midpoint'])/1e6 for t in z]); y=np.array([float(t['D']) for t in z]); ax.scatter(x,y,s=8,alpha=.65,color='#315b7d'); ax.axhline(0,color='black',lw=.6); ax.set_ylabel(lg,rotation=0,labelpad=26); ax.set_ylim(-1,1)
 if lg=='LGI': ax.axvspan(1.5,8,color='#e6a23c',alpha=.18,label='frozen introgressed interval'); ax.legend(loc='upper right',frameon=False)
axs[-1].set_xlabel('Linkage-group coordinate (Mb)'); fig.supylabel("Patterson's D (published 50-SNP windows)"); fig.suptitle('Windowed nucleotide site-pattern asymmetry'); fig.tight_layout(); save(fig,'figure_A_genome_tracks')
# B: paired LGI inside/outside contrast + block CIs
s={x['region']:x for x in rows('inside_outside_summary.tsv') if x['design']=='primary'}
u=rows('block_uncertainty.tsv'); fig,ax=plt.subplots(figsize=(6.5,4.2)); labels=['Inside','Outside']; vals=[float(s['LGI_inside']['D']),float(s['LGI_outside']['D'])]; lo=[]; hi=[]
for region in ['inside','outside']:
 q=next(x for x in u if x['block_width_bp']=='250000' and x['region']==region and x['block_start']=='SUMMARY'); lo.append(vals[len(lo)]-float(q['ci95_low'])); hi.append(float(q['ci95_high'])-vals[len(hi)])
ax.bar(labels,vals,color=['#d98e2b','#557a95']); ax.errorbar(labels,vals,yerr=[lo,hi],fmt='none',ecolor='black',capsize=4); ax.set_ylabel("Pooled Patterson's D"); ax.set_title('LGI inversion and outside regions (250-kb block jackknife)'); ax.set_ylim(0,1); fig.tight_layout(); save(fig,'figure_B_inside_outside')
# C: matched D and published q weights
r=rows('matched_nucleotide_genealogy.tsv'); x=np.array([(float(a['start'])+float(a['end']))/2e6 for a in r]); fig,ax=plt.subplots(figsize=(10,4.5)); ax.axvspan(1.5,8,color='#e6a23c',alpha=.18); ax.plot(x,[float(a['D']) for a in r],label='Nucleotide D',color='#263f5b',lw=1.5)
for k,c in [('q1','#888888'),('q2','#b04a4a'),('q3','#3a8c67')]: ax.plot(x,[float(a[k]) for a in r],label=f'TWISST {k}',alpha=.85,lw=1)
ax.set(xlabel='LGI coordinate (Mb)',ylabel='D / topology weight',ylim=(-1,1),title='Nucleotide D and related local-tree topology weights'); ax.legend(ncol=4,frameon=False); fig.tight_layout(); save(fig,'figure_C_D_and_trees')
# D: D against branch/divergence proxies
fig,axs=plt.subplots(1,3,figsize=(11,3.8));
for ax,k,label in zip(axs,['A_pair1','B_pair1','internal_pair1'],['Matched pair divergence contrast','Mean pairwise branch proxy','Pair internal branch']):
 a=np.array([float(z[k]) for z in r]); d=np.array([float(z['D']) for z in r]); ax.scatter(a,d,s=12,alpha=.55); ax.set(xlabel=label,ylabel="Patterson's D",title=label)
fig.suptitle('Descriptive spatial associations; branch metrics are substitution distances'); fig.tight_layout(); save(fig,'figure_D_D_branch_signals')
# E: whole vs paired sampling and transversion only
fig,ax=plt.subplots(figsize=(7,4)); cats=['LGI inside','LGI outside']; x=np.arange(2); w=.24
for i,(des,label) in enumerate([('primary','All verified samples'),('matched_pairs','Two documented pairs'),('primary_transversions_only','Transversions only')]):
 z={q['region']:q for q in rows('inside_outside_summary.tsv') if q['design']==des}; vals=[float(z['LGI_inside']['D']),float(z['LGI_outside']['D'])]; ax.bar(x+(i-1)*w,vals,w,label=label)
ax.set_xticks(x,cats); ax.set_ylabel("Pooled Patterson's D"); ax.set_ylim(0,1); ax.legend(frameon=False); ax.set_title('LGI estimates under sample and SNP-class sensitivity'); fig.tight_layout(); save(fig,'figure_E_sensitivity')
# F: informative SNP density/callability from primary windows
fig,axs=plt.subplots(2,1,figsize=(10,6),sharex=True);
r=rows('window_patterson_d.tsv')
for ax,lg in zip(axs,['LGI','LGII']):
 z=[q for q in r if q['linkage_group']==lg]; bins={}
 for q in z: bins.setdefault(int(float(q['midpoint'])//100000),0); bins[int(float(q['midpoint'])//100000)]+=int(q['informative_sites'])
 x=np.array(sorted(bins)); y=np.array([bins[k] for k in x]); ax.bar(x/10,y,width=.8,color='#5b8e7d'); ax.set_ylabel('Informative sites / 100 kb'); ax.set_title(lg); 
 if lg=='LGI': ax.axvspan(1.5,8,color='#e6a23c',alpha=.18)
axs[-1].set_xlabel('Coordinate (Mb)'); fig.suptitle('Informative-site density at published 50-SNP tree windows'); fig.tight_layout(); save(fig,'figure_F_informative_sites')
