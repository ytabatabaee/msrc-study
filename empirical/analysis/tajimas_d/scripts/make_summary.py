from pathlib import Path
import csv
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[4]
rows=list(csv.DictReader((root/'empirical/analysis/tajimas_d/results/dataset_feasibility.tsv').open(),delimiter='\t'))
colors={'READY FOR SNP-PANEL CLASSICAL D':'#287a51','PARTIAL/EXPLORATORY':'#d58b27','RECOVERABLE, BLOCKED THIS RUN':'#b84a4a','NOT ESTIMABLE':'#888888'}
fig,ax=plt.subplots(figsize=(9,5.7))
for y,r in zip(range(len(rows))[::-1],rows):
    ax.scatter(0,y,s=80,color=colors.get(r['status'],'#777'))
    ax.text(.12,y,r['dataset'],va='center',fontsize=9)
ax.set_xlim(-.2,3.0); ax.set_ylim(-1,len(rows)); ax.set_xticks([]); ax.set_yticks([])
for s in ax.spines.values(): s.set_visible(False)
ax.set_title("Tajima's D feasibility across empirical datasets",loc='left',weight='bold')
handles=[plt.Line2D([],[],marker='o',linestyle='',color=c,label=l) for l,c in [('Calculated (SNP-panel)',' #287a51'.strip()),('Partial / exploratory','#d58b27'),('Recoverable, blocked','#b84a4a'),('Not estimable','#888888')]]
ax.legend(handles=handles,loc='lower right',frameon=False,fontsize=8)
fig.text(.99,.01,'Neurospora is the only completed population SFS calculation in this run; no callable-base claim.',ha='right',fontsize=8,color='#555')
fig.tight_layout(rect=(0,0.04,1,1))
out=root/'empirical/analysis/tajimas_d/figures'
for ext in ('png','pdf','svg'): fig.savefig(out/f'tajimas_d_feasibility.{ext}',dpi=300,bbox_inches='tight')
