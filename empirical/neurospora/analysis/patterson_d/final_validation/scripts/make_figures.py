#!/usr/bin/env python3
"""Create four additive, publication-oriented Neurospora validation figures."""
import csv, math
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[6]
BASE=ROOT/'empirical/neurospora/analysis/patterson_d'
OUT=BASE/'final_validation'; R=OUT/'results'; FIG=OUT/'figures'; FIG.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
BLUE='#315b7d'; ORANGE='#d98e2b'; TEAL='#3a8c78'; GRAY='#70777f'; RED='#b5514a'

def read(p):
    with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def save(fig,stem):
    for ext in ['png','pdf','svg']:fig.savefig(FIG/(stem+'.'+ext),dpi=220,bbox_inches='tight')
    plt.close(fig)
def safe(value):
    try:return float(value)
    except:return float('nan')
def dval(a,b):return (a-b)/(a+b) if a+b else float('nan')
def jk(blocks):
    ab=sum(x[0] for x in blocks);ba=sum(x[1] for x in blocks);est=dval(ab,ba);n=len(blocks)
    if n<3:return est,float('nan'),float('nan'),n
    loo=[dval(ab-a,ba-b) for a,b in blocks];avg=sum(loo)/n;se=math.sqrt((n-1)/n*sum((x-avg)**2 for x in loo))
    return est,est-1.96*se,est+1.96*se,n

def fig1():
    rows=[x for x in read(R/'spatial_profile_bins.tsv') if x['linkage_group']=='LGI' and x['bin_width_bp']=='100000']
    rows.sort(key=lambda x:float(x['midpoint']))
    x=np.array([float(r['midpoint'])/1e6 for r in rows]); d=np.array([safe(r['D']) for r in rows]); info=np.array([int(r['informative_sites']) for r in rows])
    use=info>=50
    fig,(ax,ac)=plt.subplots(2,1,figsize=(10,5.4),sharex=True,gridspec_kw={'height_ratios':[3,1]})
    ax.axvspan(1.5,8.0,color=ORANGE,alpha=.16,label='Published interval (1.5–8.0 Mb)')
    for lo,hi,a in [(0.5,1.5,.08),(8,9,.08),(1,1.5,.13),(8,8.5,.13),(1.25,1.5,.18),(8,8.25,.18)]:ax.axvspan(lo,hi,color=GRAY,alpha=a,lw=0)
    sc=ax.scatter(x[use],d[use],c=info[use],cmap='viridis',s=28,edgecolors='none',zorder=3)
    ax.axhline(0,color='black',lw=.7);ax.set_ylim(-1,1);ax.set_ylabel("Pooled nucleotide D");ax.set_title('LGI D along physical coordinates (bins with ≥50 informative sites)')
    cb=fig.colorbar(sc,ax=ax,pad=.01);cb.set_label('Informative sites / 100-kb bin')
    ax.legend(frameon=False,loc='lower left')
    ac.bar(x,info,width=.085,color=TEAL,alpha=.8);ac.axhline(50,color=RED,ls='--',lw=.8,label='50-site display threshold');ac.set_ylabel('Informative\nsites');ac.set_xlabel('LGI coordinate (Mb)');ac.legend(frameon=False,loc='upper right')
    fig.text(.01,.005,'Variant-panel conditional; points are descriptive pooled bins, not independent observations or confidence intervals.',fontsize=8)
    fig.tight_layout(rect=[0,.035,1,1]);save(fig,'figure_1_LGI_D_profile')

def fig2():
    b50=[x for x in read(R/'spatial_profile_bins.tsv') if x['bin_width_bp']=='50000']
    regs=[('Focal\n1.5–8 Mb','focal'),('Left flank\n1 Mb','left1'),('Left distant\n0.5 Mb','leftd'),('Right flank\n1 Mb','right1'),('Right distant\n>9 Mb','rightd'),('LGII','LGII'),('LGIII','LGIII'),('LGIV','LGIV'),('LGV','LGV'),('LGVI','LGVI'),('LGVII','LGVII')]
    lgi=[x for x in b50 if x['linkage_group']=='LGI']
    specs={'focal':lambda m:1_500_000<=m<=8_000_000,'left1':lambda m:500_000<=m<1_500_000,'leftd':lambda m:m<500_000,'right1':lambda m:8_000_000<m<=9_000_000,'rightd':lambda m:m>9_000_000}
    estimates=[]
    for label,key in regs:
        if key in specs:
            q=[x for x in lgi if specs[key](float(x['midpoint']))]; blockdata={}
            for r in q:
                m=float(r['midpoint'])
                if key=='focal':bid=(m-1_500_000)//250_000
                elif key.startswith('left'):bid=('L',(m)//250_000)
                else:bid=('R',(m-8_000_000)//250_000)
                a=float(r['ABBA']);b=float(r['BABA']);old=blockdata.setdefault(bid,[0.,0.]);old[0]+=a;old[1]+=b
            blocks=list(blockdata.values());ab=sum(float(r['ABBA']) for r in q);ba=sum(float(r['BABA']) for r in q)
            est,lo,hi,n=jk(blocks);estimates.append((label,est,lo,hi,n,ab,ba))
        else:
            bins=read(BASE/'results/multiresolution_patterson_d.tsv')
            q=[x for x in bins if x['design']=='primary' and x['linkage_group']==key and x['bin_width_bp']=='500000'];blockdata={}
            for r in q:
                bid=int(float(r['start'])//1_000_000);z=blockdata.setdefault(bid,[0.,0.]);z[0]+=float(r['ABBA']);z[1]+=float(r['BABA'])
            blocks=list(blockdata.values());ab=sum(z[0] for z in blocks);ba=sum(z[1] for z in blocks);est,lo,hi,n=jk(blocks);estimates.append((label,est,lo,hi,n,ab,ba))
    fig,ax=plt.subplots(figsize=(10.5,4.8));y=np.arange(len(estimates)); colors=[ORANGE if x[0].startswith('Focal') else TEAL if 'flank' in x[0].lower() else BLUE if x[0].startswith('LG') else GRAY for x in estimates]
    for i,(label,est,lo,hi,n,ab,ba) in enumerate(estimates):
        ax.scatter(est,i,color=colors[i],s=38,zorder=3)
        if math.isfinite(lo) and math.isfinite(hi):ax.plot([lo,hi],[i,i],color=colors[i],lw=1.4);ax.plot([lo,lo],[i-.08,i+.08],color=colors[i]);ax.plot([hi,hi],[i-.08,i+.08],color=colors[i])
        ax.text(est+.025,i,'n=%d'%n,va='center',fontsize=7,color='#444')
    ax.axvline(0,color='black',lw=.7);ax.set_yticks(y,[x[0] for x in estimates]);ax.invert_yaxis();ax.set_xlim(-.25,1.25);ax.set_xlabel("Pooled nucleotide D (points; whiskers = 95% physical-block delete-one interval when ≥3 blocks)");ax.set_title('Focal LGI, 1-Mb flanks, distant segments, and chromosome controls')
    fig.text(.01,.005,'D pools ABBA/BABA. Flank/distant units use 250-kb blocks; other linkage groups use 1-Mb blocks. Few-block regions have no interval.',fontsize=8)
    fig.tight_layout(rect=[0,.04,1,1]);save(fig,'figure_2_regional_D_blocks')

def fig3():
    rows=read(BASE/'results/matched_nucleotide_genealogy.tsv');rows.sort(key=lambda r:float(r['start']))
    x=np.array([(float(r['start'])+float(r['end']))/2/1e6 for r in rows]);d=np.array([safe(r['D']) for r in rows]);q=np.array([safe(r['q3']) for r in rows])
    fig,ax=plt.subplots(figsize=(10,4.8));ax.axvspan(1.5,8,color=ORANGE,alpha=.15);ax.plot(x,d,color=BLUE,lw=1.5,label='Nucleotide D, N. discreta polarized');ax.plot(x,q,color=TEAL,lw=1.2,label='Published TWISST q3')
    ax.axhline(0,color='black',lw=.7);ax.set(xlabel='LGI coordinate (Mb)',ylabel='D / topology weight',ylim=(-1,1),title='Matched nucleotide D and published genealogy weight')
    ax.legend(frameon=False,ncol=2,loc='lower left');fig.text(.01,.005,'Local trees contain N. crassa and N. hispaniola but no N. discreta; q3 is not rooted D_topology. Tracks are linked and share selected SNP data.',fontsize=8)
    fig.tight_layout(rect=[0,.04,1,1]);save(fig,'figure_3_D_TWISST_q3')

def fig4():
    win=read(R/'ascertainment_diagnostics.tsv');win=[x for x in win if x['linkage_group']=='LGI']
    inside=[x for x in win if x['inside_LGI_1p5_8Mb']=='1'];outside=[x for x in win if x['inside_LGI_1p5_8Mb']=='0']
    sample=read(R/'sample_sensitivity.tsv')
    fig,axs=plt.subplots(2,2,figsize=(10,7))
    axs[0,0].boxplot([[int(x['physical_span_bp'])/1000 for x in inside],[int(x['physical_span_bp'])/1000 for x in outside]],labels=['Inside','Outside'],showfliers=False,patch_artist=True,boxprops={'facecolor':'#a8bfce'});axs[0,0].set_ylabel('50-SNP window span (kb)');axs[0,0].set_title('Variable physical span')
    axs[0,1].boxplot([[float(x['outgroup_callability']) for x in inside],[float(x['outgroup_callability']) for x in outside]],labels=['Inside','Outside'],showfliers=False,patch_artist=True,boxprops={'facecolor':'#b5d2c7'});axs[0,1].set_ylim(0,1.02);axs[0,1].set_ylabel('N. discreta call fraction');axs[0,1].set_title('Outgroup callability at SNP rows')
    for i,n in enumerate([2,4,6,8]):
        a=[float(x['D']) for x in sample if x['analysis']=='balanced_random_subset' and x['linkage_group']=='LGI' and x['n_per_mating_type']==str(n) and x['region']=='inside_LGI'];b=[float(x['D']) for x in sample if x['analysis']=='balanced_random_subset' and x['linkage_group']=='LGI' and x['n_per_mating_type']==str(n) and x['region']=='outside_LGI']
        axs[1,0].scatter(np.full(len(a),i-.12),a,s=13,color=ORANGE,alpha=.7);axs[1,0].scatter(np.full(len(b),i+.12),b,s=13,color=BLUE,alpha=.7)
    axs[1,0].set_xticks(range(4),['2/2','4/4','6/6','8/8']);axs[1,0].set_ylim(-.1,1);axs[1,0].set_ylabel('Pooled D');axs[1,0].set_xlabel('Individuals per mating type');axs[1,0].set_title('LGI balanced subsets');axs[1,0].plot([],[],'.',color=ORANGE,label='Inside');axs[1,0].plot([],[],'.',color=BLUE,label='Outside');axs[1,0].legend(frameon=False)
    for i,n in enumerate([2,4,6,8]):
        z=[float(x['D']) for x in sample if x['analysis']=='balanced_random_subset' and x['linkage_group']=='LGIV' and x['n_per_mating_type']==str(n) and x['region']=='whole'];axs[1,1].scatter(np.full(len(z),i),z,s=17,color=RED,alpha=.72)
    full=next(x for x in read(BASE/'results/inside_outside_summary.tsv') if x['design']=='primary' and x['region']=='LGIV');pair=next(x for x in read(BASE/'results/inside_outside_summary.tsv') if x['design']=='matched_pairs' and x['region']=='LGIV')
    axs[1,1].axhline(float(full['D']),color=GRAY,ls='--',label='Full 10/8 groups');axs[1,1].axhline(float(pair['D']),color=RED,ls=':',label='Two matched pairs');axs[1,1].axhline(0,color='black',lw=.6);axs[1,1].set_xticks(range(4),['2/2','4/4','6/6','8/8']);axs[1,1].set_ylabel('LGIV pooled D');axs[1,1].set_xlabel('Individuals per mating type');axs[1,1].set_title('LGIV sample sensitivity');axs[1,1].legend(frameon=False,fontsize=8)
    fig.suptitle('Ascertainment and isolate-level sensitivity');fig.tight_layout(rect=[0,0,1,.96]);save(fig,'figure_4_ascertainment_sample_sensitivity')

if __name__=='__main__':
    fig1();fig2();fig3();fig4()
