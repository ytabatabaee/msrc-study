"""Generate additive manuscript figures; no track interpolation."""
import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[6]; B=ROOT/'empirical/fire_ants_chr16/analysis/patterson_d'; R=B/'results'; F=B/'final_validation/figures'; F.mkdir(parents=True,exist_ok=True)
def rd(p):
 with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def save(fig,name):
 for ext in ('png','pdf','svg'):fig.savefig(F/f'{name}.{ext}',dpi=300,bbox_inches='tight')
 plt.close(fig)
def inv(ax):
 for s,e,label in ((17972919,28149888,'In(16)1'),(16509882,17962688,'In(16)2'),(11732910,16478110,'In(16)3')):
  ax.axvspan(s/1e6,e/1e6,facecolor='none',edgecolor='#777',hatch='///',alpha=.25,lw=.5)
 ax.text(.01,.01,'Hatched: Stolle mapped outer inversion envelopes; midpoint-assigned windows. In(16)3 has zero nucleotide windows. Coordinates are author-window axis, not SNP bases.',transform=ax.transAxes,fontsize=7,va='bottom')
def fig1():
 rows=rd(R/'window_patterson_d.tsv'); fig,axs=plt.subplots(4,1,figsize=(11,8),sharex=True,gridspec_kw={'height_ratios':[.45,.5,1,1]})
 ax=axs[0];ax.set_xlim(0,30);ax.set_ylim(-.1,3.3);ax.set_yticks([.5,1.5,2.5],['In(16)3 · 0 windows','In(16)2 · 4 windows','In(16)1 · 47 windows']);ax.set_xticks([]);ax.set_title('Stolle mapped outer envelopes (window-midpoint assignment)',loc='left',fontsize=9)
 for y,(s,e,c) in enumerate(((11732910,16478110,'#6096ba'),(16509882,17962688,'#52a675'),(17972919,28149888,'#b56576'))):ax.barh(y+.5,(e-s)/1e6,left=s/1e6,height=.48,color=c,alpha=.65,edgecolor=c)
 ax.axvspan(11.680438,27.917498,color='#e9a23b',alpha=.1);ax.set_yticklabels(['In(16)3 · 0 windows','In(16)2 · 4 windows','In(16)1 · 47 windows']);
 cov=rd(R/'window_vcf_coverage.tsv');cov=[r for r in cov if r['chromosome']=='chr16'];axs[1].scatter([float(r['window_mid'])/1e6 for r in cov],[int(r['busco_genes_with_vcf_records']) for r in cov],c=[int(r['busco_genes_with_vcf_records']) for r in cov],cmap='Greens',vmin=0,vmax=4,s=17);axs[1].axvspan(11.680438,27.917498,color='#e9a23b',alpha=.1);axs[1].set_ylim(-.25,4.25);axs[1].set_yticks([0,1,2,3,4]);axs[1].set_ylabel('BUSCOs\nwith VCF')
 for ax,q in zip(axs[2:],('Q1_primary_Argentina','Q2_primary_Argentina')):
  x=[(float(r['mid'])/1e6) for r in rows if r['quartet_id']==q and r['chrom']=='chr16']; y=[float(r['D']) if r['D']!='NA' else np.nan for r in rows if r['quartet_id']==q and r['chrom']=='chr16']
  ax.scatter(x,y,s=15,alpha=.75);ax.axvspan(11.680438,27.917498,color='#e9a23b',alpha=.12);ax.axhline(0,color='black',lw=.6);ax.set_ylim(-1,1);ax.set_ylabel(q.split('_')[0]+' D')
 axs[-1].set_xlabel('Stolle chr16 author-window midpoint (Mb)');fig.suptitle('Nucleotide D by sampled BUSCO window; gaps are not interpolated');fig.text(.11,.01,'Variant sites are on VCF scaffolds; plotted coordinates identify author windows. No SNP-base lift or interpolation.',fontsize=8);fig.tight_layout(rect=[0,.03,1,.98]);save(fig,'figure_1_D_tracks')
def fig2():
 d=rd(R/'window_patterson_d.tsv'); t=rd(R/'induced_quartet_topologies.tsv'); fig,axs=plt.subplots(2,1,figsize=(11,6),sharex=True)
 for ax,q in zip(axs,('Q1_primary_Argentina','Q2_primary_Argentina')):
  dr={int(x['window_index']):x for x in d if x['quartet_id']==q and x['chrom']=='chr16'};tr={int(x['window_index']):x for x in t if x['quartet_id']==q and x['chrom']=='chr16'}
  ix=sorted(set(dr)&set(tr));x=[float(dr[i]['mid'])/1e6 for i in ix];ax.scatter(x,[float(dr[i]['D']) if dr[i]['D']!='NA' else np.nan for i in ix],s=12,label='D_site');ax.scatter(x,[float(tr[i]['D_topology']) if tr[i]['D_topology'] not in ('NA','') else np.nan for i in ix],s=12,marker='x',label='D_topology');ax.axhline(0,color='black',lw=.6);ax.axvspan(11.680438,27.917498,color='#e9a23b',alpha=.12);ax.set_ylim(-1,1);ax.set_ylabel(q.split('_')[0]);ax.legend(frameon=False,ncol=2)
 axs[-1].set_xlabel('Stolle chr16 author-window midpoint (Mb)');fig.suptitle('Matched site and exact induced-quartet contrasts');save(fig,'figure_2_D_genealogy_tracks')
def fig3():
 r=rd(B/'final_validation/results/genotype_filter_sensitivity.tsv');fig,axs=plt.subplots(2,2,figsize=(10,7),sharex=True,sharey=True)
 for row,q in enumerate(('Q1_primary_Argentina','Q2_primary_Argentina')):
  for col,reg in enumerate(('chr16_supergene','chr16A+chr16B')):
   ax=axs[row,col];z=[x for x in r if x['quartet_id']==q and x['region']==reg]
   if reg=='chr16A+chr16B':
    z=[]
    for d in (8,10,12,15):
     for f in (.5,.6,.7,.8):
      a=next(x for x in r if x['quartet_id']==q and x['region']=='chr16A' and x['min_DP']==str(d) and x['min_ingroup_call_fraction']==str(f));b=next(x for x in r if x['quartet_id']==q and x['region']=='chr16B' and x['min_DP']==str(d) and x['min_ingroup_call_fraction']==str(f));aa=float(a['ABBA'])+float(b['ABBA']);bb=float(a['BABA'])+float(b['BABA']);ii=int(a['informative_sites'])+int(b['informative_sites']);z.append({'min_DP':str(d),'min_ingroup_call_fraction':str(f),'D':(aa-bb)/(aa+bb) if aa+bb else 'NA','informative_sites':str(ii)})
   mat=np.full((4,4),np.nan)
   for x in z:
    mat[(8,10,12,15).index(int(x['min_DP'])),(.5,.6,.7,.8).index(float(x['min_ingroup_call_fraction']))]=float(x['D']) if x['D']!='NA' else np.nan
   ax.imshow(mat,vmin=-.1,vmax=.9,cmap='coolwarm',aspect='auto');ax.set_title(q.split('_')[0]+' '+reg)
   for x in z:
    yy=(8,10,12,15).index(int(x['min_DP']));xx=(.5,.6,.7,.8).index(float(x['min_ingroup_call_fraction']));ax.text(xx,yy,('NA' if x['D']=='NA' else '%.2f'%float(x['D']))+'\n'+x['informative_sites'],ha='center',va='center',fontsize=7,color='black')
   ax.set_xticks(range(4),['50%','60%','70%','80%']);ax.set_yticks(range(4),['8','10','12','15']);ax.set_xlabel('Minimum ingroup call fraction');ax.set_ylabel('Minimum DP')
 fig.suptitle('Genotype-filter sensitivity (cell: D / informative SNP count)');fig.tight_layout();save(fig,'figure_3_filter_sensitivity')
def fig4():
 r=rd(R/'window_patterson_d.tsv'); fig,axs=plt.subplots(1,2,figsize=(10,4.6))
 for q,marker in (('Q1_primary_Argentina','o'),('Q2_primary_Argentina','s')):
  z=[x for x in r if x['quartet_id']==q and x['chrom']=='chr16' and x['region']=='chr16_supergene'];axs[0].scatter([float(x['genotype_call_rate_at_variant_sites']) for x in z],[float(x['D']) if x['D']!='NA' else np.nan for x in z],s=14,marker=marker,label=q.split('_')[0],alpha=.75)
 axs[0].set(xlabel='Call rate at variant sites',ylabel='Window D',title='Callability vs window contribution');axs[0].legend(frameon=False)
 wi=rd(B/'final_validation/results/window_influence.tsv')
 for y,q,marker in ((1,'Q1_primary_Argentina','o'),(0,'Q2_primary_Argentina','s')):
  z=[x for x in wi if x['quartet_id']==q and x['region']=='chr16_supergene'];vals=[float(x['D_supergene_without_window']) for x in z];axs[1].scatter(vals,np.full(len(vals),y)+np.linspace(-.09,.09,len(vals)),s=14,marker=marker,alpha=.55);full=.82488634196149 if q.startswith('Q1') else .8442934217014734;axs[1].axvline(full,color='#777',lw=1,ls='--')
 axs[1].set_yticks([0,1],['Q2','Q1']);axs[1].set_xlim(.75,.9);axs[1].set_xlabel('Supergene D after deleting one window');axs[1].set_title('Window-level leverage');axs[1].grid(axis='x',alpha=.2)
 fig.suptitle('Callability, window influence, and data-quality limits');fig.text(.02,.01,'Geography/outgroup subsets preserve direction. Individual LOO and gene-level influence are not computed. GT:DP lacks AD/GQ/read mapping data; mapping and copy-number bias cannot be ruled out.',fontsize=8);fig.tight_layout(rect=[0,.06,1,.95]);save(fig,'figure_4_sample_gene_reference_diagnostics')
if __name__=='__main__':fig1();fig2();fig3();fig4()
