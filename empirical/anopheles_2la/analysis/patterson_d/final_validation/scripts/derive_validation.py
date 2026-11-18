from pathlib import Path
import csv, math, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
B=Path('empirical/anopheles_2la'); V=B/'analysis/patterson_d/final_validation'; R=V/'results'; F=V/'figures'
R.mkdir(parents=True,exist_ok=True); F.mkdir(parents=True,exist_ok=True)
def read(p):
 with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rows):
 if not rows:return
 with open(R/name,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def f(x):
 try:return float(x)
 except:return float('nan')
rows=read(B/'results/patterson_d/window_results.tsv')
primary=[r for r in rows if r['filter_label']=='primary']
# robust pooled aggregation from physical-window rows
out=[]
def agg(test, subset, label, size=500000):
 a=sum(f(x['abba']) for x in subset); b=sum(f(x['baba']) for x in subset); den=a+b
 blocks={}
 for x in subset:
  k=(int(x['start'])-1)//size
  q=blocks.setdefault(k,[0.,0.,0]);q[0]+=f(x['abba']);q[1]+=f(x['baba']);q[2]+=int(float(x['n_sites_informative']))
 vals=[z for z in blocks.values() if z[0]+z[1]>0 and z[2]>0]
 D=(a-b)/den if den else float('nan'); est=[]
 for aa,bb,n in vals:
  dd=den-aa-bb
  if dd>0:est.append(((a-aa)-(b-bb))/dd)
 if len(est)>1:
  m=np.mean(est);se=math.sqrt((len(est)-1)/len(est)*sum((q-m)**2 for q in est));lo=D-1.96*se;hi=D+1.96*se
 else:se=lo=hi=float('nan')
 return dict(test_id=test,scenario=label,n_windows=len(subset),ABBA=a,BABA=b,D_pooled=D,block_bp=size,n_blocks=len(vals),jackknife_se=se,ci95_low=lo,ci95_high=hi,informative_contributions=sum(int(float(x['n_sites_informative'])) for x in subset))
# region, full tracks, exclude top influence windows, low callability quantile
for t in ['A_merus','B_merus','A_melas','B_melas','C_merus','D_merus']:
 rr=[x for x in primary if x['test_id']==t]
 for name,ss in [('inside', [x for x in rr if x['region_class']=='inside']),('outside',[x for x in rr if x['region_class']=='outside']),('left_flank',[x for x in rr if x['region_class']=='outside' and int(x['end'])<=20524057]),('right_flank',[x for x in rr if x['region_class']=='outside' and int(x['start'])>=42165533])]:out.append(agg(t,ss,name))
 # contribution influence ranked by denom; remove top 1%, 5% of inside windows
 ins=[x for x in rr if x['region_class']=='inside']
 for pct in (1,5):
  n=max(1,math.ceil(len(ins)*pct/100)); drop={x['window_id'] for x in sorted(ins,key=lambda x:f(x['abba'])+f(x['baba']),reverse=True)[:n]}
  out.append(agg(t,[x for x in ins if x['window_id'] not in drop],f'inside_drop_top_{pct}pct_denominator_windows'))
 # callability extreme low: by min average group callable chromosome metric, lower decile
 vals=np.array([min(f(x[k]) for k in ['mean_callable_P1','mean_callable_P2','mean_callable_P3','mean_callable_O']) for x in ins]);cut=np.quantile(vals,.1)
 out.append(agg(t,[x for x,v in zip(ins,vals) if v>cut], 'inside_exclude_lowest_callable_decile'))
# modest coordinate perturbations of the frozen structural interval; boundary-crossing windows excluded
for t in ['A_merus','B_merus','A_melas','B_melas']:
 rr=[x for x in primary if x['test_id']==t]
 for shift in (-250000,-100000,100000,250000):
  lo=20524058+shift; hi=42165532-shift
  if lo>=hi: continue
  inside=[x for x in rr if int(x['start'])>=lo and int(x['end'])<=hi]
  outside=[x for x in rr if int(x['end'])<lo or int(x['start'])>hi]
  out.append(agg(t,inside,'boundary_perturb_%+d_inside'%shift))
  out.append(agg(t,outside,'boundary_perturb_%+d_outside'%shift))
write('spatial_robustness.tsv',out)
# inversion quartiles by physical coordinate equal length; exact pooled effects
inv=[]
for t in ['A_merus','B_merus','A_melas','B_melas','C_merus','D_merus']:
 rr=[x for x in primary if x['test_id']==t and x['region_class']=='inside']
 for q in range(4):
  lo=20524058+(42165532-20524058+1)*q/4;hi=20524058+(42165532-20524058+1)*(q+1)/4
  ss=[x for x in rr if int(x['start'])>=lo and int(x['start'])<=hi]
  z=agg(t,ss,f'inversion_coordinate_quartile_{q+1}')
  inv.append({**z,'start_bp':int(lo),'end_bp':int(hi),'fraction_windows_same_sign':sum((f(x['D_patterson'])*z['D_pooled']>0) for x in ss)/len(ss) if ss and z['D_pooled']==z['D_pooled'] else float('nan')})
write('inversion_region_summary.tsv',inv)
# estimator/filter validation from existing results
v=[]
for t in sorted(set(x['test_id'] for x in rows)):
 for region in sorted(set(x['region_class'] for x in rows)):
  a=[x for x in rows if x['test_id']==t and x['region_class']==region and x['filter_label']=='primary']; b=[x for x in rows if x['test_id']==t and x['region_class']==region and x['filter_label']=='MAC2']
  if not a:continue
  sums=lambda s,k:sum(f(x[k]) for x in s)
  v.append(dict(test_id=t,region=region,n_windows=len(a),primary_ABBA=sums(a,'abba'),MAC2_ABBA=sums(b,'abba'),primary_BABA=sums(a,'baba'),MAC2_BABA=sums(b,'baba'),primary_pass_sites=sums(a,'n_sites_pass'),MAC2_pass_sites=sums(b,'n_sites_pass'),primary_informative_contributions=sums(a,'n_sites_informative'),MAC2_informative_contributions=sums(b,'n_sites_informative'),delta_ABBA=sums(b,'abba')-sums(a,'abba'),delta_BABA=sums(b,'baba')-sums(a,'baba'),status='derived_from_retained_full_scan'))
write('estimator_validation.tsv',v)
# coverage diagnostics summarize windows and exclusion counts
cd=[]
for t in ['A_merus','B_merus','C_merus','D_merus']:
 rr=[x for x in primary if x['test_id']==t]
 for reg in ['inside','boundary','left_flank','right_flank','outside']:
  ss=[x for x in rr if (x['region_class']==reg if reg in ['inside','boundary','outside'] else (x['region_class']=='outside' and (int(x['end'])<=20524057 if reg=='left_flank' else int(x['start'])>=42165533)))]
  if not ss:continue
  def mean(k):return float(np.mean([f(x[k]) for x in ss]))
  ex={}
  for x in ss:
   z=json.loads(x['exclusion_counts'])
   for k,n in z.items():ex[k]=ex.get(k,0)+n
  cd.append(dict(test_id=t,region=reg,n_windows=len(ss),mean_examined_variant_records=mean('n_sites_examined'),mean_pass_sites=mean('n_sites_pass'),mean_informative_contributions=mean('n_sites_informative'),mean_callable_P1=mean('mean_callable_P1'),mean_callable_P2=mean('mean_callable_P2'),mean_callable_P3=mean('mean_callable_P3'),mean_callable_outgroup=mean('mean_callable_O'),excluded_missing=sum(ex.values()) if False else ex.get('site_missing',0),excluded_outgroup_missing=ex.get('outgroup_missing',0),excluded_outgroup_polymorphic=ex.get('outgroup_polymorphic',0),excluded_multiallelic=ex.get('multiallelic',0)))
write('callability_diagnostics.tsv',cd)
# topology correspondences: window tracks include exact common windows; q_species/q_arrangement/q_third from existing topology; no exact mapping to population quartets.
top=read(B/'results/stage2r_quartet_topology_signal.tsv'); topmap={x['window_id']:x for x in top}
gene=[]
for t in ['A_merus','B_merus']:
 rr=[x for x in primary if x['test_id']==t]
 for x in rr:
  q=topmap.get(x['window_id'],{})
  if not q:continue
  gene.append(dict(window_id=x['window_id'],start=x['start'],end=x['end'],region=x['region_class'],test_id=t,D_Patterson=x['D_patterson'],ABBA=x['abba'],BABA=x['baba'],informative_contributions=x['n_sites_informative'],q_species=q.get('q_species',''),q_arrangement=q.get('q_arrangement',''),q_third=q.get('q_third',''),Delta_q_arrangement_minus_species=f(q.get('q_arrangement',''))-f(q.get('q_species','')) if q.get('q_arrangement') else '',topology_D=q.get('D',''),quartet_match='not exact: Stage2R topology uses frozen species/arrangement haplotypes, not the same population sampling as nucleotide quartet'))
write('genealogy_correspondence.tsv',gene)
# final stats and manuscript summary derived from existing region table
summ=read(B/'results/patterson_d/inside_outside_summary.tsv'); jk=read(B/'results/patterson_d/jackknife_results.tsv')
final=[]
for x in summ:
 if x['region'] not in ('inside','outside','left_flank','right_flank'):continue
 test=x['test_id'];region=x['region'];jj=[z for z in jk if z['test_id']==test and z['region']==region and int(z['block_size_bp'])==500000]
 j=jj[0] if jj else {}
 final.append(dict(test_id=test,region=region,n_windows=x['n_windows'],ABBA=x['pooled_ABBA'],BABA=x['pooled_BABA'],D_pooled=x['pooled_D'],informative_contributions=x['informative_site_sum'],block_bp=500000,n_blocks=j.get('n_blocks',''),ci95_low=j.get('ci_low',''),ci95_high=j.get('ci_high',''),status='verified_from_completed_full_scan'))
write('final_results.tsv',final)
# Create 4 coordinated figure sets
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
def savefig(fig,name):
 for ext in ['png','pdf','svg']:fig.savefig(F/f'{name}.{ext}',dpi=300,bbox_inches='tight')
 plt.close(fig)
# tracks: observed 50k D points; separate genes
fig,axs=plt.subplots(2,1,figsize=(10,5.5),sharex=True)
for ax,t in zip(axs,['A_merus','B_merus']):
 rr=sorted([x for x in primary if x['test_id']==t],key=lambda x:int(x['start']));xx=np.array([(int(x['start'])+int(x['end']))/2/1e6 for x in rr]); yy=np.array([f(x['D_patterson']) for x in rr]); ax.scatter(xx,yy,s=5,color='#31688e',alpha=.75);ax.axvspan(20.524058,42.165532,color='#f1c75b',alpha=.2);ax.axhline(0,color='black',lw=.6);ax.set_ylabel(f'{t}\n$D_{{Patterson}}$');ax.set_ylim(-1.05,1.05)
axs[-1].set_xlabel('2L coordinate (Mb; Ag3.10 reference)');fig.suptitle('Observed window Patterson’s D (50-kb source grid; points are unsmoothed)');savefig(fig,'figure_1_D_tracks')
# topology matched track inside interval
fig,axs=plt.subplots(2,1,figsize=(10,5.5),sharex=True)
for ax,t in zip(axs,['A_merus','B_merus']):
 rr=sorted([x for x in gene if x['test_id']==t and x['region']=='inside'],key=lambda x:int(x['start'])); xx=np.array([(int(x['start'])+int(x['end']))/2/1e6 for x in rr]);ax.plot(xx,[f(x['D_Patterson']) for x in rr],'.',ms=2,color='#31688e',label='$D_{Patterson}$');ax.plot(xx,[f(x['Delta_q_arrangement_minus_species']) for x in rr],'.',ms=2,color='#d1495b',label='$\\Delta_q=q_{arrangement}-q_{species}$');ax.axhline(0,color='black',lw=.5);ax.set_ylabel(t);ax.legend(frameon=False,ncol=2,loc='lower right')
axs[-1].set_xlabel('2L coordinate (Mb)');fig.suptitle('Nucleotide and Stage 2R topology tracks (different estimands; not exact quartet match)');savefig(fig,'figure_2_D_genealogy')
# region estimates with existing CI
fig,ax=plt.subplots(figsize=(9,4.5)); designs=['A_merus','A_melas','B_merus','B_melas','C_merus','D_merus']; colors=['#31688e','#5fa4a8','#d1495b','#f08a5d','#756bb1','#8c6d31'];
for i,t in enumerate(designs):
 for j,reg in enumerate(['inside','outside']):
  x=next((x for x in final if x['test_id']==t and x['region']==reg),None)
  if not x:continue
  d=f(x['D_pooled']);lo=f(x['ci95_low']);hi=f(x['ci95_high']);y=i+(j-.5)*.25;ax.errorbar(d,y,xerr=[[d-lo],[hi-d]],fmt='o',color=colors[i],capsize=2,ms=4)
ax.axvline(0,color='black',lw=.7);ax.set_yticks(range(len(designs)),designs);ax.invert_yaxis();ax.set_xlabel('Pooled nucleotide $D_{Patterson}$ (500-kb physical-block jackknife 95% CI)');ax.set_title('Inside (upper symbol) and outside (lower symbol) estimates');savefig(fig,'figure_3_region_robustness')
# diagnostics callability and n pass/informative for A/B
fig,axs=plt.subplots(1,2,figsize=(9,4))
for t,c in [('A_merus','#31688e'),('B_merus','#d1495b')]:
 ss=[x for x in cd if x['test_id']==t];regs=['left_flank','inside','right_flank','outside']; vv={x['region']:x for x in ss};axs[0].plot(regs,[f(vv[q]['mean_callable_P1']) for q in regs],'-o',label=f'{t} P1',color=c);axs[0].plot(regs,[f(vv[q]['mean_callable_outgroup']) for q in regs],'--s',label=f'{t} outgroup',color=c,alpha=.65);axs[1].plot(regs,[f(vv[q]['mean_pass_sites']) for q in regs],'-o',label=f'{t} passed',color=c);axs[1].plot(regs,[f(vv[q]['mean_informative_contributions']) for q in regs],'--s',label=f'{t} informative',color=c,alpha=.65)
axs[0].set_ylabel('Mean called chromosomes / variant record');axs[1].set_ylabel('Mean site records per 50-kb window');
for ax in axs:ax.tick_params(axis='x',rotation=25);ax.legend(fontsize=7,frameon=False)
fig.suptitle('Sampling and variant-site callability diagnostics');savefig(fig,'figure_4_sampling_diagnostics')
