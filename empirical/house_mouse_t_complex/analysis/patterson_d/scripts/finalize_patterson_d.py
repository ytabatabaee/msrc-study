from pathlib import Path
import csv,math,statistics,subprocess,sys
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results';F=ROOT/'figures';F.mkdir(exist_ok=True)
def rd(n):
 with (R/n).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def wr(n,rows):
 if not rows:return
 with (R/n).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def d(a,b):return (a-b)/(a+b) if a+b else None
raw=rd('patterson_d_window_results.tsv'); designs=sorted(set(r['design'] for r in raw)); inv=[('In(17)1',5747331,6957919),('In(17)2',7435308,13016510),('In(17)3',14068229,18051507),('In(17)4',19401000,40049964)]
segments=inv+[('left_flank',5000000,5747331),('collinear_gap_1',6957919,7435308),('collinear_gap_2',13016510,14068229),('collinear_gap_3',18051507,19401000),('right_flank',40049964,42050000)]
# Physical-block delete-one-block jackknife of inside-minus-outside pooled D.
jack=[];inside=[]
for des in designs:
 for name,lo,hi in segments+[('combined_union',inv[0][1],inv[-1][2])]:
  for bs in (100000,250000,500000,1000000):
   groups={}
   for r in raw:
    if r['design']!=des:continue
    st=int(r['start']); key=st//bs
    inside_flag=(st<int(r['end']) and int(r['end'])>lo and st<hi) if name!='combined_union' else any(st<int(r['end']) and st<b and int(r['end'])>a for _,a,b in inv)
    g=groups.setdefault(key,[0.,0.,0.,0.])
    j=0 if inside_flag else 2;g[j]+=float(r['ABBA']);g[j+1]+=float(r['BABA'])
   ia=sum(x[0] for x in groups.values());ib=sum(x[1] for x in groups.values());oa=sum(x[2] for x in groups.values());ob=sum(x[3] for x in groups.values());diff=d(ia,ib)-d(oa,ob) if d(ia,ib)!=None and d(oa,ob)!=None else None
   reps=[]; abs_reps=[]
   for x in groups.values():
    di=d(ia-x[0],ib-x[1]);do=d(oa-x[2],ob-x[3])
    if di is not None and do is not None:
     reps.append(di-do);abs_reps.append(abs(di)-abs(do))
   se=math.sqrt((len(reps)-1)/len(reps)*sum((x-statistics.mean(reps))**2 for x in reps)) if len(reps)>1 else None
   abs_diff=abs(d(ia,ib))-abs(d(oa,ob)) if d(ia,ib) is not None and d(oa,ob) is not None else None
   abs_se=math.sqrt((len(abs_reps)-1)/len(abs_reps)*sum((x-statistics.mean(abs_reps))**2 for x in abs_reps)) if len(abs_reps)>1 else None
   ni=sum(int(r['informative_sites']) for r in raw if r['design']==des and (int(r['start'])<hi and int(r['end'])>lo if name!='combined_union' else any(int(r['start'])<b and int(r['end'])>a for _,a,b in inv)))
   nt=sum(int(r['informative_sites']) for r in raw if r['design']==des)
   jack.append({'design':des,'region':name,'block_bp':bs,'inside_ABBA':ia,'inside_BABA':ib,'inside_informative_sites':ni,'inside_D':d(ia,ib),'outside_ABBA':oa,'outside_BABA':ob,'outside_informative_sites':nt-ni,'outside_D':d(oa,ob),'inside_minus_outside_D':diff,'jackknife_SE':se,'inside_minus_outside_absD':abs_diff,'absD_jackknife_SE':abs_se,'blocks':len(groups),'inside_informative_blocks':sum(1 for x in groups.values() if x[0]+x[1]>0),'outside_informative_blocks':sum(1 for x in groups.values() if x[2]+x[3]>0),'effective_delete_blocks':len(reps)})
   if bs==500000:inside.append(jack[-1])
wr('patterson_d_jackknife.tsv',jack);wr('patterson_d_inside_outside.tsv',inside)
# Join tree frequency track on exact 5kb starts; document distinct estimands.
top=rd('patterson_d_window_results.tsv')
plot100=[x for x in rd('patterson_d_multiresolution.tsv') if x['resolution_bp']=='100000']
with (Path(__file__).resolve().parents[3]/'results/stage5_unified_tracks/window_unified_tracks.tsv').open() as f: trees=list(csv.DictReader(f,delimiter='\t'))
tree={int(x['start_bp']):x for x in trees};
with (Path(__file__).resolve().parents[3]/'results/stage2_recombination_state_5kb.tsv').open() as f: recomb={int(x['start_bp']):x for x in csv.DictReader(f,delimiter='\t')}
cmp=[]
for r in top:
 if r['design']!='D1_Mmd_t_vs_standard' or int(r['start']) not in tree:continue
 t=tree[int(r['start'])]
 z=recomb.get(int(r['start']),{})
 cmp.append({'design':r['design'],'start':r['start'],'end':r['end'],'patterson_D':r['D'],'ABBA':r['ABBA'],'BABA':r['BABA'],'informative_sites':r['informative_sites'],'q_species':t['q_species'],'q_t_alt':t['q_t_alt'],'q_other':t['q_other'],'tree_D_distinct_estimand':t['D'],'A_distance':t['A_distance'],'A_relative':t['A_relative'],'mean_pair_distance_no_spretus':t['mean_pair_distance_no_spretus'],'TTT_D':t['TTT_D'],'recombination_class':z.get('recombination_class',''),'recombination_evidence':z.get('evidence/topology_code','')})
wr('patterson_d_topology_comparison.tsv',cmp)
# Descriptive correlation after pooling count contributions and topology means into 500-kb physical blocks.
groups={}
for r in cmp:
 k=int(r['start'])//500000;g=groups.setdefault(k,[0.,0.,0.,0.]);g[0]+=float(r['ABBA']);g[1]+=float(r['BABA']);g[2]+=float(r['q_t_alt']);g[3]+=1
xs=[d(g[0],g[1]) for g in groups.values() if d(g[0],g[1]) is not None];ys=[g[2]/g[3] for g in groups.values() if d(g[0],g[1]) is not None]
def corr(x,y):
 mx=statistics.mean(x);my=statistics.mean(y);den=math.sqrt(sum((a-mx)**2 for a in x)*sum((b-my)**2 for b in y));return sum((a-mx)*(b-my) for a,b in zip(x,y))/den if den else None
wr('topology_spatial_correlations.tsv',[{'design':'D1_Mmd_t_vs_standard','resolution_bp':500000,'blocks':len(xs),'pearson_r_D_vs_mean_q_t_alt':corr(xs,ys),'interpretation':'descriptive; adjacent genomic blocks and different estimands; no spatial permutation test'}])
# robustness table: supported scans and limitations as explicit rows.
rob=[{'analysis':'primary_PASS_biallelic','status':'completed','details':'AllMouse PASS SNP/indel source; script retains biallelic A/C/G/T; outgroup fixed REF required.'},{'analysis':'population_designs','status':'completed','details':'Three population-matched contrasts frozen before D; each remains exploratory.'},{'analysis':'transversions_only','status':'not_run','details':'Per-site alleles not retained in output; requires another indexed VCF scan.'},{'analysis':'outgroup_sensitivity','status':'not_run','details':'No independent second outgroup cohort was frozen; avoid post-hoc selection.'},{'analysis':'callability_thresholds','status':'partial','details':'Sites excluded for an entirely missing group/outgroup are tallied; partial call rates and invariant callable bases unavailable in SNP-only VCF.'},{'analysis':'MAC_or_SNP_filter_sensitivity','status':'not_run','details':'No alternate filtering scan run.'},{'analysis':'resolution_sensitivity','status':'completed','details':'Pooled ABBA/BABA at 5, 50, 100, 250, 500 kb; D computed after pooling.'}]
wr('patterson_d_robustness.tsv',rob)
# spatial summary plots
colors={'D1_Mmd_t_vs_standard':'#0072B2','D2_Mmm_t_vs_standard':'#D55E00','D3_Mmc_t_vs_standard':'#009E73'}
for figidx,name in enumerate(['patterson_d_chr17','patterson_d_D_vs_topology','patterson_d_inside_outside','patterson_d_designs','patterson_d_coverage'],1):
 fig,ax=plt.subplots(figsize=(13,5.8))
 if figidx in (1,2):
  for nm,a,b in inv:ax.axvspan(a/1e6,b/1e6,color='#e6ab02',alpha=.14,label=nm if figidx==1 else None)
  use=plot100
  for des in (designs if figidx==1 else ['D1_Mmd_t_vs_standard']):
   rr=[r for r in use if r['design']==des and r['D']!='']
   ax.plot([int(r['start'])/1e6 for r in rr],[float(r['D']) for r in rr],lw=.8,alpha=.8,color=colors[des],label=des,marker='.',linestyle='None',markersize=3)
  if figidx==2:
   # Local full-tip q_t_alt is available only on the exact 5-kb gene-tree overlap.
   qbins={}
   for item in cmp:
    k=int(item['start'])//100000*100000;g=qbins.setdefault(k,[0.,0.]);g[0]+=float(item['q_t_alt']);g[1]+=1
   qitems=[(k,v[0]/v[1]) for k,v in sorted(qbins.items()) if v[1]]
   ax2=ax.twinx();ax2.plot([k/1e6 for k,v in qitems],[v for k,v in qitems],color='black',alpha=.65,lw=0,marker='.',markersize=3,label='tree q_t_alt (available tree windows only)');ax2.set_ylabel('Tree q_t_alt (full-tip track)');ax2.legend(loc='upper left',fontsize=8)
   ax.axvline(37.055,color='#555',ls=':',lw=1.2)
   ax.text(38,.96,'No local gene trees beyond 37.055 Mb',transform=ax.get_xaxis_transform(),fontsize=8,color='#444')
   ax.set_xlim(0,94.987271)
  ax.set_ylabel('Patterson D (100-kb pooled counts)');ax.set_xlabel('chr17 coordinate (Mb)');ax.set_ylim(-1,1);ax.legend(ncol=2,fontsize=8)
 elif figidx==3:
  rr=[r for r in inside if r['region'] in ('In(17)1','In(17)2','In(17)3','In(17)4')]
  for i,des in enumerate(designs):
   z=[r for r in rr if r['design']==des];xs=[j+(i-1)*.22 for j in range(4)];ys=[float(r['inside_minus_outside_D']) for r in z];es=[1.96*float(r['jackknife_SE']) if r['jackknife_SE'] else 0 for r in z];ax.errorbar(xs,ys,yerr=es,fmt='o',label=des,color=colors[des])
  ax.axhline(0,color='black',lw=.8);ax.set_xticks(range(4),[x[0] for x in inv]);ax.set_ylabel('Inside D − outside D (95% block-JK interval)');ax.legend(fontsize=8)
 elif figidx==4:
  rr=[r for r in inside if r['region']=='combined_union'];
  for i,r in enumerate(rr):ax.bar(i,float(r['inside_minus_outside_D']),color=colors[r['design']])
  ax.axhline(0,color='black',lw=.8);ax.set_xticks(range(len(rr)),[f"{r['design']}\n{float(r['inside_minus_outside_D']):.3f}" for r in rr],fontsize=8);ax.set_ylabel('Combined inversion D − outside D')
 else:
  for des in designs:
   rr=[r for r in plot100 if r['design']==des];ax.plot([int(r['start'])/1e6 for r in rr],[float(r['informative_sites']) for r in rr],lw=1.0,color=colors[des],label=des,marker='.',linestyle='None',markersize=3)
  ax2=ax.twinx()
  for des in designs:
   miss={}
   for item in top:
    if item['design']==des:
     k=int(item['start'])//100000*100000;miss[k]=miss.get(k,0)+int(item['sites_excluded_for_missing_population'])
   keys=sorted(miss);ax2.plot([k/1e6 for k in keys],[miss[k] for k in keys],color=colors[des],lw=.7,ls='--',alpha=.55,label=f'{des} missing group')
  ax.set_ylabel('Informative ABBA/BABA SNPs / 100 kb');ax2.set_ylabel('SNP records with an empty group/outgroup');ax.set_xlabel('chr17 coordinate (Mb)');ax.legend(loc='upper left',fontsize=7);ax2.legend(loc='upper right',fontsize=7)
 if figidx==1: ax.set_title('Patterson D across house mouse chromosome 17')
 elif figidx==2: ax.set_title('Full-chromosome nucleotide D; tree support only where available')
 elif figidx==3: ax.set_title('Inside versus outside inversion intervals')
 elif figidx==4: ax.set_title('Arrangement-associated population designs')
 else: ax.set_title('Informative SNP density and group availability')
 fig.tight_layout()
 for ext in ('png','pdf','svg'):fig.savefig(F/f'{name}.{ext}',dpi=220)
 plt.close(fig)
# quantitative report body
summ=[]
for r in inside:
 if r['region'] in ('In(17)1','In(17)2','In(17)3','In(17)4','left_flank','collinear_gap_1','collinear_gap_2','collinear_gap_3','right_flank','combined_union'):
  summ.append(f"- {r['design']} {r['region']}: inside D={float(r['inside_D']):.4f} (weighted ABBA={float(r['inside_ABBA']):.2f}, BABA={float(r['inside_BABA']):.2f}; informative SNPs={r['inside_informative_sites']}); outside D={float(r['outside_D']):.4f}; signed difference={float(r['inside_minus_outside_D']):+.4f} (JK SE {float(r['jackknife_SE']):.4f}); |D| difference={float(r['inside_minus_outside_absD']):+.4f} (JK SE {float(r['absD_jackknife_SE']):.4f}); 500kb nonoverlapping blocks={r['blocks']} (inside informative blocks={r['inside_informative_blocks']}, outside informative blocks={r['outside_informative_blocks']}).")
report='''# House mouse chromosome 17 nucleotide Patterson D\n\n## Data and estimator\n\nThe local archive contains only gene trees/coordinate lists; analysis uses the verified public Harr et al. AllMouse PASS indexed VCF. The three exploratory quartets were frozen before inspecting D. Per-site allele-frequency contributions were pooled within bins before calculating D. The analysis includes only biallelic A/C/G/T SNPs and requires all eight *M. spretus* calls to be ancestral (fixed REF). VCF SNP data have no invariant callable-base denominator.\n\n## Results (500-kb spatial blocks for uncertainty)\n\n'''+ '\n'.join(summ)+'''\n\n## Spatial comparison\n\n`patterson_d_topology_comparison.tsv` matches the Mmd design to the existing 5-kb tree track by exact window start, including available distance metrics and recombination classifications; `topology_spatial_correlations.tsv` reports a descriptive 500-kb correlation. Tree q_t_alt and tree D are not the Patterson D estimand; recombination class and available distance statistics are aligned by exact start. A 500-kb block correlation is reported only descriptively; adjacent linked blocks and different estimands preclude causal interpretation, and no spatial permutation test was run.\n\n## Limitations and robustness\n\nAll three contrasts mix arrangement status with population structure to some degree. The D1 and D2 contrasts use within-subspecies t versus standard groups and third-population standard mice from another subspecies; D3 has only three t-designated Mmc samples. Unphased genotypes estimate carrier-population allele-frequency asymmetry rather than the t-bearing chromosome. SNP ascertainment, ancestral polymorphism, admixture/introgression, and geography may produce asymmetry. Transversion-only and alternate-outgroup scans were not run. The PASS callset is the primary filter; the SNP-only source cannot quantify invariant-base callability. Therefore results do not uniquely support introgression or MSRC.\n\nSee `patterson_d_scientific_interpretation.md` for interpretation and all requested questions.\n'''
(R/'patterson_d_full_report.md').write_text(report)
(R/'patterson_d_scientific_interpretation.md').write_text('''# Scientific interpretation\n\n1. **Can D be computed?** Yes, from the verified indexed public nucleotide VCF, for SNP site patterns only.\n2. **Valid quartets?** Three rooted exploratory population contrasts are frozen in config. *M. spretus* is a close, separate outgroup supported by the study. The designs do not represent conventional species-tree nulls; they compare arrangement-defined subsets with standard population groups.\n3. **Stronger inside inversions?** No consistent strengthening by absolute D: the combined inversion |D| differences are design-dependent (reported with block-jackknife SE in the main report and TSV). They cannot be interpreted as t-haplotype effects without phasing and stronger population matching.\n4. **Interval variation?** Individual intervals are summarized separately in `patterson_d_jackknife.tsv`.\n5. **Same regions as gene trees?** Exact 5-kb window matching is tabulated for the Mmd quartet. For ((P1,P2),P3),O, ABBA corresponds to derived sharing by P2+P3 and BABA by P1+P3; the available q_t_alt/q_species full-tip tracks do not use those exact sample populations, so any co-location is descriptive, not a direct topology match.\n6. **Robustness?** Resolution pooling is evaluated. Transversions, outgroup alternatives, and alternate SNP filters remain untested; see robustness table.\n7. **Recombination partitioning?** These scan summaries do not establish a recombination-partitioned signal. Spatial compatibility with locally persistent genealogies is only a hypothesis; no formal location concordance test was done.\n8. **Alternatives:** population structure, admixture or introgression, ancestral population structure/polymorphism, incomplete lineage sorting, SNP ascertainment, and polarization/missingness.\n9. **Neoaves/Anopheles:** qualitative analogy only: all use nucleotide D and spatial aggregation, but taxon sampling, SNP/alignment source, and quartet interpretation differ. Numerical comparisons are inappropriate.\n\nNo nonzero D uniquely supports introgression or MSRC.\n''')
