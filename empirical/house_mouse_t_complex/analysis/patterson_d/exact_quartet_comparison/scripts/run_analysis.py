#!/usr/bin/env python3
"""Score published local trees with exact frozen population sample sets."""
import io,zipfile,csv,json,sys,hashlib,math,re,statistics
from pathlib import Path
from collections import defaultdict
from Bio import Phylo
import numpy as np
BASE=Path(__file__).resolve().parents[1]; ROOT=BASE.parents[3]; REPO=BASE.parents[4]
P=BASE.parent; ORIGINAL=P/'results'; POL=P/'polarization_validation'/'results'; OUT=BASE/'results';OUT.mkdir(parents=True,exist_ok=True)
ARCH=REPO/'IST-2017-78-v1+1_Data.zip'
sys.path.insert(0,str(P/'scripts'))
from run_patterson_d import GROUPS,INV,D
MEM=P/'config'/'frozen_quartet_members.tsv'
TREE_MEMBER='Data/2-Coverage-and_AlleleRatio-Filtered_RAW_SNPs/1-Trees_for_all_5kb_windows/ML_trees.zip'
def readtsv(p):
 with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(name,rows):
 if not rows:return
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def decode_tip(tip,samples):
 for s in sorted(samples,key=len,reverse=True):
  if tip.startswith(s+'_'):
   suffix=tip[len(s)+1:]
   if suffix.startswith('tHaplSubset.fa_') or re.match(r'^\d+-\d+$',suffix):return s
 return None
def parse_coord(name):
 m=re.search(r'/(\d+)-(\d+)\.fa\.contree$',name)
 if not m:return None
 return int(m.group(1))+5_000_000,int(m.group(2))+5_000_000
# Match each edge bipartition to one of the three unrooted quartet splits.
def tree_stats(newick,groups,tipmap):
 tree=Phylo.read(io.StringIO(newick),'newick'); terminals=[t for t in tree.get_terminals() if t.name in tipmap]
 by_sample={tipmap[t.name]:t for t in terminals}; present=[[s for s in groups[g] if s in by_sample] for g in range(4)]; ns=[len(x) for x in present]; total=math.prod(ns)
 if not total:return {'n_groups':ns,'possible':0,'resolved':0,'unresolved':0,'q1':None,'q2':None,'q3':None,'D_topology':None,'mean_internal_branch':None}
 samples=list(dict.fromkeys(s for group in present for s in group)); idx={s:i for i,s in enumerate(samples)}
 dist=np.zeros((len(samples),len(samples)),dtype=float)
 for i,s in enumerate(samples):
  for j in range(i):dist[i,j]=dist[j,i]=tree.distance(by_sample[s],by_sample[samples[j]])
 arrays=[[idx[s] for s in group] for group in present]
 grid=np.indices([len(x) for x in arrays],dtype=np.int32).reshape(4,-1)
 p1,p2,p3,o=[np.asarray(arrays[i],dtype=np.int32)[grid[i]] for i in range(4)]
 sums=np.vstack((dist[p1,p2]+dist[p3,o],dist[p1,p3]+dist[p2,o],dist[p1,o]+dist[p2,p3]))
 order=np.sort(sums,axis=0); resolved=(order[1]-order[0])>1e-12
 winner=np.argmin(sums[:,resolved],axis=0);counts=np.bincount(winner,minlength=3).astype(float);nr=int(resolved.sum())
 q=(counts/nr).tolist() if nr else [None]*3
 internal=(order[1,resolved]-order[0,resolved])/2 if nr else np.asarray([])
 return {'n_groups':ns,'possible':total,'resolved':nr,'unresolved':total-nr,'q1':q[0],'q2':q[1],'q3':q[2],'D_topology':(q[2]-q[1]) if nr else None,'mean_internal_branch':float(internal.mean()) if nr else None}

def main():
 members=readtsv(MEM); designs={d:[set() for _ in range(4)] for d in GROUPS}; roleidx={'P1':0,'P2':1,'P3':2,'O':3}
 for x in members:designs[x['design']][roleidx[x['role']]].add(x['sample_id'])
 if any(designs[d]!=list(GROUPS[d]) for d in GROUPS):raise RuntimeError('Frozen member file differs from original runner membership; stop for audit.')
 samples=set().union(*(set().union(*designs[d]) for d in designs)); audit=[];
 with zipfile.ZipFile(ARCH) as outer:
  nested=zipfile.ZipFile(io.BytesIO(outer.read(TREE_MEMBER))); treefiles=[n for n in nested.namelist() if n.endswith('.contree')]
  seen=set(); treedata=[]
  for fname in treefiles:
   coord=parse_coord(fname)
   if not coord:continue
   start,end=coord; nw=nested.read(fname).decode(); tr=Phylo.read(io.StringIO(nw),'newick'); labels=[t.name for t in tr.get_terminals()];tipmap={}
   for t in labels:
    sm=decode_tip(t,samples)
    if sm:tipmap[t]=sm;seen.add(sm)
   treedata.append((start,end,nw,tipmap))
 for d,gs in designs.items():
  for i,g in enumerate(gs):
   for sm in sorted(g):audit.append({'design':d,'role':('P1','P2','P3','O')[i],'sample_id':sm,'tree_tip_found':sm in seen,'tree_representation':'t carrier uses archived tHaplSubset.fa sequence tip' if any(sm in t for _,_,_,mp in treedata for t in mp) and 'tHaplSubset' in next((t for _,_,_,mp in treedata for t in mp if mp[t]==sm),'') else 'archived individual sequence tip'})
 rows=[]
 for start,end,nw,tipmap in treedata:
  for d,gs in designs.items():
   x=tree_stats(nw,gs,tipmap);rows.append({'chromosome':'chr17','start':start,'end':end,'design':d,'n_P1':x['n_groups'][0],'n_P2':x['n_groups'][1],'n_P3':x['n_groups'][2],'n_O':x['n_groups'][3],'requested_P1':len(gs[0]),'requested_P2':len(gs[1]),'requested_P3':len(gs[2]),'requested_O':len(gs[3]),'missing_samples':sum(len(gs[i])-x['n_groups'][i] for i in range(4)),'possible_induced_quartets':x['possible'],'eligible_resolved_quartets':x['resolved'],'unresolved_quartets':x['unresolved'],'q1':x['q1'],'q2':x['q2'],'q3':x['q3'],'D_topology':x['D_topology'],'mean_internal_branch_length':x['mean_internal_branch']} )
 write('sample_mapping_audit.tsv',audit);write('matched_gene_tree_quartets.tsv',rows)
 # Use one gene tree as one locus. Aggregate resolved quartet counts within loci, then equally weight loci for q frequencies.
 pol=readtsv(POL/'corrected_window_D.tsv'); nucleotide={(x['design'],int(x['start'])):x for x in pol}; groups=defaultdict(list)
 for r in rows:groups[(r['design'],int(r['start']))].append(r)
 matched=[]
 for k,ts in groups.items():
  d,st=k; tree=ts[0]; nt=nucleotide.get(k); matched.append({'chromosome':'chr17','start':st,'end':int(tree['end']),'design':d,'ABBA':nt['ABBA'] if nt else '', 'BABA':nt['BABA'] if nt else '', 'D_site':nt['D'] if nt else '', 'informative_sites':nt['informative_sites'] if nt else '', 'q1':tree['q1'],'q2':tree['q2'],'q3':tree['q3'],'D_topology':tree['D_topology'],'eligible_quartets':tree['eligible_resolved_quartets'],'possible_quartets':tree['possible_induced_quartets'],'missing_samples':tree['missing_samples'],'mean_internal_branch_length':tree['mean_internal_branch_length'],'region':classify(st,int(tree['end']))})
 write('matched_window_dataset.tsv',matched)
 # Multi-resolution paired tracks: pool SNP counts, average locus-level q values with equal locus weight.
 multi=[]
 for res in (50000,100000,250000,500000):
  bins=defaultdict(list)
  for x in matched:bins[(x['design'],int(x['start'])//res*res)].append(x)
  for (d,st),xs in sorted(bins.items()):
   xs=[x for x in xs if x['D_site'] not in ('',None) and x['D_topology'] not in ('',None)]
   if not xs:continue
   a=sum(float(x['ABBA']) for x in xs);b=sum(float(x['BABA']) for x in xs)
   res_cov=statistics.mean(int(x['eligible_quartets'])/int(x['possible_quartets']) for x in xs if int(x['possible_quartets'])>0)
   for field in ('q1','q2','q3','D_topology','mean_internal_branch_length'):
    vals=[float(x[field]) for x in xs if x[field] not in ('',None)]
    locals()[field]=statistics.mean(vals) if vals else float('nan')
   multi.append({'design':d,'resolution_bp':res,'start':st,'end':min(st+res,37_055_000),'ABBA':a,'BABA':b,'D_site':D(a,b),'informative_sites':sum(int(x['informative_sites']) for x in xs),'q1':statistics.mean(float(x['q1']) for x in xs if x['q1'] not in ('',None)),'q2':statistics.mean(float(x['q2']) for x in xs if x['q2'] not in ('',None)),'q3':statistics.mean(float(x['q3']) for x in xs if x['q3'] not in ('',None)),'D_topology':statistics.mean(float(x['D_topology']) for x in xs if x['D_topology'] not in ('',None)),'n_loci':len(xs),'eligible_quartets_sum_diagnostic_only':sum(int(x['eligible_quartets']) for x in xs),'mean_locus_quartet_resolution_fraction':res_cov,'mean_internal_branch_length':statistics.mean(float(x['mean_internal_branch_length']) for x in xs if x['mean_internal_branch_length'] not in ('',None)),'region':classify(st,min(st+res,37_055_000))})
 write('matched_nucleotide_genealogy_tracks.tsv',multi)
 # Associations per resolution with block-bootstrap CI; circular shifts supply a spatial null, not independent Pearson p-values.
 stats=[]
 for res in (50000,100000,250000,500000):
  for d in GROUPS:
   base=[x for x in multi if x['resolution_bp']==res and x['design']==d]
   for min_sites in (1,5,10,25):
    xs=[x for x in base if int(x['informative_sites'])>=min_sites]
    for region in ('matched_all','inversion','outside_inversion'):
     sub=xs if region=='matched_all' else ([x for x in xs if any(a<=int(x['start']) and int(x['end'])<=b for _,a,b in INV)] if region=='inversion' else [x for x in xs if not any(int(x['start'])<b and int(x['end'])>a for _,a,b in INV)])
     val=association(sub)
     stats.append({'design':d,'resolution_bp':res,'minimum_informative_sites':min_sites,'minimum_quartet_resolution_fraction':'','region':region,'n_windows':len(sub),'pearson_r':val[0],'spearman_rho':val[1],'sign_agreement':val[2],'block_bootstrap_low':val[3],'block_bootstrap_high':val[4],'circular_shift_null_abs_r_p':val[5],'block_level_pearson_r':val[6]})
   for min_cov in (.8,.95,.999):
    xs=[x for x in base if float(x['mean_locus_quartet_resolution_fraction'])>=min_cov and int(x['informative_sites'])>=1]
    val=association(xs);stats.append({'design':d,'resolution_bp':res,'minimum_informative_sites':1,'minimum_quartet_resolution_fraction':min_cov,'region':'matched_all','n_windows':len(xs),'pearson_r':val[0],'spearman_rho':val[1],'sign_agreement':val[2],'block_bootstrap_low':val[3],'block_bootstrap_high':val[4],'circular_shift_null_abs_r_p':val[5],'block_level_pearson_r':val[6]})
 write('spatial_correlations.tsv',stats)
 print(json.dumps({'trees':len(treedata),'tree_range':[min(x[0] for x in treedata),max(x[1] for x in treedata)],'design_rows':len(rows),'matched_rows':len(matched),'trees_with_all_samples_by_design':{d:sum(r['design']==d and r['missing_samples']=='0' for r in rows) for d in GROUPS}},indent=2))
def classify(st,en):
 for n,a,b in INV:
  if st<b and en>a:return n
 for n,a,b in [('left_flank',5000000,5747331),('collinear_gap_1',6957919,7435308),('collinear_gap_2',13016510,14068229),('collinear_gap_3',18051507,19401000),('right_flank',40049964,42050000)]:
  if st<b and en>a:return n
 return 'matched_collinear'
def ranks(a):
 order=sorted(range(len(a)),key=lambda i:a[i]);out=[0.]*len(a)
 i=0
 while i<len(order):
  j=i+1
  while j<len(order) and a[order[j]]==a[order[i]]:j+=1
  v=(i+j-1)/2+1
  for h in order[i:j]:out[h]=v
  i=j
 return out
def pear(x,y):
 if len(x)<3:return float('nan')
 mx=statistics.mean(x);my=statistics.mean(y);dx=[a-mx for a in x];dy=[b-my for b in y];den=math.sqrt(sum(a*a for a in dx)*sum(b*b for b in dy));return sum(a*b for a,b in zip(dx,dy))/den if den else float('nan')
def association(xs):
 xs=[x for x in xs if x['informative_sites']!='' and int(x['informative_sites'])>=1 and x['D_topology'] not in ('',None)]
 x=[float(z['D_site']) for z in xs];y=[float(z['D_topology']) for z in xs];r=pear(x,y);rho=pear(ranks(x),ranks(y));sign=sum((a==0 and b==0) or a*b>0 for a,b in zip(x,y))/len(x) if x else float('nan')
 # Resample physical 500-kb blocks; keep every matched window inside each selected block together.
 import random;seed=sum((i+1)*ord(c) for i,c in enumerate(str([(z.get('design'),z.get('resolution_bp'),z.get('region')) for z in xs[:1]])));rng=random.Random(seed);by=defaultdict(list)
 for i,z in enumerate(xs):by[int(z['start'])//500000].append(i)
 keys=sorted(by);boot=[]
 if len(keys)>=3:
  for _ in range(500):
   ids=[]
   for k in (rng.choice(keys) for _ in keys):ids.extend(by[k])
   z=pear([x[i] for i in ids],[y[i] for i in ids])
   if math.isfinite(z):boot.append(z)
 lo=sorted(boot)[int(.025*len(boot))] if boot else float('nan');hi=sorted(boot)[min(len(boot)-1,int(.975*len(boot)))] if boot else float('nan')
 # Build one paired site/topology point per physical 500-kb block for a circular-shift spatial null.
 bx=[];bytop=[]
 for k in keys:
  rr=[xs[i] for i in by[k]];aa=sum(float(z['ABBA']) for z in rr);bb=sum(float(z['BABA']) for z in rr);q=[float(z['D_topology']) for z in rr]
  if aa+bb and q:bx.append((aa-bb)/(aa+bb));bytop.append(statistics.mean(q))
 rb=pear(bx,bytop);null=[]
 if len(bx)>=5:
  for sh in range(1,len(bx)):
   null.append(abs(pear(bx,bytop[sh:]+bytop[:sh])))
 p=(1+sum(z>=abs(rb) for z in null))/(1+len(null)) if null and math.isfinite(rb) else float('nan')
 return r,rho,sign,lo,hi,p,rb
if __name__=='__main__':main()
