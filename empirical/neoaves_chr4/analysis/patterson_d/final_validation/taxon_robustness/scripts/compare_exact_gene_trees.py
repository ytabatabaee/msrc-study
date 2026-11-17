#!/usr/bin/env python3
"""Compare prespecified nucleotide quartets to exact matched Stiller trees."""
import csv,gzip,io,math,random,re
from collections import Counter,defaultdict
from pathlib import Path
from Bio import Phylo
ROOT=Path(__file__).resolve().parents[7]
BASE=ROOT/'empirical/neoaves_chr4/analysis/patterson_d';PKG=BASE/'final_validation/taxon_robustness'
TREE=ROOT/'empirical/neoaves_chr4/external/stiller2024/63430.named.gene.trees.gz'
REF=ROOT/'empirical/neoaves_chr4/external/stiller2024/63K.tre'
MAN=ROOT/'empirical/neoaves_chr4/results/stage4d_locus_manifest.tsv'
LOC=PKG/'results/alternative_quartet_locus_counts.tsv'
CONFIG=PKG/'config/frozen_alternative_quartets.tsv'
INTERVALS=[('PNAS_1',25030000,32670000),('PNAS_2',33510000,34470000),('PNAS_3',44130000,56810000)]

def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 with p.open('w',newline='') as f:
  if not rows:return
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def reg(s):
 for n,a,b in INTERVALS:
  if a<=s<=b:return n
 return 'outside_all_PNAS'
def tree_splits(tree):
 # Extract exact unrooted edge bipartitions once; an induced quartet is resolved
 # iff one edge separates two selected tips from the other two.
 root=tree.root;splits=[]
 def walk(node):
  if node.is_terminal():names={node.name}
  else:
   names=set()
   for child in node.clades:names.update(walk(child))
  if node is not root:splits.append(names)
  return names
 walk(root)
 return splits

def quartet_from_splits(splits,q):
 p1,p2,p3,o=q
 pairs={frozenset((p1,p2)):'q1',frozenset((p3,o)):'q1',frozenset((p1,p3)):'q2',frozenset((p2,o)):'q2',frozenset((p2,p3)):'q3',frozenset((p1,o)):'q3'}
 for split in splits:
  side=split.intersection(q)
  if len(side)==2:
   return pairs.get(frozenset(side),pairs.get(frozenset(set(q)-side)))
 return None

def quartet(tree,q):return quartet_from_splits(tree_splits(tree),q)

def corr(x,y):
 n=len(x)
 if n<3:return None
 mx=sum(x)/n;my=sum(y)/n
 vx=sum((a-mx)**2 for a in x);vy=sum((b-my)**2 for b in y)
 if not vx or not vy:return None
 return sum((a-mx)*(b-my) for a,b in zip(x,y))/math.sqrt(vx*vy)
def ci(vals):
 vals=sorted(x for x in vals if x is not None)
 if not vals:return ('','')
 return vals[int(.025*(len(vals)-1))],vals[int(.975*(len(vals)-1))]
def main():
 cfg=read(CONFIG);qs={r['design']:(r['P1'],r['P2'],r['P3'],r['O']) for r in cfg}
 manifest={r['locus_id']:r for r in read(MAN) if r['chromosome']=='chr4'}
 nucleotide=defaultdict(dict)
 for r in read(LOC):nucleotide[r['design']][r['locus_id']]=r
 ref=Phylo.read(str(REF),'newick');refres={}
 for name,q in qs.items():
  missing=sorted(set(q)-{t.name for t in ref.get_terminals()})
  if missing:res='unavailable'
  else:
   t=Phylo.read(str(REF),'newick')
   for tip in list(t.get_terminals()):
    if tip.name not in q:t.prune(tip)
   res=quartet(t,q) or 'unresolved'
  refres[name]=res
 tree_counts={d:Counter() for d in qs};whole_taxa=Counter();matched_taxa=Counter();topologies={d:{} for d in qs}
 nchr4=0;matched_bins=[]
 with gzip.open(TREE,'rt') as f:
  for line in f:
   line=line.rstrip('\n')
   if '\t' not in line:continue
   locus,nwk=line.split('\t',1)
   tips=set(re.findall(r'(?:\(|,)([^():,]+)(?=[:),])',nwk))
   whole_taxa.update(tips)
   if locus not in manifest:continue
   nchr4+=1;m=manifest[locus];st=int(m['start']);rgn=reg(st)
   active=[]
   for d,q in qs.items():
    if set(q)<=tips:
     tree_counts[d]['taxa_present']+=1;active.append((d,q))
    else:tree_counts[d]['missing_taxa']+=1
   if active:
    tr=Phylo.read(io.StringIO(nwk),'newick');splits=tree_splits(tr)
    for d,q in active:
     top=quartet_from_splits(splits,q)
     if top:topologies[d][locus]=(top,st,rgn);tree_counts[d][top]+=1
     else:tree_counts[d]['unresolved']+=1
   if nchr4%1000==0:print('matched chr4 gene trees',nchr4,flush=True)
 # Candidate availability in alignment archive is inferred from exact joint-coverage count
 # for each design, while per-species coverage comes from the archive taxon inventory.
 qrows=[]
 for c in cfg:
  d=c['design'];ct=tree_counts[d];resolved=sum(ct[x] for x in ('q1','q2','q3'))
  qrows.append(dict(design=d,P1=c['P1'],P2=c['P2'],P3=c['P3'],O=c['O'],source_reference_topology=refres[d],chr4_exact_quartet_trees=resolved,q1=ct['q1'],q2=ct['q2'],q3=ct['q3'],unresolved_exact_quartets=ct['unresolved'],trees_missing_one_or_more_taxa=ct['missing_taxa'],D_topology_q3_minus_q2=((ct['q3']-ct['q2'])/resolved if resolved else ''),interpretation='q1 is intended rooted null under exact Stiller 363-taxon reference; local tree counts are conditional on four exact taxa.'))
 write(PKG/'results/quartet_null_validation.tsv',qrows)
 # Taxon audit: species/order identities are explicit; no broad clade is assumed monophyletic.
 taxmeta={
 'Tauraco_erythrolophus':('Otidimorphae','Musophagiformes; turaco','Primary P1 representative; one of three order-level Otidimorphae lineages.'),
 'Cuculus_canorus':('Otidimorphae','Cuculiformes; cuckoo','Alternative P1 from another Otidimorphae order.'),
 'Ardeotis_kori':('Otidimorphae','Otidiformes; bustard','Alternative P1 from a third Otidimorphae order.'),
 'Columba_livia':('Columbimorphae','Columbiformes; pigeon','Primary P2 representative.'),
 'Pterocles_gutturalis':('Columbimorphae','Pterocliformes; sandgrouse','Alternative P2 from a distinct order.'),
 'Mesitornis_unicolor':('Columbimorphae (Jarvis/Columbea framework)','Mesitornithiformes; mesite','Alternative P2; mesite affinity/order-level clade placement is topology-sensitive.'),
 'Ciconia_maguari':('P3 side of N61 quartet','Aequornithes; Ciconiiformes; stork','Primary P3. It is not itself Telluraves; P3 denotes the other side of the focal sister relationship in the frozen quartet.'),
 'Balaeniceps_rex':('P3 broad-complement candidate','Pelecaniformes; shoebill; Aequornithes','P3 waterbird alternative; not a close Ciconia family replicate.'),
 'Nipponia_nippon':('P3 broad-complement candidate','Ciconiiformes; ibis','P3 waterbird alternative; not a close Ciconia family replicate.'),
 'Aquila_chrysaetos':('P3 Telluraves control','Accipitriformes; eagle; Telluraves','Explicit broad-side control, not a same-clade Ciconia substitute.'),
 'Gallus_gallus':('Outgroup','Galliformes','Primary external outgroup.'),
 'Meleagris_gallopavo':('Outgroup sensitivity','Galliformes','Alternative external galliform; correlated outgroup sensitivity.')}
 # Full tree archive tip counts (exact names) and reference membership are computed below in a separate pass.
 alltree=whole_taxa
 refnames={t.name for t in ref.get_terminals()}
 # Alignment taxon coverage is calculated from the in-scan inventory file once written.
 acov={}
 acpath=PKG/'results/candidate_alignment_coverage.tsv'
 if acpath.exists():acov={r['taxon']:r for r in read(acpath)}
 trows=[]
 for tax,(group,assignment,relationship) in taxmeta.items():
  ar=acov.get(tax,{})
  designs_for=[d for d,q in qs.items() if tax in q]
  qstates=[f'{d}:{refres[d]}' for d in designs_for]
  trows.append(dict(scientific_name=tax.replace('_',' '),taxon_id=tax,published_lineage_assignment=group,order_assignment=assignment,present_in_nucleotide_alignments=ar.get('present_in_chr4_alignments','pending'),usable_chr4_locus_coverage=ar.get('present_in_chr4_alignments','pending'),percent_of_chr4_loci=ar.get('percent_of_6247_chr4_loci','pending'),present_in_gene_tree_archive=('yes' if alltree[tax] else 'no'),gene_tree_loci_with_taxon=alltree[tax],present_in_Stiller_363_reference=('yes' if tax in refnames else 'no'),published_species_history_relationship=relationship,designs_in_panel=';'.join(designs_for),Stiller_reference_induced_quartets=';'.join(qstates),independent_48_taxon_exact_relationship='not established; Jarvis broad tree lacks Ciconia maguari',chr4_exact_quartet_trees_for_any_panel=len(set().union(*(set(topologies[d]) for d in qs if tax in qs[d])))))
 write(PKG/'results/candidate_taxa.tsv',trows)
 # Bin matched Dsite and exact-tree Dtopology by exact locus IDs; classify spatial bins using their start coordinate, matching the frozen final-validation method.
 robust=[]
 for cfgrow in cfg:
  d=cfgrow['design'];q=qs[d];tops=topologies[d]
  for bs in (250000,500000,1000000):
   bins=defaultdict(lambda:Counter())
   for locus,(top,st,rgn) in tops.items():
    nr=nucleotide[d].get(locus)
    if nr is None:continue
    k=st//bs;c=bins[k];c['n']+=1;c[top]+=1;c['ABBA']+=int(nr['ABBA']);c['BABA']+=int(nr['BABA']);c['tv_ABBA']+=int(nr['transversion_ABBA']);c['tv_BABA']+=int(nr['transversion_BABA']);c['inside']+=int(rgn!='outside_all_PNAS');c['outside']+=int(rgn=='outside_all_PNAS')
   for filt,ak,bk in [('all','ABBA','BABA'),('transversions_only','tv_ABBA','tv_BABA')]:
    data=[]
    for k,c in sorted(bins.items()):
     nq=sum(c[z] for z in ('q1','q2','q3'));ns=c[ak]+c[bk]
     if not nq or not ns:continue
     row=dict(design=d,filter=filt,resolution_bp=bs,window_start=k*bs,window_end=(k+1)*bs,region='inside' if any(a<=k*bs<=b for _,a,b in INTERVALS) else 'outside',n_exact_gene_trees=c['n'],q1=c['q1'],q2=c['q2'],q3=c['q3'],D_topology=(c['q3']-c['q2'])/nq,D_site=(c[ak]-c[bk])/ns,ABBA=c[ak],BABA=c[bk],informative_sites=ns)
     matched_bins.append(row)
     bin_start=k*bs
     bin_region='inside' if any(a<=bin_start<=b for _,a,b in INTERVALS) else 'outside'
     data.append(dict(block=k,mid=(k+.5)*bs,region=bin_region,n=c['n'],q1=c['q1']/nq,q2=c['q2']/nq,q3=c['q3']/nq,dt=row['D_topology'],ds=row['D_site'],ABBA=c[ak],BABA=c[bk]))
    def add_corr(label,dd):
     x=[z['ds'] for z in dd];y=[z['dt'] for z in dd]
     whole=corr(x,y)
     if label=='whole':
      centered_data=[z for z in data if z['region'] in ('inside','outside')]
      means={r:(sum(z['ds'] for z in centered_data if z['region']==r)/max(1,sum(z['region']==r for z in centered_data)),sum(z['dt'] for z in centered_data if z['region']==r)/max(1,sum(z['region']==r for z in centered_data))) for r in ('inside','outside')}
      xc=[z['ds']-means[z['region']][0] for z in centered_data];yc=[z['dt']-means[z['region']][1] for z in centered_data]
      centered=corr(xc,yc)
     else:centered=''
     nwt=sum(z['n'] for z in dd)
     weighted=lambda key:(sum(z[key]*z['n'] for z in dd)/nwt if nwt else '')
     robust.append(dict(design=d,filter=filt,resolution_bp=bs,scope=label,n_bins=len(dd),n_exact_gene_trees=nwt,mean_q1=weighted('q1'),mean_q2=weighted('q2'),mean_q3=weighted('q3'),mean_D_topology=weighted('dt'),mean_D_site=weighted('ds'),pearson_Dsite_Dtopology=whole,within_region_centered_pearson=centered,block_bootstrap_CI_low='',block_bootstrap_CI_high='',centered_bootstrap_CI_low='',centered_bootstrap_CI_high='',bootstrap_block_bp=2000000,bootstrap_replicates=1000,seed=20261118,interpretation='Descriptive locus-matched bin correlations; linked blocks resampled, quartets are not independent.'))
    for scope in ('whole','inside','outside'):
     dd=data if scope=='whole' else [z for z in data if z['region']==scope]
     add_corr(scope,dd)
    # Bootstrap whole-chromosome and region-centered correlations by resampling 2-Mb blocks.
    blocks=defaultdict(list)
    for z in data:blocks[int(z['mid']//2000000)].append(z)
    keys=list(blocks);rng=random.Random(20261118+bs+sum(map(ord,d))+len(filt))
    for target in ('whole','inside','outside','within_region_centered'):
     boots=[]
     if len(keys)>=2:
      for _ in range(1000):
       sample=[]
       for k in (rng.choice(keys) for j in range(len(keys))):sample.extend(blocks[k])
       x=[z['ds'] for z in sample];y=[z['dt'] for z in sample]
       if target=='within_region_centered':
        means={r:(sum(z['ds'] for z in sample if z['region']==r)/max(1,sum(z['region']==r for z in sample)),sum(z['dt'] for z in sample if z['region']==r)/max(1,sum(z['region']==r for z in sample))) for r in ('inside','outside')}
        sample=[z for z in sample if z['region'] in ('inside','outside')]
        if len(sample)<3:continue
        means={r:(sum(z['ds'] for z in sample if z['region']==r)/max(1,sum(z['region']==r for z in sample)),sum(z['dt'] for z in sample if z['region']==r)/max(1,sum(z['region']==r for z in sample))) for r in ('inside','outside')}
        x=[z['ds']-means[z['region']][0] for z in sample];y=[z['dt']-means[z['region']][1] for z in sample]
       boots.append(corr(x,y))
     low,high=ci(boots)
     for z in robust:
      if z['design']==d and z['filter']==filt and z['resolution_bp']==bs and z['scope']==('whole' if target=='within_region_centered' else target):
       if target=='within_region_centered':z['centered_bootstrap_CI_low']=low;z['centered_bootstrap_CI_high']=high
       else:z['block_bootstrap_CI_low']=low;z['block_bootstrap_CI_high']=high
 write(PKG/'results/nucleotide_topology_robustness.tsv',robust)
 write(PKG/'results/nucleotide_topology_matched_bins.tsv',matched_bins)
 print('tree designs',len(qrows),'correlation rows',len(robust),'taxon audit rows',len(trows),flush=True)
if __name__=='__main__':main()
