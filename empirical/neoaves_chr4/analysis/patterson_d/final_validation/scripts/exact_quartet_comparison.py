#!/usr/bin/env python3
"""Match frozen nucleotide D counts to exact induced quartets in Stiller trees."""
import csv, gzip, io, math
from collections import Counter, defaultdict
from pathlib import Path
from Bio import Phylo

ROOT=Path(__file__).resolve().parents[6]
BASE=ROOT/'empirical/neoaves_chr4/analysis/patterson_d'
OUT=BASE/'final_validation/results'
TREE=ROOT/'empirical/neoaves_chr4/external/stiller2024/63430.named.gene.trees.gz'
MAN=ROOT/'empirical/neoaves_chr4/results/stage4d_locus_manifest.tsv'
LOC=BASE/'results/patterson_d_locus_results.tsv'
REF=ROOT/'empirical/neoaves_chr4/external/stiller2024/63K.tre'
QUARTETS={
 'Columbea_primary':('Columba_livia','Phoenicopterus_ruber','Ciconia_maguari','Gallus_gallus'),
 'N61_primary':('Tauraco_erythrolophus','Columba_livia','Ciconia_maguari','Gallus_gallus'),
 'N62_primary':('Columba_livia','Pterocles_gutturalis','Tauraco_erythrolophus','Gallus_gallus')}

def read_tsv(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows,fields):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def induced(tree,q):
 names={x.name for x in tree.get_terminals()}
 if not set(q)<=names:return None
 d=lambda a,b:tree.distance(a,b)
 sums=(d(q[0],q[1])+d(q[2],q[3]),d(q[0],q[2])+d(q[1],q[3]),d(q[0],q[3])+d(q[1],q[2]))
 m=min(sums); ix=[i for i,v in enumerate(sums) if math.isclose(v,m,rel_tol=1e-9,abs_tol=1e-9)]
 if len(ix)!=1:return None,None
 other=sorted(sums)
 return ('q1','q2','q3')[ix[0]],max(0.0,(other[1]-other[0])/2.0)

def main():
 manifest={r['locus_id']:r for r in read_tsv(MAN) if r['chromosome']=='chr4'}
 sites=defaultdict(tuple)
 for r in read_tsv(LOC):
  if r['filter']=='all' and r['design'] in QUARTETS:sites[(r['design'],r['locus_id'])]=(int(r['ABBA']),int(r['BABA']))
 exact={d:{} for d in QUARTETS}; branches={d:defaultdict(list) for d in QUARTETS}; present=Counter()
 with gzip.open(TREE,'rt') as fh:
  for line in fh:
   locus,nwk=line.rstrip('\n').split('\t',1)
   if locus not in manifest:continue
   tree=Phylo.read(io.StringIO(nwk),'newick')
   for d,q in QUARTETS.items():
    present[d]+=int(all(x in nwk for x in q))
    found=induced(tree,q)
    if found and found[0]:
     top,edge=found; exact[d][locus]=top; branches[d][top].append(edge)
 print('resolved trees', {d:len(exact[d]) for d in QUARTETS},flush=True)
 # 1-Mb bins on the Stiller locus-coordinate axis; only exact quartet trees contribute to both statistics.
 bins=defaultdict(lambda:Counter())
 for d,loci in exact.items():
  for locus,top in loci.items():
   b=int(manifest[locus]['start'])//1_000_000; a,bb=sites.get((d,locus),(0,0)); c=bins[(d,b)]
   c['n']+=1;c[top]+=1;c['ABBA']+=a;c['BABA']+=bb
 rows=[]
 for (d,b),c in sorted(bins.items()):
  qden=sum(c[x] for x in ('q1','q2','q3'));den=c['ABBA']+c['BABA']
  rows.append(dict(design=d,window_start=b*1_000_000,window_end=(b+1)*1_000_000,n_exact_quartet_trees=c['n'],n_q1=c['q1'],n_q2=c['q2'],n_q3=c['q3'],D_topology=(c['q3']-c['q2'])/qden if qden else '',nucleotide_ABBA=c['ABBA'],nucleotide_BABA=c['BABA'],D_site_pooled=(c['ABBA']-c['BABA'])/den if den else '',informative_sites=den))
 fields=list(rows[0]) if rows else []
 write_tsv(OUT/'nucleotide_topology_comparison.tsv',rows,fields)
 # Full chr4 matched quartet frequencies; and induced species-history topology from published 363-taxon reference.
 ref=Phylo.read(str(REF),'newick'); out=[]
 for d,q in QUARTETS.items():
  y=Phylo.read(str(REF),'newick')
  for tip in list(y.get_terminals()):
   if tip.name not in q:y.prune(tip)
  ref_top,_=induced(y,q); n=len(exact[d]); ct=Counter(exact[d].values())
  out.append(dict(design=d,source_tree='63K.tre',source_kind='published Stiller 363-taxon reference; genome-wide tree, not chr4-only',n_chr4_loci=len(manifest),n_exact_quartet_trees=n,q1_frequency=ct['q1']/n if n else '',q2_frequency=ct['q2']/n if n else '',q3_frequency=ct['q3']/n if n else '',D_topology_q3_minus_q2=(ct['q3']-ct['q2'])/n if n else '',topology_for_original_order=ref_top,interpretation='Reference tree supplies a species-history hypothesis, not proof of a neutral MSC D=0 null; discordance, introgression, ancestral structure, and finite-sites effects remain possible.'))
 write_tsv(OUT/'quartet_topology_validation.tsv',out,list(out[0]))
 branch_rows=[]
 for d in QUARTETS:
  for top in ('q1','q2','q3'):
   vals=branches[d][top]; vals=sorted(vals)
   med=vals[len(vals)//2] if vals else ''
   branch_rows.append(dict(design=d,topology=top,n_gene_trees=len(vals),median_internal_branch_substitutions_per_site=med,mean_internal_branch_substitutions_per_site=(sum(vals)/len(vals) if vals else ''),coordinate_scope='all exact quartet-bearing chr4 trees; branch length is induced unrooted quartet internal edge'))
 write_tsv(OUT/'quartet_branch_lengths.tsv',branch_rows,list(branch_rows[0]))
 print('reference topology',[(r['design'],r['topology_for_original_order']) for r in out],flush=True)
if __name__=='__main__':main()
