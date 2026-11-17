#!/usr/bin/env python3
"""Add independent Jarvis 48-taxon evidence and explicit quartet identities."""
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[6]
P=ROOT/'empirical/neoaves_chr4/analysis/patterson_d/final_validation/results/quartet_topology_validation.tsv'
Q={'Columbea_primary':('Columba_livia','Phoenicopterus_ruber','Ciconia_maguari','Gallus_gallus'),'N61_primary':('Tauraco_erythrolophus','Columba_livia','Ciconia_maguari','Gallus_gallus'),'N62_primary':('Columba_livia','Pterocles_gutturalis','Tauraco_erythrolophus','Gallus_gallus')}
def read():
 with P.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def main():
 base=read(); rows=[]
 for r in base:
  d=r['design'];q=Q[d]
  r.update(taxon_P1=q[0],taxon_P2=q[1],taxon_P3=q[2],outgroup=q[3],induced_topology=r['topology_for_original_order'],support='published tree label 1; not a posterior probability',independent_evidence='Stiller 2024 published 363-taxon genome-wide reference tree; contains all four labels. Not chr4-only, but not independent of the underlying Stiller data.',uncertainty='Single resolved reference cannot represent species-tree uncertainty, ILS, introgression/ancestral structure, within-species histories, or finite-sites bias.',interpretability='D=0 null is conditional on this topology and MSC assumptions; not established by this tree alone.')
  rows.append(r)
 for d,q in Q.items():
  if d=='Columbea_primary':
   top='exact quartet unavailable; broad implication q1 only if Ciconia is in Passerea'
   ev='Jarvis et al. 2014 TENT supports Columbea and describes pigeons and flamingos in sister clades; Ciconia maguari is absent from the 48-taxon panel, so exact four-tip topology is unavailable.'
  elif d=='N61_primary':
   top='not exactly evaluable; Ciconia maguari absent'
   ev='Jarvis et al. 2014 48-taxon TENT supplies broad Otidimorphae/Columbimorphae/Passerea placements, but Ciconia maguari is absent; the quartet cannot be induced from the 48 sampled tips.'
  else:
   top='broad implication q1; exact induced quartet not independently recalculated'
   ev='Jarvis et al. 2014 TENT supports broad Columbea membership for Columbimorphae (pigeon/sandgrouse) and Passerea for Tauraco; an exact Newick tree is not available in the local source package for recomputation.'
  rows.append(dict(design=d,source_tree='Jarvis et al. 2014 TENT; published 48-taxon tree (Fig. 1/S1)',source_kind='independent historical genome-scale species-tree analysis; Newick not locally available',n_chr4_loci='',n_exact_quartet_trees='',q1_frequency='',q2_frequency='',q3_frequency='',D_topology_q3_minus_q2='',topology_for_original_order=top,interpretation='Broad clade evidence only; do not treat as exact quartet validation.',taxon_P1=q[0],taxon_P2=q[1],taxon_P3=q[2],outgroup=q[3],induced_topology=top,support='published clade/tree evidence; exact quartet support unavailable',independent_evidence=ev,uncertainty='Rapid early Neoaves radiation and disagreement among genome-scale backbones; exact four-species sampling incomplete for Ciconia.',interpretability='Qualitative clade evidence; insufficient alone to certify a quartet-specific D=0 null.'))
 fields=[]
 for r in rows:
  for k in r:
   if k not in fields:fields.append(k)
 with P.open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
 print('wrote 3 Stiller exact rows plus 3 Jarvis broad-evidence rows')
if __name__=='__main__':main()
