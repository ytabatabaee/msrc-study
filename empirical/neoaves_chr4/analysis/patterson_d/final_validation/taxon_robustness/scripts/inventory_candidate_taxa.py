#!/usr/bin/env python3
"""Inventory taxa in local chr4 alignments without computing site-pattern statistics."""
import csv,re,tarfile
from collections import Counter
from pathlib import Path
ROOT=Path(__file__).resolve().parents[7]
PKG=ROOT/'empirical/neoaves_chr4/analysis/patterson_d/final_validation/taxon_robustness'
ARCH=ROOT/'empirical/neoaves_chr4/external/stiller2024/63430.alns.tar.gz'
MAN=ROOT/'empirical/neoaves_chr4/results/stage4d_locus_manifest.tsv'

def stem(s):
 b=Path(s).name;b=re.sub(r'_fasta\.treefile$','',b);b=re.sub(r'\.(?:fasta|fas|fa|aln)$','',b,flags=re.I);return b.replace('.1k.start','_1k_start')
def main():
 m={}
 with MAN.open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['chromosome']=='chr4':m[stem(r['locus_id'])]=r
 with (PKG/'config/frozen_alternative_quartets.tsv').open() as f:cfg=list(csv.DictReader(f,delimiter='\t'))
 taxa=sorted({c[k] for c in cfg for k in ('P1','P2','P3','O')});counts=Counter();n=0
 with tarfile.open(ARCH,'r:gz') as tf:
  for z in tf:
   if not z.isfile() or '._' in z.name or not z.name.lower().endswith(('.fasta','.fas','.fa')):continue
   key=stem(z.name)
   if key not in m:continue
   n+=1;present=set()
   for line in tf.extractfile(z):
    if line.startswith(b'>'):
     x=line[1:].split()[0].decode('ascii','replace')
     if x in taxa:present.add(x)
   counts.update(present)
   if n%1000==0:print('chr4 alignments inventoried',n,flush=True)
 rows=[]
 for t in taxa:rows.append(dict(taxon=t,scientific_name=t.replace('_',' '),present_in_chr4_alignments=counts[t],usable_chr4_locus_coverage=counts[t],percent_of_6247_chr4_loci=round(100*counts[t]/n,3) if n else '',candidate_set='frozen_alternative_quartets.tsv'))
 with (PKG/'results/candidate_alignment_coverage.tsv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 print('matched',n,'taxa',len(rows),flush=True)
if __name__=='__main__':main()
