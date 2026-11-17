#!/usr/bin/env python3
"""Recompute matched and full-span delete-block jackknife D contrasts."""
from pathlib import Path
import csv, math
P=Path(__file__).resolve().parents[2]; OUT=P/'final_validation/results'
def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rows):
 with open(p,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def pat(A,B):return (A-B)/(A+B) if A+B else float('nan')
def calc(ins,out):return pat(ins[0],ins[1])-pat(out[0],out[1])
inv=[(5747331,6957919),(7435308,13016510),(14068229,18051507),(19401000,40049964)]
full=read(P/'polarization_validation/results/corrected_window_D.tsv')
matched=read(P/'exact_quartet_comparison/results/matched_window_dataset.tsv')
allrows=[]
for design in sorted(set(x['design'] for x in full)):
 base={int(x['start']):(float(x['ABBA']),float(x['BABA'])) for x in full if x['design']==design and int(x['start'])<94_987_271}
 tr={int(x['start']):x for x in matched if x['design']==design}
 for scope in ('full_contig','exact_tree_matched_span'):
  obs=base if scope=='full_contig' else {s:v for s,v in base.items() if s in tr}
  for bp in (250000,500000,1000000):
   block={};inside=[0.,0.];outside=[0.,0.]
   for st,(a,b) in obs.items():
    if scope=='exact_tree_matched_span': isin=tr[st]['recombination_inversion_class'].startswith('In(17)')
    else:isin=any(st < end and st+5000 > start for start,end in inv)
    zone='inside' if isin else 'outside';key=st//bp
    q=block.setdefault(key,{'inside':[0.,0.],'outside':[0.,0.]})[zone];q[0]+=a;q[1]+=b
    dest=inside if isin else outside;dest[0]+=a;dest[1]+=b
   theta=calc(inside,outside);keys=list(block);jack=[]
   for k in keys:
    i=[inside[0]-block[k]['inside'][0],inside[1]-block[k]['inside'][1]]
    o=[outside[0]-block[k]['outside'][0],outside[1]-block[k]['outside'][1]]
    jack.append(calc(i,o))
   mean=sum(jack)/len(jack);se=math.sqrt((len(jack)-1)/len(jack)*sum((x-mean)**2 for x in jack))
   allrows.append({'design':design,'scope':scope,'block_bp':bp,'estimate_inside_minus_outside_D':theta,'jackknife_SE':se,'normal_95_low':theta-1.96*se,'normal_95_high':theta+1.96*se,'inside_ABBA':inside[0],'inside_BABA':inside[1],'outside_ABBA':outside[0],'outside_BABA':outside[1],'delete_blocks':len(keys),'note':'Delete physical genomic blocks; D is recomputed after deletion; exact_tree_matched_span restricts to shared tree/SNP start coordinates'})
write(OUT/'block_uncertainty_validation.tsv',allrows)
print('wrote',len(allrows),'block jackknife summaries')
