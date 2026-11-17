#!/usr/bin/env python3
"""Summarize optional scan only at exact tree/SNP matched 5-kb starts."""
from pathlib import Path
import csv
P=Path(__file__).resolve().parents[2]
OUT=P/'final_validation/results'
def read(path):return list(csv.DictReader(open(path),delimiter='\t'))
def write(path,rows):
 with open(path,'w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
raw=read(P/'polarization_validation/results/robustness_window_counts.tsv')
trees=read(P/'exact_quartet_comparison/results/matched_window_dataset.tsv')
labels={(r['design'],int(r['start'])):r['recombination_inversion_class'] for r in trees}
regions=['all_exact_matched_tree_windows','inversion_union','matched_tree_complement','In(17)1','In(17)2','In(17)3','In(17)4','left_flank','right_flank','collinear_gap_1','collinear_gap_2','collinear_gap_3']
out=[]
for design in sorted(set(r['design'] for r in raw)):
 for threshold in (.8,.9,1.0):
  for filt in ('all_biallelic_SNPs','transversions_only'):
   rr=[r for r in raw if r['design']==design and float(r['callrate_threshold'])==threshold and r['filter']==filt and (design,int(r['start'])) in labels]
   for region in regions:
    z=[]
    for r in rr:
     label=labels[(design,int(r['start']))]
     inside=label.startswith('In(17)')
     take=(region=='all_exact_matched_tree_windows' or
       region=='inversion_union' and inside or
       region=='matched_tree_complement' and not inside or
       region==label)
     if take:z.append(r)
    a=sum(float(r['ABBA']) for r in z);b=sum(float(r['BABA']) for r in z);n=sum(int(r['informative_sites']) for r in z)
    out.append({'design':design,'callrate_threshold':threshold,'filter':filt,'region':region,'matched_5kb_windows':len(z),'ABBA':a,'BABA':b,'D':(a-b)/(a+b) if a+b else '', 'informative_sites':n,'coverage_note':'Only exact 5-kb start matches present in matched_window_dataset.tsv; labels use frozen published intervals'})
write(OUT/'robustness_exact_tree_matched.tsv',out)
# Append clearly identified exact-match sensitivity rows to final table (idempotent).
fpath=OUT/'final_results.tsv'; final=read(fpath)
final=[r for r in final if r.get('estimand')!='robustness_exact_tree_matched_D_site']
for r in out:
 if r['region'] in ('all_exact_matched_tree_windows','inversion_union','matched_tree_complement'):
  final.append({'design':r['design'],'estimand':'robustness_exact_tree_matched_D_site','region':r['region'],'resolution_bp':'5kb exact tree windows','value':r['D'],'comparison_value':'','difference':'','ABBA':r['ABBA'],'BABA':r['BABA'],'informative_sites':r['informative_sites'],'status':f"{r['filter']};callrate>={r['callrate_threshold']};complete",'note':'Only exact matched tree starts; not all windows in physical span'})
write(fpath,final)
