#!/usr/bin/env python3
"""Inventory downloaded Matschiner et al. VCF sample headers and callability."""
from __future__ import annotations
import csv, gzip, re
from pathlib import Path

REPO = Path(__file__).resolve().parents[5]
RAW = REPO / 'data/atlantic_cod/raw/zenodo/local_downloads'
OUT = REPO / 'empirical/atlantic_cod/analysis/patterson_d/results/vcf_sample_inventory.tsv'
FILES = sorted(RAW.glob('*.vcf.gz'))


def population(sample: str) -> str:
    if sample.startswith('Gadmor_'):
        return re.sub(r'\d+$', '', sample)
    return sample

rows=[]
for path in FILES:
    samples=[]; call={}; sites=0; snps=0; multi=0; indel=0
    contigs={}; refmeta=[]; header_reference=''
    with gzip.open(path,'rt') as fh:
        for line in fh:
            if line.startswith('##'):
                if 'gadMor2.fasta' in line: header_reference='gadMor2.fasta (GATK command line)'
                if line.startswith('##reference='): header_reference=line.strip().split('=',1)[1]
                continue
            if line.startswith('#CHROM'):
                samples=line.rstrip('\n').split('\t')[9:]
                call={s:0 for s in samples}
                continue
            v=line.rstrip('\n').split('\t')
            sites+=1; chrom=v[0]; pos=int(v[1]); contigs.setdefault(chrom,[pos,pos]); contigs[chrom][0]=min(contigs[chrom][0],pos);contigs[chrom][1]=max(contigs[chrom][1],pos)
            if len(v[3])==1 and len(v[4])==1 and v[3] in 'ACGT' and v[4] in 'ACGT' and ',' not in v[4]: snps+=1
            elif ',' in v[4]: multi+=1
            else: indel+=1
            for s,g in zip(samples,v[9:]):
                gt=g.split(':',1)[0]
                if '.' not in gt: call[s]+=1
    summary=f"{';'.join(f'{c}:{r[0]}-{r[1]}' for c,r in contigs.items())}"
    for s in samples:
        rows.append(['sample',path.name,header_reference,sites,snps,multi,indel,len(samples),summary,s,population(s),call[s],sites-call[s],f'{call[s]/sites:.6f}' if sites else 'NA'])
    print(f'{path.name}: {sites} records; {snps} biallelic SNPs; {len(samples)} samples; {len(contigs)} contigs',flush=True)
with OUT.open('w') as f:
    w=csv.writer(f,delimiter='\t',lineterminator='\n')
    w.writerow(['row_type','file','reference_header','records','biallelic_snp_records','multiallelic_records','non_snp_records','sample_count','contig_coordinate_ranges','sample_id','population_group','called_records','missing_records','sample_call_rate'])
    w.writerows(rows)
print('wrote',OUT)
