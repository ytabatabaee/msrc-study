#!/usr/bin/env python3
"""Freeze additional cod D designs using baseline/arrangement metadata only."""
import csv, gzip
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[5]
AN=ROOT/'empirical/atlantic_cod/analysis/patterson_d'
DATA=ROOT/'data/atlantic_cod/processed'
RAW=ROOT/'data/atlantic_cod/raw/zenodo/local_downloads'
OUT=AN/'validation/config'
LGS=('LG01','LG02','LG07','LG12')

def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rs):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open('w') as f:
        w=csv.DictWriter(f,list(rs[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rs)
def vcf_samples(lg):
    p=RAW/f'gadus_morhua_supergenes_{lg.lower()}_snapp.vcf.gz'
    with gzip.open(str(p),'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):return line.rstrip('\n').split('\t')[9:]
    raise RuntimeError('VCF header missing: '+str(p))

def main():
    arrangements={r['tree_label'][:-4]:r for r in rows(DATA/'population_arrangements.tsv')}
    primary={r['lg']:r for r in rows(AN/'config/frozen_quartets.tsv')}
    stage3=rows(DATA/'stage3_baseline_quartets.tsv')
    frozen=[];membership=[]
    for lg in LGS:
        samples=vcf_samples(lg);p0=primary[lg]
        primary_pops=[p0['P1'],p0['P2'],p0['P3']]
        primary_states=[arrangements[p][lg+'_state'] for p in primary_pops]
        candidates={'same_arrangement_q3':[],'alternative_arrangement_q2':[],'baseline_pair_control':[]}
        for r in stage3:
            if r['lg']!=lg:continue
            pairs=[x.split(',') for x in r['baseline_split'].split('|')]
            for pp,other in ((pairs[0],pairs[1]),(pairs[1],pairs[0])):
                for p3tree in other:
                    names=[pp[0][:-4],pp[1][:-4],p3tree[:-4]]
                    if len(set(names))!=3 or set(names)==set(primary_pops):continue
                    ids={p:sorted(s for s in samples if s.startswith(p)) for p in names}
                    if any(len(ids[p])<2 for p in names):continue
                    st=[arrangements[p][lg+'_state'] for p in names]
                    if st[1]==st[2] and st[0]!=st[1]:category='same_arrangement_q3'
                    elif st[0]==st[2] and st[1]!=st[0]:category='alternative_arrangement_q2'
                    elif st[0]==st[1] and st[2]!=st[0]:category='baseline_pair_control'
                    else:continue
                    reuse=len(set(names)&set(primary_pops));pattern=sum(a==b for a,b in zip(st,primary_states))
                    source_order=int(r['quartet_id'].split('SQ')[-1])
                    eco_same=arrangements[names[0]]['ecotype']==arrangements[names[2]]['ecotype']
                    rank=(reuse,-pattern,-int(eco_same),source_order,tuple(names))
                    candidates[category].append((rank,names,ids,r,st))
        for c in candidates:candidates[c].sort(key=lambda z:z[0])
        classes=[
          ('same_arrangement_q3','ALTQ3','Alternative population triple with P2/P3 sharing arrangement, as in primary design, and an outside-inversion baseline P1/P2 sister pair.'),
          ('alternative_arrangement_q2','ALTQ2','Alternative arrangement prediction: P1/P3 share arrangement, mapping to BABA (q2) with this frozen ordering.'),
          ('baseline_pair_control','BASELINE_CONTROL','P1/P2 baseline sisters share arrangement while P3 differs; arrangement matches the baseline pair.'),
        ]
        for category,tag,rationale in classes:
            if not candidates[category]:raise RuntimeError('No candidate for %s %s'%(lg,category))
            _,names,ids,r,st=candidates[category][0]
            row={'lg':lg,'design_id':lg+'_VALID_'+tag,'design_class':category,'P1':names[0],'P2':names[1],'P3':names[2],
                 'O':'Borsai','O_sample_id':'Borsai','outgroup_species':'Boreogadus saida',
                 'P1_sample_ids':','.join(ids[names[0]][:2]),'P2_sample_ids':','.join(ids[names[1]][:2]),'P3_sample_ids':','.join(ids[names[2]][:2]),
                 'O_sample_ids':'Borsai',
                 'n_P1':2,'n_P2':2,'n_P3':2,'n_O':1,
                 'P1_arrangement':st[0],'P2_arrangement':st[1],'P3_arrangement':st[2],
                 'arrangement_state_pattern':''.join('A' if s=='ancestral' else 'D' for s in st),
                 'stage3_candidate_id':r['quartet_id'],'baseline_split':r['baseline_split'],'arrangement_split':r['arrangement_split'],
                 'population_locations':';'.join(arrangements[n]['location'] for n in names),
                 'population_ecotypes':';'.join(arrangements[n]['ecotype'] for n in names),
                 'rationale_and_confounders':rationale+' Geography/ecotype remain correlated with arrangement and ancestry; this is not a randomized control.'}
            frozen.append(row)
            for role in ('P1','P2','P3','O'):
                for sample in row[role+'_sample_ids'].split(','):
                    membership.append({'lg':lg,'design_id':row['design_id'],'role':role,'population':row[role],
                                       'sample_id':sample,'sample_count_for_role':row['n_'+role]})
    write(OUT/'frozen_validation_quartets.tsv',frozen)
    write(OUT/'frozen_validation_sample_membership.tsv',membership)
    print('Frozen',len(frozen),'additional rooted quartet designs from Stage-3 splits and arrangement metadata; no D values were read.')
    for r in frozen:print(r['design_id'],r['P1'],r['P2'],r['P3'],r['arrangement_state_pattern'],r['stage3_candidate_id'])

if __name__=='__main__':main()
