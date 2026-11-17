#!/usr/bin/env python3
"""Compute prespecified alternative N61 quartet D and raw-locus validation."""
import csv, hashlib, math, random, re, tarfile, sys
from collections import defaultdict
from pathlib import Path


ROOT=Path(__file__).resolve().parents[7]
BASE=ROOT/'empirical/neoaves_chr4/analysis/patterson_d'
PKG=BASE/'final_validation/taxon_robustness'
ARCH=ROOT/'empirical/neoaves_chr4/external/stiller2024/63430.alns.tar.gz'
MAN=ROOT/'empirical/neoaves_chr4/results/stage4d_locus_manifest.tsv'
FROZEN=PKG/'config/frozen_alternative_quartets.tsv'
OLD=BASE/'results/patterson_d_locus_results.tsv'
sys.path.insert(0,str(BASE/'scripts'))
from patterson_d import count_sites  # reuse the validated production site counter for new quartet estimates
INTERVALS=[('PNAS_1',25030000,32670000),('PNAS_2',33510000,34470000),('PNAS_3',44130000,56810000)]
BASES=set('ACGT'); TRANS={frozenset(('A','G')),frozenset(('C','T'))}

def read_tsv(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        if not rows:return
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def stem(name):
    b=Path(name).name
    b=re.sub(r'_fasta\.treefile$','',b,flags=re.I)
    b=re.sub(r'\.(?:fasta|fas|fa|aln)$','',b,flags=re.I)
    return b.replace('.1k.start','_1k_start')
def parse_fasta(handle):
    seqs={};duplicates=[];name=None;chunks=[]
    for raw in handle:
        line=raw.decode('ascii','replace') if isinstance(raw,bytes) else raw
        line=line.strip()
        if not line:continue
        if line.startswith('>'):
            if name is not None:
                if name in seqs:duplicates.append(name)
                seqs[name]=''.join(chunks).upper()
            name=line[1:].split()[0];chunks=[]
        else:chunks.append(''.join(line.split()))
    if name is not None:
        if name in seqs:duplicates.append(name)
        seqs[name]=''.join(chunks).upper()
    return seqs,duplicates

def count_independent(seq, tv_only=False):
    a=b=amb=multi=other=0;length=len(seq['P1'])
    for i in range(length):
        vals=[seq[k][i] for k in ('P1','P2','P3','O')]
        if any(x not in BASES for x in vals):amb+=1;continue
        p1,p2,p3,o=vals
        if p1==o and p2==p3 and p2!=o:pat='ABBA'
        elif p2==o and p1==p3 and p1!=o:pat='BABA'
        else:
            states=set(vals)
            if len(states)>2:multi+=1
            else:other+=1
            continue
        if tv_only and frozenset((o,p3)) in TRANS:continue
        if pat=='ABBA':a+=1
        else:b+=1
    return dict(ABBA=a,BABA=b,informative_sites=a+b,ambiguous_sites=amb,multistate_nonpattern_sites=multi,other_biallelic_sites=other,alignment_length=length)

def d(a,b):return (a-b)/(a+b) if a+b else None
def jk(vals, bs, pred_inside):
    # vals tuples: start, ABBA, BABA, region name. Delete the physical block in
    # both inside/outside groups, preserving one matched jackknife replicate.
    g=defaultdict(lambda:[0,0,0,0])
    for st,a,b,reg in vals:
        v=g[st//bs]; idx=0 if pred_inside(reg) else 2;v[idx]+=a;v[idx+1]+=b
    ia=sum(v[0] for v in g.values());ib=sum(v[1] for v in g.values());oa=sum(v[2] for v in g.values());ob=sum(v[3] for v in g.values())
    di=d(ia,ib);do=d(oa,ob);delta=di-do if di is not None and do is not None else None
    reps={'inside':[],'outside':[],'delta':[]}
    for x,y,z,w in g.values():
        xi=d(ia-x,ib-y);yo=d(oa-z,ob-w)
        if xi is not None:reps['inside'].append(xi)
        if yo is not None:reps['outside'].append(yo)
        if xi is not None and yo is not None:reps['delta'].append(xi-yo)
    def se(x):
        if len(x)<2:return None
        m=sum(x)/len(x);return math.sqrt((len(x)-1)/len(x)*sum((v-m)**2 for v in x))
    counts={'n_inside_loci':sum(pred_inside(r) for _,_,_,r in vals),'n_outside_loci':sum(not pred_inside(r) for _,_,_,r in vals),
            'inside_informative_sites':ia+ib,'outside_informative_sites':oa+ob,'inside_ABBA':ia,'inside_BABA':ib,'outside_ABBA':oa,'outside_BABA':ob,'n_blocks':len(g)}
    out=[]
    for scope,theta in [('inside',di),('outside',do),('delta',delta)]:
        s=se(reps[scope]);out.append((theta,s,theta-1.96*s if s is not None and theta is not None else None,theta+1.96*s if s is not None and theta is not None else None))
    return out,counts

def region(st):
    for n,a,b in INTERVALS:
        if a<=st<=b:return n
    return 'outside_all_PNAS'

def sample_validation(oldrows):
    rng=random.Random(20261118)
    base=[r for r in oldrows if r['design']=='N61_primary' and r['filter']=='all']
    info=sorted(int(r['informative']) for r in base)
    lo=info[max(0,int(.2*len(info))-1)];hi=info[min(len(info)-1,int(.8*len(info)))]
    def pick(rows,n):
        rows=list(rows);rng.shuffle(rows);return rows[:min(n,len(rows))]
    chosen={}
    for name,_,_ in INTERVALS:
        here=[r for r in base if r['region']==name]
        groups=[pick([r for r in here if int(r['informative'])>=hi and int(r['excluded'])>0],2),
                pick([r for r in here if int(r['informative'])<=lo],2),
                pick([r for r in here if int(r['excluded'])>0],2),pick(here,4)]
        for g in groups:
            for r in g:chosen[r['locus_id']]=set(chosen.get(r['locus_id'],set()))|{name}
    outside=[r for r in base if r['region']=='outside_all_PNAS']
    for label,subset in [('outside_high_info',[r for r in outside if int(r['informative'])>=hi]),
                         ('outside_low_info',[r for r in outside if int(r['informative'])<=lo]),
                         ('outside_ambiguity',[r for r in outside if int(r['excluded'])>0])]:
        for r in pick(subset,6):chosen[r['locus_id']]=set(chosen.get(r['locus_id'],set()))|{label}
    return chosen,lo,hi

def main():
    designs=read_tsv(FROZEN); old=read_tsv(OLD)
    oldmap={(r['design'],r['locus_id'],r['filter']):r for r in old}
    chosen,lo,hi=sample_validation(old)
    with (PKG/'results/raw_validation_sample_manifest.tsv').open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['locus_id','strata','seed','N61_informative_lower_quintile_cut','N61_informative_upper_quintile_cut'])
        for locus,labels in sorted(chosen.items()):w.writerow([locus,','.join(sorted(labels)),20261118,lo,hi])
    manifest={}
    for r in read_tsv(MAN):
        if r['chromosome']=='chr4':manifest[stem(r['locus_id'])]=r
    records=[];raw=[];seen=set();nmatch=0;noarch=0;archive_dupes=0
    with tarfile.open(ARCH,'r:gz') as tf:
      for m in tf:
        if not m.isfile() or '._' in m.name or not m.name.lower().endswith(('.fasta','.fas','.fa')):continue
        key=stem(m.name)
        if key not in manifest:continue
        if key in seen:archive_dupes+=1;continue
        seen.add(key);nmatch+=1;mr=manifest[key];loc=mr['locus_id'];st=int(mr['start']);reg=region(st)
        seqs,dups=parse_fasta(tf.extractfile(m));lens={len(v) for v in seqs.values()}
        alignment_equal=len(lens)<=1
        for cfg in designs:
            taxa=[cfg[k] for k in ('P1','P2','P3','O')]
            if not all(t in seqs for t in taxa):continue
            s={k:seqs[cfg[k]] for k in ('P1','P2','P3','O')}
            if len({len(v) for v in s.values()})!=1:continue
            a,b,ex,inf=count_sites(s,False);ta,tb,tex,tinf=count_sites(s,True)
            records.append(dict(design=cfg['design'],locus_id=loc,start=st,end=mr['end'],midpoint=mr['midpoint'],region=reg,region_class='inside' if reg!='outside_all_PNAS' else 'outside',ABBA=a,BABA=b,D=d(a,b) if a+b else '',informative=a+b,ambiguous_or_gap_sites=ex,transversion_ABBA=ta,transversion_BABA=tb,transversion_D=d(ta,tb) if ta+tb else '',transversion_informative=ta+tb,transition_ABBA=a-ta,transition_BABA=b-tb,alignment_length=len(s['P1'])))
        if loc in chosen:
            if dups:archive_dupes+=0
            for cfg in read_tsv(BASE/'config/frozen_quartets.tsv'):
                taxa=[cfg[k] for k in ('P1','P2','P3','O')]
                if not all(t in seqs for t in taxa):
                    raw.append(dict(locus_id=loc,design=cfg['design'],strata=','.join(sorted(chosen[loc])),presence='missing_taxon',missing_taxa=','.join(t for t in taxa if t not in seqs),duplicate_sequence_ids=';'.join(sorted(set(dups))),all_taxa_alignment_lengths_equal=alignment_equal,independent_ABBA='',independent_BABA='',saved_ABBA='',saved_BABA='',ABBA_match='NA',BABA_match='NA',status='UNAVAILABLE'));continue
                s={k:seqs[cfg[k]] for k in ('P1','P2','P3','O')}
                same4=len({len(v) for v in s.values()})==1
                if not same4:
                    raw.append(dict(locus_id=loc,design=cfg['design'],strata=','.join(sorted(chosen[loc])),presence='present',missing_taxa='',duplicate_sequence_ids=';'.join(sorted(set(dups))),all_taxa_alignment_lengths_equal=alignment_equal,independent_ABBA='',independent_BABA='',saved_ABBA='',saved_BABA='',ABBA_match='NA',BABA_match='NA',status='ALIGNMENT_LENGTH_MISMATCH'));continue
                one=count_independent(s,False);tv=count_independent(s,True)
                oldall=oldmap.get((cfg['design'],loc,'all'));oldtv=oldmap.get((cfg['design'],loc,'transversions_only'))
                for filt,got,saved in [('all',one,oldall),('transversions_only',tv,oldtv)]:
                    sa=int(saved['ABBA']) if saved else None;sb=int(saved['BABA']) if saved else None
                    am=(got['ABBA']==sa) if sa is not None else None;bm=(got['BABA']==sb) if sb is not None else None
                    raw.append(dict(locus_id=loc,design=cfg['design'],filter=filt,strata=','.join(sorted(chosen[loc])),presence='present',missing_taxa='',duplicate_sequence_ids=';'.join(sorted(set(dups))),all_taxa_alignment_lengths_equal=alignment_equal,alignment_length=one['alignment_length'],independent_ABBA=got['ABBA'],independent_BABA=got['BABA'],independent_informative=got['informative_sites'],independent_ambiguous=got['ambiguous_sites'],independent_multistate=got['multistate_nonpattern_sites'],saved_ABBA=sa if sa is not None else '',saved_BABA=sb if sb is not None else '',ABBA_match=am if am is not None else 'NA',BABA_match=bm if bm is not None else 'NA',status='PASS' if am is True and bm is True and not dups and alignment_equal else 'FAIL' if am is False or bm is False or dups or not alignment_equal else 'NO_SAVED_ROW'))
        if nmatch%1000==0:print('matched chr4 alignments',nmatch,flush=True)
    if nmatch!=len(manifest):raise RuntimeError(f'archive mapping mismatch: saw {nmatch} expected {len(manifest)} chr4 loci')
    write_tsv(PKG/'results/alternative_quartet_locus_counts.tsv',records)
    write_tsv(PKG/'results/raw_alignment_validation.tsv',raw)
    # Pooled inside/outside and named interval summaries, preserving count pooling.
    sums=[]
    for cfg in designs:
      rr=[r for r in records if r['design']==cfg['design']]
      for filt,ak,bk in [('all','ABBA','BABA'),('transversions_only','transversion_ABBA','transversion_BABA')]:
       for scope in ['PNAS_1','PNAS_2','PNAS_3','all_PNAS','outside_all_PNAS']:
        pred=lambda rg,s=scope: rg!='outside_all_PNAS' if s=='all_PNAS' else rg==s if s.startswith('PNAS_') else rg=='outside_all_PNAS'
        val=[(int(r['start']),int(r[ak]),int(r[bk]),r['region']) for r in rr]
        if scope=='outside_all_PNAS':
          aa=sum(a for _,a,b,reg in val if pred(reg));bb=sum(b for _,a,b,reg in val if pred(reg));nl=sum(pred(reg) for _,a,b,reg in val);inf=aa+bb;dd=d(aa,bb)
          sums.append(dict(design=cfg['design'],filter=filt,scope=scope,inside_D='',outside_D=dd,delta_D='',inside_SE='',outside_SE='',delta_SE='',inside_CI_low='',inside_CI_high='',outside_CI_low='',outside_CI_high='',delta_CI_low='',delta_CI_high='',inside_loci='',outside_loci=nl,inside_informative_sites='',outside_informative_sites=inf,inside_ABBA='',inside_BABA='',outside_ABBA=aa,outside_BABA=bb,block_size_bp='NA',contrast_status='outside_descriptor_only'));continue
        for bs in (500000,1000000,2000000):
          # For each named interval compare that interval with all other chromosome-4 loci.
          is_inside=lambda rg,s=scope: (rg!='outside_all_PNAS') if s=='all_PNAS' else (rg==s)
          stats,ct=jk(val,bs,is_inside)
          (di,si,li,ui),(do,so,lo_ci,up_o),(delta,sd,ld,ud)=stats
          sums.append(dict(design=cfg['design'],filter=filt,scope=scope,inside_D=di,outside_D=do,delta_D=delta,inside_SE=si,outside_SE=so,delta_SE=sd,inside_CI_low=li,inside_CI_high=ui,outside_CI_low=lo_ci,outside_CI_high=up_o,delta_CI_low=ld,delta_CI_high=ud,inside_loci=ct['n_inside_loci'],outside_loci=ct['n_outside_loci'],inside_informative_sites=ct['inside_informative_sites'],outside_informative_sites=ct['outside_informative_sites'],inside_ABBA=ct['inside_ABBA'],inside_BABA=ct['inside_BABA'],outside_ABBA=ct['outside_ABBA'],outside_BABA=ct['outside_BABA'],block_size_bp=bs,contrast_status='defined' if delta is not None else 'undefined'))
    write_tsv(PKG/'results/taxon_robustness_block_summary.tsv',sums)
    # Report scan details and raw validation verdicts.
    passed=sum(r.get('status')=='PASS' for r in raw);failed=sum(str(r.get('status','')).startswith('FAIL') for r in raw);unavail=sum(r.get('status')=='UNAVAILABLE' for r in raw)
    report=['# Independent raw-alignment validation','',f'Deterministic stratified sample seed: 20261118. Unique loci sampled: {len(chosen)} (inside PNAS_1/2/3, outside high-information, outside low-information, and outside with ambiguity/missing calls). N61-primary informative-count lower/upper quintile cutoffs were {lo}/{hi}.','',f'The checksum-verified local archive was streamed once; {nmatch} chr4 alignments matched the manifest. Alternative quartet per-locus all-site and transversion counts were computed in this pass. Raw validation covered the three original primary and their three predeclared Meleagris-outgroup designs; each independently decoded row was compared to the saved per-locus all-site and transversion ABBA/BABA integers.','',f'Validation rows PASS: {passed}; discrepancy/format failures: {failed}; unavailable quartet/locus rows: {unavail}. Duplicated chr4 archive keys: {archive_dupes}.','', 'The reference implementation uses a separate FASTA parser and direct four-base pattern comparisons; it does not import or invoke the production estimator. It records duplicate IDs, full-alignment length equality, four-taxon alignment length, ambiguous/gap site count, multistate non-pattern sites, and exact ABBA/BABA counts. Outgroup polarization is the observed O base; only patterns with one derived P3/O state and the exact ABBA or BABA configuration count. Transversions-only excludes A↔G and C↔T differences between P3 and O.','', '## Interpretation','', 'Exact integer agreement is required for both ABBA and BABA, separately for all informative patterns and transversions-only. The per-row comparison is in `raw_alignment_validation.tsv`; the deterministic sample list and strata are in `raw_validation_sample_manifest.tsv`. No production filters were modified to force agreement.','']
    (PKG/'results/raw_alignment_validation_report.md').write_text('\n'.join(report))
    print('n_records',len(records),'raw_rows',len(raw),'PASS',passed,'FAIL',failed,'UNAVAILABLE',unavail,flush=True)
if __name__=='__main__':main()
