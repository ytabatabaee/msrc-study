"""Core helpers for the cod nucleotide Patterson D analysis."""
from __future__ import annotations
import csv, gzip, math, re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

WINDOW = 250_000
TRANSITIONS = {('A','G'),('G','A'),('C','T'),('T','C')}


def parse_gt(gt: str):
    """Return a diploid biallelic genotype as an allele tuple, else None."""
    gt = gt.split(':', 1)[0]
    sep = '|' if '|' in gt else ('/' if '/' in gt else None)
    if sep is None:
        return None
    vals = gt.split(sep)
    if len(vals) != 2 or any(x not in ('0','1') for x in vals):
        return None
    return tuple(int(x) for x in vals)


def alt_frequency(genotypes: Iterable[str], min_called_individuals: int = 1):
    called=0; alt=0; alleles=0
    for gt in genotypes:
        g=parse_gt(gt)
        if g is None: continue
        called += 1; alt += sum(g); alleles += len(g)
    if called < min_called_individuals or alleles == 0:
        return None, called
    return alt / alleles, called


def fixed_outgroup_state(gt: str):
    g=parse_gt(gt)
    if g is None: return None
    if g[0] == g[1]: return g[0]
    return None


def site_contributions(p1_alt: float, p2_alt: float, p3_alt: float, outgroup_allele: int):
    """Return ABBA/BABA frequency contributions after outgroup polarization."""
    if outgroup_allele not in (0,1):
        return None
    if outgroup_allele == 0:  # REF ancestral, ALT derived
        p1,p2,p3=p1_alt,p2_alt,p3_alt
    else:                     # ALT ancestral, REF derived
        p1,p2,p3=1-p1_alt,1-p2_alt,1-p3_alt
    return (1-p1)*p2*p3, p1*(1-p2)*p3


def d_from_sums(abba: float, baba: float):
    den=abba+baba
    return (abba-baba)/den if den > 0 else None


def window_id(lg: str, pos: int, size: int = WINDOW):
    start=((pos-1)//size)*size+1
    return f'{lg}_{start:09d}_{start+size-1:09d}',start,start+size-1


def pop_group(sample: str):
    return re.sub(r'\d+$','',sample) if sample.startswith('Gadmor_') else sample


def read_designs(path: Path):
    return list(csv.DictReader(path.open(), delimiter='\t'))


def read_tree_windows(path: Path):
    out=defaultdict(dict)
    with path.open() as f:
        for r in csv.DictReader(f,delimiter='\t'):
            out[r['lg']][r['window_id']]={**r,'start':int(r['start']),'end':int(r['end']),'midpoint':int(r['midpoint'])}
    return out


def scan_vcf(path: Path, designs: dict, tree_windows: dict, outgroups=None,
             min_call_fraction: float=1.0, sample_mode: str='frozen',
             transversions_only: bool=False, keep_sites: bool=True):
    """Scan a VCF for designs on selected LGs. Each design has frozen sample IDs."""
    outgroups=outgroups or {lg:d['O_sample_id'] for lg,d in designs.items()}
    results={}
    with gzip.open(path,'rt') as fh:
        samples=[]; sample_index={}; format_gt_i=0
        for line in fh:
            if line.startswith('##'): continue
            if line.startswith('#CHROM'):
                samples=line.rstrip('\n').split('\t')[9:]
                sample_index={s:i for i,s in enumerate(samples)}
                continue
            v=line.rstrip('\n').split('\t')
            lg=v[0]
            if lg not in designs: continue
            pos=int(v[1]); wid,wstart,wend=window_id(lg,pos)
            if wid not in tree_windows.get(lg,{}): continue
            d=designs[lg]; key=d['design_id']; r=results.setdefault(key,{'lg':lg,'design_id':key,'source':path.name,'outgroup':outgroups[lg],
                'records_in_tree_windows':0,'biallelic_snps':0,'multi_or_other_alt':0,'non_snp_records':0,'nonpass_records':0,
                'ingroup_callability_exclusions':0,'outgroup_missing':0,'outgroup_heterozygous':0,'outgroup_hom_ref':0,'outgroup_hom_alt':0,
                'polarized_sites':0,'informative_sites':0,'sum_abba':0.0,'sum_baba':0.0,'windows':defaultdict(lambda:defaultdict(float)),'sites':[]})
            ws=r['windows'][wid]
            r['records_in_tree_windows']+=1;ws['records']+=1
            ref,alt=v[3],v[4]
            if v[6] not in ('PASS','.'): r['nonpass_records']+=1;ws['nonpass']+=1;continue
            if ',' in alt: r['multi_or_other_alt']+=1;ws['multiallelic']+=1;continue
            if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT':
                r['non_snp_records']+=1;ws['non_snp']+=1;continue
            if transversions_only and (ref,alt) in TRANSITIONS: ws['transition_filtered']+=1;continue
            r['biallelic_snps']+=1;ws['biallelic_snps']+=1
            format_fields=v[8].split(':')
            try: gt_i=format_fields.index('GT')
            except ValueError: r['ingroup_callability_exclusions']+=1;ws['missingness_exclusions']+=1;continue
            def getgt(sample):
                ix=sample_index.get(sample)
                if ix is None: return None
                vals=v[9+ix].split(':')
                return vals[gt_i] if gt_i < len(vals) else None
            def choose_samples(pop, role):
                if sample_mode=='all_available':
                    return [s for s in samples if pop_group(s)==pop]
                if sample_mode=='first_available':
                    arr=sorted(s for s in samples if pop_group(s)==pop)
                    return arr[:1]
                frozen=d[f'{role}_sample_ids'].split(',')
                return [s for s in frozen if s in sample_index]
            popvals=[]; enough=True; missing_cells=0
            for role in ('P1','P2','P3'):
                pop=d[role]
                ids=choose_samples(pop,role)
                if not ids:
                    enough=False; popvals.append((None,0));continue
                min_called=max(1,math.ceil(min_call_fraction*len(ids)-1e-12))
                freq,ncall=alt_frequency((getgt(s) or '.' for s in ids),min_called)
                missing_cells += len(ids)-ncall
                if freq is None: enough=False
                popvals.append((freq,ncall))
            outid=outgroups[lg]; og_gt=getgt(outid)
            og_state=fixed_outgroup_state(og_gt or '.')
            if og_gt is None or '.' in og_gt: r['outgroup_missing']+=1;ws['outgroup_missing']+=1
            elif og_state is None: r['outgroup_heterozygous']+=1;ws['outgroup_heterozygous']+=1
            elif og_state==0: r['outgroup_hom_ref']+=1;ws['outgroup_hom_ref']+=1
            else: r['outgroup_hom_alt']+=1;ws['outgroup_hom_alt']+=1
            if not enough or og_state is None:
                r['ingroup_callability_exclusions']+=1;ws['missingness_exclusions']+=1
                if keep_sites: r['sites'].append((pos,wid,0.0,0.0,False,False,og_state))
                continue
            r['polarized_sites']+=1;ws['polarized_sites']+=1
            contrib=site_contributions(popvals[0][0],popvals[1][0],popvals[2][0],og_state)
            abba,baba=contrib
            r['sum_abba']+=abba;r['sum_baba']+=baba
            ws['sum_abba']+=abba;ws['sum_baba']+=baba
            info=(abba+baba)>0
            if info:r['informative_sites']+=1;ws['informative_sites']+=1
            if keep_sites:r['sites'].append((pos,wid,abba,baba,True,info,og_state))
    for r in results.values():
        r['D']=d_from_sums(r['sum_abba'],r['sum_baba'])
        for w in r['windows'].values(): w['D']=d_from_sums(w.get('sum_abba',0.0),w.get('sum_baba',0.0))
    return results


def add_window_stats(target, source):
    for wid,s in source.items():
        t=target[wid]
        for k,v in s.items():
            if k=='D': continue
            t[k]=t.get(k,0)+v
    for s in target.values(): s['D']=d_from_sums(s.get('sum_abba',0.0),s.get('sum_baba',0.0))
    return target


def pool_windows(windows: dict, selected_ids: set[str]):
    vals=[windows.get(wid,{}) for wid in selected_ids]
    keys=('records','biallelic_snps','polarized_sites','informative_sites','sum_abba','sum_baba','outgroup_missing','outgroup_heterozygous','outgroup_hom_ref','outgroup_hom_alt','missingness_exclusions')
    out={k:sum(float(v.get(k,0)) for v in vals) for k in keys}
    out['D']=d_from_sums(out['sum_abba'],out['sum_baba'])
    return out


def block_jackknife(windows: dict, window_regions: dict, block_size: int):
    """Return inside/outside/contrast block-delete jackknife summaries."""
    blockvals=defaultdict(lambda:defaultdict(lambda:[0.0,0.0,0]))
    for wid,s in windows.items():
        reg=window_regions.get(wid)
        if reg not in ('inside','outside'): continue
        start=int(wid.split('_')[1]); block=(start-1)//block_size
        x=blockvals[block][reg];x[0]+=s.get('sum_abba',0.0);x[1]+=s.get('sum_baba',0.0);x[2]+=s.get('informative_sites',0)
    all_blocks=sorted(blockvals)
    totals={r:[sum(blockvals[b][r][i] for b in all_blocks) for i in (0,1)] for r in ('inside','outside')}
    def se(estimates):
        n=len(estimates)
        if n<2:return None
        m=sum(estimates)/n
        return math.sqrt((n-1)/n*sum((x-m)**2 for x in estimates))
    out=[]
    for reg in ('inside','outside'):
        bs=[b for b in all_blocks if blockvals[b][reg][2]>0]
        est=[]
        for b in bs:
            a=totals[reg][0]-blockvals[b][reg][0]; bb=totals[reg][1]-blockvals[b][reg][1]
            d=d_from_sums(a,bb)
            if d is not None: est.append(d)
        d0=d_from_sums(*totals[reg]); s=se(est)
        out.append({'region':reg,'block_size_bp':block_size,'n_blocks':len(bs),'D':d0,'SE':s,'CI_low':d0-1.96*s if s is not None and d0 is not None else None,'CI_high':d0+1.96*s if s is not None and d0 is not None else None})
    bs=all_blocks; est=[]
    for b in bs:
        a_in=totals['inside'][0]-blockvals[b]['inside'][0];b_in=totals['inside'][1]-blockvals[b]['inside'][1]
        a_out=totals['outside'][0]-blockvals[b]['outside'][0];b_out=totals['outside'][1]-blockvals[b]['outside'][1]
        di=d_from_sums(a_in,b_in);do=d_from_sums(a_out,b_out)
        if di is not None and do is not None:est.append(di-do)
    di=d_from_sums(*totals['inside']);do=d_from_sums(*totals['outside']);delta=di-do if di is not None and do is not None else None;s=se(est)
    out.append({'region':'inside_minus_outside','block_size_bp':block_size,'n_blocks':len(bs),'D':delta,'SE':s,'CI_low':delta-1.96*s if s is not None and delta is not None else None,'CI_high':delta+1.96*s if s is not None and delta is not None else None})
    return out


def triple_topology(newick: str, taxa: tuple[str,str,str]):
    """Induce a rooted three-tip topology using the rooted published SNAPP tree."""
    from Bio import Phylo
    from io import StringIO
    tree=Phylo.read(StringIO(newick),'newick')
    tips={t.name for t in tree.get_terminals()}
    if not set(taxa)<=tips:return None
    pairs=[(taxa[0],taxa[1],'q1'),(taxa[0],taxa[2],'q2'),(taxa[1],taxa[2],'q3')]
    scores=[]
    for a,b,q in pairs:
        lca=tree.common_ancestor(a,b)
        scores.append((len(tree.get_path(lca)),q))
    mx=max(x[0] for x in scores); winners=[q for depth,q in scores if depth==mx]
    if len(winners)!=1:return None
    return winners[0]


def d_correlation(xs,ys):
    from scipy.stats import pearsonr,spearmanr
    pairs=[(float(x),float(y)) for x,y in zip(xs,ys) if x is not None and y is not None and math.isfinite(float(x)) and math.isfinite(float(y))]
    if len(pairs)<3:return len(pairs),None,None
    x=[p[0] for p in pairs];y=[p[1] for p in pairs]
    if len(set(x))<2 or len(set(y))<2:return len(pairs),None,None
    return len(pairs),float(pearsonr(x,y).statistic),float(spearmanr(x,y).statistic)
