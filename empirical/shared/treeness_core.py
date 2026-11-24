"""Shared treeness calculations and spatial summaries for empirical local trees."""
from __future__ import annotations
import csv, hashlib, io, math, random, statistics
from collections import Counter
from pathlib import Path
import numpy as np
from Bio import Phylo


def read_tsv(path):
    with Path(path).open(newline='', errors='replace') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write_tsv(path, rows, fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if fields is None: fields=list(rows[0]) if rows else []
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader()
        for r in rows:w.writerow({k:fmt(r.get(k)) for k in fields})


def fmt(x):
    if x is None:return 'NA'
    if isinstance(x,(float,np.floating)):
        return 'NA' if not math.isfinite(float(x)) else f'{float(x):.10g}'
    return str(x)


def calculate(newick, outgroups=(), rooted=False):
    """Calculate branch sums; root edge/stem is excluded, missing lengths stay missing.

    Internal edge means both sides of the edge contain >=2 retained tips. This
    edge-partition rule works for rooted and unrooted Newick representations.
    For the source tree, use all taxa. Optional outgroup sensitivity recomputes
    edge classes after restricting to the explicitly supplied ingroup set.
    """
    tree=Phylo.read(io.StringIO(newick.strip()),'newick')
    tips=[str(x.name) for x in tree.get_terminals()]
    if len(tips)!=len(set(tips)):raise ValueError('duplicate_tip_labels')
    if not tips:raise ValueError('no_tips')
    outgroups=set(outgroups)&set(tips)
    descendants={}
    for clade in tree.find_clades(order='postorder'):
        if clade.is_terminal():descendants[clade]={str(clade.name)}
        else:
            d=set()
            for child in clade.clades:d.update(descendants[child])
            descendants[clade]=d
    nonroot=[c for c in tree.find_clades(order='preorder') if c is not tree.root]
    def edge_sums(keep):
        n_keep=len(keep);internal=terminal=0.0;zero_i=zero_t=missing=invalid=valid=0;terminal_edges=[];ninternal=nterminal=0
        for clade in nonroot:
            desc=descendants[clade]&keep;nd=len(desc);other=n_keep-nd
            if nd==0 or other==0:continue
            is_internal=(nd>=2 if rooted else (nd>=2 and other>=2))
            is_terminal=(nd==1 if rooted else (nd==1 or other==1))
            if not is_internal and not is_terminal:continue
            x=clade.branch_length
            if x is None:missing+=1;continue
            try:x=float(x)
            except (TypeError,ValueError):invalid+=1;continue
            if not math.isfinite(x) or x<0:invalid+=1;continue
            valid+=1
            if is_internal:
                ninternal+=1;internal+=x
                if x==0:zero_i+=1
            else:
                nterminal+=1;terminal+=x;terminal_edges.append(x)
                if x==0:zero_t+=1
        total=internal+terminal
        complete=missing==0 and invalid==0
        treeness=(internal/total if complete and total>0 else math.nan)
        return {'internal_branch_length_sum':internal,
                'terminal_branch_length_sum':terminal,
                'total_tree_length':total if complete else math.nan,
                'treeness':treeness,'n_valid_nonroot_branches':valid,
                'n_internal_branches':ninternal,'n_terminal_branches':nterminal,
                'n_zero_length_internal_branches':zero_i,'n_zero_length_terminal_branches':zero_t,
                'max_terminal_branch_length':max(terminal_edges) if terminal_edges else math.nan,
                'max_terminal_branch_fraction':max(terminal_edges)/total if terminal_edges and complete and total>0 else math.nan,
                'n_missing_branch_lengths':missing,'n_invalid_branch_lengths':invalid,'n_tips':len(keep),
                'root_stem_length_excluded':tree.root.branch_length is not None}
    m=edge_sums(set(tips))
    m['taxa']='|'.join(sorted(tips));m['taxon_set_sha256']=hashlib.sha256('\n'.join(sorted(tips)).encode()).hexdigest()
    m['n_outgroup_tips']=len(outgroups)
    if outgroups and len(set(tips)-outgroups)>=3:
        p=edge_sums(set(tips)-outgroups)
        m['treeness_no_outgroups']=p['treeness'];m['delta_treeness_all_minus_no_outgroups']=(m['treeness']-p['treeness'] if math.isfinite(m['treeness']) and math.isfinite(p['treeness']) else math.nan)
        m['internal_sum_no_outgroups']=p['internal_branch_length_sum'];m['terminal_sum_no_outgroups']=p['terminal_branch_length_sum']
    else:
        m.update({'treeness_no_outgroups':math.nan,'delta_treeness_all_minus_no_outgroups':math.nan,'internal_sum_no_outgroups':math.nan,'terminal_sum_no_outgroups':math.nan})
    return m


def mark_region(start,end,intervals):
    if not intervals:return 'unclassified',''
    overlaps=[];full=[]
    for region_id,a,b in intervals:
        if start<=b and end>=a:
            overlaps.append(region_id)
            if start>=a and end<=b:full.append(region_id)
    if full:return 'inside',';'.join(full)
    if overlaps:return 'boundary',';'.join(overlaps)
    return 'outside',''


def _runs(rows, value_key):
    rows=sorted(rows,key=lambda r:(r['chromosome'],float(r['start'])))
    if not rows:return []
    widths=[max(1,float(r['end'])-float(r['start'])+1) for r in rows]
    gaplim=max(statistics.median(widths)*2,1)
    runs=[];cur=[];last=None;chrom=None
    for r in rows:
        c=r['chromosome'];st=float(r['start']);en=float(r['end'])
        if cur and (c!=chrom or st-last>gaplim):runs.append(cur);cur=[]
        cur.append(float(r[value_key]));last=en;chrom=c
    if cur:runs.append(cur)
    return runs


def _sample_runs(runs, block_len, rng):
    vals=[]
    for run in runs:
        effective=max(1,min(block_len,max(1,len(run)//2)))
        if len(run)<3:vals.extend(run);continue
        sample=[]
        while len(sample)<len(run):
            start=rng.randrange(len(run));sample.extend(run[(start+j)%len(run)] for j in range(effective))
        vals.extend(sample[:len(run)])
    return vals


def compare_inside_outside(rows, reps=2000, seed=20261124, block_bp=2_000_000):
    """Block-bootstrap effects within each chromosome/stratum, no window iid test."""
    eligible=[r for r in rows if r.get('region_class') in ('inside','outside') and math.isfinite(float(r['treeness']))]
    if not eligible:return []
    boundary_counts=Counter((r['tree_set'],r['chromosome']) for r in rows if r.get('region_class')=='boundary' and r.get('valid')=='true')
    groups={}
    for r in eligible:groups.setdefault((r['tree_set'],r['chromosome']),[]).append(r)
    out=[];rng=random.Random(seed)
    for (tree_set,chrom),group in groups.items():
        region_id=';'.join(sorted({r.get('region_id','') for r in group if r['region_class']=='inside' and r.get('region_id')}))
        sigs={}
        for r in group:sigs.setdefault(r['taxon_set_sha256'],[]).append(r)
        modal=max(sigs,key=lambda x:len(sigs[x]));comparable=sigs[modal]
        ins=[r for r in comparable if r['region_class']=='inside'];outs=[r for r in comparable if r['region_class']=='outside']
        if not ins or not outs:continue
        iw=[float(r['end'])-float(r['start'])+1 for r in comparable]
        bw=max(1,int(round(block_bp/statistics.median(iw))))
        in_t=[float(r['treeness']) for r in ins];out_t=[float(r['treeness']) for r in outs]
        boot=[]
        if min(len(ins),len(outs))>=3:
            ir=_runs(ins,'treeness');orr=_runs(outs,'treeness')
            for _ in range(reps):boot.append(statistics.mean(_sample_runs(ir,bw,rng))-statistics.mean(_sample_runs(orr,bw,rng)))
        row={'tree_set':tree_set,'chromosome':chrom,'region_id':region_id,'n_inside':len(ins),'n_outside':len(outs),'n_boundary':boundary_counts.get((tree_set,chrom),0),'n_taxa_comparison':int(comparable[0]['n_taxa']),'n_excluded_nonmodal_taxa':len(group)-len(comparable),'median_treeness_inside':statistics.median(in_t),'median_treeness_outside':statistics.median(out_t),'delta_median_treeness':statistics.median(in_t)-statistics.median(out_t),'mean_treeness_inside':statistics.mean(in_t),'mean_treeness_outside':statistics.mean(out_t),'delta_mean_treeness':statistics.mean(in_t)-statistics.mean(out_t),'bootstrap_block_bp':block_bp,'bootstrap_replicates':len(boot),'delta_mean_ci_low':float(np.quantile(boot,.025)) if boot else math.nan,'delta_mean_ci_high':float(np.quantile(boot,.975)) if boot else math.nan,'inference_note':'moving-block bootstrap over genomic runs; one chromosome is a genomic series, not biological replication'}
        for metric,outcol,incol in [('internal_branch_length_sum','median_internal_inside','median_internal_outside'),('terminal_branch_length_sum','median_terminal_inside','median_terminal_outside')]:
            row[outcol]=statistics.median(float(r[metric]) for r in ins);row[incol]=statistics.median(float(r[metric]) for r in outs);row['delta_'+outcol.replace('median_','')]=row[outcol]-row[incol]
        out.append(row)
    return out
