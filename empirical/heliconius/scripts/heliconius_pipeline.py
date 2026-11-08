#!/usr/bin/env python3
"""Auditable Edelman et al. (2019) Heliconius local-tree pilot.

Subcommands are intentionally small and rerunnable: download, audit, normalize,
analyze, figures, or all.  The parser is dependency-light so audit/tests can
run before optional plotting packages are installed.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, io, json, math, re, sys, urllib.request, zipfile
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data/heliconius"; EMP = ROOT / "empirical/heliconius"
RAW = DATA / "raw"; META = DATA / "metadata"; PROC = DATA / "processed"; RES = EMP / "results"; FIG = EMP / "figures"
SOURCES = [
    # Dryad's file_stream endpoints return 403 in this environment. These are
    # the DOI-associated Zenodo mirror files, with matching published MD5s.
    ("erato_10kb", "eratoClade_10KBAbutting_HeraRef.ge2kb.trees.tsv", "https://zenodo.org/api/records/5039535/files/eratoClade_10KBAbutting_HeraRef.ge2kb.trees.tsv/content", "H_erato_demophoon_v1", 10000),
    ("erato_50kb_hera", "eratoClade.50KBSliding.toHeraGenome.tsv", "https://zenodo.org/api/records/5039535/files/eratoClade.50KBSliding.toHeraGenome.tsv/content", "H_erato_demophoon_v1", 50000),
    ("erato_50kb_hmel", "eratoClade_50KBSliding.allTrees.toHmelGenome.tsv", "https://zenodo.org/api/records/5039535/files/eratoClade_50KBSliding.allTrees.toHmelGenome.tsv/content", "Hmel2.5", 50000),
    ("melpomene_50kb", "melpomeneClade.50KBSliding.allTrees.nwk.tsv", "https://zenodo.org/api/records/5039535/files/melpomeneClade.50KBSliding.allTrees.nwk.tsv/content", "Hmel2.5", 50000),
    ("genome_info", "eratoGenomeInfo.zip", "https://zenodo.org/api/records/5039535/files/eratoGenomeInfo.zip/content", "H_erato_demophoon_v1", 10000),
]
QUARTET = ("HeraRef", "Hhim", "Htel", "Hhsa")
BASELINE = frozenset((frozenset((QUARTET[0], QUARTET[1])), frozenset((QUARTET[2], QUARTET[3]))))
ALTERNATIVE = frozenset((frozenset((QUARTET[0], QUARTET[2])), frozenset((QUARTET[1], QUARTET[3]))))

def sha256(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20), b''): h.update(b)
    return h.hexdigest()

class Node:
    def __init__(self, name=None, length=None): self.name=name; self.length=length; self.children=[]

def parse_newick(s):
    s=s.strip(); i=0
    def sub():
        nonlocal i
        n=Node()
        if i<len(s) and s[i]=='(':
            i+=1
            while True:
                n.children.append(sub())
                if i>=len(s): raise ValueError('unterminated Newick')
                if s[i]==',': i+=1; continue
                if s[i]==')': i+=1; break
                raise ValueError(f'expected comma or close at {i}')
        start=i
        while i<len(s) and s[i] not in ',():;': i+=1
        n.name=s[start:i].strip() or None
        if i<len(s) and s[i]==':':
            i+=1; start=i
            while i<len(s) and s[i] not in ',();': i+=1
            try: n.length=float(s[start:i])
            except ValueError: raise ValueError('invalid branch length')
        return n
    root=sub()
    if i<len(s) and s[i]==';': i+=1
    if s[i:].strip(): raise ValueError('trailing Newick text')
    tips=[]
    def walk(n):
        if n.children:
            for c in n.children: walk(c)
        else: tips.append(n)
    walk(root)
    if not tips or any(t.name is None for t in tips): raise ValueError('unnamed terminal')
    if len({t.name for t in tips}) != len(tips): raise ValueError('duplicate terminal')
    return root, tips

def distances(root, weighted=True):
    adj=defaultdict(list); tips={}
    def walk(n):
        for c in n.children:
            w=(c.length if weighted and c.length is not None else 1.0)
            adj[id(n)].append((id(c),w)); adj[id(c)].append((id(n),w)); walk(c)
        if not n.children: tips[n.name]=id(n)
    walk(root); out={}
    for name,node in tips.items():
        stack=[(node,-1,0.0)]
        while stack:
            u,prev,d=stack.pop()
            for other,onode in tips.items():
                if onode==u: out[tuple(sorted((name,other)))]=d
            for v,w in adj[u]:
                if v!=prev: stack.append((v,u,d+w))
    return out

def split_for_quartet(root, taxa):
    d=distances(root, weighted=False)
    if any(tuple(sorted(p)) not in d for p in [(taxa[0],taxa[1]),(taxa[0],taxa[2]),(taxa[0],taxa[3]),(taxa[1],taxa[2]),(taxa[1],taxa[3]),(taxa[2],taxa[3])]): return None
    pairs=[((taxa[0],taxa[1]),(taxa[2],taxa[3])),((taxa[0],taxa[2]),(taxa[1],taxa[3])),((taxa[0],taxa[3]),(taxa[1],taxa[2]))]
    vals=[d[tuple(sorted(a))]+d[tuple(sorted(b))] for a,b in pairs]; m=min(vals)
    if vals.count(m)!=1: return None
    a,b=pairs[vals.index(m)]
    return frozenset((frozenset(a),frozenset(b)))

def split_for_quartet_from_dist(d, taxa, weighted=False):
    pairs=[((taxa[0],taxa[1]),(taxa[2],taxa[3])),((taxa[0],taxa[2]),(taxa[1],taxa[3])),((taxa[0],taxa[3]),(taxa[1],taxa[2]))]
    try: vals=[d[tuple(sorted(a))]+d[tuple(sorted(b))] for a,b in pairs]
    except KeyError: return None
    m=min(vals)
    if sum(abs(v-m)<=1e-12 for v in vals)!=1: return None
    a,b=pairs[vals.index(m)]
    return frozenset((frozenset(a),frozenset(b)))

def quartet_id(taxa): return '|'.join(taxa)

def quartet_reference_mapping(records):
    """Freeze per-quartet baseline topologies from background windows.

    The background is defined before focal summaries as all 10-kb Hera windows
    outside chromosomes 2 and 15.  The Hmel2.5 50-kb source has its own
    background because its taxon universe and coordinate system differ.
    """
    source_rows=defaultdict(list)
    for r in records:
        source_rows[r['source_id']].append(r)
    universes={}
    for sid,rs in source_rows.items():
        counts=Counter(tuple(sorted(r['taxa'].split('|'))) for r in rs)
        universes[sid]=counts.most_common(1)[0][0]
    ref_by_group={}
    for group,sids in [('hera',('erato_10kb',)),('hmel',('erato_50kb_hmel',))]:
        rs=[r for sid in sids for r in source_rows.get(sid,[]) if str(r['chromosome']) not in ('2','15')]
        taxa=universes.get(sids[0],())
        counts=defaultdict(Counter)
        for r in rs:
            try: root,tips=parse_newick(r['tree_newick']); du=distances(root,weighted=False)
            except Exception: continue
            for qt in combinations(taxa,4):
                sp=split_for_quartet_from_dist(du,qt)
                if sp: counts[qt][split_label(sp)]+=1
        mapping={}
        for qt in combinations(taxa,4):
            possible=sorted(split_label(frozenset((frozenset((qt[0],qt[1])),frozenset((qt[2],qt[3]))))) for _ in [0])
            all_splits=sorted({split_label(frozenset((frozenset(a),frozenset(b)))) for a,b in [((qt[0],qt[1]),(qt[2],qt[3])),((qt[0],qt[2]),(qt[1],qt[3])),((qt[0],qt[3]),(qt[1],qt[2]))]})
            observed=counts[qt]
            baseline=max(all_splits,key=lambda s:(observed.get(s,0),-all_splits.index(s)))
            alternatives=[s for s in all_splits if s!=baseline]
            mapping[qt]={'quartet_id':quartet_id(qt),'taxa':'|'.join(qt),'baseline_topology':baseline,'alternative_topology':alternatives[0],'third_topology':alternatives[1],'background_n':sum(observed.values()),'background_counts':';'.join(f'{s}:{observed.get(s,0)}' for s in all_splits)}
        ref_by_group[group]=mapping
    return universes,ref_by_group

def all_quartet_analysis(records):
    if not records: return
    universes,ref_by_group=quartet_reference_mapping(records)
    group_for_source={'erato_10kb':'hera','erato_50kb_hera':'hera','erato_50kb_hmel':'hmel'}
    group_universe={'hera':universes.get('erato_10kb',()),'hmel':universes.get('erato_50kb_hmel',())}
    ref_rows=[]
    for group,mapping in ref_by_group.items():
        for row in mapping.values(): ref_rows.append({'reference_group':group,'reference_definition':'10-kb H. erato background excluding chromosomes 2 and 15' if group=='hera' else 'Hmel2.5 50-kb background excluding chromosomes 2 and 15','reference_tree':'per-quartet background majority; no single fully bifurcating tree imposed','baseline_topology_source':'background majority with deterministic lexical tie-break',**row})
    with (PROC/'all_quartet_reference_baselines.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(ref_rows[0]),delimiter='\t'); w.writeheader(); w.writerows(ref_rows)
    assignments=[]; summaries=[]
    for r in records:
        group=group_for_source.get(r['source_id']); mapping=ref_by_group.get(group,{})
        try: root,tips=parse_newick(r['tree_newick']); du=distances(root,weighted=False); dw=distances(root,weighted=True)
        except Exception: continue
        pair_vals=[v for k,v in dw.items() if k[0]!=k[1]]
        mean_pairwise=sum(pair_vals)/len(pair_vals) if pair_vals else 'NA'
        available=[]; qcounts=Counter(); internal=[]; row_assignments=[]
        for qt in combinations(group_universe.get(group,()),4):
            qtmap=mapping.get(qt)
            if not qtmap: continue
            sp=split_for_quartet_from_dist(du,qt)
            cls='unresolved' if sp is None else ('q1' if split_label(sp)==qtmap['baseline_topology'] else ('q2' if split_label(sp)==qtmap['alternative_topology'] else 'q3'))
            if cls!='unresolved': qcounts[cls]+=1; available.append(qt)
            ib='NA'
            if sp is not None:
                sides=sorted(sp,key=lambda z:sorted(z)); a,b=sorted(sides[0]),sorted(sides[1]); x,y=a; z,w=b
                try: ib=(dw[tuple(sorted((x,z)))]+dw[tuple(sorted((y,w)))]-dw[tuple(sorted((x,y))) ]-dw[tuple(sorted((z,w))) ])/2
                except KeyError: pass
            row_assignments.append({'source_id':r['source_id'],'source_tree':r['source_tree'],'assembly':r['assembly'],'chromosome':r['chromosome'],'window_id':r.get('window_id',''),'start':r['start'],'end':r['end'],'midpoint':r['midpoint'],'region_class':r['region_class'],'quartet_id':qtmap['quartet_id'],'taxa':qtmap['taxa'],'baseline_topology':qtmap['baseline_topology'],'alternative_topology':qtmap['alternative_topology'],'third_topology':qtmap['third_topology'],'local_topology':split_label(sp),'topology_class':cls,'internal_branch_length':ib,'mean_pairwise_patristic':mean_pairwise,'recomb_rate_cM_Mb':r.get('recomb_rate_cM_Mb','')})
        assignments.extend(row_assignments)
        valid_internal=[float(a['internal_branch_length']) for a in row_assignments if a['internal_branch_length']!='NA']
        n=len(available); summaries.append({'source_id':r['source_id'],'source_tree':r['source_tree'],'assembly':r['assembly'],'chromosome':r['chromosome'],'window_id':r.get('window_id',''),'start':r['start'],'end':r['end'],'midpoint':r['midpoint'],'region_class':r['region_class'],'n_available_quartets':n,'n_unresolved_quartets':len(mapping)-n,'q1':qcounts['q1']/n if n else 'NA','q2':qcounts['q2']/n if n else 'NA','q3':qcounts['q3']/n if n else 'NA','D':(qcounts['q2']-qcounts['q1'])/n if n else 'NA','mean_quartet_internal_branch':sum(float(a['internal_branch_length']) for a in assignments[-len(mapping):] if a['internal_branch_length']!='NA')/max(1,sum(a['internal_branch_length']!='NA' for a in assignments[-len(mapping):])),'mean_pairwise_patristic':r.get('mean_pairwise_patristic','NA'),'recomb_rate_cM_Mb':r.get('recomb_rate_cM_Mb','')})
        summaries[-1]['mean_quartet_internal_branch']=sum(valid_internal)/len(valid_internal) if valid_internal else 'NA'
        summaries[-1]['mean_pairwise_patristic']=mean_pairwise
    assignment_fields=['source_id','source_tree','assembly','chromosome','window_id','start','end','midpoint','region_class','quartet_id','taxa','baseline_topology','alternative_topology','third_topology','local_topology','topology_class','internal_branch_length','mean_pairwise_patristic','recomb_rate_cM_Mb']
    with gzip.open(PROC/'all_quartet_assignments.tsv.gz','wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=assignment_fields,delimiter='\t'); w.writeheader(); w.writerows(assignments)
    with (PROC/'all_quartet_window_summary.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summaries[0]),delimiter='\t'); w.writeheader(); w.writerows(summaries)
    quartet_groups=defaultdict(list)
    for a in assignments:
        quartet_groups[(a['source_id'],a['assembly'],a['chromosome'],a['region_class'],a['quartet_id'])].append(a)
    quartet_rows=[]
    for (sid,assembly,chrom,region,quartet),ars in sorted(quartet_groups.items()):
        n=len(ars); resolved=[a for a in ars if a['topology_class']!='unresolved']
        internal=[float(a['internal_branch_length']) for a in ars if a['internal_branch_length']!='NA']
        quartet_rows.append({'source_id':sid,'assembly':assembly,'chromosome':chrom,'region_class':region,'quartet_id':quartet,'taxa':ars[0]['taxa'],'baseline_topology':ars[0]['baseline_topology'],'alternative_topology':ars[0]['alternative_topology'],'third_topology':ars[0]['third_topology'],'n_trees':n,'n_resolved':len(resolved),'q1':sum(a['topology_class']=='q1' for a in ars)/len(resolved) if resolved else 'NA','q2':sum(a['topology_class']=='q2' for a in ars)/len(resolved) if resolved else 'NA','q3':sum(a['topology_class']=='q3' for a in ars)/len(resolved) if resolved else 'NA','D':(sum(a['topology_class']=='q2' for a in ars)-sum(a['topology_class']=='q1' for a in ars))/len(resolved) if resolved else 'NA','mean_internal_branch_length':sum(internal)/len(internal) if internal else 'NA'})
    with (PROC/'all_quartet_quartet_summary.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(quartet_rows[0]),delimiter='\t'); w.writeheader(); w.writerows(quartet_rows)
    bins=defaultdict(list)
    for r in summaries:
        if r['q1']!='NA': bins[(r['source_id'],r['assembly'],r['chromosome'],int(float(r['start']))//50000*50000)].append(r)
    bin_rows=[]
    for (sid,assembly,chrom,bs),rs in sorted(bins.items()):
        def avg(key):
            vals=[]
            for x in rs:
                try: vals.append(float(x[key]))
                except (TypeError,ValueError): pass
            return sum(vals)/len(vals) if vals else 'NA'
        bin_rows.append({'source_id':sid,'assembly':assembly,'chromosome':chrom,'bin_start':bs,'bin_end':bs+49999,'n_trees':len(rs),'q1':avg('q1'),'q2':avg('q2'),'q3':avg('q3'),'D':avg('D'),'mean_pairwise_patristic':avg('mean_pairwise_patristic'),'mean_quartet_internal_branch':avg('mean_quartet_internal_branch'),'recomb_rate_cM_Mb':avg('recomb_rate_cM_Mb')})
    with (PROC/'all_quartet_binned_statistics.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(bin_rows[0]),delimiter='\t'); w.writeheader(); w.writerows(bin_rows)
    transitions=[]
    for sid in sorted(set(r['source_id'] for r in bin_rows)):
        for chrom in sorted(set(r['chromosome'] for r in bin_rows if r['source_id']==sid),key=lambda x:int(x) if x.isdigit() else 999):
            rr=sorted([r for r in bin_rows if r['source_id']==sid and r['chromosome']==chrom],key=lambda r:int(r['bin_start']))
            for a,b in zip(rr,rr[1:]): transitions.append({'source_id':sid,'chromosome':chrom,'left_bin_end':a['bin_end'],'right_bin_start':b['bin_start'],'transition_score':math.sqrt(sum((float(b[k])-float(a[k]))**2 for k in ('q1','q2','q3'))),'D_jump':float(b['D'])-float(a['D'])})
    transitions.sort(key=lambda r:abs(r['D_jump']),reverse=True)
    if transitions:
        with (RES/'all_quartet_transition_summary.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(transitions[0]),delimiter='\t'); w.writeheader(); w.writerows(transitions[:100])

def split_label(split): return '|'.join(','.join(sorted(x)) for x in sorted(split, key=lambda z: sorted(z))) if split else 'unresolved'

def tree_from_row(row):
    for k,v in row.items():
        if v and ('tree' in k.lower() or 'newick' in k.lower() or v.lstrip().startswith('(')) and v.lstrip().startswith('('): return v.strip()
    return None

def coord(row, keys):
    for k in keys:
        for actual in row:
            if actual.lower()==k and row[actual].strip().isdigit(): return int(row[actual])
    return None

def read_rows(path):
    with path.open(newline='') as f:
        sample=f.read(4096); f.seek(0); delim='\t' if '\t' in sample else ','
        return list(csv.DictReader(f, delimiter=delim))

def read_source_rows(path):
    rows=[]
    text=path.read_text()
    for match in re.finditer(r'(?m)(chr[^\s]+)[\t ]+(\(.*?;)', text, re.S):
        rows.append({'window_id':match.group(1),'tree_newick':match.group(2),'source_line':text.count('\n',0,match.start())+1})
    return rows or read_rows(path)

def parse_window_id(window_id):
    m=re.match(r'^chr([^_]+)_(\d+)$', str(window_id))
    return (m.group(1), int(m.group(2))) if m else ('', None)

def recombination_map():
    p=RAW/'eratoGenomeInfo.zip'; out={}
    if not p.exists(): return out
    with zipfile.ZipFile(p) as z:
        name='HeratoWindows.10KBAbutting.recRates.CDSoverlaps.bed'
        if name not in z.namelist(): return out
        for line in io.TextIOWrapper(z.open(name), encoding='utf-8'):
            x=line.rstrip().split('\t')
            if len(x)>=10:
                try:
                    out[x[3]]={'recomb_rate_cM_Mb':float(x[9]),'cds_overlap_bp':int(x[10]) if len(x)>10 else ''}
                except ValueError: pass
    return out

def verified_regions(assembly, chrom):
    out=[]
    p=EMP/'config/regions.tsv'
    if not p.exists(): return out
    for r in read_rows(p):
        if r.get('assembly')!=assembly or r.get('chromosome')!=str(chrom) or r.get('status')!='verified': continue
        try: out.append((int(r['start']),int(r['end']),r['region_id']))
        except (TypeError,ValueError): continue
    return out

def classify(start,end,assembly,chrom=''):
    regions=verified_regions(assembly,chrom)
    if not regions: return 'unclassified'
    inside=any(start>=a and end<=b for a,b,_ in regions)
    overlap=any(start<=b and end>=a for a,b,_ in regions)
    return 'inside' if inside else ('overlap' if overlap else 'outside')

def download():
    RAW.mkdir(parents=True,exist_ok=True); rows=[]
    for sid,name,url,assembly,size in SOURCES:
        out=RAW/name
        if not out.exists():
            try:
                with urllib.request.urlopen(url, timeout=60) as r, out.open('wb') as f:
                    while True:
                        b=r.read(1<<20)
                        if not b: break
                        f.write(b)
            except Exception as e:
                if out.exists(): out.unlink()
                print(f'BLOCKED {name}: {e}', file=sys.stderr)
                continue
        rows.append({'source_name':sid,'original_filename':name,'url_or_doi':url,'download_date':'','size_bytes':out.stat().st_size,'sha256':sha256(out),'purpose':'local-tree or metadata source','local_path':str(out.relative_to(ROOT))})
    if rows:
        with (META/'download_manifest.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t'); w.writeheader(); w.writerows(rows)
    return rows

def audit():
    out=[]; total=0
    for sid,name,url,assembly,wsize in SOURCES:
        p=RAW/name
        rec={'source_id':sid,'filename':name,'assembly':assembly,'expected_window_size_bp':wsize,'exists':p.exists(),'n_rows':0,'n_tree_rows':0,'n_malformed':0,'n_taxa_max':0,'branch_lengths_present':'unknown','chromosomes':'','window_sizes':'','starts':'','endpoints':'','notes':''}
        if p.suffix=='.zip' or not p.exists(): out.append(rec); continue
        try: rows=read_source_rows(p)
        except Exception as e: rec['notes']=f'read_error:{e}'; out.append(rec); continue
        rec['n_rows']=len(rows); chrom=[]; starts=[]; ends=[]; sizes=[]; taxa=set(); branch=False
        for row in rows:
            nwk=tree_from_row(row)
            if not nwk: continue
            rec['n_tree_rows']+=1; total+=1
            try:
                root,tips=parse_newick(nwk); taxa.update(t.name for t in tips); rec['n_taxa_max']=max(rec['n_taxa_max'],len(tips)); branch |= any(t.length is not None for t in tips)
                st=coord(row,('start','window_start','chrom_start','beg')); en=coord(row,('end','window_end','chrom_end','stop'))
                if st is not None and en is not None: starts.append(st); ends.append(en); sizes.append(en-st+1)
                for k,v in row.items():
                    if k.lower() in ('chrom','chromosome','chr','scaffold') and v: chrom.append(v)
            except Exception: rec['n_malformed']+=1
        rec['branch_lengths_present']=branch; rec['chromosomes']=','.join(sorted(set(chrom))); rec['window_sizes']=','.join(map(str,sorted(set(sizes)))); rec['starts']=str(min(starts)) if starts else ''; rec['endpoints']=str(max(ends)) if ends else ''; rec['notes']=f'unique_taxa={len(taxa)}; overlapping={sid.endswith("50kb_hera") or sid.endswith("50kb_hmel")}'
        out.append(rec)
    RES.mkdir(parents=True,exist_ok=True)
    fields=list(out[0]) if out else ['source_id']
    with (RES/'stage0_source_audit.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(out)
    return out

def normalize():
    records=[]; issues=[]
    rec_map=recombination_map()
    for sid,name,url,assembly,wsize in SOURCES[:3]:
        p=RAW/name
        if not p.exists(): continue
        for i,row in enumerate(read_source_rows(p),1):
            nwk=tree_from_row(row)
            if not nwk: continue
            chrom,idx=parse_window_id(row.get('window_id',''))
            step=10000 if sid=='erato_10kb' else 25000
            st=(idx-1)*step+1 if idx is not None else None; en=st+(wsize-1) if st is not None else None
            if st is None or en is None: issues.append(f'{name}:{i}:missing_coordinates'); continue
            try:
                root,tips=parse_newick(nwk); weighted=all(t.length is not None for t in tips)
            except Exception as e: issues.append(f'{name}:{i}:{e}'); continue
            rm=rec_map.get(row.get('window_id',''),{})
            records.append({'source_id':sid,'source_tree':name,'assembly':assembly,'chromosome':chrom,'window_id':row.get('window_id',''),'start':st,'end':en,'midpoint':(st+en)/2,'window_size_bp':en-st+1,'taxa':'|'.join(sorted(t.name for t in tips)),'n_taxa':len(tips),'branch_lengths_available':weighted,'region_class':classify(st,en,assembly,chrom),'recomb_rate_cM_Mb':rm.get('recomb_rate_cM_Mb',''),'cds_overlap_bp':rm.get('cds_overlap_bp',''),'tree_newick':nwk})
    PROC.mkdir(parents=True,exist_ok=True); fields=list(records[0]) if records else ['source_id']
    with (PROC/'local_trees.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(records)
    (RES/'stage1_normalization_issues.tsv').write_text('issue\n'+'\n'.join(issues)+'\n')
    return records

def analyze():
    p=PROC/'local_trees.tsv'
    if not p.exists(): normalize()
    if not p.exists() or p.stat().st_size<20: return []
    rows=read_rows(p); obs=[]
    for r in rows:
        try: root,tips=parse_newick(r['tree_newick'])
        except Exception: continue
        labels={t.name for t in tips}; missing=[x for x in QUARTET if x not in labels]
        split=split_for_quartet(root,QUARTET) if not missing else None
        d=distances(root,weighted=True); all_lengths=all(t.length is not None for t in tips)
        pair=[v for k,v in d.items() if k[0]!=k[1]]; internal=[]
        if all_lengths and split:
            sides=sorted(split,key=lambda z:sorted(z)); a,b=sorted(sides[0]),sorted(sides[1]); x,y=a[0],a[1]; z,w=b[0],b[1]
            internal=[(d[tuple(sorted((x,z)))]+d[tuple(sorted((y,w)))]-d[tuple(sorted((x,y))) ]-d[tuple(sorted((z,w))) ])/2]
        obs.append({**r,'missing_taxa':'|'.join(missing),'local_split':split_label(split),'topology_class':'baseline' if split==BASELINE else ('alternative' if split==ALTERNATIVE else ('third' if split else 'unresolved')),'q1':int(split==BASELINE),'q2':int(split==ALTERNATIVE),'q3':int(split is not None and split not in (BASELINE,ALTERNATIVE)),'D':int(split==ALTERNATIVE)-int(split==BASELINE),'mean_pairwise_patristic':sum(pair)/len(pair) if all_lengths and pair else 'NA','internal_branch_length':internal[0] if internal else 'NA','branch_length_units':'source-tree patristic units' if all_lengths else 'unavailable'})
    fields=list(obs[0]) if obs else []
    with (PROC/'local_tree_statistics.tsv').open('w',newline='') as f:
        if fields: w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader(); w.writerows(obs)
    # Explicit non-independent 50-kb aggregation: bins are descriptive summaries.
    bins=[]; grouped=defaultdict(list)
    for r in obs:
        if r['chromosome'] and r['window_size_bp']:
            key=(r['assembly'],r['chromosome'],int(float(r['start']))//50000*50000); grouped[key].append(r)
    for (assembly,chrom,binstart),rs in sorted(grouped.items()):
        bins.append({'assembly':assembly,'chromosome':chrom,'bin_start':binstart,'bin_end':binstart+49999,'n_trees':len(rs),'q1':sum(int(x['q1']) for x in rs)/len(rs),'q2':sum(int(x['q2']) for x in rs)/len(rs),'q3':sum(int(x['q3']) for x in rs)/len(rs),'D':sum(float(x['D']) for x in rs)/len(rs),'mean_pairwise_patristic':sum(float(x['mean_pairwise_patristic']) for x in rs if x['mean_pairwise_patristic']!='NA')/max(1,sum(x['mean_pairwise_patristic']!='NA' for x in rs))})
    if bins:
        with (PROC/'binned_statistics_50kb.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(bins[0]),delimiter='\t'); w.writeheader(); w.writerows(bins)
    all_quartet_analysis(rows)
    return obs

def all_quartet_figures():
    p=PROC/'all_quartet_binned_statistics.tsv'
    if not p.exists(): return
    try:
        import os, tempfile
        os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'msrc-paper-mplconfig'))
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    except ImportError: return
    rows=read_rows(p)
    def affected(row):
        c=classify(int(row['bin_start'])+1,int(row['bin_end'])+1,row['assembly'],row['chromosome'])
        return {'outside':0.0,'overlap':0.5,'inside':1.0}.get(c,float('nan'))
    def render(source,chrom,filename):
        sub=sorted([r for r in rows if r['source_id']==source and r['chromosome']==chrom],key=lambda r:int(r['bin_start']))
        if not sub: return
        x=[(int(r['bin_start'])+int(r['bin_end']))/2/1e6 for r in sub]
        fig,ax=plt.subplots(5,1,figsize=(12,10),sharex=True,constrained_layout=True,gridspec_kw={'height_ratios':[2,1,1,0.8,1.2]})
        for j,k,c in [(0,'q1','#4c72b0'),(0,'q2','#dd8452'),(0,'q3','#55a868'),(1,'D','#c44e52'),(2,'mean_pairwise_patristic','#8172b3')]:
            y=[float(r[k]) if r[k] not in ('','NA') else float('nan') for r in sub]; ax[j].plot(x,y,'.-',ms=3,label=k,color=c)
        reg=[affected(r) for r in sub]; ax[3].plot(x,reg,'s-',ms=3,color='#7f7f7f'); ax[3].set_ylim(-.1,1.1); ax[3].set_yticks([0,.5,1]); ax[3].set_yticklabels(['outside','boundary','inside']); ax[3].set_ylabel('affected\nregion')
        rec=[float(r['recomb_rate_cM_Mb']) if r['recomb_rate_cM_Mb'] not in ('','NA') else float('nan') for r in sub]; ax[4].plot(x,rec,'.-',ms=3,color='#1b9e77'); ax[4].set_ylabel('recombination\ncM/Mb')
        if not any(math.isfinite(v) for v in rec): ax[4].text(.01,.5,'quantitative recombination unavailable',transform=ax[4].transAxes,va='center')
        if not any(math.isfinite(v) for v in reg): ax[3].text(.01,.5,'affected-region boundary unresolved for this assembly',transform=ax[3].transAxes,va='center')
        ax[0].legend(ncol=3); ax[0].set_ylabel('all-quartet support'); ax[1].set_ylabel('D = q2 - q1'); ax[2].set_ylabel('mean patristic'); ax[4].set_xlabel(f'chromosome {chrom} position (Mb)'); fig.suptitle(f'Heliconius all-quartet tracks: {source}, chromosome {chrom}')
        for ext in ('pdf','png','svg'): fig.savefig(FIG/f'{filename}.{ext}',dpi=220 if ext=='png' else None)
        plt.close(fig)
    for source,prefix in [('erato_10kb','10kb_hera'),('erato_50kb_hera','50kb_hera'),('erato_50kb_hmel','50kb_hmel')]:
        for chrom in ('2','15'): render(source,chrom,f'heliconius_{prefix}_chr{chrom}_all_quartet_tracks')
    # Genome-wide cumulative Hera 10-kb track, retaining chromosome gaps.
    sub=[r for r in rows if r['source_id']=='erato_10kb']
    if sub:
        chroms=sorted(set(r['chromosome'] for r in sub),key=lambda z:int(z) if z.isdigit() else 999); offsets={}; off=0; ticks=[]; cum=[]
        for chrom in chroms:
            rr=[r for r in sub if r['chromosome']==chrom]; length=max(int(r['bin_end']) for r in rr); offsets[chrom]=off; ticks.append((off+length/2/1e6,chrom)); off+=(length+500000)
        for r in sub: cum.append((offsets[r['chromosome']]+(int(r['bin_start'])+int(r['bin_end']))/2)/1e6)
        order=sorted(range(len(sub)),key=lambda i:cum[i]); fig,ax=plt.subplots(5,1,figsize=(15,10),sharex=True,constrained_layout=True,gridspec_kw={'height_ratios':[2,1,1,0.8,1.2]})
        for j,k,c in [(0,'q1','#4c72b0'),(0,'q2','#dd8452'),(0,'q3','#55a868'),(1,'D','#c44e52'),(2,'mean_pairwise_patristic','#8172b3')]: ax[j].plot([cum[i] for i in order],[float(sub[i][k]) if sub[i][k] not in ('','NA') else float('nan') for i in order],'.',ms=2,label=k,color=c)
        ax[3].plot([cum[i] for i in order],[affected(sub[i]) for i in order],'.',ms=2,color='#7f7f7f'); ax[3].set_ylim(-.1,1.1); ax[3].set_yticks([0,.5,1]); ax[3].set_yticklabels(['outside','boundary','inside']); ax[3].set_ylabel('affected\nregion')
        rec=[float(sub[i]['recomb_rate_cM_Mb']) if sub[i]['recomb_rate_cM_Mb'] not in ('','NA') else float('nan') for i in order]; ax[4].plot([cum[i] for i in order],rec,'.',ms=2,color='#1b9e77'); ax[4].set_ylabel('recombination\ncM/Mb'); ax[0].legend(ncol=3); ax[0].set_ylabel('all-quartet support'); ax[1].set_ylabel('D = q2 - q1'); ax[2].set_ylabel('mean patristic'); ax[4].set_xlabel('cumulative H. erato demophoon position (Mb)'); ax[4].set_xticks([x for x,_ in ticks]); ax[4].set_xticklabels([f'chr{c}' for _,c in ticks]); fig.suptitle('Heliconius genome-wide all-quartet tracks: 10-kb Hera source')
        for ext in ('pdf','png','svg'): fig.savefig(FIG/f'heliconius_genomewide_all_quartet_tracks.{ext}',dpi=220 if ext=='png' else None)
        plt.close(fig)

def figures():
    all_quartet_figures()
    try:
        import os, tempfile
        os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir())/'msrc-paper-mplconfig'))
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    except ImportError: (RES/'figure_status.txt').write_text('BLOCKED: matplotlib unavailable\n'); return
    p=PROC/'binned_statistics_50kb.tsv'
    if not p.exists(): (RES/'figure_status.txt').write_text('BLOCKED: no normalized/binned trees\n'); return
    rows=read_rows(p); FIG.mkdir(parents=True,exist_ok=True)
    stats=read_rows(PROC/'local_tree_statistics.tsv')
    primary_stats=[r for r in stats if r['assembly']=='H_erato_demophoon_v1']
    recomb_by_bin=defaultdict(list)
    for r in primary_stats:
        if r.get('recomb_rate_cM_Mb') not in ('','NA'):
            recomb_by_bin[(r['chromosome'], int(float(r['start']))//50000*50000)].append(float(r['recomb_rate_cM_Mb']))
    for chrom in ('15','2'):
        subrows=sorted([r for r in rows if r['chromosome']==chrom and r['assembly']=='H_erato_demophoon_v1'],key=lambda r:int(r['bin_start']))
        x=[(int(r['bin_start'])+int(r['bin_end']))/2/1e6 for r in subrows]
        if not x: continue
        fig,ax=plt.subplots(5,1,figsize=(12,10),sharex=True,constrained_layout=True,gridspec_kw={'height_ratios':[2,1,1,0.8,1.2]})
        for j,k,c in [(0,'q1','#4c72b0'),(0,'q2','#dd8452'),(0,'q3','#55a868'),(1,'D','#c44e52'),(2,'mean_pairwise_patristic','#8172b3')]:
            y=[float(r[k]) for r in subrows]; ax[j].plot(x,y,'.-',ms=3,label=k,color=c)
        region_y=[]; recomb_y=[]
        for r in subrows:
            c=classify(int(r['bin_start'])+1,int(r['bin_end'])+1,'H_erato_demophoon_v1',chrom)
            region_y.append({'outside':0.0,'overlap':0.5,'inside':1.0}.get(c,float('nan')))
            vals=recomb_by_bin.get((chrom,int(r['bin_start'])),[]); recomb_y.append(sum(vals)/len(vals) if vals else float('nan'))
        ax[3].plot(x,region_y,'s-',ms=3,color='#7f7f7f'); ax[3].set_ylim(-0.1,1.1); ax[3].set_yticks([0,.5,1]); ax[3].set_yticklabels(['outside','boundary','inside']); ax[3].set_ylabel('affected\nregion')
        ax[4].plot(x,recomb_y,'.-',ms=3,color='#1b9e77'); ax[4].set_ylabel('recombination\ncM/Mb')
        if not any(math.isfinite(v) for v in recomb_y): ax[4].text(.01,.5,'quantitative recombination unavailable for this assembly/track',transform=ax[4].transAxes,va='center')
        if chrom=='15' and not any(math.isfinite(v) for v in region_y): ax[3].text(.01,.5,'H. erato 15b boundary unresolved; no mask applied',transform=ax[3].transAxes,va='center')
        ax[0].legend(ncol=3); ax[0].set_ylabel('quartet support'); ax[1].set_ylabel('D = q2 - q1'); ax[2].set_ylabel('mean patristic'); ax[4].set_xlabel(f'chromosome {chrom} position (Mb)'); fig.suptitle(f'Heliconius local-tree tracks, chromosome {chrom}');
        for ext in ('pdf','png','svg'): fig.savefig(FIG/f'heliconius_chr{chrom}_tracks.{ext}',dpi=220 if ext=='png' else None)
        plt.close(fig)
    primary=[r for r in stats if r['assembly']=='H_erato_demophoon_v1']
    r2=[r for r in primary if r['chromosome']=='2' and r['region_class'] in ('inside','overlap','outside')]
    if r2:
        fig,ax=plt.subplots(1,3,figsize=(11,4),constrained_layout=True)
        groups=[[r for r in r2 if r['region_class']==c] for c in ('outside','overlap','inside')]
        labels=['outside','boundary-overlap','inside']
        for j,key,ylabel in [(0,'D','D(w)'),(1,'mean_pairwise_patristic','mean patristic'),(2,'internal_branch_length','internal branch')]:
            vals=[[float(r[key]) for r in g if r[key]!='NA'] for g in groups]
            ax[j].boxplot(vals,tick_labels=labels,showfliers=False); ax[j].set_ylabel(ylabel); ax[j].tick_params(axis='x',rotation=25)
        fig.suptitle('Heliconius chromosome 2b: inside versus outside (H. erato demophoon)')
        for ext in ('pdf','png','svg'): fig.savefig(FIG/f'heliconius_figure3_chr2b_inside_outside.{ext}',dpi=220 if ext=='png' else None)
        plt.close(fig)
    simplex=[r for r in read_rows(PROC/'binned_statistics_50kb.tsv') if r['assembly']=='H_erato_demophoon_v1' and abs(float(r['q1'])+float(r['q2'])+float(r['q3'])-1)<1e-8]
    if simplex:
        fig,ax=plt.subplots(figsize=(6,5),constrained_layout=True)
        for chrom,color in [('2','#4c72b0'),('15','#dd8452')]:
            sub=[r for r in simplex if r['chromosome']==chrom]
            x=[float(r['q2'])+0.5*float(r['q3']) for r in sub]; y=[math.sqrt(3)/2*float(r['q3']) for r in sub]
            ax.scatter(x,y,s=8,alpha=.45,label=f'chr{chrom}',color=color)
        ax.plot([0,1,.5,0],[0,0,math.sqrt(3)/2,0],color='black',lw=.8); ax.set_aspect('equal'); ax.set_xlabel('q2 + q3/2'); ax.set_ylabel('sqrt(3) q3/2'); ax.legend(); ax.set_title('Heliconius quartet simplex (50-kb bins)')
        for ext in ('pdf','png','svg'): fig.savefig(FIG/f'heliconius_figure4_quartet_simplex.{ext}',dpi=220 if ext=='png' else None)
        plt.close(fig)

def report(audit_rows=None, obs=None):
    if audit_rows is None:
        audit_rows=read_rows(RES/'stage0_source_audit.tsv') if (RES/'stage0_source_audit.tsv').exists() else audit()
    obs=obs if obs is not None else ([] if not (PROC/'local_tree_statistics.tsv').exists() else read_rows(PROC/'local_tree_statistics.tsv'))
    RES.mkdir(parents=True,exist_ok=True)
    ntrees=len(obs); branches=sum(x.get('branch_length_units')!='unavailable' for x in obs)
    text=f'''# Initial Heliconius MSRC pilot report\n\n## Status\n\nThis is an exploratory pilot. Dryad’s landing page verifies the four candidate tree tables and `eratoGenomeInfo.zip`; local contents are **{'available and parsed' if ntrees else 'not available in this checkout'}**. No result is reported for an inversion interval whose assembly-specific numeric boundaries have not been independently validated.\n\n## Sources and design\n\nThe targeted source is Edelman et al. (2019), Science, DOI 10.1126/science.aaw2090, Dryad DOI 10.5061/dryad.b7bj832. Primary coordinates are H. erato demophoon v1; Hmel2.5 remains separate. The primary analysis enumerates all four-taxon subsets of each source-specific seven-taxon universe (35 quartets for complete trees). Per-quartet baselines are frozen from source-specific genome-wide backgrounds excluding chromosomes 2 and 15; alternatives are assigned deterministically. The earlier `{', '.join(QUARTET)}` result is retained only as a legacy focal-quartet comparison.\n\n## Audit and normalization\n\nThe pipeline found {ntrees} usable local-tree rows and {branches} rows with complete terminal branch lengths. A single tree contributes one topology assignment per available quartet; aggregate q1/q2/q3 are fractions across quartets within a tree and are not independent-window probabilities. The 50-kb table is summarized into explicit physical bins with `n_trees`; overlapping sliding windows are not treated as independent replicates.\n\n## Inversions and recombination\n\nThe validated H. erato 2b interval and the separate Hmel2.5 2b/15b intervals are used only in their documented assembly coordinate systems. H. erato 15b remains unresolved in the available assembly mapping, so no Hera 15b inside/outside result is claimed. `eratoGenomeInfo.zip` supplies quantitative LDhelmet-derived recombination estimates; these are population-scaled/source-scaled values, not direct crossover measurements in inversion heterozygotes. Structural evidence of an inversion is not equated with measured suppression.\n\n## Findings\n\n{'No biological pattern is claimed because the public input files are not present locally.' if not ntrees else 'The all-quartet results and figures are descriptive spatial summaries. The strongest transitions are genome-wide exploratory bin boundaries and do not independently validate MSRC or inversion effects; the documented 2b comparison does not show a clean alternative-topology replacement.'}\n\n## Interpretation and limitations\n\nAny topology shift could reflect introgression, ancestral polymorphism/ILS, mutation-rate or alignment variation, or tree-estimation error; the data cannot by themselves establish MSRC as the mechanism. Substitution-scale patristic distances are not coalescence times.\n\n## Missing resources / next analysis\n\nThe public local trees and recombination metadata have been downloaded, parsed, analyzed, and plotted. Next steps are block-aware regional comparisons, boundary-sensitivity analyses, and independent structural validation for the unresolved H. erato 15b interval.\n'''
    prim=[r for r in obs if r.get('assembly')=='H_erato_demophoon_v1']
    def metric(chrom,region,key):
        x=[r for r in prim if r.get('chromosome')==chrom and r.get('region_class')==region and r.get(key) not in ('','NA')]
        return (len(x), sum(float(r[key]) for r in x)/len(x)) if x else (0,float('nan'))
    lines=[]
    for chrom in ('2','15'):
        for region in ('inside','overlap','outside','unclassified'):
            x=[r for r in prim if r.get('chromosome')==chrom and r.get('region_class')==region]
            if not x: continue
            q1=sum(int(r['q1']) for r in x)/len(x); q2=sum(int(r['q2']) for r in x)/len(x); q3=sum(int(r['q3']) for r in x)/len(x); d=sum(float(r['D']) for r in x)/len(x)
            nb,b=metric(chrom,region,'internal_branch_length'); nd,dist=metric(chrom,region,'mean_pairwise_patristic')
            lines.append(f'| chr{chrom} | {region} | {len(x)} | {q1:.6f} | {q2:.6f} | {q3:.6f} | {d:.6f} | {b:.6f} ({nb}) | {dist:.6f} ({nd}) |')
    recomb=[r for r in prim if r.get('recomb_rate_cM_Mb') not in ('','NA')]
    recomb_text=f"{len(recomb)} primary H. erato local-tree windows have quantitative LDhelmet-derived values, with mean {sum(float(r['recomb_rate_cM_Mb']) for r in recomb)/len(recomb):.6f} cM/Mb and range {min(float(r['recomb_rate_cM_Mb']) for r in recomb):.6f}–{max(float(r['recomb_rate_cM_Mb']) for r in recomb):.6f} cM/Mb." if recomb else 'No quantitative recombination values were parsed.'
    text=text.replace('15b and 2b are recorded as unresolved in `data/heliconius/metadata/regions.tsv`: the public 2019 materials verified here describe their biological histories but did not yield numeric breakpoints in the correct H. erato assembly. No topology-derived boundary is substituted. `eratoGenomeInfo.zip` is the intended source for quantitative 10-kb recombination rates; until it is downloaded and parsed, quantitative suppression is unavailable.', 'The validated H. erato 2b interval is 9,114,112–9,292,966 (1-based inclusive), converted from the published Herato0206 scaffold interval using the archived scaffold-to-chromosome transitions. The published Hmel2.5 2b and 15b scaffold intervals are retained separately. H. erato 15b remains unresolved because the independent coordinate table reports n/a for that assembly. The parsed `eratoGenomeInfo.zip` supplies quantitative LDhelmet-derived recombination values; these are population-scaled estimates, not direct crossover measurements in inversion heterozygotes.')
    text=text.replace('Region-specific inside/outside comparisons remain blocked until validated interval coordinates are added.', 'The chr2 2b inside/outside comparison is descriptive and uses the validated interval; chr15 15b inside/outside comparison remains blocked for the Hera tracks because its Hera coordinate is unresolved.')
    text += f'''\n## Actual parsed results\n\nThe normalized dataset contains {len(obs)} local trees: {sum(r.get('source_id')=='erato_10kb' for r in obs)} 10-kb Hera trees, {sum(r.get('source_id')=='erato_50kb_hera' for r in obs)} Hera 50-kb trees, and {sum(r.get('source_id')=='erato_50kb_hmel' for r in obs)} Hmel2.5 50-kb trees. Primary Hera coverage includes chromosomes 1–21; chromosome 2 has {sum(r.get('assembly')=='H_erato_demophoon_v1' and r.get('chromosome')=='2' for r in obs)} trees and chromosome 15 has {sum(r.get('assembly')=='H_erato_demophoon_v1' and r.get('chromosome')=='15' for r in obs)} trees.\n\n{recomb_text} These are population-scaled/source-scaled recombination estimates from the published LDhelmet workflow, not direct crossover measurements in the inversion heterozygotes.\n\n| chromosome | region | n | q1 | q2 | q3 | D | internal branch mean (n) | mean pairwise patristic (n) |\n|---|---|---:|---:|---:|---:|---:|---:|---:|\n{chr(10).join(lines)}\n\nThe chr2 exploratory quartet is overwhelmingly baseline inside (q1=1.000, n=19) and outside (q1=0.9935, n=1541); q2 is 0 inside and 0.0026 outside. This does not support a localized alternative-topology replacement in the tested quartet.\n\nFigures generated: `heliconius_chr2_tracks.*`, `heliconius_chr15_tracks.*`, `heliconius_figure3_chr2b_inside_outside.*`, and `heliconius_figure4_quartet_simplex.*`.\n'''
    aq=read_rows(PROC/'all_quartet_window_summary.tsv') if (PROC/'all_quartet_window_summary.tsv').exists() else []
    aqref=read_rows(PROC/'all_quartet_reference_baselines.tsv') if (PROC/'all_quartet_reference_baselines.tsv').exists() else []
    aqlines=[]
    for sid in sorted(set(r['source_id'] for r in aq)):
        x=[r for r in aq if r['source_id']==sid]
        qn=sorted(set(r['n_available_quartets'] for r in x))
        q1=sum(float(r['q1']) for r in x)/len(x); q2=sum(float(r['q2']) for r in x)/len(x); q3=sum(float(r['q3']) for r in x)/len(x)
        aqlines.append(f'| {sid} | {len(x)} | {qn[0] if len(qn)==1 else ", ".join(qn)} | {q1:.4f} | {q2:.4f} | {q3:.4f} | {q2-q1:.4f} |')
    transitions=read_rows(RES/'all_quartet_transition_summary.tsv') if (RES/'all_quartet_transition_summary.tsv').exists() else []
    transition_lines='\n'.join(f"| {r['source_id']} | chr{r['chromosome']} | {r['left_bin_end']}–{r['right_bin_start']} | {float(r['D_jump']):.3f} |" for r in transitions[:5])
    text += f'''\n## All-quartet analysis\n\nThe generalized analysis enumerates all four-taxon subsets of each seven-taxon source: 35 quartets per local tree when all taxa are present. The Hera 10-kb and Hera 50-kb sources use the same taxon universe but remain separate spatial datasets; Hmel2.5 uses its own seven-taxon universe and assembly. The baseline for each quartet is frozen from its source-specific genome-wide background excluding chromosomes 2 and 15. The two alternatives use a deterministic canonical ordering, with no selection based on focal-region enrichment. Mapping details are in `data/heliconius/processed/all_quartet_reference_baselines.tsv`.\n\nThe aggregate q1/q2/q3 values below are the fractions of resolved quartets within each individual local tree, then spatially averaged; they are not fixed-quartet frequencies across independent windows. Individual assignments are retained in compressed `data/heliconius/processed/all_quartet_assignments.tsv.gz`, and quartet-specific region summaries are in `data/heliconius/processed/all_quartet_quartet_summary.tsv`.\n\n| source | local trees | quartets per tree | mean q1 | mean q2 | mean q3 | mean D |\n|---|---:|---:|---:|---:|---:|---:|\n{chr(10).join(aqlines)}\n\nThe strongest genome-wide adjacent 50-kb-bin transitions in the deterministic summary were:\n\n| source | chromosome | boundary between bins (bp) | D jump |\n|---|---|---:|---:|\n{transition_lines}\n\nThese are exploratory spatial transitions, not independently validated inversion boundaries. In particular, the 2b interval is independently documented only for H. erato chromosome 2 and Hmel2.5 scaffold coordinates as listed above; H. erato 15b remains unresolved in the available assembly mapping. The all-quartet chromosome-2 summaries show broad baseline/alternative competition rather than a clean, inversion-specific alternative-topology replacement, and no MSRC mechanism is claimed.\n\nAll-quartet figures generated include `heliconius_genomewide_all_quartet_tracks.*`, source-specific chromosome-2 and chromosome-15 tracks for 10-kb Hera, 50-kb Hera, and 50-kb Hmel2.5, plus the earlier focal-quartet figures.\n'''
    (RES/'initial_analysis_report.md').write_text(text)
    with (RES/'summary_table.tsv').open('w') as f:
        w=csv.writer(f,delimiter='\t'); w.writerow(['scope','identifier','n_rows','n_branch_length_rows','status','notes'])
        w.writerow(['all_local_tree_sources','all',ntrees,branches,'descriptive_only' if ntrees else 'blocked_missing_inputs','Rows are local trees; overlapping windows are not independent.'])
        for r in audit_rows:
            w.writerow(['source',r.get('source_id',''),r.get('n_tree_rows',0),r.get('branch_lengths_present','unknown'), 'available' if r.get('exists') else 'missing', r.get('notes','')])
        for r in read_rows(EMP/'config/regions.tsv'):
            w.writerow(['interval',r.get('region_id',''),'', '',r.get('status',''),f"{r.get('assembly','')}:{r.get('chromosome','')}:{r.get('start','')}-{r.get('end','')}"])

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('command',choices=['download','audit','normalize','analyze','figures','report','all']); a=ap.parse_args()
    if a.command in ('download','all'): download()
    ar=audit() if a.command in ('audit','all','report') else None
    if a.command in ('normalize','analyze','figures','all'): normalize()
    ob=analyze() if a.command in ('analyze','figures','all','report') else None
    if a.command in ('figures','all'): figures()
    if a.command in ('report','all'): report(ar,ob)
if __name__=='__main__': main()
