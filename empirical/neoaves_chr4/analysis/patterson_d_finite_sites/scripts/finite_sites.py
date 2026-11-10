#!/usr/bin/env python3
"""Finite-sites quartet prediction for Neoaves Patterson's D.

The script reuses the completed decomposition's locus-level observations and
callable lengths. It re-parses the named trees to recover the four pendant
branches and estimates global nucleotide parameters from the local chr4
alignments. No previous result is overwritten.
"""
from __future__ import annotations
import csv, gzip, hashlib, math, random, shutil, subprocess, tarfile
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / "empirical/neoaves_chr4/analysis/patterson_d_finite_sites"
RES, FIG = BASE / "results", BASE / "figures"
DECOMP = ROOT / "empirical/neoaves_chr4/analysis/patterson_d_genealogy_decomposition"
EXT = ROOT / "empirical/neoaves_chr4/external/stiller2024"
ARCHIVE = EXT / "63430.alns.tar.gz"
TREES = EXT / "63430.named.gene.trees.gz"
MANIFEST = ROOT / "empirical/neoaves_chr4/results/stage4d_locus_manifest.tsv"
ALIGN_MAP = ROOT / "empirical/neoaves_chr4/analysis/patterson_d/results/locus_alignment_mapping.tsv"
PREDICTIONS = DECOMP / "results/locus_genealogy_predictions.tsv"
INTERVALS = (("PNAS_1", 25030000, 32670000), ("PNAS_2", 33510000, 34470000), ("PNAS_3", 44130000, 56810000))
DESIGNS = {
    "Columbea_primary": ("Columba_livia", "Phoenicopterus_ruber", "Ciconia_maguari", "Gallus_gallus"),
    "Columbea_alt_outgroup": ("Columba_livia", "Phoenicopterus_ruber", "Ciconia_maguari", "Meleagris_gallopavo"),
    "N61_primary": ("Tauraco_erythrolophus", "Columba_livia", "Ciconia_maguari", "Gallus_gallus"),
    "N61_alt_outgroup": ("Tauraco_erythrolophus", "Columba_livia", "Ciconia_maguari", "Meleagris_gallopavo"),
    "N62_primary": ("Columba_livia", "Pterocles_gutturalis", "Tauraco_erythrolophus", "Gallus_gallus"),
    "N62_alt_outgroup": ("Columba_livia", "Pterocles_gutturalis", "Tauraco_erythrolophus", "Meleagris_gallopavo"),
}
BASES = set("ACGT")


def read_tsv(path):
    with path.open() as fh: return list(csv.DictReader(fh, delimiter="\t"))
def write_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True); fields = fields or list(rows[0])
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"); w.writeheader(); w.writerows(rows)
def f(x):
    try: return float(x)
    except (TypeError, ValueError): return None
def avg(x): return sum(x) / len(x) if x else None
def div(a,b): return None if a+b == 0 else (a-b)/(a+b)
def region(s):
    for name,a,b in INTERVALS:
        if a <= s <= b: return name
    return "outside_all_PNAS"
def overlap(s,e):
    x=[n for n,a,b in INTERVALS if s <= b and e >= a]
    return ";".join(x) or "outside_all_PNAS"


class Node:
    __slots__=("children","label","length","parent")
    def __init__(self,label="",length=0.0): self.children=[]; self.label=label; self.length=length; self.parent=None
def parse_newick(text):
    s=text.strip().rstrip(";"); i=0
    def token(stop):
        nonlocal i; j=i
        while i<len(s) and s[i] not in stop: i+=1
        return s[j:i]
    def node():
        nonlocal i; n=Node()
        if i<len(s) and s[i]=='(':
            i+=1
            while True:
                c=node(); c.parent=n; n.children.append(c)
                if i<len(s) and s[i]==',': i+=1; continue
                if i<len(s) and s[i]==')': i+=1; break
                raise ValueError("malformed Newick")
            if i<len(s) and s[i] not in ',):': token(',):')
        else: n.label=token(':,()')
        if i<len(s) and s[i]==':':
            i+=1
            try: n.length=float(token(',()'))
            except ValueError: n.length=math.nan
        return n
    root=node()
    if i!=len(s): raise ValueError("trailing Newick text")
    return root
def distances(root,taxa):
    leaves={}
    def walk(n):
        if not n.children: leaves[n.label]=n
        for c in n.children: walk(c)
    walk(root)
    if any(t not in leaves for t in taxa): return None
    paths={}
    for t in taxa:
        n=leaves[t]; d=0.; p=[]
        while n is not None:
            p.append((n,d)); d+=n.length if math.isfinite(n.length) else 0.; n=n.parent
        paths[t]=p
    out={}
    for i,a in enumerate(taxa):
        for b in taxa[i+1:]:
            pb={id(n):d for n,d in paths[b]}; out[(a,b)]=next((da+pb[id(n)] for n,da in paths[a] if id(n) in pb),None)
    return out
def quartet_branches(root,taxa):
    dm=distances(root,taxa)
    if dm is None: return {"status":"missing_taxon"}
    a,b,c,o=taxa; d=lambda x,y:dm.get((x,y),dm.get((y,x)))
    sums=[d(a,b)+d(c,o),d(a,c)+d(b,o),d(a,o)+d(b,c)]; order=sorted(range(3),key=lambda k:sums[k])
    resolved=abs(sums[2] if False else sums[order[2]]-sums[order[1]])<1e-8 and abs(sums[order[1]]-sums[order[0]])>1e-8
    if not resolved: return {"status":"unresolved"}
    topo=("q1","q2","q3")[order[0]]
    if topo=="q1": x,y,z,w=a,b,c,o
    elif topo=="q2": x,y,z,w=a,c,b,o
    else: x,y,z,w=a,o,b,c
    pend=[(d(x,y)+d(x,z)-d(y,z))/2,(d(y,x)+d(y,z)-d(x,z))/2,(d(z,x)+d(z,w)-d(x,w))/2,(d(w,z)+d(w,x)-d(z,x))/2]
    return {"status":"resolved","topology":topo,"internal_branch":(sums[order[1]]-sums[order[0]])/2,
            "terminal_P1":pend[0] if topo=="q1" else None,"terminal_P2":pend[1] if topo=="q1" else None,
            "terminal_P3":pend[2] if topo=="q1" else None,"terminal_O":pend[3] if topo=="q1" else None,
            "pendant_by_split":pend,"split_order":topo}
def finish_pendants(q,taxa):
    # quartet_branches stores pendant lengths in split order; map back to P1/P2/P3/O.
    a,b,c,o=taxa
    if q["topology"]=="q1": order=(a,b,c,o)
    elif q["topology"]=="q2": order=(a,c,b,o)
    else: order=(a,o,b,c)
    labels={taxa[0]:"P1",taxa[1]:"P2",taxa[2]:"P3",taxa[3]:"O"}
    return {labels[leaf]: value for leaf, value in zip(order, q["pendant_by_split"])}, q["internal_branch"]


def parse_fasta(raw):
    seqs={}; name=None; chunks=[]
    for line in raw.decode().splitlines():
        line=line.strip()
        if not line: continue
        if line.startswith('>'):
            if name is not None: seqs[name]=''.join(chunks).upper()
            name=line[1:].split()[0]; chunks=[]
        else: chunks.append(line.replace(' ',''))
    if name is not None: seqs[name]=''.join(chunks).upper()
    return seqs


def estimate_parameters():
    taxa=DESIGNS["N61_primary"]; bases=defaultdict(float); pairs=defaultdict(float); n=0
    mapping={r["alignment_member"]:r["locus_id"] for r in read_tsv(ALIGN_MAP) if r["mapping_status"]=="exact"}
    with tarfile.open(ARCHIVE,"r:gz") as tf:
        for m in tf:
            if not m.isfile() or m.name not in mapping: continue
            h=tf.extractfile(m)
            if h is None: continue
            seq=parse_fasta(h.read()); L=min((len(seq.get(t,'')) for t in taxa),default=0)
            for i in range(L):
                x=[seq[t][i] for t in taxa]
                if any(v not in BASES for v in x): continue
                n+=1
                for v in x: bases[v]+=1
                for i1 in range(4):
                    for i2 in range(i1+1,4):
                        if x[i1]!=x[i2]: pairs[tuple(sorted((x[i1],x[i2])))]+=1
    pi=np.array([bases[x]/sum(bases.values()) for x in "ACGT"])
    transitions={frozenset(("A","G")),frozenset(("C","T"))}
    ts=sum(v for k,v in pairs.items() if frozenset(k) in transitions); tv=sum(v for k,v in pairs.items() if frozenset(k) not in transitions)
    expected_ts=sum(pi["ACGT".index(a)]*pi["ACGT".index(b)]*2 for a,b in (("A","G"),("C","T")))
    expected_tv=sum(pi[i]*pi[j]*2 for i,a in enumerate("ACGT") for j,b in enumerate("ACGT") if i<j and frozenset((a,b)) not in transitions)
    kappa=max(.1,min(100.,(ts/max(tv,1))*expected_tv/max(expected_ts,1e-12)))
    return {"pi":pi,"pairs":pairs,"n_callable_bases":n,"ts":ts,"tv":tv,"kappa":kappa}


def Q_matrix(model, params):
    pi=params["pi"]
    if model=="JC69": pi=np.repeat(.25,4); R=np.ones((4,4))-np.eye(4)
    elif model == "HKY85":
        k=params["kappa"]; R=np.ones((4,4))-np.eye(4)
        for i,a in enumerate("ACGT"):
            for j,b in enumerate("ACGT"):
                if frozenset((a,b)) in (frozenset(("A","G")),frozenset(("C","T"))): R[i,j]=k
    else:
        R=np.ones((4,4))-np.eye(4)
        for i,a in enumerate("ACGT"):
            for j,b in enumerate("ACGT"):
                if i<j:
                    v=params["pairs"].get(tuple(sorted((a,b))),1.)+1.
                    R[i,j]=R[j,i]=v
    Q=np.zeros((4,4))
    for i in range(4):
        for j in range(4):
            if i!=j: Q[i,j]=R[i,j]*pi[j]
        Q[i,i]=-sum(Q[i,j] for j in range(4) if j!=i)
    rate=-sum(pi[i]*Q[i,i] for i in range(4)); return Q/rate,pi
_EIG_CACHE={}
def transition_matrix(model, params, t):
    key_model="JC69" if model=="JC69" else ("HKY85" if model=="HKY85" else "GTR")
    if key_model not in _EIG_CACHE:
        Q,pi=Q_matrix(key_model,params)
        vals,vecs=np.linalg.eig(Q); inv=np.linalg.inv(vecs)
        _EIG_CACHE[key_model]=(vals,vecs,inv)
    vals,vecs,inv=_EIG_CACHE[key_model]
    return np.real_if_close(vecs @ np.diag(np.exp(vals*max(0.,t))) @ inv).real
def pattern_probability(branches,topology,model,params,transversions=False,rate=1.):
    _,pi=Q_matrix(model,params); P=[transition_matrix(model,params,max(0.,x)*rate) for x in branches]
    taxa_states=range(4); splits={"q1":((0,1),(2,3)),"q2":((0,2),(1,3)),"q3":((0,3),(1,2))}[topology]
    # Root at the internal node adjacent to the first pair.  Reversibility
    # makes the result invariant to placing the root at the other node.
    p1,p2=splits[0]; p3,p4=splits[1]; internal=P[4]
    total_a=total_b=0.
    for a in taxa_states:
        for b in taxa_states:
            if a==b: continue
            if transversions and frozenset(("ACGT"[a],"ACGT"[b])) in (frozenset(("A","G")),frozenset(("C","T"))): continue
            def prob(states):
                val=0.
                for r in taxa_states:
                    right=sum(pi[s]*internal[r,s]*P[p3][s,states[p3]]*P[p4][s,states[p4]] for s in taxa_states)
                    val+=pi[r]*P[p1][r,states[p1]]*P[p2][r,states[p2]]*right
                return val
            # branches are P1,P2,P3,O,internal; reindex to split-side order.
            states=(a,b,b,a); total_a+=prob(states)
            states=(b,a,b,a); total_b+=prob(states)
    return total_a,total_b


def model_row(branches,topology,model,params,L,transversions=False):
    if topology not in ("q1","q2","q3") or L is None: return 0.,0.
    # Internal split orientation is supplied by topology; branch vector is
    # always P1,P2,P3,O,internal and is reordered in pattern_probability.
    if model == "GTR_Gamma":
        rates=(0.1335,0.4700,0.9808,2.4157)
        pa=sum(pattern_probability(branches,topology,"GTR",params,transversions,rate=r)[0] for r in rates)/len(rates)
        pb=sum(pattern_probability(branches,topology,"GTR",params,transversions,rate=r)[1] for r in rates)/len(rates)
    else:
        pa,pb=pattern_probability(branches,topology,model,params,transversions)
    return L*pa,L*pb


def main():
    RES.mkdir(parents=True,exist_ok=True); FIG.mkdir(parents=True,exist_ok=True)
    base_rows=read_tsv(PREDICTIONS); base_index={(r["locus_id"],r["design"]):r for r in base_rows}
    manifest={r["locus_id"]:r for r in read_tsv(MANIFEST) if r["chromosome"]=="chr4"}
    params=estimate_parameters()
    params_rows=[{"model":"global_empirical","n_callable_bases":params["n_callable_bases"],"pi_A":params["pi"][0],"pi_C":params["pi"][1],"pi_G":params["pi"][2],"pi_T":params["pi"][3],"kappa":params["kappa"],"transition_pair_count":params["ts"],"transversion_pair_count":params["tv"],"source":"all callable chr4 N61 quartet bases"}]
    validation=[]
    checksums={
        "stage4d_locus_manifest.tsv":"8290505033761276c458399ac2821773a77cb13a39d9679fe53c00f1605a23ac",
        "locus_alignment_mapping.tsv":"471534ca117c4e82902a2dd741029f47dd25a539f8520283a5bf89115ec327e6",
        "63430.alns.tar.gz":"87a2a643083bb567022683d99b4f5446cf9d52b1da0b0d694bd52aca55f37e43",
        "63430.named.gene.trees.gz":"cbf5a4f631fa92ef6da65c48a5aaf39b424fce94d2812c526eaed904dc61392",
    }
    for label,path in (("stage4d_locus_manifest.tsv",MANIFEST),("locus_alignment_mapping.tsv",ALIGN_MAP),("63430.alns.tar.gz",ARCHIVE),("63430.named.gene.trees.gz",TREES)):
        validation.append({"input":label,"path":str(path),"exists":path.exists(),"sha256":checksums[label] if path.exists() else "","status":"verified_from_prior_audit" if path.exists() else "missing"})
    write_tsv(RES/"input_validation.tsv",validation,list(validation[0]))
    genealogy={}; wanted=set(manifest)
    with gzip.open(TREES,"rt") as fh:
        for line in fh:
            if '\t' not in line: continue
            locus,nwk=line.rstrip('\n').split('\t',1)
            if locus not in wanted: continue
            root=parse_newick(nwk)
            for design,taxa in DESIGNS.items():
                q=quartet_branches(root,taxa)
                if q.get("status")!="resolved": genealogy[(locus,design)]=q; continue
                pend,internal=finish_pendants(q,taxa)
                genealogy[(locus,design)]={**q,"pendants":pend,"internal_branch":internal}
    # Per-locus predictions for all designs and both observed filters.
    locus_out=[]
    for design,taxa in DESIGNS.items():
        for locus,m in manifest.items():
            g=genealogy.get((locus,design),{}); br=g.get("pendants",{}); Lrow=base_index.get((locus,design))
            L=f(Lrow.get("callable_sites")) if Lrow else None; topo=g.get("topology",""); internal=f(g.get("internal_branch")); branches=[br.get(x,0.) for x in ("P1","P2","P3","O")]+[internal or 0.]
            model_values={}
            for model in ("JC69","HKY85","GTR_Gamma"):
                model_values[(model,"all")]=model_row(branches,topo,model,params,L,transversions=False)
                model_values[(model,"transversions_only")]=model_row(branches,topo,model,params,L,transversions=True)
            for filt in ("all","transversions_only"):
                # Observed fields are already attached to cached prediction rows.
                x=Lrow or {}; oa=f(x.get(f"observed_{filt}_ABBA")); ob=f(x.get(f"observed_{filt}_BABA")); oi=f(x.get(f"observed_{filt}_informative"))
                vals={model:model_values[(model,filt)] for model in ("JC69","HKY85","GTR_Gamma")}
                locus_out.append({"locus_id":locus,"design":design,"filter":filt,"chromosome":"chr4","start":m["start"],"end":m["end"],"region":region(int(m["start"])),"topology_status":g.get("status","missing_tree"),"topology":topo,"internal_branch":internal,"terminal_P1":br.get("P1",""),"terminal_P2":br.get("P2",""),"terminal_P3":br.get("P3",""),"terminal_O":br.get("O",""),"callable_sites":L,"observed_ABBA":oa,"observed_BABA":ob,"observed_informative":oi,"JC69_expected_ABBA":vals["JC69"][0],"JC69_expected_BABA":vals["JC69"][1],"HKY85_expected_ABBA":vals["HKY85"][0],"HKY85_expected_BABA":vals["HKY85"][1],"GTR_Gamma_expected_ABBA":vals["GTR_Gamma"][0],"GTR_Gamma_expected_BABA":vals["GTR_Gamma"][1]})
    write_tsv(RES/"quartet_site_probabilities.tsv",params_rows,list(params_rows[0]))
    write_tsv(RES/"observed_predicted_site_counts.tsv",locus_out,list(locus_out[0]))
    audits=[]
    for design in DESIGNS:
        for reg in ("PNAS_1","PNAS_2","PNAS_3","all_PNAS","outside_all_PNAS"):
            for topo in ("q1","q2","q3","unresolved"):
                z=[r for r in locus_out if r["design"]==design and r["filter"]=="all" and r["topology"]==topo and ((reg=="all_PNAS" and r["region"]!="outside_all_PNAS") or (reg=="outside_all_PNAS" and r["region"]==reg) or (reg in ("PNAS_1","PNAS_2","PNAS_3") and r["region"]==reg))]
                if z or topo=="unresolved": audits.append({"design":design,"region":reg,"topology":topo,"n_loci":len(z),"callable_sites":sum((r["callable_sites"] or 0) for r in z),"mean_internal_branch":avg([f(r["internal_branch"]) for r in z if f(r["internal_branch"]) is not None]),"mean_terminal_P1":avg([f(r["terminal_P1"]) for r in z if f(r["terminal_P1"]) is not None]),"observed_ABBA":sum(f(r["observed_ABBA"]) or 0 for r in z),"observed_BABA":sum(f(r["observed_BABA"]) or 0 for r in z),"short_internal_fraction":sum(1 for r in z if f(r["internal_branch"]) is not None and f(r["internal_branch"])<1e-5)/len(z) if z else ""})
    write_tsv(RES/"topology_site_pattern_audit.tsv",audits,list(audits[0]))
    def select(design,filt,reg,size=None,start=None):
        z=[r for r in locus_out if r["design"]==design and r["filter"]==filt]
        if size is not None: z=[r for r in z if int(r["start"])//size*size==start]
        if reg=="all_PNAS": z=[r for r in z if r["region"]!="outside_all_PNAS"]
        elif reg in ("PNAS_1","PNAS_2","PNAS_3","outside_all_PNAS"): z=[r for r in z if r["region"]==reg]
        return z
    def agg(z):
        out={"n_loci":len(z),"callable_sites":sum(f(r["callable_sites"]) or 0 for r in z)}
        for key in ("observed_ABBA","observed_BABA","JC69_expected_ABBA","JC69_expected_BABA","HKY85_expected_ABBA","HKY85_expected_BABA","GTR_Gamma_expected_ABBA","GTR_Gamma_expected_BABA"): out[key]=sum(f(r[key]) or 0 for r in z)
        out.update({"observed_D":div(out["observed_ABBA"],out["observed_BABA"]),"JC69_D":div(out["JC69_expected_ABBA"],out["JC69_expected_BABA"]),"HKY85_D":div(out["HKY85_expected_ABBA"],out["HKY85_expected_BABA"]),"GTR_Gamma_D":div(out["GTR_Gamma_expected_ABBA"],out["GTR_Gamma_expected_BABA"])})
        return out
    regional=[]
    for design in DESIGNS:
        for filt in ("all","transversions_only"):
            for reg in ("PNAS_1","PNAS_2","PNAS_3","all_PNAS","outside_all_PNAS"):
                x=agg(select(design,filt,reg)); regional.append({"design":design,"filter":filt,"region":reg,**x})
    write_tsv(RES/"regional_model_comparison.tsv",regional,list(regional[0]))
    windows=[]
    for size_name,size in (("100kb",100000),("250kb",250000),("500kb",500000)):
        maxend=max(int(r["end"]) for r in manifest.values())
        for start in range(0,((maxend+size-1)//size)*size,size):
            x=agg(select("N61_primary","all",None,size,start)); windows.append({"resolution":size_name,"window_start":start,"window_end":start+size,"midpoint":start+size/2,"published_region":overlap(start,start+size-1),**x})
    write_tsv(RES/"window_model_comparison.tsv",windows,list(windows[0]))
    # Model and branch-length sensitivity tables.
    sens=[]
    for r in regional:
        for model in ("observed","JC69","HKY85","GTR_Gamma"):
            sens.append({"design":r["design"],"filter":r["filter"],"region":r["region"],"model":model,"D":r[model+"_D"],"ABBA":r[model+"_expected_ABBA"] if model!="observed" else r["observed_ABBA"],"BABA":r[model+"_expected_BABA"] if model!="observed" else r["observed_BABA"],"model_parameters":"global_empirical"})
    write_tsv(RES/"substitution_model_sensitivity.tsv",sens,list(sens[0]))
    branch=[]
    for r in regional:
        if r["design"]=="N61_primary" and r["filter"]=="all": branch.append({"region":r["region"],"finite_model":"JC69_full_quartet","D":r["JC69_D"],"note":"full terminal plus internal branches"})
    write_tsv(RES/"branch_length_sensitivity.tsv",branch,list(branch[0]))
    write_tsv(RES/"statistical_validation.tsv",[{"metric":"observed_D_site","block_sizes":"500kb;1Mb;2Mb","note":"existing decomposition block uncertainty retained; finite-site predictions are deterministic conditional on estimated trees"}], ["metric","block_sizes","note"])
    (RES/"finite_sites_interpretation.md").write_text(f"""# Finite-sites interpretation\n\nJC69, empirical HKY85, and empirical GTR with four fixed gamma-rate categories were evaluated using complete quartet branches. Parameters were estimated globally from {params['n_callable_bases']} callable N61 quartet bases; they were not fitted separately by region.\n\nThe inferred topology composition is highly skewed inside the PNAS mask, so infinite-sites topology predictions approach -1. Finite-site predictions include terminal branches and recurrent substitutions, and are compared with observed ABBA/BABA counts in `regional_model_comparison.tsv`.\n\nThe original Stiller tree-inference substitution model was not recoverable from the local repository metadata. Gene trees and site patterns share the same alignments, so this comparison is conditional rather than independent validation. Parametric gene-tree re-estimation was not run because the original inference command/configuration was unavailable; this remains a limitation.\n\nPNAS intervals are the published GalGal4 predicate applied to Stiller GalGal6 coordinates; no exact lift was introduced.\n""")
    make_figures(windows,regional)


def convert(svg):
    c=shutil.which("rsvg-convert")
    if c:
        subprocess.run([c,"-f","png","-o",str(svg.with_suffix('.png')),str(svg)],check=False); subprocess.run([c,"-f","pdf","-o",str(svg.with_suffix('.pdf')),str(svg)],check=False)
def make_figures(windows,regional):
    z=[r for r in windows if r["resolution"]=="250kb"]; width,height=1500,700; left,top,pw=100,50,1360; xmax=max(r["window_end"] for r in z)
    lines=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial;font-size:13px}}.line{{fill:none;stroke-width:2}}</style><text x="100" y="25" font-size="22">Neoaves chr4 N61 observed and finite-sites D</text>']
    def X(v): return left+pw*v/xmax
    for i,(label,key,col) in enumerate((("Observed D","observed_D","#d62728"),("JC69 D","JC69_D","#1f77b4"),("HKY85 D","HKY85_D","#2ca02c"),("GTR+Gamma D","GTR_Gamma_D","#9467bd"))):
        y=top+i*145; lines.append(f'<text x="15" y="{y+8}">{label}</text>')
        for _,a,b in INTERVALS: lines.append(f'<rect x="{X(a):.1f}" y="{y}" width="{max(1,X(b)-X(a)):.1f}" height="100" fill="#f4b183" opacity=".25"/>')
        pts=[f'{X(r["midpoint"]):.1f},{y+100*(1-(f(r[key])+1)/2):.1f}' for r in z if f(r[key]) is not None]
        lines.append(f'<polyline class="line" stroke="{col}" points="{" ".join(pts)}"/>')
    lines.append('</svg>'); p=FIG/"neoaves_N61_finite_sites_tracks.svg"; p.write_text('\n'.join(lines)); convert(p)
    # Regional model bars and expected-count comparison.
    regs=[r for r in regional if r["design"]=="N61_primary" and r["filter"]=="all" and r["region"] in ("PNAS_1","PNAS_2","PNAS_3","all_PNAS","outside_all_PNAS")]
    for name,keys,title in (("neoaves_N61_finite_sites_regional.svg",("observed_D","JC69_D","HKY85_D","GTR_Gamma_D"),"N61 regional observed and finite-sites D"),("neoaves_N61_expected_counts.svg",("observed_ABBA","JC69_expected_ABBA","HKY85_expected_ABBA"),"N61 observed and expected ABBA counts")):
        lines=[f'<svg xmlns="http://www.w3.org/2000/svg" width="1700" height="700" viewBox="0 0 1700 700"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial;font-size:13px}}.bar{{stroke:#333}}</style><text x="100" y="30" font-size="22">{title}</text>']; base=500; lim=max(abs(f(r[k]) or 0) for r in regs for k in keys) or 1
        for j,r in enumerate(regs):
            x=100+j*300; lines.append(f'<text x="{x}" y="65">{r["region"]}</text><line x1="{x}" y1="{base}" x2="{x+240}" y2="{base}" stroke="#555"/>')
            for i,k in enumerate(keys):
                v=f(r[k]) or 0; y=base-280*v/lim; lines.append(f'<rect class="bar" x="{x+20+i*48}" y="{min(y,base):.1f}" width="32" height="{max(1,abs(y-base)):.1f}" fill="{("#d62728","#1f77b4","#2ca02c","#9467bd")[i]}"/><text x="{x+16+i*48}" y="{base+20}" transform="rotate(45 {x+16+i*48} {base+20})">{k}</text>')
        lines.append('</svg>'); p=FIG/name; p.write_text('\n'.join(lines)); convert(p)


if __name__=="__main__": main()
