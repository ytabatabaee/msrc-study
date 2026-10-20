#!/usr/bin/env python3
"""Prepare and run Phlag on deterministic matched-taxon representative sets.

Representative choices are made from frozen mappings before any inversion
coordinates are read.  This script never collapses a tree to a consensus.
"""
from __future__ import annotations

import argparse, csv, gzip, hashlib, io, json, random, re, shutil, subprocess, time
from pathlib import Path
from collections import defaultdict

import dendropy
from Bio import Phylo

ROOT=Path(__file__).resolve().parents[4]
OUT=ROOT/"empirical/cross_system_boundary_detection/stage13b/phlag"
RAW=OUT/"raw_results"
PHLAG=Path("/private/tmp/phlag-venv/bin/phlag")
PHLAG_COMMIT="76597d0e4ebd0e507fdafebf7f5163ef634a36d6"
SEED=20261020

def read(p):
    with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def write(p,rows,fields):
    p.parent.mkdir(parents=True,exist_ok=True)
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows({k:r.get(k,"") for k in fields} for r in rows)
def sha(p):
    h=hashlib.sha256();
    with p.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()

def label_tree(src,out,leaf_rename=None):
    tree=Phylo.read(str(src),"newick"); labels=[]; counter=1
    for c in tree.find_clades(order="preorder"):
        if c.is_terminal():
            if leaf_rename and c.name in leaf_rename: c.name=leaf_rename[c.name]
            continue
        if c is tree.root: c.name="ROOT"; continue
        c.name=f"N{counter:03d}"; labels.append(c.name); counter+=1
    Phylo.write(tree,str(out),"newick"); out.write_text("[&R] "+out.read_text().strip()+"\n")
    return labels

def tips(text):
    t=Phylo.read(io.StringIO(text),"newick")
    return [n.name for n in t.get_terminals()]

def prune_rename(text, selected):
    t=Phylo.read(io.StringIO(text),"newick"); keep=set(selected)
    for node in list(t.get_terminals()):
        if node.name not in keep:t.prune(node)
    for node in t.get_terminals(): node.name=selected[node.name]
    return t.format("newick").strip()

def fire_choices(rep):
    m=read(ROOT/"data/fire_ants_chr16/processed/stage6c_tree_tip_supplement_match_audit.tsv"); groups=defaultdict(list)
    for r in m:groups[r["group"]].append(r["tree_tip"])
    rng=random.Random(SEED+rep); return {g:rng.choice(sorted(v)) for g,v in sorted(groups.items())}

def anopheles_choices(rep):
    m=read(ROOT/"data/anopheles_2la/processed/stage4r_haplotype_species_map.tsv"); groups=defaultdict(list)
    for r in m:groups[r["species"]].append(r["haplotype_id"])
    rng=random.Random(SEED+1000+rep); return {g:rng.choice(sorted(v)) for g,v in sorted(groups.items())}

def mouse_choices(rep):
    m=read(ROOT/"data/house_mouse_t_complex/metadata/tip_mapping.tsv"); groups=defaultdict(list)
    for r in m:
        if r["t_status"]=="standard_noncarrier":groups[r["population"]].append(r["tree_tip"])
        if r["t_status"]=="outgroup_not_t_haplotype":groups["SPRE"].append(r["tree_tip"])
    rng=random.Random(SEED+2000+rep); return {g:rng.choice(sorted(v)) for g,v in sorted(groups.items())}

def prepare(dataset,rep):
    base=RAW/"matched_inputs"/dataset/f"replicate_{rep:03d}"; base.mkdir(parents=True,exist_ok=True)
    if dataset=="fire_ants_chr16":
        tree_src=ROOT/"empirical/fire_ants_chr16/results/astral4/stage6_astral4_background_reference.nwk"; local=ROOT/"data/fire_ants_chr16/processed/stage1_window_trees.tsv"; choices=fire_choices(rep); group_labels={"invicta/macdonaghi_SB":"inv_mac_SB","invicta/macdonaghi_Sb":"inv_mac_Sb","richteri_SB":"richteri_SB","richteri_Sb":"richteri_Sb","geminata":"geminata","saevissima":"saevissima","pusillignis":"pusillignis"}; rows=read(local); trees=[(r["window_index"],r["start"],r["end"],r["tree_newick"]) for r in rows];
    elif dataset=="anopheles_2la":
        tree_src=ROOT/"empirical/anopheles_2la/results/astral6/T_all6_species.nwk"; local=ROOT/"data/anopheles_2la/processed/stage4r_all72_haplotype_local_trees.nwk.gz"; choices=anopheles_choices(rep); rows=[]; 
        with gzip.open(local,"rt") as f:
            for line in f:
                parts=line.strip().split("\t",1)
                if len(parts)==2:rows.append((parts[0],parts[1]))
        trees=[(x[0],"","",x[1]) for x in rows]
    else:
        tree_src=ROOT/"empirical/house_mouse_t_complex/results/stage2_population_astral/P1_ALL.nwk"; local=ROOT/"data/house_mouse_t_complex/processed/house_mouse_t_complex_ml_5kb.tre"; choices=mouse_choices(rep); rows=[]; 
        with local.open() as f:
            for i,line in enumerate(f):trees=[(str(i+1),"","",line.strip()) for i,line in enumerate(f)] if False else None
        trees=[(str(i+1),"","",line.strip()) for i,line in enumerate(local.open())]
    safe_leaf_names={}
    if dataset=="fire_ants_chr16":
        # DendroPy's default taxon namespace is case-insensitive.  The
        # frozen group names SB and Sb therefore collide; use stable, neutral
        # lineage labels for the Phlag input while retaining the biological
        # group mapping in representatives.tsv.
        safe_leaf_names={"inv_mac_SB":"inv_mac_group1","inv_mac_Sb":"inv_mac_group2",
                         "richteri_SB":"richteri_group1","richteri_Sb":"richteri_group2"}
    species=base/"species_tree.nwk"; labels=label_tree(tree_src,species,safe_leaf_names)
    if dataset=="fire_ants_chr16":
        label_by_old={tip:safe_leaf_names.get(group_labels[group],group_labels[group]) for group,tip in choices.items()}
    elif dataset=="anopheles_2la":
        # Local labels include a short species suffix after '|'.
        first_tips=tips(trees[0][3]); label_by_old={next(x for x in first_tips if x.startswith(hap+"|")):species_name for species_name,hap in choices.items()}
    else:
        label_by_old={tip:("Czech_Republic" if population=="Czech Republic" else population) for population,tip in choices.items()}
    selected=set(label_by_old); gene=base/"gene_trees.nwk"; mapping=[]
    with gene.open("w") as out:
        for idx,(wid,start,end,text) in enumerate(trees):
            present=set(tips(text)); chosen={k:v for k,v in label_by_old.items() if k in present}
            if len(chosen)!=len(label_by_old): continue
            out.write(prune_rename(text,chosen)+"\n"); mapping.append({"tree_index":len(mapping),"window_id":wid,"start":start,"end":end})
    write(base/"window_mapping.tsv",mapping,["tree_index","window_id","start","end"]); write(base/"representatives.tsv",[{"dataset":dataset,"replicate":rep,"reference_taxon":v,"source_tip":k} for k,v in label_by_old.items()],["dataset","replicate","reference_taxon","source_tip"])
    # Verify exact set equality using the first written tree.
    local_tips=set(tips(gene.read_text().splitlines()[0])); species_tips=set(tips(species.read_text()));
    if local_tips!=species_tips: raise RuntimeError(f"matched taxon failure {dataset} rep {rep}: {local_tips ^ species_tips}")
    # Phlag collapses the basal bifurcation before resolving labels.  The
    # first non-root label is the basal child in these rooted reference trees
    # and is therefore intentionally excluded from the eligible focal-edge
    # list, matching the upstream preprocessing rather than inventing an edge.
    eligible=labels[1:] if labels else labels
    return {"dataset":dataset,"replicate":rep,"species_tree":species,"gene_trees":gene,"mapping":base/"window_mapping.tsv","n_windows":len(mapping),"focal_edges":eligible,"representatives":base/"representatives.tsv"}

def parse_output(path,mapping):
    data=[x.strip() for x in path.read_text(errors="replace").splitlines() if x.strip() and not x.startswith("#") and "," in x]
    state=next((x for x in data if all(y in {"0","1","nan","-1"} for y in x.split(","))),""); probs=data[data.index(state)+1].split(",") if state and data.index(state)+1<len(data) else []
    states=state.split(",") if state else []; finite=all(x.lower()!="nan" for x in probs) and len(probs)==len(mapping)
    cls="numerically_degenerate" if not finite else ("valid_no_detection" if len(set(states))<=1 else "valid_informative")
    return cls,states,probs

def run_one(info,dataset,rep,focal):
    raw=RAW/"runs"/dataset/f"replicate_{rep:03d}";raw.mkdir(parents=True,exist_ok=True); out=raw/f"{focal}.txt"; cmd=[str(PHLAG),"-s",str(info["species_tree"]),"-g",str(info["gene_trees"]),"-e",focal,"-o",str(out)]
    t=time.time(); cp=subprocess.run(cmd,text=True,capture_output=True); (raw/f"{focal}.stdout").write_text(cp.stdout);(raw/f"{focal}.stderr").write_text(cp.stderr); runtime=time.time()-t
    if cp.returncode:return {"dataset":dataset,"replicate":rep,"focal_branch":focal,"status":"execution_failure","returncode":cp.returncode,"runtime_seconds":runtime,"n_windows":info["n_windows"]}
    mapping=read(info["mapping"]); cls,states,probs=parse_output(out,mapping); return {"dataset":dataset,"replicate":rep,"focal_branch":focal,"status":cls,"returncode":0,"runtime_seconds":runtime,"n_windows":len(mapping),"n_state1":sum(x=="1" for x in states),"posterior_finite":str(all(x.lower()!="nan" for x in probs)).lower(),"posterior_max":max((float(x) for x in probs if x.lower()!="nan"),default=float("nan")),"output":str(out)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--datasets",default="anopheles_2la,fire_ants_chr16,house_mouse_t_complex");ap.add_argument("--replicates",type=int,default=1);ap.add_argument("--run",action="store_true");ap.add_argument("--max-focal-edges",type=int,default=1);args=ap.parse_args()
    rows=[]; runs=[]
    for dataset in args.datasets.split(","):
        for rep in range(args.replicates):
            info=prepare(dataset,rep); rows.append({"dataset":dataset,"replicate":rep,"n_windows":info["n_windows"],"n_focal_edges":len(info["focal_edges"]),"focal_edges":",".join(info["focal_edges"]),"representatives":str(info["representatives"]),"taxon_set_valid":"true"})
            if args.run:
                for focal in info["focal_edges"][:args.max_focal_edges]: runs.append(run_one(info,dataset,rep,focal))
    write(OUT/"phlag_sampling_manifest.tsv",rows,["dataset","replicate","n_windows","n_focal_edges","focal_edges","representatives","taxon_set_valid"])
    write(OUT/"phlag_run_manifest.tsv",runs,list(runs[0]) if runs else ["dataset","replicate","focal_branch","status"])
    write(OUT/"phlag_numerical_diagnostics.tsv",runs,list(runs[0]) if runs else ["dataset","replicate","focal_branch","status"])
    (OUT/"phlag_manifest.json").write_text(json.dumps({"stage":"13b","phlag_commit":PHLAG_COMMIT,"seed":SEED,"datasets":args.datasets.split(","),"replicates":args.replicates,"run":args.run,"max_focal_edges":args.max_focal_edges,"representative_policy":"random eligible tip per reference taxon, fixed across windows; no inversion coordinates used"},indent=2)+"\n")

if __name__=="__main__":main()
