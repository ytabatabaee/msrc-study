#!/usr/bin/env python3
"""Extract and match Matschiner et al. ED Fig. 3 source-data tracks."""
import csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[5]
SRC=ROOT/"data/atlantic_cod/raw/nature_source_data/41559_2022_1661_MOESM10_ESM_source_data_EDFig3.txt"
OUT=ROOT/"empirical/atlantic_cod/analysis/joint_genomic_tracks/results"
TREE=ROOT/"data/atlantic_cod/processed/cod_window_trees.tsv"
section=None; rows=[]
for line in SRC.read_text().splitlines():
    if line.startswith("### Fst and dxy (gadMor2)"): section="gadMor2"; continue
    if line.startswith("### Fst and dxy (gadMor3)"): section="gadMor3"; continue
    if not section or not line.strip() or line.startswith(("lg ","chr ")): continue
    p=line.split()
    if len(p)!=16: continue
    chrom=p[0].lower()
    if section=="gadMor2": lg=chrom.upper()
    else: lg=f"LG{int(chrom[3:]):02d}"
    rows.append(dict(reference=section,lg=lg,window_center=int(float(p[1])),num_snps=int(p[2]),gc_content=p[3],prop_var1=p[4],prop_var2=p[5],prop_var1_2=p[6],prop_het1=p[7],prop_het2=p[8],prop_het1_2=p[9],pi_snp_ascertained=p[10],f_st=p[11],d_xy_physical_window_snp_ascertained=p[12],d_xy_per_snp_record=p[13],d_f=p[14],afd=p[15]))
OUT.mkdir(parents=True,exist_ok=True)
with (OUT/"atlantic_cod_published_edfig3_stats.tsv").open("w") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
with TREE.open() as f: tree=list(csv.DictReader(f,delimiter="\t"))
by={(r["lg"],int(r["midpoint"])):r for r in tree}; matched=[]
for r in rows:
    if r["reference"]!="gadMor2":continue
    t=by.get((r["lg"],r["window_center"]))
    matched.append(dict(lg=r["lg"],window_center=r["window_center"],source_start=r["window_center"]-124999,source_end=r["window_center"]+125000,tree_window_id=t["window_id"] if t else "",tree_start=t["start"] if t else "",tree_end=t["end"] if t else "",matched=bool(t)))
with (OUT/"atlantic_cod_published_edfig3_gadMor2_window_match.tsv").open("w") as f:
    w=csv.DictWriter(f,fieldnames=list(matched[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(matched)
print(f"rows={len(rows)} gadMor2={sum(r['reference']=='gadMor2' for r in rows)} gadMor3={sum(r['reference']=='gadMor3' for r in rows)} gadMor2_matched={sum(r['matched'] for r in matched)} gadMor2_unmatched={sum(not r['matched'] for r in matched)}")
