#!/usr/bin/env python3
"""Build coordinate and haploblock inventories from the published Stage 0 files."""
from __future__ import annotations
import csv
import re
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data/sunflowers/raw"
RESULTS = ROOT / "empirical/sunflowers/results"
COORD_RE = re.compile(r"Ha412HOChr(\d+)-(\d+)-(\d+)")
TIP_RE = re.compile(r"(?:^|[,(])([^(),:]+):")
BRANCH_RE = re.compile(r":(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?")

def read_trees():
    trees = []
    for line in (RAW / "10.2kbphysical.trees.txt").read_text().splitlines():
        key, newick = line.split("\t", 1)
        m = COORD_RE.fullmatch(key)
        if not m:
            raise ValueError(f"Unrecognized tree coordinate key: {key}")
        chrom, start, end = map(int, m.groups())
        if end <= start or not newick.endswith(";") or newick.count("(") != newick.count(")"):
            raise ValueError(f"Invalid coordinate or Newick record: {key}")
        tips = TIP_RE.findall(newick)
        branches = BRANCH_RE.findall(newick)
        trees.append({"tree_id": key, "chrom": chrom, "start": start, "end": end,
                      "length_bp": end-start, "tip_count": len(tips),
                      "branch_length_count": len(branches),
                      "all_branch_lengths_present": bool(branches) and len(branches) >= len(tips)})
    return trees

def write_tree_inventory(trees):
    path = RESULTS / "tree_coordinate_inventory.tsv"
    cols = ["tree_id", "chromosome", "start_bp", "end_bp", "window_length_bp", "tip_count", "branch_length_count", "all_branch_lengths_present"]
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(cols)
        for t in trees:
            w.writerow([t["tree_id"], f"Ha412HOChr{t['chrom']:02d}", t["start"], t["end"], t["length_bp"], t["tip_count"], t["branch_length_count"], str(t["all_branch_lengths_present"]).lower()])

def write_haploblock_inventory(trees):
    source = RESULTS / "haploblocks_reported.tsv"
    out = RESULTS / "haploblock_inventory.tsv"
    with source.open() as f: rows = list(csv.DictReader(f, delimiter="\t"))
    cols = list(rows[0]) + ["suppression_evidence", "structural_evidence", "evolutionary_age", "interval_precision", "overlapping_tree_count"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n"); w.writeheader()
        for row in rows:
            regions = []
            for part in row["reported_regions_mbp"].split(","):
                left, right = map(float, part.strip().split("-"))
                regions.append((int(left*1_000_000), int(right*1_000_000)))
            overlap = sum(1 for t in trees if t["chrom"] == int(row["chromosome"]) and any(t["start"] < b and t["end"] > a for a,b in regions))
            codes = {x.strip() for x in row["support_codes"].split(",") if x.strip()}
            row.update(suppression_evidence="genetic-map suppression reported" if "g" in codes else "not marked g in Extended Data Table 1",
                       structural_evidence=",".join(sorted(codes - {"g"})) or "none recorded in table",
                       evolutionary_age="not resolved per block in table; paper reports 32/37 older than annuus-argophyllus split",
                       interval_precision="rounded Mbp; some blocks comprise discontiguous pieces",
                       overlapping_tree_count=overlap)
            w.writerow(row)

def main():
    RESULTS.mkdir(parents=True, exist_ok=True)
    trees = read_trees()
    write_tree_inventory(trees)
    write_haploblock_inventory(trees)
    with tarfile.open(RAW / "ind_trees.tar.gz", "r:gz") as tf:
        archive_ids = {m.name.rsplit("/", 1)[-1][:-len(".treefile")]
                       for m in tf.getmembers() if m.name.endswith(".treefile")}
    tree_ids = {t["tree_id"] for t in trees}
    if archive_ids != tree_ids:
        raise ValueError(f"Tree archive/table ID mismatch: archive-only={len(archive_ids-tree_ids)}, table-only={len(tree_ids-archive_ids)}")
    archived = len(archive_ids)
    manifest = []
    for line in (RAW / "10kb_window_sizes.txt").read_text().splitlines():
        name, called_bases = line.split()
        manifest.append((name[:-6] if name.endswith(".fasta") else name, int(called_bases)))
    retained_manifest = sum(n >= 2000 for _, n in manifest)
    no_tree = sum(n >= 2000 and key not in tree_ids for key, n in manifest)
    print(f"records={len(trees)} archive_treefiles={archived} unique_windows={len(tree_ids)}")
    print(f"tips={min(t['tip_count'] for t in trees)}..{max(t['tip_count'] for t in trees)}; branch_lengths_present={all(t['all_branch_lengths_present'] for t in trees)}")
    print(f"manifest_windows={len(manifest)}; manifest_called_bases_ge_2000={retained_manifest}; qualifying_manifest_windows_without_tree={no_tree}")
    print("outputs=tree_coordinate_inventory.tsv,haploblock_inventory.tsv")

if __name__ == "__main__": main()
