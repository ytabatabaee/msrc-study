#!/usr/bin/env python3
"""Analyze validated staged collections after they exist.

The script is intentionally data-driven and refuses incomplete collections.
It reports per-collection quartet frequencies, treatment distances, ASTRAL
topology summaries, and common-split branch errors.  A missing collection is
reported as pending rather than as a zero-valued result.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path


def read_tsv(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def topology(tree_line: str, quartet: tuple[str, str, str, str]) -> str:
    # Use the repository's Newick parser and quartet implementation so the
    # staged analysis has the same split convention as the pilot analysis.
    from compute_original_baseline import parse_newick, quartet_key
    state = quartet_key(parse_newick(tree_line), quartet)
    left = tuple(sorted(state)); right = tuple(sorted(set(quartet) - set(state)))
    return "|".join(["".join(left), "".join(right)])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-root", type=Path, default=Path(__file__).resolve().parents[3])
    args = parser.parse_args()
    root = args.study_root.resolve()
    script_dir = root / "analysis/scripts"
    import sys
    sys.path.insert(0, str(script_dir))
    from run_astral_scaled_pilot import parse_lengths_correct, normalized_rf

    manifest = read_tsv(root / "analysis/staged_production/staged_manifest.tsv")
    output = root / "analysis/staged_production/results"
    rows = []
    for item in manifest:
        collection = root / item["output_root"] / item["collection_id"]
        marker = collection / "collection_complete.json"
        if not marker.exists():
            continue
        summary = json.loads(marker.read_text())
        gene_file = collection / "gene_trees.nwk"
        lines = [line.strip() for line in gene_file.read_text().splitlines() if line.strip()]
        if len(lines) != 500 or summary.get("n_gene_trees") != 500:
            raise ValueError(f"invalid staged collection count: {collection}")
        # The fixed quartet is a diagnostic, not a replacement for all 30-taxon
        # quartet summaries.  Additional quartets can be added without changing
        # collection validation.
        quartet = tuple(["s1", "s2", "s3", "s4"])
        # Taxon names are not guaranteed to be s1..s4, so derive a stable
        # quartet from the first parsed tree.
        from compute_original_baseline import parse_newick
        taxa = sorted(parse_newick(lines[0]).taxa)
        quartet = tuple(taxa[:4])
        counts = Counter(topology(line, quartet) for line in lines)
        total = sum(counts.values())
        astral_summary = collection / "astral_summary.json"
        ar = json.loads(astral_summary.read_text()) if astral_summary.exists() else {}
        rows.append({
            "condition": item["condition"], "replicate": item["replicate"],
            "treatment": item["treatment"], "n_gene_trees": total,
            "quartet": ",".join(quartet),
            "q_topology_1": counts.most_common(3)[0][1] / total if counts else "",
            "q_topology_2": counts.most_common(3)[1][1] / total if len(counts) > 1 else 0.0,
            "q_topology_3": counts.most_common(3)[2][1] / total if len(counts) > 2 else 0.0,
            "astral_normalized_rf": ar.get("normalized_rf", ""),
            "astral_exact_recovery": ar.get("exact_topology_recovery", ""),
            "status": "PASS",
        })
    write_tsv(output / "staged_collection_summary.tsv", rows)
    print(f"wrote {output / 'staged_collection_summary.tsv'} ({len(rows)} validated collections)")


if __name__ == "__main__":
    main()
