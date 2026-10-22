#!/usr/bin/env python3
"""Recompute branch-length diagnostics from scaled species trees and ASTRAL outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.study_root.resolve()
    out = args.output or root / "analysis/staged_production/results/branch_length_comparison.tsv"
    # Delegate parsing to the already tested project parser.  This script is
    # intentionally a thin, reproducible wrapper around the validated source
    # and ASTRAL formats; it does not compare raw SimPhy lengths directly.
    import sys
    sys.path.insert(0, str(root / "analysis/scripts"))
    from run_astral_scaled_pilot import parse_lengths_correct

    manifest = read_tsv(root / "analysis/staged_production/staged_manifest.tsv")
    rows = []
    for item in manifest:
        collection = root / item["output_root"] / item["collection_id"]
        astral = collection / "astral_species_tree.nwk"
        if not astral.exists():
            continue
        source = root / item["species_tree_path"]
        taxa, true_splits, true_lengths = parse_lengths_correct(source.read_text())
        inferred_taxa, inferred_splits, inferred_lengths = parse_lengths_correct(astral.read_text())
        if taxa != inferred_taxa:
            raise ValueError(f"taxon mismatch in {collection}")
        n = int(item["scaled_effective_population_size"])
        for split in sorted(true_splits & inferred_splits, key=lambda s: tuple(sorted(s))):
            true_cu = true_lengths[split] / (2.0 * n)
            est = inferred_lengths[split]
            rows.append({
                "condition": item["condition"], "replicate": item["replicate"],
                "treatment": item["treatment"], "split": ",".join(sorted(split)),
                "true_culength": true_cu, "estimated_culength": est,
                "signed_error_cu": est - true_cu, "absolute_error_cu": abs(est - true_cu),
                "relative_error": (est - true_cu) / true_cu if true_cu else "",
            })
    if rows:
        write_tsv(out, rows)
    print(f"wrote {out} ({len(rows)} matching splits)")


if __name__ == "__main__":
    main()
