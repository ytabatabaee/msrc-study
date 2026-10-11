#!/usr/bin/env python3
"""Estimate quartet margins and candidate epsilon reversals from source pools."""

from __future__ import annotations

import csv
import random
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from compute_original_baseline import parse_newick, quartet_key  # noqa: E402

OUT = ROOT / "analysis/results/mai30_epsilon_quartets.tsv"
SELECTED = [("007", "outgroup.0.species.1.5.genes.1.5"), ("042", "outgroup.0.species.1.5.genes.1.5")]
TREATMENTS = ("msrc_no_suppression", "msrc_moderate_suppression", "msrc_strong_suppression")


def keys_for(path: Path, quartets: list[tuple[str, str, str, str]]) -> list[list[int]]:
    trees = [parse_newick(x) for x in path.read_text().splitlines() if x.strip()]
    return [[hash(quartet_key(tree, q)) for q in quartets] for tree in trees]


def main() -> None:
    rng = random.Random(20261011)
    rows = []
    for rep, condition in SELECTED:
        species_path = ROOT / "datasets/raw" / condition / rep / "s_tree.trees"
        species = parse_newick(species_path.read_text())
        taxa = sorted(species.taxa)
        all_quartets = []
        for i in range(len(taxa)):
            for j in range(i + 1, len(taxa)):
                for k in range(j + 1, len(taxa)):
                    for l in range(k + 1, len(taxa)):
                        all_quartets.append((taxa[i], taxa[j], taxa[k], taxa[l]))
        quartets = rng.sample(all_quartets, min(2000, len(all_quartets)))
        msc_path = ROOT / "datasets/processed/mai30_scaled_msc_controls" / condition / rep / "msc_gene_trees.nwk"
        msc_keys = keys_for(msc_path, quartets)
        msc_history = ["MSC_direct"] * len(msc_keys)
        species_keys = [hash(quartet_key(species, q)) for q in quartets]
        for treatment in TREATMENTS:
            group = f"{condition}__replicate_{rep}__{treatment}__freq_0.5"
            msrc_dir = ROOT / "datasets/processed/mai30_scaled_wf_pilot" / group
            if not (msrc_dir / "collection_summary.json").exists():
                continue
            if len([x for x in (msrc_dir / "gene_trees.nwk").read_text().splitlines() if x.strip()]) != 500:
                continue
            msrc_keys = keys_for(msrc_dir / "gene_trees.nwk", quartets)
            with (msrc_dir / "locus_history_mapping.tsv").open() as handle:
                msrc_history = [r["history_id"] for r in csv.DictReader(handle, delimiter="\t")]
            for q_index, quartet in enumerate(quartets):
                species_key = species_keys[q_index]
                def counts(matrix, history):
                    observed = sorted({row[q_index] for row in matrix if row[q_index] != species_key})
                    if len(observed) > 2:
                        raise ValueError("more than three quartet topologies")
                    grouped = defaultdict(lambda: [0, 0, 0])
                    flat = [0, 0, 0]
                    for row, hid in zip(matrix, history):
                        key = row[q_index]
                        if key == species_key: idx = 0
                        elif key == observed[0]: idx = 1
                        else: idx = 2
                        grouped[hid][idx] += 1; flat[idx] += 1
                    return np.asarray(flat) / sum(flat), grouped
                p0, _ = counts(msc_keys, msc_history)
                p1, grouped = counts(msrc_keys, msrc_history)
                delta0 = p0[0] - max(p0[1:])
                delta1 = p1[0] - max(p1[1:])
                threshold = ""
                if delta0 * delta1 < 0:
                    threshold = -delta0 / (delta1 - delta0)
                rows.append({
                    "replicate": rep, "condition": condition, "treatment": treatment,
                    "quartet": ",".join(quartet), "delta_msc": delta0, "delta_msrc": delta1,
                    "candidate_epsilon_threshold": threshold,
                    "msc_species_dominant": int(delta0 > 0), "msrc_species_dominant": int(delta1 > 0),
                    "endpoint_reversal": int(delta0 * delta1 < 0),
                    "n_quartets_sampled": len(quartets), "n_msrc_histories": len(set(msrc_history)),
                })
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    print(f"wrote {OUT} ({len(rows)} quartet-treatment rows)")


if __name__ == "__main__":
    main()
