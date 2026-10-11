#!/usr/bin/env python3
"""History-level uncertainty for candidate quartet epsilon thresholds."""

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

IN = ROOT / "analysis/results/mai30_epsilon_quartets.tsv"
OUT = ROOT / "analysis/results/mai30_epsilon_threshold_bootstrap.tsv"


def main() -> None:
    all_rows = list(csv.DictReader(IN.open(), delimiter="\t"))
    output = []
    rng = random.Random(20261011)
    for rep, treatment in sorted({(r["replicate"], r["treatment"]) for r in all_rows}):
        candidates = [r for r in all_rows if r["replicate"] == rep and r["treatment"] == treatment and r["endpoint_reversal"] == "1"][:20]
        if not candidates:
            output.append({"replicate": rep, "treatment": treatment, "candidate_quartets": 0, "bootstrap_reversal_fraction_median": 0, "threshold_q05": "", "threshold_median": "", "threshold_q95": ""})
            continue
        condition = "outgroup.0.species.1.5.genes.1.5"
        msrc_group = f"{condition}__replicate_{rep}__{treatment}__freq_0.5"
        msrc_path = ROOT / "datasets/processed/mai30_scaled_wf_pilot" / msrc_group / "gene_trees.nwk"
        msc_path = ROOT / "datasets/processed/mai30_scaled_msc_controls" / condition / rep / "msc_gene_trees.nwk"
        msrc_trees = [parse_newick(x) for x in msrc_path.read_text().splitlines() if x.strip()]
        msc_trees = [parse_newick(x) for x in msc_path.read_text().splitlines() if x.strip()]
        species = parse_newick((ROOT / "datasets/raw" / condition / rep / "s_tree.trees").read_text())
        with (msrc_path.parent / "locus_history_mapping.tsv").open() as handle:
            histories = [r["history_id"] for r in csv.DictReader(handle, delimiter="\t")]
        by_history = defaultdict(list)
        for tree, history in zip(msrc_trees, histories): by_history[history].append(tree)
        history_ids = sorted(by_history)
        prepared = []
        for row in candidates:
            q = tuple(row["quartet"].split(",")); key = quartet_key(species, q)
            def tally(trees):
                counts = defaultdict(int)
                for tree in trees: counts[quartet_key(tree, q)] += 1
                return key, counts
            prepared.append((q, key, tally(msc_trees), {h: tally(by_history[h])[1] for h in history_ids}))
        boot_thresholds = []; boot_reversals = []
        for b in range(100):
            sampled_histories = [rng.choice(history_ids) for _ in history_ids]
            sampled_msrc = [tree for history in sampled_histories for tree in by_history[history]]
            sampled_msc = [rng.choice(msc_trees) for _ in range(500)]
            thresholds = []; reversals = 0
            for q, key, msc_counts, history_counts in prepared:
                d0_counts = msc_counts[1]
                d1_counts = defaultdict(int)
                for history in sampled_histories:
                    for topology, count in history_counts[history].items(): d1_counts[topology] += count
                def margin(counts, total):
                    vals = sorted((count / total for k, count in counts.items() if k != key), reverse=True)
                    return counts[key] / total - (vals[0] if vals else 0)
                d0 = margin(d0_counts, 500); d1 = margin(d1_counts, 500)
                if d0 * d1 < 0:
                    reversals += 1; thresholds.append(-d0 / (d1 - d0))
            boot_reversals.append(reversals / len(candidates)); boot_thresholds.extend(thresholds)
        output.append({"replicate": rep, "treatment": treatment, "candidate_quartets": len(candidates), "bootstrap_reversal_fraction_median": float(np.median(boot_reversals)), "bootstrap_reversal_fraction_q05": float(np.quantile(boot_reversals, .05)), "bootstrap_reversal_fraction_q95": float(np.quantile(boot_reversals, .95)), "threshold_q05": float(np.quantile(boot_thresholds, .05)) if boot_thresholds else "", "threshold_median": float(np.median(boot_thresholds)) if boot_thresholds else "", "threshold_q95": float(np.quantile(boot_thresholds, .95)) if boot_thresholds else ""})
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t"); writer.writeheader(); writer.writerows(output)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
