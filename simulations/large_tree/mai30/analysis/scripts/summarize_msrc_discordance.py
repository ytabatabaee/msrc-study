#!/usr/bin/env python3
"""Summarize the two primary Mai30 MSRC comparisons.

The first output measures MSRC gene-tree discordance from the fixed Mai
species tree.  The second compares the MSRC gene-tree distribution with the
original Mai gene-tree distribution for the same replicate.  Independent
gene-tree realizations are compared as distributions; they are never paired
by locus.

The script includes complete 500-tree collections and the currently available
20-tree smoke collections, with an explicit scope column.  Downstream
comparisons should normally filter ``analysis_scope == complete_500``.
"""

from __future__ import annotations

import csv
import math
import random
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "datasets" / "processed"
RAW = ROOT / "datasets" / "raw"
OUT = ROOT / "analysis" / "results"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from compute_original_baseline import ParsedTree, parse_newick, quartet_key, rf  # noqa: E402


N_QUARTETS = 1000
SEED = 20261011


def read_trees(path: Path) -> list[ParsedTree]:
    trees: list[ParsedTree] = []
    with path.open() as handle:
        for line_no, line in enumerate(handle, 1):
            if line.strip():
                try:
                    trees.append(parse_newick(line))
                except Exception as exc:
                    raise ValueError(f"{path}:{line_no}: malformed Newick") from exc
    if not trees:
        raise ValueError(f"{path}: no trees")
    taxa = trees[0].taxa
    for index, tree in enumerate(trees, 1):
        if tree.taxa != taxa:
            raise ValueError(f"{path}: tree {index} has an inconsistent taxon set")
        if tree.n_internal != len(taxa) - 1:
            raise ValueError(f"{path}: tree {index} is not fully bifurcating")
    return trees


def quantile(values: list[float], p: float) -> float:
    values = sorted(values)
    if len(values) == 1:
        return values[0]
    x = p * (len(values) - 1)
    lo = math.floor(x)
    hi = math.ceil(x)
    return values[lo] + (values[hi] - values[lo]) * (x - lo)


def empirical_wasserstein(a: list[float], b: list[float]) -> float:
    """One-dimensional equal-weight empirical Wasserstein distance."""
    return sum(abs(quantile(a, (i + 0.5) / max(len(a), len(b))) -
                   quantile(b, (i + 0.5) / max(len(a), len(b))))
               for i in range(max(len(a), len(b)))) / max(len(a), len(b))


def ks_distance(a: list[float], b: list[float]) -> float:
    values = sorted(set(a + b))
    return max(abs(sum(x <= value for x in a) / len(a) -
                   sum(x <= value for x in b) / len(b)) for value in values)


def sampled_quartets(taxa: list[str], count: int) -> list[tuple[str, str, str, str]]:
    rng = random.Random(SEED)
    return [tuple(rng.sample(taxa, 4)) for _ in range(count)]


def species_metrics(trees: list[ParsedTree], species: ParsedTree, quartets: list[tuple[str, str, str, str]]) -> dict[str, float]:
    rf_values = [rf(tree, species)[1] for tree in trees]
    split_values = [sum(split in tree.splits for split in species.splits) / len(species.splits) for tree in trees]
    quartet_values = [sum(quartet_key(tree, q) == quartet_key(species, q) for tree in trees) / len(trees)
                      for q in quartets]
    return {
        "mean_normalized_rf": sum(rf_values) / len(rf_values),
        "median_normalized_rf": quantile(rf_values, 0.5),
        "exact_species_topology_fraction": sum(value == 0 for value in rf_values) / len(rf_values),
        "mean_species_split_concordance": sum(split_values) / len(split_values),
        "mean_sampled_quartet_concordance": sum(quartet_values) / len(quartet_values),
    }


def topology_distribution_metrics(msrc: list[ParsedTree], original: list[ParsedTree], species: ParsedTree, quartets: list[tuple[str, str, str, str]]) -> dict[str, float]:
    msrc_rf = [rf(tree, species)[1] for tree in msrc]
    original_rf = [rf(tree, species)[1] for tree in original]
    msrc_split_counts = Counter(split for tree in msrc for split in tree.splits)
    original_split_counts = Counter(split for tree in original for split in tree.splits)
    all_splits = set(msrc_split_counts) | set(original_split_counts)
    split_l1 = sum(abs(msrc_split_counts[s] / len(msrc) - original_split_counts[s] / len(original))
                   for s in all_splits) / 2.0

    quartet_l1 = []
    quartet_margin_difference = []
    for quartet in quartets:
        msrc_counts = Counter(quartet_key(tree, quartet) for tree in msrc)
        original_counts = Counter(quartet_key(tree, quartet) for tree in original)
        states = set(msrc_counts) | set(original_counts)
        msrc_freq = {state: msrc_counts[state] / len(msrc) for state in states}
        original_freq = {state: original_counts[state] / len(original) for state in states}
        quartet_l1.append(sum(abs(msrc_freq[state] - original_freq[state]) for state in states) / 2.0)
        msrc_values = sorted(msrc_freq.values(), reverse=True) + [0.0, 0.0]
        original_values = sorted(original_freq.values(), reverse=True) + [0.0, 0.0]
        msrc_margin = msrc_values[0] - msrc_values[1]
        original_margin = original_values[0] - original_values[1]
        quartet_margin_difference.append(msrc_margin - original_margin)

    return {
        "msrc_mean_rf_minus_original_mean_rf": sum(msrc_rf) / len(msrc_rf) - sum(original_rf) / len(original_rf),
        "rf_distribution_wasserstein": empirical_wasserstein(msrc_rf, original_rf),
        "rf_distribution_ks": ks_distance(msrc_rf, original_rf),
        "msrc_exact_topology_minus_original": sum(x == 0 for x in msrc_rf) / len(msrc_rf) - sum(x == 0 for x in original_rf) / len(original_rf),
        "msrc_species_split_concordance_minus_original": (
            sum(sum(split in tree.splits for split in species.splits) / len(species.splits) for tree in msrc) / len(msrc)
            - sum(sum(split in tree.splits for split in species.splits) / len(species.splits) for tree in original) / len(original)
        ),
        "mean_quartet_frequency_l1": sum(quartet_l1) / len(quartet_l1),
        "mean_quartet_margin_difference": sum(quartet_margin_difference) / len(quartet_margin_difference),
    }


def parse_group(group: str) -> tuple[str, str, str]:
    parts = group.split("__")
    condition = parts[0]
    rep = next(part.removeprefix("replicate_") for part in parts if part.startswith("replicate_"))
    treatment = next(part for part in parts if part.startswith("msrc_"))
    return condition, rep, treatment


def main() -> None:
    astral = {}
    astral_path = OUT / "mai30_scaled_wf_astral.tsv"
    if astral_path.exists():
        with astral_path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                astral[row["group"]] = row

    species_rows: list[dict[str, object]] = []
    distribution_rows: list[dict[str, object]] = []
    quartets_by_taxa: dict[frozenset[str], list[tuple[str, str, str, str]]] = {}
    pilot = PROCESSED / "mai30_scaled_wf_pilot"
    for summary_path in sorted(pilot.glob("*/collection_summary.json")):
        import json
        summary = json.loads(summary_path.read_text())
        group = summary["group"]
        condition, rep, treatment = parse_group(group)
        msrc_path = Path(summary["gene_tree_path"])
        # The summaries predate this workspace layout on some machines.
        if not msrc_path.exists():
            msrc_path = summary_path.parent / "gene_trees.nwk"
        species_path = RAW / condition / rep / "s_tree.trees"
        original_path = RAW / condition / rep / "truegenetrees"
        species = parse_newick(species_path.read_text())
        msrc = read_trees(msrc_path)
        original = read_trees(original_path)
        if species.taxa not in quartets_by_taxa:
            quartets_by_taxa[species.taxa] = sampled_quartets(sorted(species.taxa), N_QUARTETS)
        quartets = quartets_by_taxa[species.taxa]
        sm = species_metrics(msrc, species, quartets)
        dm = topology_distribution_metrics(msrc, original, species, quartets)
        ar = astral.get(group, {})
        scope = "complete_500" if len(msrc) == 500 and summary["n_independent_histories"] >= 25 else "partial_smoke"
        base = {
            "group": group, "condition": condition, "replicate": rep, "treatment": treatment,
            "n_gene_trees": len(msrc), "n_independent_histories": summary["n_independent_histories"],
            "analysis_scope": scope, "original_n_gene_trees": len(original),
            "original_mean_normalized_rf": species_metrics(original, species, quartets)["mean_normalized_rf"],
            "astral_normalized_rf": ar.get("normalized_rf", ""),
            "astral_exact_topology_recovery": ar.get("exact_topology_recovery", ""),
        }
        species_rows.append({**base, **sm})
        distribution_rows.append({**base, **dm})

    OUT.mkdir(parents=True, exist_ok=True)
    for path, rows in [(OUT / "mai30_msrc_vs_species.tsv", species_rows),
                       (OUT / "mai30_msrc_vs_original_distribution.tsv", distribution_rows)]:
        with path.open("w", newline="") as handle:
            fields = list(rows[0])
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"wrote {path} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
