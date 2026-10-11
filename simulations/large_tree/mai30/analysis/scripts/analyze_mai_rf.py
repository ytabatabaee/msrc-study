#!/usr/bin/env python3
"""Corrected RF analysis and strict matching for the Mai30 experiment."""

from __future__ import annotations

import csv
import json
import math
import random
import re
import statistics
import sys
import zlib
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "datasets/raw"
PROCESSED = ROOT / "datasets/processed"
RESULTS = ROOT / "analysis/results"
FIGURES = ROOT / "analysis/figures"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rf_metrics import TreeValidationError, Tree, read_tree_file, rf_distance  # noqa: E402


MAIN_CONDITIONS = (
    "outgroup.0.species.0.15.genes.0.15",
    "outgroup.0.species.1.5.genes.1.5",
    "outgroup.0.species.5.genes.5",
)
N_REPLICATES = 100
EXPECTED_GENES = 500
WITHIN_PAIRS = 5000
CROSS_PAIRS = 5000
SEED = 20261011


def stable_seed(*values: str) -> int:
    return SEED + zlib.crc32("|".join(values).encode())


def mean_rf(a: list[Tree], b: list[Tree], pairs: list[tuple[int, int]]) -> float:
    values = [rf_distance(a[i], b[j]).normalized for i, j in pairs]
    return statistics.mean(values)


def random_pairs(n_a: int, n_b: int, count: int, seed: int, distinct_within: bool = False) -> list[tuple[int, int]]:
    if distinct_within:
        possible = n_a * (n_a - 1) // 2
        if count >= possible:
            return [(i, j) for i in range(n_a) for j in range(i + 1, n_a)]
    rng = random.Random(seed)
    pairs: list[tuple[int, int]] = []
    seen: set[tuple[int, int]] = set()
    target = min(count, n_a * (n_a - 1) // 2) if distinct_within else count
    while len(pairs) < target:
        i = rng.randrange(n_a)
        j = rng.randrange(n_a if distinct_within else n_b)
        if distinct_within:
            if i == j:
                continue
            pair = (min(i, j), max(i, j))
            if pair in seen:
                continue
        else:
            pair = (i, j)
        seen.add(pair)
        pairs.append(pair)
    return pairs


def species_mean(trees: list[Tree], species: Tree) -> float:
    return statistics.mean(rf_distance(tree, species).normalized for tree in trees)


def get_msrc_collections() -> dict[tuple[str, str, str], tuple[Path, dict]]:
    result = {}
    root = PROCESSED / "mai30_scaled_wf_pilot"
    for summary_path in sorted(root.glob("*/collection_summary.json")):
        summary = json.loads(summary_path.read_text())
        group = summary["group"]
        parts = group.split("__")
        try:
            condition = parts[0]
            replicate = next(x.removeprefix("replicate_") for x in parts if x.startswith("replicate_"))
            treatment = next(x for x in parts if x.startswith("msrc_"))
        except StopIteration as exc:
            raise RuntimeError(f"cannot parse MSRC collection group: {group}") from exc
        key = (condition, replicate, treatment)
        if key in result:
            raise RuntimeError(f"duplicate MSRC collection for {key}")
        tree_path = summary_path.parent / "gene_trees.nwk"
        result[key] = (tree_path, summary)
    return result


def write_table(path: Path, rows: list[dict], fields: list[str], delimiter: str = "\t") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def metric_rows() -> tuple[list[dict], list[dict]]:
    collections = get_msrc_collections()
    treatments = sorted({key[2] for key in collections})
    if not treatments:
        raise RuntimeError("no MSRC collections found")
    rows: list[dict] = []
    errors: list[dict] = []
    for condition in MAIN_CONDITIONS:
        for rep_number in range(1, N_REPLICATES + 1):
            replicate = f"{rep_number:03d}"
            species_path = RAW / condition / replicate / "s_tree.trees"
            original_path = RAW / condition / replicate / "truegenetrees"
            row_base = {"model_condition": condition, "replicate_id": replicate,
                        "original_gene_tree_count": "", "msrc_gene_tree_count": "",
                        "taxon_count": "", "validation_status": "", "validation_message": ""}
            try:
                if not species_path.exists() or not original_path.exists():
                    raise TreeValidationError("missing original species tree or gene-tree file")
                species = read_tree_file(species_path, 1)[0]
                original = read_tree_file(original_path, EXPECTED_GENES)
                if len(species.taxa) != 30:
                    raise TreeValidationError(f"expected 30 taxa, found {len(species.taxa)}")
                if any(tree.taxa != species.taxa for tree in original):
                    raise TreeValidationError("original gene-tree taxon set differs from species tree")
                original_species_rf = species_mean(original, species)
                original_pairs = random_pairs(len(original), len(original), WITHIN_PAIRS,
                                               stable_seed(condition, replicate, "original"), True)
                original_within = mean_rf(original, original, original_pairs)
                for treatment in treatments:
                    key = (condition, replicate, treatment)
                    msrc_entry = collections.get(key)
                    row = dict(row_base)
                    row.update({"treatment": treatment, "original_gene_tree_count": len(original),
                                "taxon_count": len(species.taxa), "original_species_rf": original_species_rf,
                                "original_within_rf": original_within})
                    if msrc_entry is None:
                        row.update({"validation_status": "MISSING_MSRC", "validation_message": "no matched MSRC collection"})
                        rows.append(row)
                        continue
                    msrc_path, summary = msrc_entry
                    try:
                        msrc = read_tree_file(msrc_path)
                        if any(tree.taxa != species.taxa for tree in msrc):
                            raise TreeValidationError("MSRC gene-tree taxon set differs from species tree")
                        n_msrc = len(msrc)
                        status = "PASS" if n_msrc == EXPECTED_GENES else "PARTIAL"
                        msg = "" if status == "PASS" else f"MSRC collection has {n_msrc}, expected {EXPECTED_GENES}"
                        cross_pairs = [(random.Random(stable_seed(condition, replicate, treatment, "cross") + i).randrange(len(original)),
                                        random.Random(stable_seed(condition, replicate, treatment, "cross") + i + 1).randrange(n_msrc))
                                       for i in range(CROSS_PAIRS)]
                        index_pairs = [(i, i) for i in range(min(len(original), n_msrc))]
                        msrc_pairs = random_pairs(n_msrc, n_msrc, WITHIN_PAIRS,
                                                  stable_seed(condition, replicate, treatment, "msrc"), True)
                        row.update({
                            "msrc_gene_tree_count": n_msrc,
                            "cross_distribution_rf": mean_rf(original, msrc, cross_pairs),
                            "index_paired_rf_diagnostic": mean_rf(original, msrc, index_pairs),
                            "msrc_within_rf": mean_rf(msrc, msrc, msrc_pairs),
                            "msrc_species_rf": species_mean(msrc, species),
                            "validation_status": status,
                            "validation_message": msg,
                        })
                    except (OSError, TreeValidationError) as exc:
                        row.update({"validation_status": "INVALID_MSRC", "validation_message": str(exc)})
                        errors.append({**row, "error_source": "MSRC"})
                    rows.append(row)
            except (OSError, TreeValidationError) as exc:
                errors.append({**row_base, "error_source": "ORIGINAL", "validation_status": "INVALID_ORIGINAL", "validation_message": str(exc)})
                for treatment in treatments:
                    rows.append({**row_base, "treatment": treatment, "validation_status": "INVALID_ORIGINAL", "validation_message": str(exc)})
    return rows, errors


METRICS = ("cross_distribution_rf", "original_within_rf", "msrc_within_rf", "original_species_rf", "msrc_species_rf")


def aggregate(rows: list[dict]) -> list[dict]:
    output = []
    for condition in MAIN_CONDITIONS:
        for treatment in sorted({row["treatment"] for row in rows}):
            subset = [row for row in rows if row["model_condition"] == condition and row["treatment"] == treatment]
            valid = [row for row in subset if row["validation_status"] == "PASS"]
            original_valid = [row for row in subset if row["validation_status"] != "INVALID_ORIGINAL" and row.get("original_species_rf", "") != ""]
            out = {"model_condition": condition, "treatment": treatment, "expected_replicates": N_REPLICATES,
                   "valid_complete_replicates": len(valid),
                   "partial_replicates": sum(row["validation_status"] == "PARTIAL" for row in subset),
                   "missing_msrc_replicates": sum(row["validation_status"] == "MISSING_MSRC" for row in subset),
                   "invalid_replicates": sum(row["validation_status"].startswith("INVALID") for row in subset),
                   "aggregate_status": "COMPLETE_100" if len(valid) == N_REPLICATES else "INCOMPLETE"}
            for metric in METRICS:
                source_rows = original_valid if metric.startswith("original_") else valid
                values = [float(row[metric]) for row in source_rows if row.get(metric, "") != ""]
                out[f"{metric}_mean"] = statistics.mean(values) if values else ""
                out[f"{metric}_sd"] = statistics.stdev(values) if len(values) > 1 else ""
                out[f"{metric}_ci95_low"] = (statistics.mean(values) - 1.96 * statistics.stdev(values) / math.sqrt(len(values))) if len(values) > 1 else (values[0] if values else "")
                out[f"{metric}_ci95_high"] = (statistics.mean(values) + 1.96 * statistics.stdev(values) / math.sqrt(len(values))) if len(values) > 1 else (values[0] if values else "")
            output.append(out)
    return output


def plot(rows: list[dict]) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    metrics = [("cross_distribution_rf", "cross-distribution nRF"),
               ("original_within_rf", "original within-pool nRF"),
               ("msrc_within_rf", "MSRC within-pool nRF"),
               ("original_species_rf", "original gene-tree/species nRF"),
               ("msrc_species_rf", "MSRC gene-tree/species nRF")]
    treatments = sorted({row["treatment"] for row in rows})
    fig, axes = plt.subplots(1, len(metrics), figsize=(16, 4), sharey=False)
    for ax, (metric, label) in zip(axes, metrics):
        for treatment in treatments:
            vals_by_condition = []
            for condition in MAIN_CONDITIONS:
                vals = [float(row[metric]) for row in rows if row["model_condition"] == condition and row["treatment"] == treatment and row["validation_status"] == "PASS" and row.get(metric, "") != ""]
                vals_by_condition.append(vals)
            positions = list(range(len(MAIN_CONDITIONS)))
            for x, vals in zip(positions, vals_by_condition):
                if vals:
                    ax.scatter([x] * len(vals), vals, s=12, alpha=0.65, label=treatment if x == 0 else None)
        ax.set_title(label, fontsize=9)
        ax.set_xticks(range(len(MAIN_CONDITIONS)), ["0.15", "1.5", "5"])
        ax.set_xlabel("Mai model condition", fontsize=8)
        ax.set_ylim(0, 1.02)
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("normalized unrooted RF", fontsize=9)
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIGURES / "mai30_rf_corrected.png", dpi=220)


def main() -> None:
    rows, errors = metric_rows()
    fields = list(rows[0])
    write_table(RESULTS / "mai30_rf_replicate_results.tsv", rows, fields)
    write_table(RESULTS / "mai30_rf_replicate_results.csv", rows, fields, ",")
    aggregate_rows = aggregate(rows)
    aggregate_fields = list(aggregate_rows[0])
    write_table(RESULTS / "mai30_rf_aggregate_results.tsv", aggregate_rows, aggregate_fields)
    write_table(RESULTS / "mai30_rf_aggregate_results.csv", aggregate_rows, aggregate_fields, ",")
    manifest_fields = ["model_condition", "replicate_id", "treatment", "original_gene_tree_count",
                       "msrc_gene_tree_count", "taxon_count", "validation_status", "validation_message"]
    write_table(RESULTS / "mai30_rf_validation_manifest.tsv", rows, manifest_fields)
    write_table(RESULTS / "mai30_rf_validation_manifest.csv", rows, manifest_fields, ",")
    plot(rows)
    print(f"wrote {len(rows)} per-replicate rows")
    print(f"complete MSRC rows: {sum(row['validation_status'] == 'PASS' for row in rows)}")
    print(f"partial MSRC rows: {sum(row['validation_status'] == 'PARTIAL' for row in rows)}")
    print(f"missing MSRC rows: {sum(row['validation_status'] == 'MISSING_MSRC' for row in rows)}")
    print(f"wrote {len(aggregate_rows)} aggregate rows")


if __name__ == "__main__":
    main()
