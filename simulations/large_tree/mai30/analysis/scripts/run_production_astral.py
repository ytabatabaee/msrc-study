#!/usr/bin/env python3
"""Run one validation-gated ASTRAL task for a production collection."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve()
STUDY_ROOT = SCRIPT.parents[2]
import sys
sys.path.insert(0, str(SCRIPT.parent))
from run_astral_scaled_pilot import parse_lengths_correct, normalized_rf  # noqa: E402


def row_for(manifest: Path, task_id: int) -> dict[str, str]:
    with manifest.open() as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if int(row["task_id"]) == task_id]
    if len(rows) != 1:
        raise ValueError(f"task_id {task_id} is not unique")
    return rows[0]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", type=Path, default=STUDY_ROOT / "datasets/processed/mai30_production_manifest.tsv")
    p.add_argument("--task-id", type=int, required=True)
    p.add_argument("--astral", type=Path, required=True)
    args = p.parse_args()
    row = row_for(args.manifest, args.task_id)
    collection = STUDY_ROOT / row["output_root"] / row["collection_id"]
    marker = collection / "collection_complete.json"
    if not marker.exists():
        raise RuntimeError(f"missing validated collection marker: {marker}")
    gene_trees = collection / "gene_trees.nwk"
    source_species = STUDY_ROOT / row["species_tree_path"]
    output = collection / "astral_species_tree.nwk"
    log = collection / "astral.log"
    command = [str(args.astral), "-i", str(gene_trees), "-o", str(output), "-t", "2", "-u", "1", "--length", "CULength"]
    completed = subprocess.run(command, cwd=STUDY_ROOT, capture_output=True, text=True, check=False)
    log.write_text(completed.stdout + completed.stderr)
    if completed.returncode != 0:
        raise RuntimeError(f"ASTRAL failed for task {args.task_id}; see {log}")
    inferred_taxa, inferred_splits, inferred_lengths = parse_lengths_correct(output.read_text())
    true_taxa, true_splits, true_lengths = parse_lengths_correct(source_species.read_text())
    if inferred_taxa != true_taxa:
        raise ValueError(f"ASTRAL taxon mismatch for task {args.task_id}")
    n = len(true_taxa)
    n_scaled = int(row["scaled_effective_population_size"])
    true_cu = {split: length / (2.0 * n_scaled) for split, length in true_lengths.items() if split in true_splits}
    common = sorted(set(inferred_lengths) & set(true_cu), key=lambda x: tuple(sorted(x)))
    errors = [inferred_lengths[s] - true_cu[s] for s in common]
    summary = {"status": "PASS", "task_id": args.task_id, "collection_id": row["collection_id"], "n_taxa": n, "normalized_rf": normalized_rf(inferred_splits, true_splits, n), "exact_topology_recovery": int(inferred_splits == true_splits), "common_internal_splits": len(common), "mean_absolute_cu_error_common_splits": sum(abs(x) for x in errors) / len(errors) if errors else None, "mean_cu_bias_common_splits": sum(errors) / len(errors) if errors else None, "astral_output": str(output.relative_to(STUDY_ROOT)), "astral_command": " ".join(command)}
    (collection / "astral_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
