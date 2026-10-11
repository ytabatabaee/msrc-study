#!/usr/bin/env python3
"""Run ASTRAL-IV on complete epsilon mixtures and summarize inference."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
MIX_ROOT = ROOT / "datasets/processed/mai30_epsilon_mixtures"
OUT = ROOT / "analysis/results/mai30_epsilon_astral.tsv"
ASTRAL = Path("/Users/ytabatabaee/Desktop/ASTER/bin/astral")
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_astral_scaled_pilot import parse_lengths_correct, normalized_rf  # noqa: E402


def main() -> None:
    rows = []
    for cfg_path in sorted(MIX_ROOT.glob("replicate_*/mixture_config.json")):
        cfg = json.loads(cfg_path.read_text())
        group = cfg["msrc_group"]
        out = cfg_path.parent
        trees = [x for x in (out / "gene_trees.nwk").read_text().splitlines() if x.strip()]
        if len(trees) != 500:
            raise ValueError(f"mixture count is not 500: {out}")
        astral_out = out / "astral_species_tree.nwk"
        command = [str(ASTRAL), "-i", str(out / "gene_trees.nwk"), "-o", str(astral_out), "-t", "2", "-u", "1", "--length", "CULength"]
        result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=False)
        (out / "astral.log").write_text(result.stdout + result.stderr)
        if result.returncode != 0:
            raise RuntimeError(f"ASTRAL failed for {out}: {result.stderr[-1000:]}")
        inferred_taxa, inferred_splits, inferred_lengths = parse_lengths_correct(astral_out.read_text())
        run_dir = ROOT / "datasets/processed/mai30_scaled_wf_runs" / group / "history_000"
        cfg_history = json.loads((run_dir / "simulation_config.json").read_text())
        true_taxa, true_splits, true_lengths = parse_lengths_correct((run_dir / "species_tree.nwk").read_text())
        if inferred_taxa != true_taxa:
            raise ValueError(f"taxon set mismatch: {out}")
        n_scaled = int(cfg_history["species_tree"]["default_effective_population_size"])
        true_cu = {s: true_lengths[s] / (2 * n_scaled) for s in true_splits}
        common = set(true_cu) & set(inferred_lengths)
        errors = [inferred_lengths[s] - true_cu[s] for s in common]
        abs_errors = [abs(x) for x in errors]
        rows.append({
            "replicate": cfg["replicate"],
            "treatment": group.split("__")[2],
            "epsilon": cfg["epsilon"],
            "mixture_replicate": cfg["mixture_replicate"],
            "n_msrc": cfg["n_msrc"],
            "n_msc": cfg["n_msc"],
            "n_gene_trees": len(trees),
            "normalized_rf": normalized_rf(inferred_splits, true_splits, len(true_taxa)),
            "exact_topology_recovery": int(inferred_splits == true_splits),
            "correct_internal_splits": len(inferred_splits & true_splits),
            "common_internal_splits": len(common),
            "mean_signed_cu_error_common": float(np.mean(errors)) if errors else "",
            "mean_absolute_cu_error_common": float(np.mean(abs_errors)) if abs_errors else "",
            "median_absolute_cu_error_common": float(np.median(abs_errors)) if abs_errors else "",
            "astral_output": str(astral_out),
        })
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    print(f"wrote {OUT} ({len(rows)} datasets)")


if __name__ == "__main__":
    main()
