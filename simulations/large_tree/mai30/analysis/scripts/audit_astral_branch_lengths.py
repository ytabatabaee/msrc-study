#!/usr/bin/env python3
"""Audit ASTER CULength values against scaled Mai species-tree edges."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[2]
sys.path.insert(0, str(SCRIPT.parent))
from run_astral_scaled_pilot import parse_lengths_correct  # noqa: E402

RESULT = ROOT / "analysis/results/mai30_branch_length_edges.tsv"
FIG = ROOT / "analysis/figures/mai30_branch_length_audit.png"
PROCESSED = ROOT / "datasets/processed"


def find_source(group: str) -> tuple[Path, int]:
    if group.endswith("__msc_control"):
        condition, rep = group.rsplit("__replicate_", 1)
        rep = rep.replace("__msc_control", "")
        pilot = next(PROCESSED.joinpath("mai30_scaled_wf_runs").glob(f"{condition}__replicate_{rep}__*__freq_0.5/history_000"))
        cfg = json.loads((pilot / "simulation_config.json").read_text())
        return pilot / "species_tree.nwk", int(cfg["species_tree"]["default_effective_population_size"])
    run = PROCESSED / "mai30_scaled_wf_runs" / group / "history_000"
    cfg = json.loads((run / "simulation_config.json").read_text())
    return run / "species_tree.nwk", int(cfg["species_tree"]["default_effective_population_size"])


def main() -> None:
    astral = ROOT / "analysis/results/mai30_scaled_wf_astral.tsv"
    rows = []
    with astral.open() as handle:
        for result in csv.DictReader(handle, delimiter="\t"):
            group = result["group"]
            source, n_scaled = find_source(group)
            taxa, true_splits, true_lengths = parse_lengths_correct(source.read_text())
            _, inferred_splits, inferred_lengths = parse_lengths_correct(Path(result["astral_output"]).read_text())
            for split in sorted(true_splits, key=lambda s: (len(s), tuple(sorted(s)))):
                true_cu = true_lengths[split] / (2.0 * n_scaled)
                estimate = inferred_lengths.get(split, np.nan)
                rows.append({
                    "group": group,
                    "source_type": "MSC" if group.endswith("__msc_control") else "MSRC",
                    "replicate": group.split("__replicate_")[1].split("__")[0],
                    "n_gene_trees": result["n_gene_trees"],
                    "split_size": len(split),
                    "true_cu": true_cu,
                    "estimated_cu": estimate,
                    "error_cu": estimate - true_cu if np.isfinite(estimate) else np.nan,
                    "absolute_error_cu": abs(estimate - true_cu) if np.isfinite(estimate) else np.nan,
                    "common_split": int(split in inferred_splits),
                    "topology_exact": result["exact_topology_recovery"],
                })
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    with RESULT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    selected = [r for r in rows if r["replicate"] == "042" and r["source_type"] == "MSC" and r["common_split"] == 1]
    x = np.array([float(r["true_cu"]) for r in selected])
    y = np.array([float(r["estimated_cu"]) for r in selected])
    axes[0].scatter(x, y, s=22, color="#377eb8")
    lim = max(np.max(x), np.max(y)) * 1.05
    axes[0].plot([0, lim], [0, lim], "k--", lw=0.8)
    axes[0].set(xlabel="true internal branch length (CU)", ylabel="ASTRAL CULength estimate", title="Replicate 042 MSC: common splits")
    axes[1].scatter(x, y - x, s=22, color="#d73027")
    axes[1].axhline(0, color="black", lw=0.8)
    axes[1].set(xlabel="true internal branch length (CU)", ylabel="estimate − truth (CU)", title="Error grows with long true branches")
    fig.tight_layout()
    FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG, dpi=220)
    print(f"wrote {RESULT}\nwrote {FIG}")


if __name__ == "__main__":
    main()
