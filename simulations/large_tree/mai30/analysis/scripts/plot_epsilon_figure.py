#!/usr/bin/env python3
"""Plot the preliminary epsilon main figure from complete source pools."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis/figures/mai30_epsilon_main_figure.pdf"


def main() -> None:
    rows = list(csv.DictReader((ROOT / "analysis/results/mai30_epsilon_astral.tsv").open(), delimiter="\t"))
    for r in rows:
        r["epsilon"] = float(r["epsilon"]); r["rf"] = float(r["normalized_rf"]); r["cu"] = float(r["mean_absolute_cu_error_common"])
    qrows = list(csv.DictReader((ROOT / "analysis/results/mai30_epsilon_quartets.tsv").open(), delimiter="\t"))
    reversal = defaultdict(float)
    for key in {(r["replicate"], r["treatment"]) for r in qrows}:
        sub = [r for r in qrows if (r["replicate"], r["treatment"]) == key]
        reversal[key] = sum(r["endpoint_reversal"] == "1" for r in sub) / len(sub)
    colors = {"msrc_no_suppression": "#377eb8", "msrc_moderate_suppression": "#ff8c42", "msrc_strong_suppression": "#d73027"}
    labels = {"msrc_no_suppression": "suppression 1.0", "msrc_moderate_suppression": "suppression 0.5", "msrc_strong_suppression": "suppression 0.1"}
    fig, axes = plt.subplots(2, 2, figsize=(10, 7.2), constrained_layout=True)
    ax = axes[0, 0]
    source = next((ROOT / "datasets/processed/mai30_scaled_wf_runs").glob("*replicate_042*/history_000/species_tree.nwk"))
    ax.text(0.02, 0.82, "Mai 30-taxon fixed backbone", transform=ax.transAxes, fontsize=14)
    ax.text(0.02, 0.64, "MSC loci  ───────────────┐", transform=ax.transAxes, family="monospace")
    ax.text(0.02, 0.52, "MSRC loci ── rearrangement ┘", transform=ax.transAxes, family="monospace")
    ax.text(0.02, 0.30, "ε = fraction of loci sampled\nfrom the MSRC source pool", transform=ax.transAxes, fontsize=11)
    ax.text(0.02, 0.08, "scaled WF: c = 0.001; complete pools only", transform=ax.transAxes, fontsize=9, color="#555")
    ax.set_title("A  Evolutionary scenario and locus mixture", loc="left", fontsize=10)
    ax.axis("off")
    for ax, metric, ylabel, title in [(axes[0, 1], "rf", "normalized RF error", "B  Topology error versus ε"), (axes[1, 0], "cu", "mean absolute CULength error", "C  Branch-length error versus ε")]:
        for rep, style in [("042", "-"), ("007", "--")]:
            for treatment, color in colors.items():
                sub = [r for r in rows if r["replicate"] == rep and r["treatment"] == treatment]
                sub.sort(key=lambda r:r["epsilon"])
                if not sub: continue
                ax.plot([r["epsilon"] for r in sub], [r[metric] for r in sub], marker="o", color=color, linestyle=style, label=f"{rep} {labels[treatment]}")
        ax.set_xlabel("ε"); ax.set_ylabel(ylabel); ax.set_title(title, loc="left", fontsize=10); ax.grid(alpha=.25)
    axes[0, 1].legend(fontsize=6, frameon=False, ncol=2)
    axes[1, 1].set_title("D  Candidate population quartet reversals", loc="left", fontsize=10)
    x = np.arange(2); width = .24
    for i, treatment in enumerate(colors):
        vals = [reversal.get((rep, treatment), 0) for rep in ("042", "007")]
        axes[1, 1].bar(x + (i-1)*width, vals, width, color=colors[treatment], label=labels[treatment])
    axes[1, 1].set_xticks(x, ["low ILS\n042", "extreme ILS\n007"]); axes[1, 1].set_ylabel("fraction of sampled quartets with endpoint sign change"); axes[1, 1].set_ylim(0, 1); axes[1, 1].grid(axis="y", alpha=.25)
    fig.suptitle("Mai30 epsilon-mixture experiment (complete matched source pools)", fontsize=13)
    OUT.parent.mkdir(parents=True, exist_ok=True); fig.savefig(OUT); print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
