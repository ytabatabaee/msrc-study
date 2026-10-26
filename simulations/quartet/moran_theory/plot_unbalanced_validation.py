"""Plot analytical versus simulated pectinate-quartet validation results."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


COLORS = {"q1": "#0072B2", "q2": "#D55E00", "q3": "#009E73"}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def condition_key(row: dict[str, str]) -> str:
    return (
        f"M={row['M']}, f={row['fraction']}, rho={row['rho']}, "
        f"S={row['sister_length']}, I={row['intermediate_length']}, R={row['root_length']}"
    )


def build_figure(rows: list[dict[str, str]], output: Path) -> None:
    conditions: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        conditions[condition_key(row)].append(row)
    condition_names = list(conditions)
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.1), gridspec_kw={"width_ratios": (1.15, 1.35, 1.2)})

    ax = axes[0]
    for topology in ("q1", "q2", "q3"):
        subset = [r for r in rows if r["topology"] == topology]
        x = np.array([float(r["analytic"]) for r in subset])
        y = np.array([float(r["observed"]) for r in subset])
        se = np.array([float(r["mc_se"]) for r in subset])
        ax.errorbar(x, y, yerr=se, fmt="o", ms=4.5, alpha=0.75, color=COLORS[topology],
                    capsize=2, lw=0.7, label=topology)
    ax.plot([0.25, 0.40], [0.25, 0.40], color="#555555", ls="--", lw=0.8, zorder=0)
    ax.set(xlim=(0.25, 0.40), ylim=(0.25, 0.40), xlabel="Analytical probability", ylabel="Simulated frequency",
           title="A  Analytical vs simulated")
    ax.legend(frameon=False, fontsize=8)
    ax.set_aspect("equal", adjustable="box")

    ax = axes[1]
    x_positions = np.arange(len(rows))
    for i, row in enumerate(rows):
        ax.scatter(i, float(row["standardized_error"]), s=20, alpha=0.8,
                   color=COLORS[row["topology"]], edgecolors="none")
    ax.axhline(0.0, color="#555555", lw=0.8)
    ax.axhline(2.0, color="#999999", lw=0.6, ls=":")
    ax.axhline(-2.0, color="#999999", lw=0.6, ls=":")
    ax.set(xlim=(-1, len(rows)), xlabel="Condition index (three topology points per condition)", ylabel="Standardized residual",
           title="B  Monte Carlo residuals")
    ax.set_xticks([3 * i + 1 for i in range(len(condition_names))])
    ax.set_xticklabels([str(i + 1) for i in range(len(condition_names))], fontsize=7)
    ax.grid(axis="y", color="#dddddd", lw=0.5)

    ax = axes[2]
    max_errors = [max(float(r["absolute_error"]) for r in conditions[name]) for name in condition_names]
    ax.scatter(np.arange(len(condition_names)), max_errors, s=28, color="#4C78A8", alpha=0.8, edgecolors="none")
    ax.axhline(0.0, color="#555555", lw=0.8)
    ax.set(xlabel="Condition index", ylabel="Maximum absolute error", title="C  Error by condition")
    ax.set_xticks(np.arange(len(condition_names)))
    ax.set_xticklabels([str(i + 1) for i in range(len(condition_names))], fontsize=7)
    ax.grid(axis="y", color="#dddddd", lw=0.5)

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=8)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path(__file__).parent / "results" / "unbalanced_species_tree_validation.tsv")
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results" / "unbalanced_species_tree_validation")
    args = parser.parse_args()
    build_figure(read_rows(args.input), args.output)


if __name__ == "__main__":
    main()
