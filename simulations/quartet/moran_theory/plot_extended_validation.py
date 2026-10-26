"""Plot the targeted extended pectinate-theory validation."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COLORS = {"q1": "#0072B2", "q2": "#D55E00", "q3": "#009E73"}


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_figure(rows: list[dict[str, str]], conditions: list[dict[str, str]], output: Path) -> None:
    by_index = {r["condition_index"]: r for r in conditions}
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.2), constrained_layout=True)

    ax = axes[0, 0]
    for topology in ("q1", "q2", "q3"):
        subset = [r for r in rows if r["topology"] == topology]
        x = np.asarray([float(r["analytic"]) for r in subset])
        y = np.asarray([float(r["observed"]) for r in subset])
        se = np.asarray([float(r["mc_se"]) for r in subset])
        ax.errorbar(x, y, yerr=se, fmt="o", ms=4.5, alpha=.78, color=COLORS[topology],
                    capsize=2, lw=.7, label=topology)
    ax.plot([0.20, .50], [.20, .50], ls="--", lw=.8, color="#555555")
    ax.set(xlim=(.20, .50), ylim=(.20, .50), xlabel="Analytical probability",
           ylabel="Simulated frequency", title="A  Predicted versus simulated")
    ax.legend(frameon=False, fontsize=8)
    ax.set_aspect("equal", adjustable="box")

    ax = axes[0, 1]
    for topology in ("q1", "q2", "q3"):
        subset = [r for r in rows if r["topology"] == topology]
        x = np.asarray([int(r["condition_index"]) for r in subset])
        y = np.asarray([float(r["standardized_error"]) for r in subset])
        ax.scatter(x, y, s=20, color=COLORS[topology], alpha=.78, edgecolors="none", label=topology)
    ax.axhline(0, color="#555555", lw=.8)
    ax.axhline(2, color="#999999", ls=":", lw=.7)
    ax.axhline(-2, color="#999999", ls=":", lw=.7)
    ax.set(xlabel="Condition index", ylabel="Standardized residual", title="B  Monte Carlo residuals")
    ax.set_xticks(range(1, len(conditions) + 1))
    ax.tick_params(axis="x", labelsize=7)
    ax.grid(axis="y", color="#dddddd", lw=.5)

    ax = axes[1, 0]
    x = np.asarray([float(r["population_size_ratio"]) for r in conditions])
    y = np.asarray([float(r["max_absolute_error"]) for r in conditions])
    c = np.asarray([float(r["intermediate_length"]) for r in conditions])
    sc = ax.scatter(x, y, c=c, cmap="viridis", s=30, alpha=.82, edgecolors="none")
    ax.set(xlabel="Max/min chromosome-copy count", ylabel="Maximum absolute error",
           title="C  Error versus population-size asymmetry")
    cb = fig.colorbar(sc, ax=ax, pad=.02)
    cb.set_label("Intermediate duration", fontsize=8)
    ax.grid(color="#dddddd", lw=.5)

    ax = axes[1, 1]
    for topology in ("q1", "q2", "q3"):
        subset = [r for r in rows if r["topology"] == topology]
        x = np.asarray([float(r["intermediate_length"]) for r in subset])
        theory = np.asarray([float(r["analytic"]) for r in subset])
        observed = np.asarray([float(r["observed"]) for r in subset])
        ax.plot(x, theory, linestyle="", marker="_", ms=8, mew=1.3, color=COLORS[topology], label=f"{topology} theory")
        ax.scatter(x, observed, s=18, color=COLORS[topology], alpha=.62, edgecolors="none", marker="o")
    ax.set(xlabel="Intermediate duration $T_I$", ylabel="Quartet probability",
           title="D  Effect of intermediate duration")
    ax.legend(frameon=False, fontsize=7, ncol=2)
    ax.grid(color="#dddddd", lw=.5)

    for ax in axes.flat:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=8)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    root = Path(__file__).parent / "results"
    parser.add_argument("--input", type=Path, default=root / "unbalanced_species_tree_extended_validation.tsv")
    parser.add_argument("--conditions", type=Path, default=root / "unbalanced_species_tree_extended_conditions.tsv")
    parser.add_argument("--output", type=Path, default=root / "unbalanced_species_tree_extended_validation")
    args = parser.parse_args()
    build_figure(read(args.input), read(args.conditions), args.output)


if __name__ == "__main__":
    main()
