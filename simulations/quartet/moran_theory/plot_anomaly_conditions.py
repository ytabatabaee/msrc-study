"""Plot exact operator anomaly-condition calculations."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def read(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_figure(region, thresholds, simplex, output: Path):
    rho = sorted({float(r["rho"]) for r in region})
    I = sorted({float(r["intermediate_length"]) for r in region})
    grid_rho, grid_I = np.meshgrid(rho, I, indexing="ij")
    gamma_u = np.array([[float(next(r["gamma_unbalanced"] for r in region if float(r["rho"]) == rr and float(r["intermediate_length"]) == ii)) for ii in I] for rr in rho])
    gamma_b = np.array([[float(next(r["gamma_balanced"] for r in region if float(r["rho"]) == rr and float(r["intermediate_length"]) == ii)) for ii in I] for rr in rho])

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.0), constrained_layout=True)
    ax = axes[0, 0]
    cf = ax.contourf(grid_I, grid_rho, gamma_u, levels=15, cmap="RdBu_r")
    ax.contour(grid_I, grid_rho, gamma_u, levels=[0], colors="black", linewidths=1.2)
    fig.colorbar(cf, ax=ax, label=r"$\Gamma(0101)$")
    ax.set(xlabel=r"Intermediate duration $T_I$", ylabel=r"Switching $\rho$",
           title="A  Unbalanced anomaly region")
    ax.text(.98, .04, r"conservative bound: $T_I>10(1-p)\log2=3.47$\nfor $p=1/2$",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7)

    ax = axes[0, 1]
    vals = [(float(r["intermediate_length"]), float(r["root_q2_equals_q1"]))
            for r in thresholds if r["root_q2_equals_q1"] != "NaN"]
    if vals:
        ax.plot([x for x, _ in vals], [y for _, y in vals], marker="o", color="#0072B2", lw=1.2)
    ax.set(xlabel=r"Intermediate duration $T_I$", ylabel=r"Numerical $\rho^*$",
           title="B  q2 = q1 crossing for Z=0101")
    ax.text(.98, .06, "operator root; uniqueness not proved", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7)
    ax.grid(color="#dddddd", lw=.5)

    ax = axes[1, 0]
    selected_I = min(I, key=lambda value: abs(value - .2))
    rows = [r for r in region if abs(float(r["intermediate_length"]) - selected_I) < 1e-10]
    rows.sort(key=lambda r: float(r["rho"]))
    x = [float(r["rho"]) for r in rows]
    ax.plot(x, [float(r["gamma_unbalanced"]) for r in rows], "o-", color="#0072B2", label="unbalanced $\Gamma$")
    ax.plot(x, [float(r["gamma_balanced"]) for r in rows], "o--", color="#D55E00", label="balanced $\Gamma$")
    ax.plot(x, [float(r["delta_unbalanced"]) for r in rows], "s-", color="#009E73", label="unbalanced $\Delta$")
    ax.plot(x, [float(r["delta_balanced"]) for r in rows], "s--", color="#CC79A7", label="balanced $\Delta$")
    ax.axhline(0, color="#555555", lw=.8)
    ax.set(xlabel=r"Switching $\rho$", ylabel="Margin", title="C  Balanced versus unbalanced")
    ax.legend(frameon=False, fontsize=7, ncol=2)
    ax.grid(color="#dddddd", lw=.5)

    ax = axes[1, 1]
    q = np.array([[float(r["q1"]), float(r["q2"]), float(r["q3"])] for r in simplex])
    xy = np.column_stack((.5 * q[:, 0] + q[:, 2], np.sqrt(3) / 2 * q[:, 0]))
    vertices = np.array([[.5, np.sqrt(3)/2], [0, 0], [1, 0], [.5, np.sqrt(3)/2]])
    ax.plot(vertices[:, 0], vertices[:, 1], color="#333333", lw=.8)
    ax.plot([.5, .5], [0, np.sqrt(3)/2], ls="--", color="#777777", lw=.8)
    sc = ax.scatter(xy[:, 0], xy[:, 1], c=[float(r["intermediate_length"]) for r in simplex],
                    cmap="viridis", s=34, edgecolors="white", linewidths=.3, zorder=3)
    ax.text(.5, np.sqrt(3)/2 + .06, "q1", ha="center", fontsize=8)
    ax.text(-.03, -.04, "q2", ha="right", fontsize=8)
    ax.text(1.03, -.04, "q3", ha="left", fontsize=8)
    ax.set(xlim=(-.08,1.08), ylim=(-.08,.98), aspect="equal", title="D  Exact conditional path, Z=0101")
    ax.axis("off")
    fig.colorbar(sc, ax=ax, fraction=.046, pad=.04, label=r"$T_I$")

    for axis in axes.flat:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.tick_params(labelsize=8)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    root = Path(__file__).parent / "results"
    parser.add_argument("--region", type=Path, default=root / "anomaly_region_grid.tsv")
    parser.add_argument("--thresholds", type=Path, default=root / "anomaly_switching_thresholds.tsv")
    parser.add_argument("--simplex", type=Path, default=root / "anomaly_simplex_path.tsv")
    parser.add_argument("--output", type=Path, default=root / "anomaly_conditions_validation")
    args = parser.parse_args()
    build_figure(read(args.region), read(args.thresholds), read(args.simplex), args.output)


if __name__ == "__main__":
    main()
