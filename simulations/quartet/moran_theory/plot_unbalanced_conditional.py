"""Plot exact arrangement-conditional pectinate quartet probabilities."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

PATTERN_COLORS = {
    "all_same": "#4C78A8", "q1_pair": "#F58518", "q2_pair": "#54A24B",
    "q3_pair": "#E45756", "singleton_or_triple": "#B279A2",
}
TOPOLOGY_COLORS = {"q1": "#0072B2", "q2": "#D55E00", "q3": "#009E73"}


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pattern_class(pattern: str) -> str:
    if pattern in {"0000", "1111"}:
        return "all_same"
    if pattern in {"0011", "1100"}:
        return "q1_pair"
    if pattern in {"0101", "1010"}:
        return "q2_pair"
    if pattern in {"0110", "1001"}:
        return "q3_pair"
    return "singleton_or_triple"


def simplex_xy(q: np.ndarray) -> tuple[float, float]:
    # Top q1=1, lower-left q2=1, lower-right q3=1.
    return float(q[0] * .5 + q[2]), float(q[0] * np.sqrt(3) / 2)


def draw_simplex(ax, rows: list[dict[str, str]], value_prefix: str = "") -> None:
    vertices = np.array([[.5, np.sqrt(3) / 2], [0, 0], [1, 0], [.5, np.sqrt(3) / 2]])
    ax.plot(vertices[:, 0], vertices[:, 1], color="#333333", lw=.8)
    # q2=q3 is the vertical median in this orientation.
    ax.plot([.5, .5], [0, np.sqrt(3) / 2], ls="--", lw=.8, color="#777777")
    for row in rows:
        q = np.array([float(row[f"{value_prefix}q{i}"]) for i in (1, 2, 3)])
        x, y = simplex_xy(q)
        pz = float(row["pattern_probability"])
        ax.scatter(x, y, s=30 + 600 * pz, color=PATTERN_COLORS[pattern_class(row["pattern"])],
                   alpha=.78, edgecolors="white", linewidths=.25, zorder=3)
    ax.text(.5, np.sqrt(3) / 2 + .065, "q1", ha="center", fontsize=8)
    ax.text(-.035, -.045, "q2", ha="right", fontsize=8)
    ax.text(1.035, -.045, "q3", ha="left", fontsize=8)
    ax.set(xlim=(-.08, 1.08), ylim=(-.08, .98), aspect="equal")
    ax.axis("off")


def build_figure(rows: list[dict[str, str]], parameters: list[dict[str, str]], output: Path) -> None:
    base = [row for row in rows if row["condition"] == "base"]
    patterns = [row["pattern"] for row in base]
    x = np.arange(len(patterns))
    fig = plt.figure(figsize=(11.2, 8.5), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=(1.1, 1.35))

    ax = fig.add_subplot(grid[0, 0])
    draw_simplex(ax, base)
    ax.set_title("A  Conditional quartet simplex", loc="left", fontsize=11, pad=8)
    handles = [plt.Line2D([], [], marker="o", ls="", color=color, label=label.replace("_", " "), ms=6)
               for label, color in PATTERN_COLORS.items()]
    ax.legend(handles=handles, frameon=False, fontsize=7, loc="lower center", bbox_to_anchor=(.5, -.17), ncol=2)

    ax = fig.add_subplot(grid[0, 1])
    delta = np.asarray([float(row["delta_q2_minus_q3"]) for row in base])
    gamma = np.asarray([float(row["gamma_anomaly_margin"]) for row in base])
    ax.axhline(0, color="#555555", lw=.8)
    ax.bar(x - .18, delta, width=.36, color="#0072B2", alpha=.82, label=r"$\Delta=q_2-q_3$")
    ax.bar(x + .18, gamma, width=.36, color="#D55E00", alpha=.82, label=r"$\Gamma=\max(q_2,q_3)-q_1$")
    ax.set_xticks(x, patterns, rotation=45, ha="right")
    ax.set(ylabel="Conditional probability margin", title="B  Discordance and anomaly")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(axis="y", color="#dddddd", lw=.5)

    ax = fig.add_subplot(grid[1, 0])
    for parameter, title in (("rho", r"Switching $\rho$"), ("intermediate_length", r"Intermediate duration $T_I$")):
        pass
    for topology in ("q1", "q2", "q3"):
        sub = [r for r in parameters if r["parameter"] == "rho" and r["topology"] == topology]
        ax.plot([float(r["value"]) for r in sub], [float(r["probability"]) for r in sub], marker="o", ms=3.5,
                lw=1.1, color=TOPOLOGY_COLORS[topology], label=topology)
    ax.set(xlabel=r"Switching $\rho$", ylabel="P(Q | Z=0101)", title="C  Parameter dependence")
    ax.grid(color="#dddddd", lw=.5)
    ax.legend(frameon=False, fontsize=7, ncol=3)

    ax2 = fig.add_subplot(grid[1, 1])
    for topology in ("q1", "q2", "q3"):
        sub = [r for r in parameters if r["parameter"] == "intermediate_length" and r["topology"] == topology]
        ax2.plot([float(r["value"]) for r in sub], [float(r["probability"]) for r in sub], marker="o", ms=3.5,
                 lw=1.1, color=TOPOLOGY_COLORS[topology], label=topology)
    ax2.set(xlabel=r"Intermediate duration $T_I$", ylabel="P(Q | Z=0101)", title="D  Matched shape comparison")
    # Add the matched balanced values at the same T_I values.
    for topology in ("q1", "q2", "q3"):
        sub = [r for r in rows if r["condition"] == "base" and r["pattern"] == "0101"]
        ax2.axhline(float(sub[0][f"balanced_{topology}"]), color=TOPOLOGY_COLORS[topology], ls=":", lw=.9)
    ax2.text(.98, .03, "solid: unbalanced; dotted: balanced\nmatched tip depths", transform=ax2.transAxes,
             ha="right", va="bottom", fontsize=7)
    ax2.grid(color="#dddddd", lw=.5)
    for axis in fig.axes:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.tick_params(labelsize=8)
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(output.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    root = Path(__file__).parent / "results"
    parser.add_argument("--input", type=Path, default=root / "unbalanced_conditional_analytic.tsv")
    parser.add_argument("--parameters", type=Path, default=root / "unbalanced_conditional_parameters.tsv")
    parser.add_argument("--output", type=Path, default=root / "unbalanced_conditional_validation")
    args = parser.parse_args()
    build_figure(read(args.input), read(args.parameters), args.output)


if __name__ == "__main__":
    main()
