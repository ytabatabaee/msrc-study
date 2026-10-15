#!/usr/bin/env python3
"""Create the three-panel literature-plus-quartet suppression figure."""
from __future__ import annotations
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))
from house_mouse_stage2_utils import FIGURES, PROCESSED, RESULTS, read_tsv

OUT_PNG = FIGURES / "house_mouse_t_complex_suppression.png"
OUT_PDF = FIGURES / "house_mouse_t_complex_suppression.pdf"
Q_COLORS = {"q_species": "#4c78a8", "q_t_alt": "#f58518", "q_other": "#54a24b"}
INVERSION_COLORS = ["#6a3d9a", "#1b9e77", "#d95f02", "#7570b3"]
T_BP, TF_BP = 8_438_000, 27_300_000
COVERAGE_START_BP, COVERAGE_END_BP = 5_000_000, 37_055_000

def quartets(rows: list[dict[str, str]], width: int = 250_000):
    groups = {}
    for row in rows:
        if row["treatment"] != "ALL_TIPS":
            continue
        start = COVERAGE_START_BP + ((int(float(row["start_bp"])) - COVERAGE_START_BP) // width) * width
        groups.setdefault(start, []).append(row)
    return [(start + width / 2, {key: float(np.mean([float(r[key]) for r in block])) for key in Q_COLORS})
            for start, block in sorted(groups.items())]

def retained_gaps():
    metadata = sorted(read_tsv(PROCESSED / "house_mouse_t_complex_ml_5kb_metadata.tsv"), key=lambda row: int(row["start_bp"]))
    return [(int(left["end_bp"]) + 1, int(right["start_bp"]) - 1) for left, right in zip(metadata, metadata[1:])
            if int(right["start_bp"]) > int(left["end_bp"]) + 1]

def shade_gaps(ax):
    for start, end in retained_gaps():
        ax.axvspan(start / 1e6, end / 1e6, color="#eeeeee", zorder=-5)

def chromosome_axis(ax):
    ax.set_xlim(5, 40.5)
    ax.set_xticks(np.arange(5, 41, 5))
    ax.grid(axis="x", alpha=.16)

def panel_a(ax):
    intervals = read_tsv(RESULTS / "t_complex_inversion_intervals.tsv")
    for index, row in enumerate(intervals):
        start, end = int(row["start_bp"]) / 1e6, int(row["end_bp"]) / 1e6
        ax.plot([start, end], [index, index], color=INVERSION_COLORS[index], lw=10, solid_capstyle="butt", zorder=2)
        ax.text((start + end) / 2, index, row["inversion_name"], ha="center", va="center", fontsize=8, color="white", weight="bold", zorder=3)
    ax.axvspan(COVERAGE_START_BP / 1e6, COVERAGE_END_BP / 1e6, color="#bdbdbd", alpha=.12, zorder=0)
    ax.text((COVERAGE_START_BP + COVERAGE_END_BP) / 2e6, 3.52, "available 5-kb gene-tree windows: chr17:5.00–37.05 Mb", ha="center", va="bottom", fontsize=8, color="#555555")
    ax.plot([T_BP / 1e6, TF_BP / 1e6], [4.03, 4.03], color="#222222", lw=1.1)
    for position in (T_BP / 1e6, TF_BP / 1e6):
        ax.plot([position, position], [3.94, 4.12], color="#222222", lw=1.1)
    ax.text(T_BP / 1e6, 4.22, "T", ha="center", va="bottom", fontsize=9, weight="bold")
    ax.text(TF_BP / 1e6, 4.22, "tf", ha="center", va="bottom", fontsize=9, weight="bold")
    ax.text((T_BP + TF_BP) / 2e6, 4.10, "T–tf marker interval", ha="center", va="bottom", fontsize=8)
    ax.axvline(T_BP / 1e6, color="#333333", ls=":", lw=.8)
    ax.axvline(TF_BP / 1e6, color="#333333", ls=":", lw=.8)
    ax.annotate("In(17)4 extends beyond\ngene-tree coverage", xy=(40.0, 3), xytext=(34.0, 3.65),
                arrowprops={"arrowstyle": "-", "color": "#555555", "lw": .8}, fontsize=8, color="#555555", ha="center")
    ax.set_ylim(-.6, 4.55)
    ax.set_yticks(np.arange(len(intervals)), [row["inversion_name"] for row in intervals])
    ax.invert_yaxis()
    ax.set_ylabel("Published inversion")
    ax.set_title("A. Chromosome 17 structure", loc="left", fontsize=11, weight="bold")
    ax.set_xlabel("GRCm38/mm10 chr17 position (Mb)")
    chromosome_axis(ax)

def panel_b(ax):
    ax.plot([7, 12], [1, 1], color="#b2182b", lw=10, solid_capstyle="butt")
    ax.plot([.1, .5], [0, 0], color="#2166ac", lw=10, solid_capstyle="butt")
    for value, y, color in [(7, 1, "#b2182b"), (12, 1, "#b2182b"), (.1, 0, "#2166ac"), (.5, 0, "#2166ac")]:
        ax.plot([value, value], [y - .15, y + .15], color=color, lw=1.2)
    ax.text(12.25, 1, "7–12 cM", va="center", fontsize=10, color="#b2182b", weight="bold")
    ax.text(.75, 0, "0.1–0.5 cM", va="center", fontsize=10, color="#2166ac", weight="bold")
    ax.set_xlim(0, 13.5); ax.set_ylim(-.65, 1.65)
    ax.set_yticks([1, 0], ["Conventional map", "Naturally occurring t/+ heterozygotes"])
    ax.set_xlabel("T–tf genetic distance (cM)")
    ax.set_title("B. Recombination suppression between T and tf", loc="left", fontsize=11, weight="bold")
    ax.text(.02, .04, "Historical linkage-map comparison (Silver & Artzt, 1981).\nPublished ranges; not point estimates.",
            transform=ax.transAxes, fontsize=8.5, color="#444444", va="bottom")
    ax.grid(axis="x", alpha=.16)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

def panel_c(ax_q, ax_d, rows):
    raw = [row for row in rows if row["treatment"] == "ALL_TIPS"]
    x_raw = np.array([float(row["midpoint_bp"]) / 1e6 for row in raw])
    shade_gaps(ax_q); shade_gaps(ax_d)
    for key, color in Q_COLORS.items():
        ax_q.scatter(x_raw, [float(row[key]) for row in raw], s=2, alpha=.08, color=color)
    summary = quartets(rows)
    x_summary = np.array([midpoint / 1e6 for midpoint, _ in summary])
    labels = {"q_species": r"$q_{\mathrm{species}}$", "q_t_alt": r"$q_{t\mathrm{-alt}}$", "q_other": r"$q_{\mathrm{other}}$"}
    for key, color in Q_COLORS.items():
        ax_q.plot(x_summary, [values[key] for _, values in summary], color=color, lw=1.5, label=labels[key])
    ax_q.axhline(1 / 3, color="#777777", ls="--", lw=.7)
    ax_q.set_ylim(0, 1); ax_q.set_ylabel("Quartet fraction")
    ax_q.set_title("C. Local quartet support varies across the rearranged region", loc="left", fontsize=11, weight="bold")
    ax_q.legend(frameon=False, fontsize=9, ncol=3, loc="upper right")
    ax_q.text(.01, .04, "points = original 5-kb windows; lines = documented 250-kb means", transform=ax_q.transAxes, fontsize=8, color="#555555")
    delta = np.array([values["q_t_alt"] - values["q_species"] for _, values in summary])
    ax_d.plot(x_summary, delta, color="#333333", lw=1.5)
    ax_d.axhline(0, color="#777777", lw=.8)
    ax_d.set_ylim(-.8, .8); ax_d.set_ylabel(r"$D(w)$"); ax_d.set_xlabel("GRCm38/mm10 chr17 position (Mb)")
    ax_d.set_title(r"$D(w)=q_{t\mathrm{-alt}}(w)-q_{\mathrm{species}}(w)$", loc="left", fontsize=10)
    ax_d.text(.01, .06, r"$D(w)>0$: t-alt support exceeds species support; $D(w)<0$: the reverse", transform=ax_d.transAxes, fontsize=8, color="#555555")
    chromosome_axis(ax_q); chromosome_axis(ax_d)

def main():
    rows = read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv")
    fig = plt.figure(figsize=(12.5, 11.5))
    grid = fig.add_gridspec(3, 1, height_ratios=[1.25, .9, 2.0], hspace=.34)
    ax_a = fig.add_subplot(grid[0]); ax_b = fig.add_subplot(grid[1])
    c_grid = grid[2].subgridspec(2, 1, height_ratios=[1.8, .9], hspace=.16)
    ax_q = fig.add_subplot(c_grid[0]); ax_d = fig.add_subplot(c_grid[1], sharex=ax_q)
    panel_a(ax_a); panel_b(ax_b); panel_c(ax_q, ax_d, rows)
    ax_q.tick_params(axis="x", labelbottom=False)
    fig.suptitle("House-mouse t-complex: independent recombination suppression and heterogeneous local quartet support", fontsize=14, weight="bold", y=.985)
    fig.text(.01, .01, "Structural evidence and historical linkage measurements are independent of the local gene-tree quartet profiles.", fontsize=8, color="#555555")
    fig.subplots_adjust(left=.13, right=.98, top=.93, bottom=.07)
    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    fig.savefig(OUT_PDF, bbox_inches="tight")
    plt.close(fig)

if __name__ == "__main__":
    main()
