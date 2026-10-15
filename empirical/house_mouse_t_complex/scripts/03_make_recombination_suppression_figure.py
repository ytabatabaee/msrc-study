#!/usr/bin/env python3
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
INTERVAL_COLORS = ["#6a3d9a", "#1b9e77", "#d95f02", "#7570b3"]

def bins(rows, width=250_000):
    groups = {}
    for row in rows:
        if row["treatment"] != "ALL_TIPS":
            continue
        start = 5_000_000 + ((int(row["start_bp"]) - 5_000_000) // width) * width
        groups.setdefault(start, []).append(row)
    return [(start + width / 2, {key: float(np.mean([float(r[key]) for r in block])) for key in Q_COLORS})
            for start, block in sorted(groups.items())]

def gaps():
    meta = sorted(read_tsv(PROCESSED / "house_mouse_t_complex_ml_5kb_metadata.tsv"), key=lambda r: int(r["start_bp"]))
    return [(int(a["end_bp"]) + 1, int(b["start_bp"]) - 1) for a, b in zip(meta, meta[1:])
            if int(b["start_bp"]) > int(a["end_bp"]) + 1]

def shade_gaps(ax):
    for start, end in gaps():
        ax.axvspan(start / 1e6, end / 1e6, color="#eeeeee", zorder=-5)

def axis(ax):
    ax.set_xlim(5, 40.5)
    ax.grid(axis="x", alpha=0.15)

def panel_a(ax):
    rows = read_tsv(RESULTS / "t_complex_inversion_intervals.tsv")
    for i, row in enumerate(rows):
        start, end = int(row["start_bp"]) / 1e6, int(row["end_bp"]) / 1e6
        ax.add_patch(plt.Rectangle((start, i + .18), end - start, .64, color=INTERVAL_COLORS[i], alpha=.72))
        ax.text((start + end) / 2, i + .5, row["inversion_name"], ha="center", va="center", fontsize=8)
    ax.set_ylim(0, 4); ax.set_yticks(np.arange(4) + .5, [r["inversion_name"] for r in rows]); ax.invert_yaxis()
    ax.set_ylabel("Published inversion")
    ax.set_title("A. Four t-complex inversions (sequence-validated mm10 intervals)", loc="left", fontsize=10)
    ax.text(40.35, 3.9, "junction-region uncertainty\nreported by source", ha="right", va="bottom", fontsize=7, color="#555555")
    axis(ax)

def panel_b(ax):
    shade_gaps(ax)
    t_mid, tf_mid = 8.438, 27.309
    ax.plot([t_mid, tf_mid], [6.5, 6.5], color="#555555", lw=2)
    ax.scatter([t_mid, tf_mid], [6.5, 6.5], color="#222222", s=18, zorder=3)
    ax.text(t_mid, 6.65, "T", ha="center", va="top", fontsize=8)
    ax.text(tf_mid, 6.65, "tf/Itpr3", ha="center", va="top", fontsize=8)
    ax.text((t_mid + tf_mid) / 2, 6.82, "historical T–tf pair; approximate modern mm10 anchors", ha="center", va="top", fontsize=7)
    ax.errorbar([t_mid], [9.5], yerr=[[2.5], [2.5]], fmt="none", ecolor="#b2182b", capsize=4, lw=2)
    ax.errorbar([t_mid], [0.3], yerr=[[0.2], [0.2]], fmt="none", ecolor="#2166ac", capsize=4, lw=2)
    ax.scatter([t_mid, t_mid], [9.5, 0.3], color=["#b2182b", "#2166ac"], zorder=3)
    ax.text(t_mid + .65, 9.5, "+/+ conventional map: 7–12 cM", va="center", fontsize=8, color="#b2182b")
    ax.text(t_mid + .65, 0.3, "+/t heterozygotes: 0.1–0.5 cM", va="center", fontsize=8, color="#2166ac")
    ax.text(5.2, 2.0, "Region-wide summary: four inversions ≈20 cM;\nreported suppression ≈100–200-fold (not a coordinate-resolved rate)", fontsize=8, va="top", bbox={"facecolor": "white", "edgecolor": "#999999"})
    ax.set_ylim(0, 11.5); ax.set_ylabel("Reported map distance (cM)")
    ax.set_title("B. Direct suppression measurements are marker-pair / region-wide", loc="left", fontsize=10)
    ax.text(40.35, 10.9, "No CI or window-level map recovered", ha="right", va="top", fontsize=7, color="#555555")
    axis(ax)

def panel_c(ax, rows):
    raw = [r for r in rows if r["treatment"] == "ALL_TIPS"]
    x = np.array([float(r["midpoint_bp"]) / 1e6 for r in raw])
    shade_gaps(ax)
    for key, color in Q_COLORS.items():
        ax.scatter(x, [float(r[key]) for r in raw], s=2, alpha=.07, color=color)
    b = bins(rows); xb = np.array([p / 1e6 for p, _ in b])
    for key, color in Q_COLORS.items():
        ax.plot(xb, [d[key] for _, d in b], lw=1.5, color=color, label=key)
    ax.axhline(1 / 3, color="#777777", ls="--", lw=.7)
    ax.set_ylim(0, 1); ax.set_ylabel("Quartet fraction")
    ax.set_title("C. Frozen ALL_TIPS quartet support; raw 5-kb windows and 250-kb means", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper right")
    axis(ax)

def panel_d(ax, rows):
    shade_gaps(ax)
    b = bins(rows); x = np.array([p / 1e6 for p, _ in b]); delta = np.array([d["q_t_alt"] - d["q_species"] for _, d in b])
    ax.plot(x, delta, color="#333333", lw=1.5)
    ax.axhline(0, color="#777777", lw=.8)
    ax.set_ylim(-.8, .8); ax.set_ylabel("D(w)"); ax.set_xlabel("chr17 position (Mb)")
    ax.set_title("D. Quartet contrast: D(w) = q_t-alt(w) − q_species(w)", loc="left", fontsize=10)
    ax.text(40.35, -.72, "Quartet association with the\ncoarse literature estimates is\nnot tested", ha="right", va="bottom", fontsize=7, color="#555555")
    axis(ax)

def main():
    rows = read_tsv(RESULTS / "stage2_fixed_quartet_scan.tsv")
    fig, axes = plt.subplots(4, 1, figsize=(13, 12), sharex=True, gridspec_kw={"height_ratios": [1.0, 1.8, 1.5, 1.0], "hspace": .22})
    panel_a(axes[0]); panel_b(axes[1]); panel_c(axes[2], rows); panel_d(axes[3], rows)
    fig.suptitle("House-mouse t-complex: independent recombination-suppression evidence and quartet support", fontsize=14)
    fig.savefig(OUT_PNG, dpi=220, bbox_inches="tight"); fig.savefig(OUT_PDF, bbox_inches="tight"); plt.close(fig)

if __name__ == "__main__":
    main()
