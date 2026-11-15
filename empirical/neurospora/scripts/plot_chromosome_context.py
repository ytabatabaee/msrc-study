#!/usr/bin/env python3
"""Plot published LGI intervals and locally inferred topology runs."""
from __future__ import annotations
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "empirical/neurospora/results"
FIGURES = ROOT / "empirical/neurospora/figures"
SUMMARY = RESULTS / "region_annotation_summary.tsv"
XMAX = 9.5

def read_summary():
    with SUMMARY.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def plot():
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.labelsize": 9,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "svg.fonttype": "none"})
    rows = read_summary()
    by = {(r["feature"], r["lineage_or_comparison"]): r for r in rows}
    fig, ax = plt.subplots(figsize=(8.6, 5.3))
    ax.set_xlim(0, XMAX)
    ax.set_ylim(-0.75, 8.65)
    ax.set_yticks([8, 7, 6, 5, 4, 3, 2, 1, 0])
    ax.set_yticklabels(["L10: recombination suppression", "L10: introgression-associated",
                        "L10: topology shift inferred here", "L10 CJ01 mat-A inversion",
                        "L9 965A inversion", "L8 mat-a inversion", "L7 mat-A inversion",
                        "L5 mat-A / L6 reference inversion", "L4 7585 inversion"])
    ax.xaxis.set_ticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    ax.set_xlabel("Linkage group I position (Mb; 2509a reference axis)")
    ax.grid(axis="x", color="#d5d8dc", lw=.6, zorder=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#687078")
    ax.tick_params(axis="y", length=0, pad=8)
    ax.tick_params(axis="x", colors="#39424e")

    def interval(feature, lineage, y, color, alpha=1.0, edge=None, hatch=None, height=.50):
        r = by[(feature, lineage)]
        left, right = int(r["start_bp"])/1e6, int(r["end_bp"])/1e6
        ax.barh(y, right-left, left=left, height=height, color=color, alpha=alpha,
                edgecolor=edge or color, linewidth=.8, hatch=hatch, zorder=3)
        return left, right

    interval("recombination_suppressed_region", "L10 mating-type chromosome", 8, "#b9d9e8", .85, "#5f98b2", height=.58)
    interval("introgression_impacted_region", "L10 mat-A; introgression from N. hispaniola", 7, "#f2bd72", .85, "#c98228", height=.58)
    # The own-analysis topology signal consists of two distinct q3-dominant runs.
    for start, end in ((1.2, 5.7), (5.8, 8.2)):
        ax.barh(6, end-start, left=start, height=.48, color="#7868a6", edgecolor="#594982", linewidth=.8, zorder=3)
    rearrangements = [
        ("L10 strain CJ01 mat-A", 5),
        ("L9 strain 965A", 4),
        ("L8 all mat-a strains", 3),
        ("L7 all mat-A strains", 2),
        ("L5 all mat-A strains; corresponding breakpoint pattern in L6 reference genome", 1),
        ("L4 strain 7585", 0),
    ]
    for lineage, y in rearrangements:
        interval("rearrangement_inversion", lineage, y, "#d9795b", .95, "#a8412b", height=.43)
    ax.text(1.38, 8, "1.38–7.58 Mb", ha="left", va="center", fontsize=7, color="#315d72")
    ax.text(1.38, 7, "L10 mat-A, 6.20 Mb", ha="left", va="center", fontsize=7, color="#80551c")
    ax.text(1.2, 6.33, "1.2–5.7", ha="left", va="bottom", fontsize=6.7, color="#594982")
    ax.text(8.2, 6.33, "5.8–8.2 Mb", ha="right", va="bottom", fontsize=6.7, color="#594982")
    for lineage, y in rearrangements:
        r = by[("rearrangement_inversion", lineage)]
        size = (int(r["end_bp"])-int(r["start_bp"]))/1e6
        ax.text(9.45, y, "{:.2f} Mb".format(size), ha="right", va="center", fontsize=7, color="#71301f")
    ax.legend(handles=[Patch(facecolor="#b9d9e8", edgecolor="#5f98b2", label="Published L10 suppression interval"),
                       Patch(facecolor="#f2bd72", edgecolor="#c98228", label="Published introgression-associated interval"),
                       Patch(facecolor="#7868a6", edgecolor="#594982", label="Topology runs inferred from local trees"),
                       Patch(facecolor="#d9795b", edgecolor="#a8412b", label="Sun et al. inversion interval")],
              loc="upper center", bbox_to_anchor=(.5, 1.13), ncol=2, frameon=False, fontsize=7.5,
              columnspacing=1.6, handlelength=1.7)
    fig.subplots_adjust(left=.28, right=.98, top=.84, bottom=.21)
    fig.text(.012, .025, "Lineage-specific rearrangements are plotted separately; L10's 3.63-Mb CJ01 inversion is not fixed across L10.\n"
             "Axis rounded to 9.5 Mb; local analysis reaches 9,464,227 bp. Suppression extent varies among lineages.",
             ha="left", va="bottom", fontsize=7, color="#4c5661")
    FIGURES.mkdir(parents=True, exist_ok=True)
    out = FIGURES / "figure0_chromosome_context"
    for suffix, kwargs in (("png", {"dpi": 300}), ("pdf", {}), ("svg", {})):
        fig.savefig(out.with_suffix("." + suffix), bbox_inches="tight", **kwargs)
    plt.close(fig)

if __name__ == "__main__":
    plot()
