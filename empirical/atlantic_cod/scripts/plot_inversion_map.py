#!/usr/bin/env python3
"""Plot Atlantic cod inversion spans from published breakpoint intervals.

Breakpoint coordinates and assembly provenance are recorded in
 data/atlantic_cod/metadata/inversion_breakpoints.tsv. For each inversion,
 the two published breakpoint intervals bound the rearranged span. Track
 extents use the published 250-kb tree-window data where available; they are
 observed spans, not assembly chromosome lengths.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MultipleLocator

ROOT = Path(__file__).resolve().parents[3]
BREAKPOINT_FILE = ROOT / "data/atlantic_cod/metadata/inversion_breakpoints.tsv"
WINDOW_FILE = ROOT / "data/atlantic_cod/processed/cod_window_trees.tsv"
OUTPUT_DIR = ROOT / "empirical/atlantic_cod/results"
FOCAL_LGS = ("LG01", "LG02", "LG07", "LG12")
COLORS = {"LG01_inv1": "#D1495B", "LG01_inv2": "#8F4A8E", "LG02_inv1": "#D1495B", "LG07_inv1": "#D1495B", "LG12_inv1": "#D1495B"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_data():
    inversions: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(BREAKPOINT_FILE):
        if row["lg"] in FOCAL_LGS:
            for field in ("left_break_start", "left_break_end", "right_break_start", "right_break_end"):
                row[field] = int(row[field])
            inversions[row["lg"]].append(row)

    observed_ends = {lg: 0 for lg in FOCAL_LGS}
    for row in read_tsv(WINDOW_FILE):
        if row["lg"] in observed_ends:
            observed_ends[row["lg"]] = max(observed_ends[row["lg"]], int(row["end"]))

    for lg in FOCAL_LGS:
        expected = 2 if lg == "LG01" else 1
        if len(inversions[lg]) != expected:
            raise ValueError(f"Expected {expected} inversion(s) for {lg}; found {len(inversions[lg])}")
        if observed_ends[lg] == 0:
            raise ValueError(f"No published window span found for {lg}")
        farthest_break = max(x["right_break_end"] for x in inversions[lg])
        observed_ends[lg] = max(observed_ends[lg], farthest_break)
    return inversions, observed_ends


def format_boundary(start: int, end: int) -> str:
    return f"{start:,}" if start == end else f"{start:,}–{end:,}"


def make_figure(output_dir: Path = OUTPUT_DIR) -> tuple[Path, Path]:
    inversions, observed_ends = load_data()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, ax = plt.subplots(figsize=(11.5, 5.4))
    y_positions = dict(zip(FOCAL_LGS, range(len(FOCAL_LGS) - 1, -1, -1)))
    max_span = max(observed_ends.values()) / 1e6

    for lg in FOCAL_LGS:
        y = y_positions[lg]
        span = observed_ends[lg] / 1e6
        ax.plot([0, span], [y, y], color="#A5AFB9", lw=8, solid_capstyle="round", zorder=1)
        for j, row in enumerate(inversions[lg]):
            left_lo, left_hi = row["left_break_start"], row["left_break_end"]
            right_lo, right_hi = row["right_break_start"], row["right_break_end"]
            color = COLORS[row["inversion_id"]]
            ybar = y + (0.10 if j else 0)
            # Fill between the outer edges of the two published breakpoint intervals.
            ax.plot([left_lo / 1e6, right_hi / 1e6], [ybar, ybar], color=color, lw=10,
                    solid_capstyle="butt", zorder=2)
            # Show the breakpoint uncertainty intervals explicitly as narrow bands.
            for lo, hi in ((left_lo, left_hi), (right_lo, right_hi)):
                ax.fill_betweenx([ybar - 0.085, ybar + 0.085], lo / 1e6, hi / 1e6,
                                 color=color, alpha=0.24, zorder=2.5)
                ax.plot([lo / 1e6, hi / 1e6], [ybar, ybar], color="#602D39", lw=2.1, zorder=3)

            label = (f"{row['inversion_id'].replace(lg + '_', '').upper()}: "
                     f"{format_boundary(left_lo, left_hi)} – {format_boundary(right_lo, right_hi)} bp")
            # Place the label above its own interval; on LG01 stagger the two labels.
            midpoint = (left_lo + right_hi) / 2e6
            label_y = y + 0.30 + 0.26 * j
            ax.annotate(label, xy=(midpoint, ybar), xytext=(midpoint, label_y),
                        ha="center", va="bottom", fontsize=7.8, color=color,
                        arrowprops={"arrowstyle": "-", "color": color, "lw": 0.7}, zorder=4)

        assembly = inversions[lg][0]["assembly"]
        ax.text(-0.35, y, f"{lg}\n{assembly}", va="center", ha="right", fontweight="bold", fontsize=8.3)
        ax.text(span + 0.35, y, f"{span:.2f} Mb data span", va="center", ha="left",
                color="#59636E", fontsize=8)

    ax.set_xlim(0, max_span + 5.2)
    ax.set_ylim(-0.63, len(FOCAL_LGS) - 0.04)
    ax.set_yticks([])
    ax.set_xlabel("Linkage-group coordinate (Mb)", labelpad=10)
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _pos: f"{x:g}"))
    ax.grid(axis="x", color="#E5E9ED", lw=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#B8C0C8")
    ax.tick_params(axis="x", colors="#59636E")
    fig.suptitle("Atlantic cod inversion locations", x=0.13, y=0.98, ha="left", fontsize=15, fontweight="bold")
    fig.text(0.13, 0.925,
             "Colored spans are bounded by the published breakpoint intervals; narrow shaded bands show breakpoint uncertainty.",
             ha="left", va="top", fontsize=8.6, color="#59636E")
    fig.text(0.13, 0.025,
             "LG01 breakpoint coordinates are from gadMor_Celtic (Kirubakaran et al. 2020); LG02/LG07/LG12 are gadMor2 (Matschiner et al. 2022). "
             "Track ends use observed 250-kb window spans, not assembly chromosome lengths.",
             ha="left", va="bottom", fontsize=8, color="#59636E", wrap=True)
    fig.subplots_adjust(left=0.17, right=0.88, top=0.82, bottom=0.18)

    output_dir.mkdir(parents=True, exist_ok=True)
    png = output_dir / "cod_inversion_map.png"
    pdf = output_dir / "cod_inversion_map.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png, pdf


if __name__ == "__main__":
    for path in make_figure():
        print(path)
