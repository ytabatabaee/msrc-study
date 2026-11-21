#!/usr/bin/env python3
"""Regenerate Atlantic cod quartet-support and assembly-consistent interval figures.

The gene-tree windows and frozen high-LD masks use gadMor2 coordinates.
Breakpoints from other assemblies are only plotted if converted coordinates
are present and their audit status verifies the target assembly.
"""
from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MultipleLocator

ROOT = Path(__file__).resolve().parents[3]
AUDIT = ROOT / "empirical/atlantic_cod/results/coordinate_assembly_audit.tsv"
WINDOWS = ROOT / "data/atlantic_cod/processed/stage6_local_quartet_support.tsv"
WINDOW_TREES = ROOT / "data/atlantic_cod/processed/cod_window_trees.tsv"
REGIONS = ROOT / "data/atlantic_cod/metadata/region_manifest.tsv"
FIGURES = ROOT / "empirical/atlantic_cod/figures"
RESULTS = ROOT / "empirical/atlantic_cod/results"
LGS = ("LG01", "LG02", "LG07", "LG12")
COLORS = {"q_arrangement": "#2ca02c", "q_baseline": "#9467bd", "q_third": "#8c564b"}
LABELS = {"q_arrangement": "arrangement", "q_baseline": "baseline", "q_third": "third"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def verified_inversions():
    out: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
    unresolved: dict[str, list[str]] = defaultdict(list)
    for row in read_tsv(AUDIT):
        if not row["feature_id"].endswith("_inv1") and not row["feature_id"].endswith("_inv2"):
            continue
        if (row["conversion_status"] == "same_assembly_no_conversion"
                and row["original_assembly"] == "gadMor2"
                and row["target_assembly"] == "gadMor2"):
            # Build a full inversion span only after both endpoints verify on the target build.
            continue
        unresolved[row["lg"]].append(row["feature_id"])
    audit = read_tsv(AUDIT)
    by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in audit:
        if "_inv" in row["feature_id"]:
            by_id[row["feature_id"]].append(row)
    for inv_id, rows in by_id.items():
        left = next((r for r in rows if r["boundary_id"] == "left"), None)
        right = next((r for r in rows if r["boundary_id"] == "right"), None)
        if not left or not right:
            continue
        if any(r["conversion_status"] != "same_assembly_no_conversion"
               or r["original_assembly"] != "gadMor2"
               or r["target_assembly"] != "gadMor2" for r in (left, right)):
            continue
        lg = left["lg"]
        out[lg].append((int(left["converted_start"]), int(right["converted_end"]), inv_id))
    return out, unresolved


def load_regions():
    return {r["lg"]: (int(r["start"]), int(r["end"])) for r in read_tsv(REGIONS) if r["lg"] in LGS}


def make_support_figure() -> tuple[Path, Path]:
    rows = read_tsv(WINDOWS)
    regions = load_regions()
    span: dict[str, int] = {lg: 0 for lg in LGS}
    for row in read_tsv(WINDOW_TREES):
        if row["lg"] in span:
            span[row["lg"]] = max(span[row["lg"]], int(row["end"]))
    valid, unresolved = verified_inversions()

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(4, 1, figsize=(7.5, 8.0), sharey=True)
    for ax, lg in zip(axes, LGS):
        panel = sorted((r for r in rows if r["lg"] == lg), key=lambda r: int(r["midpoint"]))
        x = [int(r["midpoint"]) / 1e6 for r in panel]
        # This interval comes from the same gadMor2 window/coordinate workflow.
        lo, hi = regions[lg]
        ax.axvspan(lo / 1e6, hi / 1e6, color="#d9d9d9", alpha=0.58, label="cod_phylogenomics supergene mask", zorder=0)
        for start, end, inv_id in valid.get(lg, []):
            ax.axvspan(start / 1e6, end / 1e6, color="#e76f51", alpha=0.16,
                       label="breakpoint-bounded inversion" if lg == "LG02" else None, zorder=0.2)
        for key in COLORS:
            ax.plot(x, [float(r[key]) for r in panel], lw=1.0, marker="o", ms=2,
                    color=COLORS[key], label=LABELS[key], zorder=2)
        ax.set_ylim(-0.03, 1.03)
        ax.set_ylabel(lg)
        ax.axhline(0.5, color="black", lw=0.4, alpha=0.5)
        ax.set_xlim(0, span[lg] / 1e6)
    axes[0].legend(loc="upper right", fontsize=7.5, ncol=4, framealpha=0.92)
    axes[-1].set_xlabel("gadMor2 genomic position (Mb)")
    fig.suptitle("Atlantic cod local quartet-support fractions", y=0.995)
    fig.text(0.5, 0.006,
             "Gray = frozen cod_phylogenomics supergene masks (gadMor2). LG01 contains two adjacent inversions; "
             "their gadMor_Celtic breakpoints are not plotted because no reliable lift-over to gadMor2 is available.",
             ha="center", va="bottom", fontsize=7.1, color="#4d5964")
    fig.tight_layout(rect=(0, 0.035, 1, 0.98))
    FIGURES.mkdir(parents=True, exist_ok=True)
    png = FIGURES / "atlantic_cod_quartet_support_tracks.png"
    pdf = FIGURES / "atlantic_cod_quartet_support_tracks.pdf"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def make_interval_map() -> tuple[Path, Path]:
    regions = load_regions()
    valid, unresolved = verified_inversions()
    ends = {lg: max(int(r["end"]) for r in read_tsv(WINDOW_TREES) if r["lg"] == lg) for lg in LGS}
    fig, ax = plt.subplots(figsize=(10.5, 4.5))
    ys = dict(zip(LGS, (3, 2, 1, 0)))
    for lg in LGS:
        y = ys[lg]
        ax.plot([0, ends[lg] / 1e6], [y, y], color="#a5afb9", lw=8, solid_capstyle="round", zorder=1)
        a, b = regions[lg]
        ax.plot([a / 1e6, b / 1e6], [y, y], color="#b7c0c8", lw=10, solid_capstyle="butt", zorder=2)
        for start, end, inv_id in valid.get(lg, []):
            ax.plot([start / 1e6, end / 1e6], [y, y], color="#d1495b", lw=7,
                    solid_capstyle="butt", zorder=3)
        if lg == "LG01":
            note = "2 adjacent inversions; gadMor_Celtic breakpoints unresolved on gadMor2"
            ax.text((a + b) / 2e6, y + 0.24, note, ha="center", va="bottom", fontsize=7.6, color="#8a3442")
        else:
            ax.text((a + b) / 2e6, y + 0.22, "inversion breakpoints", ha="center", va="bottom", fontsize=7.5, color="#8a3442")
        ax.text(-0.35, y, lg, ha="right", va="center", fontweight="bold")
        ax.text(ends[lg] / 1e6 + 0.3, y, f"{ends[lg]/1e6:.2f} Mb data span", va="center", fontsize=7.5, color="#59636e")
    ax.set_xlim(0, max(ends.values()) / 1e6 + 5.2)
    ax.set_ylim(-0.65, 3.65)
    ax.set_yticks([])
    ax.set_xlabel("gadMor2 genomic position (Mb)")
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _p: f"{x:g}"))
    ax.grid(axis="x", color="#e5e9ed", lw=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_title("Atlantic cod supergene envelopes and verified inversion spans", loc="left", fontweight="bold")
    fig.text(0.08, 0.015, "Gray bars: cod_phylogenomics supergene masks. Red bars: same-build breakpoint-bounded inversions. "
             "LG01 inversion bars are omitted pending lift-over; no coordinates are inferred.", fontsize=7.5, color="#59636e")
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    RESULTS.mkdir(parents=True, exist_ok=True)
    png = RESULTS / "cod_inversion_map.png"
    pdf = RESULTS / "cod_inversion_map.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png, pdf


if __name__ == "__main__":
    for path in (*make_support_figure(), *make_interval_map()):
        print(path)
