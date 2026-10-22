#!/usr/bin/env python3
"""Render publication figures for the Neoaves MSRC threshold analysis.

The figure is deliberately limited to precomputed tables: no tree inference or
quartet enumeration is performed here.
"""
from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parents[1] / "figures" / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "empirical" / "neoaves_chr4"
RESULTS = BASE / "results" / "stage4e"
FIGURES = BASE / "figures"
LOCUS_TABLE = BASE / "results" / "locus_table_chr4.tsv"
MANIFEST = BASE / "results" / "stage4d_locus_manifest.tsv"

QCOLORS = {"q1": "#222222", "q2": "#2166ac", "q3": "#b2182b"}
TREATMENT_LABELS = {
    "J48_T0": "T0",
    "J48_T_PNAS": "PNAS",
    "J48_T_STRUCT_PRIMARY": "struct.",
    "J48_T_BLOCK_500K": "500 kb",
    "J48_NONCHR4": "non-chr4",
}
TREATMENTS = tuple(TREATMENT_LABELS)
PNAS_INTERVALS = ((25.03, 32.67), (33.51, 34.47), (44.13, 56.81))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def rolling(values: list[float], window: int = 101) -> list[float]:
    """Centered rolling mean, with bounded windows at the ends."""
    if not values:
        return []
    half = window // 2
    prefix = [0.0]
    for value in values:
        prefix.append(prefix[-1] + value)
    out = []
    for i in range(len(values)):
        lo, hi = max(0, i - half), min(len(values), i + half + 1)
        out.append((prefix[hi] - prefix[lo]) / (hi - lo))
    return out


def merge_spans(rows: list[dict[str, str]], flag: str) -> list[tuple[float, float]]:
    intervals = sorted(
        (float(row["start"]) / 1e6, float(row["end"]) / 1e6)
        for row in rows
        if row["chromosome"] == "chr4" and row[flag].lower() == "true"
    )
    merged: list[list[float]] = []
    for start, end in intervals:
        if merged and start <= merged[-1][1] + 0.01:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def load_tracks() -> dict[str, dict[str, list[float]]]:
    groups: dict[str, list[dict[str, str]]] = {"N61": [], "N62": []}
    for row in read_tsv(LOCUS_TABLE):
        if row["clade"] in groups and row["Chromosome"] == "chr4":
            groups[row["clade"]].append(row)
    tracks = {}
    for clade, rows in groups.items():
        rows.sort(key=lambda row: float(row["midpoint"]))
        tracks[clade] = {
            "x": [float(row["midpoint"]) / 1e6 for row in rows],
            **{q: rolling([float(row[q]) for row in rows]) for q in QCOLORS},
        }
    return tracks


def load_mixtures() -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    for row in read_tsv(RESULTS / "genomewide_quartet_thresholds.tsv"):
        if row["comparison"] != "chr4_vs_nonchr4" or row["pair"] != "q1=q2":
            continue
        result[row["clade"]] = {
            "background": [float(row[f"background_q{i}"]) for i in (1, 2, 3)],
            "affected": [float(row[f"affected_q{i}"]) for i in (1, 2, 3)],
            "observed": float(row["observed_weight"]),
            "epsilon": float(row["epsilon_star"]) if row["status"] == "valid_interior" else None,
            "status": row["status"],
        }
    return result


def load_j48() -> dict[str, dict[str, object]]:
    rows = read_tsv(RESULTS / "j48_topology_comparison.tsv")
    rf = {t: 0 for t in TREATMENTS}
    for row in rows:
        if row["left"] == "J48_T0":
            rf[row["right"]] = int(row["rf_symmetric_difference"])
        elif row["right"] == "J48_T0":
            rf[row["left"]] = int(row["rf_symmetric_difference"])
    return {
        treatment: {"rf": rf[treatment], "class": "baseline" if rf[treatment] == 0 else "alternative"}
        for treatment in TREATMENTS
    }


def style_axis(axis):
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.tick_params(labelsize=8, length=3)
    axis.grid(axis="y", color="#d9d9d9", linewidth=0.5, alpha=0.7)


def plot_tracks(axis, tracks, pnas_spans, structural_spans):
    for start, end in pnas_spans:
        axis.axvspan(start, end, color="#e69f00", alpha=0.18, lw=0, zorder=0)
    for start, end in structural_spans:
        axis.axvspan(start, end, color="#7b3294", alpha=0.12, lw=0, zorder=0)
    offsets = {"N61": 1.05, "N62": 0.0}
    for clade in ("N61", "N62"):
        offset = offsets[clade]
        track = tracks[clade]
        for q, color in QCOLORS.items():
            axis.plot(track["x"], [offset + value for value in track[q]], color=color, lw=0.75)
        axis.axhline(offset + 1 / 3, color="#888888", lw=0.5, ls=":")
        axis.text(0.25, offset + 0.91, clade, fontsize=8, fontweight="bold", va="top")
    axis.set_xlim(0, 60)
    axis.set_ylim(-0.02, 2.05)
    axis.set_yticks([0.0, 0.33, 0.67, 1.05, 1.38, 1.72])
    axis.set_yticklabels(["0", "1/3", "2/3", "0", "1/3", "2/3"])
    axis.set_ylabel("QQS support", fontsize=8)
    axis.set_xlabel("GalGal4 chr4 position (Mb)", fontsize=8)
    style_axis(axis)


def plot_mixture(axis, clade: str, mixture: dict[str, object], show_ylabel=True):
    background = mixture["background"]
    affected = mixture["affected"]
    epsilon = [i / 200 for i in range(201)]
    for index, q in enumerate(("q1", "q2", "q3")):
        values = [(1 - e) * background[index] + e * affected[index] for e in epsilon]
        axis.plot(epsilon, values, color=QCOLORS[q], lw=1.4, label=q)
    observed = float(mixture["observed"])
    axis.axvline(observed, color="#777777", lw=0.9, ls="--")
    axis.text(observed + 0.015, 0.98, "observed", transform=axis.get_xaxis_transform(), rotation=90,
              va="top", ha="left", fontsize=7, color="#666666")
    crossover = mixture["epsilon"]
    if crossover is not None:
        axis.axvline(float(crossover), color="#7b3294", lw=1.0, ls=":")
        axis.text(float(crossover) + 0.015, 0.06, rf"$\epsilon^*={float(crossover):.3f}$",
                  transform=axis.get_xaxis_transform(), rotation=90, va="bottom", fontsize=7,
                  color="#7b3294")
    else:
        axis.text(0.03, 0.08, "no interior crossover", transform=axis.transAxes, fontsize=7,
                  color="#555555")
    axis.set_title(clade, loc="left", fontsize=9, fontweight="bold", pad=3)
    axis.set_xlim(0, 1)
    axis.set_ylim(0.20, 0.52)
    axis.set_xticks([0, 0.25, 0.5, 0.75, 1])
    axis.set_xlabel(r"Affected-region mixture weight $\epsilon$", fontsize=8)
    if show_ylabel:
        axis.set_ylabel("Quartet support", fontsize=8)
    style_axis(axis)


def plot_outcomes(axis, j48):
    baseline = "#4d9221"
    alternative = "#c51b7d"
    for i, treatment in enumerate(TREATMENTS):
        entry = j48[treatment]
        color = baseline if entry["class"] == "baseline" else alternative
        axis.bar(i, 1, width=0.86, color=color, edgecolor="white", linewidth=0.8)
        axis.text(i, 0.5, "baseline" if entry["class"] == "baseline" else "alternative",
                  ha="center", va="center", fontsize=7, color="white", rotation=90)
        axis.text(i, -0.13, f"RF {entry['rf']}", ha="center", va="top", fontsize=7, color="#444444")
    axis.set_xticks(range(len(TREATMENTS)), [TREATMENT_LABELS[t] for t in TREATMENTS])
    axis.set_ylim(-0.25, 1.05)
    axis.set_yticks([])
    axis.set_ylabel("Inferred topology class", fontsize=8)
    axis.set_title("C  Completed 48-taxon ASTRAL treatments", loc="left", fontsize=10, fontweight="bold")
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.tick_params(axis="x", labelsize=8, length=0)


def save_figure(fig, stem: str, legacy: bool = False):
    names = [stem] + (["neoaves_msrc_threshold"] if legacy else [])
    for name in names:
        fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight")
        fig.savefig(FIGURES / f"{name}.png", dpi=300, bbox_inches="tight")
        fig.savefig(FIGURES / f"{name}.svg", bbox_inches="tight")


def make_main(tracks, mixtures, j48, pnas_spans, structural_spans):
    fig = plt.figure(figsize=(7.25, 6.8), constrained_layout=True)
    grid = fig.add_gridspec(3, 1, height_ratios=(1.45, 1.25, 0.82), hspace=0.25)
    axis_a = fig.add_subplot(grid[0])
    axis_a.set_title("A  Spatial quartet signal on chromosome 4", loc="left", fontsize=10, fontweight="bold")
    plot_tracks(axis_a, tracks, pnas_spans, structural_spans)
    subgrid = grid[1].subgridspec(1, 2, wspace=0.12)
    axis_b1, axis_b2 = fig.add_subplot(subgrid[0]), fig.add_subplot(subgrid[1])
    plot_mixture(axis_b1, "N61", mixtures["N61"])
    plot_mixture(axis_b2, "N62", mixtures["N62"], show_ylabel=False)
    axis_b1.set_title("B  Empirical quartet-mixture curves | N61", loc="left", fontsize=9, fontweight="bold", pad=3)
    axis_b2.set_title("N62", loc="left", fontsize=9, fontweight="bold", pad=3)
    axis_b2.tick_params(labelleft=False)
    axis_b1.legend(handles=[Line2D([0], [0], color=QCOLORS[q], lw=1.4, label=q) for q in QCOLORS],
                   loc="upper right", frameon=False, fontsize=7, ncol=3, handlelength=1.5)
    axis_c = fig.add_subplot(grid[2])
    plot_outcomes(axis_c, j48)
    axis_a.legend(handles=[Patch(facecolor="#e69f00", alpha=0.3, label="PNAS outlier intervals"),
                           Patch(facecolor="#7b3294", alpha=0.22, label="structural-primary mask")],
                  loc="upper right", frameon=False, fontsize=7)
    save_figure(fig, "neoaves_msrc_threshold_main", legacy=True)
    plt.close(fig)


def make_supplement(tracks, mixtures, pnas_spans, structural_spans):
    fig, axes = plt.subplots(2, 2, figsize=(7.25, 6.0), constrained_layout=True)
    axes[0, 0].set_title("A  Columbea mixture", loc="left", fontsize=9, fontweight="bold")
    plot_mixture(axes[0, 0], "Columbea", mixtures["Columbea"])
    axes[0, 0].set_title("A  Columbea mixture", loc="left", fontsize=9, fontweight="bold", pad=3)
    axes[0, 0].legend(handles=[Line2D([0], [0], color=QCOLORS[q], lw=1.4, label=q) for q in QCOLORS],
                      loc="upper right", frameon=False, fontsize=7, ncol=3, handlelength=1.5)
    axes[0, 1].set_title("B  N61 chr4 track", loc="left", fontsize=9, fontweight="bold")
    plot_single_track(axes[0, 1], tracks["N61"], pnas_spans, structural_spans)
    axes[1, 0].set_title("C  N62 chr4 track", loc="left", fontsize=9, fontweight="bold")
    plot_single_track(axes[1, 0], tracks["N62"], pnas_spans, structural_spans)
    axes[1, 1].axis("off")
    axes[1, 1].text(0.04, 0.92, "Interpretation", fontsize=9, fontweight="bold", va="top")
    axes[1, 1].text(0.04, 0.78, "N61 and N62 have interior\nq1=q2 crossovers; Columbea\ndoes not.", fontsize=8, va="top")
    axes[1, 1].text(0.04, 0.48, "QQS endpoints are recovered\npublished full-source statistics,\nnot exact 48-taxon quartet\nestimates.", fontsize=8, va="top", color="#555555")
    save_figure(fig, "neoaves_msrc_threshold_supplement")
    plt.close(fig)


def plot_single_track(axis, track, pnas_spans, structural_spans):
    for start, end in pnas_spans:
        axis.axvspan(start, end, color="#e69f00", alpha=0.18, lw=0)
    for start, end in structural_spans:
        axis.axvspan(start, end, color="#7b3294", alpha=0.12, lw=0)
    for q, color in QCOLORS.items():
        axis.plot(track["x"], track[q], color=color, lw=0.8, label=q)
    axis.axhline(1 / 3, color="#888888", lw=0.5, ls=":")
    axis.set_xlim(0, 60)
    axis.set_ylim(0.0, 0.9)
    axis.set_xlabel("GalGal4 chr4 position (Mb)", fontsize=8)
    axis.set_ylabel("QQS support", fontsize=8)
    style_axis(axis)


def main():
    FIGURES.mkdir(parents=True, exist_ok=True)
    manifest = read_tsv(MANIFEST)
    pnas_spans = PNAS_INTERVALS
    structural_spans = merge_spans(manifest, "frozen_structural_primary")
    tracks = load_tracks()
    mixtures = load_mixtures()
    j48 = load_j48()
    make_main(tracks, mixtures, j48, pnas_spans, structural_spans)
    make_supplement(tracks, mixtures, pnas_spans, structural_spans)


if __name__ == "__main__":
    main()
