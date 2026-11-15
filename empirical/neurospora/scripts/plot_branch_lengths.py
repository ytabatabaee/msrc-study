#!/usr/bin/env python3
"""Plot binned Neurospora LGI quartet distances from the frozen tree inventory."""
from __future__ import annotations
import csv
import math
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "empirical/neurospora/results"
FIGURES = ROOT / "empirical/neurospora/figures"
INTERVAL = (1_500_000, 8_000_000)
VALUE_FIELDS = {"pair1_contrast": "pair1_A_hisp_minus_a_hisp", "pair2_contrast": "pair2_A_hisp_minus_a_hisp", "B": "pair1_B_mean_six_pairwise"}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def finite(value: str) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def aggregate_bins(rows: list[dict[str, str]], width: int, chromosome: str = "LGI") -> list[dict[str, float | int]]:
    """Average available per-window tree estimates on a complete coordinate grid.

    Empty bins are retained as NaN, so a line plot breaks at missing intervals.
    Each observed window contributes one estimate to its bin (window mean).
    """
    if width <= 0:
        raise ValueError("width must be positive")
    source = [r for r in rows if r.get("linkage_group") == chromosome]
    if not source:
        return []
    starts = [int(r["start"]) for r in source]
    first = min(starts) // width * width
    last = max(int(r["end"]) for r in source) // width * width
    bins: dict[int, dict[str, list[float]]] = {}
    for row in source:
        start = int(row["start"])
        key = start // width * width
        bucket = bins.setdefault(key, {name: [] for name in VALUE_FIELDS})
        for name, field in VALUE_FIELDS.items():
            value = finite(row.get(field, ""))
            pair = "pair2" if name == "pair2_contrast" else "pair1"
            if value is not None and row.get(f"{pair}_available") == "1":
                bucket[name].append(value)
    output = []
    for start in range(first, last + 1, width):
        bucket = bins.get(start, {})
        record: dict[str, float | int] = {"start": start, "end": start + width,
                                         "mid": start + width / 2,
                                         "n_windows_pair1": len(bucket.get("pair1_contrast", [])),
                                         "n_windows_pair2": len(bucket.get("pair2_contrast", [])),
                                         "n_windows_B": len(bucket.get("B", []))}
        for name in VALUE_FIELDS:
            values = bucket.get(name, [])
            record[name] = float(np.mean(values)) if values else float("nan")
        contrast_values = [float(record[name]) for name in ("pair1_contrast", "pair2_contrast") if math.isfinite(float(record[name]))]
        record["lineage_contrast"] = float(np.mean(contrast_values)) if contrast_values else float("nan")
        output.append(record)
    return output


def write_bins(path: Path, rows: list[dict[str, float | int]], width: int) -> None:
    fields = ["bin_width_bp", "start", "end", "mid", "n_windows_pair1", "pair1_contrast", "n_windows_pair2", "pair2_contrast", "lineage_contrast", "n_windows_B", "B"]
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({"bin_width_bp": width, **{k: ("NA" if isinstance(v, float) and not math.isfinite(v) else v) for k, v in row.items()}})


def plot(rows: list[dict[str, float | int]], output: Path) -> None:
    plt.rcParams.update({"font.family": "Arial", "font.size": 8, "axes.labelsize": 9,
                         "xtick.labelsize": 8, "ytick.labelsize": 8, "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 4.5), sharex=True,
                             gridspec_kw={"height_ratios": [1, 1]}, constrained_layout=True)
    x = np.array([float(r["mid"]) / 1e6 for r in rows])
    vals = [("lineage_contrast", "Mean mat-A/a contrast (subs/site)", "#762a83"),
            ("B", "B(w) (subs/site)", "#2166ac")]
    for ax, (field, ylabel, color), panel in zip(axes, vals, ("a", "b")):
        ax.axvspan(INTERVAL[0]/1e6, INTERVAL[1]/1e6, color="#fdae61", alpha=.18, lw=0, zorder=0)
        if field == "lineage_contrast":
            for pair, label in (("pair1_contrast", "CJ03/CJ04"), ("pair2_contrast", "CJ15/CJ16")):
                y_pair = np.array([float(r[pair]) for r in rows])
                ax.plot(x, y_pair, color=color, lw=.65, alpha=.34, zorder=2, label=label)
            y = np.array([float(r["lineage_contrast"]) for r in rows])
            ax.plot(x, y, color=color, lw=1.25, solid_capstyle="round", zorder=3, label="pair mean")
            ax.axhline(0, color="#555555", lw=.65, zorder=1)
        else:
            y = np.array([float(r[field]) for r in rows])
            ax.plot(x, y, color=color, lw=1.0, solid_capstyle="round", zorder=2)
        ax.set_ylabel(ylabel)
        ax.text(-.09, .98, panel, transform=ax.transAxes, ha="left", va="top", fontweight="bold", fontsize=10)
        ax.grid(axis="y", color="#bdbdbd", alpha=.35, lw=.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[-1].set_xlabel("LGI genomic position (Mb)")
    axes[-1].set_xlim(float(x[0]-.05), float(x[-1]+.05))
    fig.text(.5, 1.015, "Published introgressed / recombination-suppressed interval: 1.5–8.0 Mb",
             ha="center", va="bottom", fontsize=8, color="#555555")
    for suffix, kwargs in (("png", {"dpi": 300}), ("pdf", {}), ("svg", {})):
        fig.savefig(output.with_suffix(f".{suffix}"), bbox_inches="tight", **kwargs)
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    source = read_rows(RESULTS / "branch_length_inventory.tsv")
    if not source:
        raise ValueError("No frozen per-window branch-length estimates found")
    for width in (50_000, 100_000, 200_000):
        rows = aggregate_bins(source, width)
        write_bins(RESULTS / f"lgi_{width//1000}kb_branch_bins.tsv", rows, width)
        if width == 100_000:
            plot(rows, FIGURES / "figure2_spatial_branch_lengths")
    caption = """**Figure 2. Binned local-tree distances across Neurospora linkage group I.** (a) Faint lines show the two independently documented matched pairs, CJ03 (mat A)/CJ04 (mat a) and CJ15 (mat A)/CJ16 (mat a); the prominent line is their equal-weight mean within each 100-kb bin. For each pair the plotted contrast is d(mat A, *N. hispaniola*) − d(mat a, *N. hispaniola*), where d is patristic distance in substitutions per site. This is not the manuscript's arrangement-informed A(w): the Neurospora design has mating-type pairs but no independent same- versus opposite-arrangement pair classes. (b) B(w) is the mean of the six unordered pairwise patristic distances among CJ03, CJ04, *N. crassa*, and *N. hispaniola*. Values are arithmetic means of available 50-variable-site window estimates in non-overlapping 100-kb bins; empty bins remain gaps. The shaded interval is the independently reported 1.5–8.0 Mb introgressed/recombination-suppressed region. All branch lengths are local BIONJ substitution distances per site, not CU or calibrated time. The two matched pairs are each documented mat-A/mat-a homokaryons."""
    (RESULTS / "figure2_caption.md").write_text(caption + "\n")

if __name__ == "__main__":
    main()
