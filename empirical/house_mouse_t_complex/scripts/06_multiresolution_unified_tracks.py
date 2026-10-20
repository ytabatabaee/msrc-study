#!/usr/bin/env python3
"""Generate multi-resolution unified house-mouse topology/distance tracks.

This is a visualization sensitivity analysis.  It does not rerun boundary
detection.  Quartet support is pooled from the frozen ALL_TIPS counts and
distance features are arithmetic means of valid local-tree observations.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from collections import defaultdict
from pathlib import Path

import numpy as np

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc-paper-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical/house_mouse_t_complex"
RES = EMP / "results"
OUT = RES / "stage5_unified_tracks"
FIG = EMP / "figures"
Q = RES / "stage2_fixed_quartet_scan.tsv"
DIST = RES / "stage4_distance/stage4_window_distance.tsv"
META = ROOT / "data/house_mouse_t_complex/processed/house_mouse_t_complex_ml_5kb_metadata.tsv"
INT = RES / "t_complex_inversion_intervals.tsv"
SEED = 20261020
RESOLUTIONS = (25_000, 50_000, 100_000, 250_000)
COLORS = {"q_species": "#4c78a8", "q_t_alt": "#f58518", "q_other": "#54a24b"}
INT_COLORS = ["#6a3d9a", "#1b9e77", "#d95f02", "#7570b3"]


def read(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path: Path, records, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        out = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        out.writeheader()
        out.writerows({key: row.get(key, "") for key in fields} for row in records)


def number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interval_list():
    return [(int(row["start_bp"]), int(row["end_bp"]), row["inversion_name"]) for row in read(INT)]


def observed_bp(rows):
    """Union of observed half-open window intervals within a bin."""
    spans = sorted((int(row["start_bp"]), int(row["end_bp"]) + 1) for row in rows)
    total = 0
    current = None
    for start, end in spans:
        if current is None:
            current = [start, end]
        elif start <= current[1]:
            current[1] = max(current[1], end)
        else:
            total += current[1] - current[0]
            current = [start, end]
    if current is not None:
        total += current[1] - current[0]
    return total


def aggregate(rows, size):
    groups = defaultdict(list)
    anchor = 5_000_000
    for row in rows:
        start = int(row["start_bp"])
        groups[anchor + ((start - anchor) // size) * size].append(row)
    output = []
    for start, group in sorted(groups.items()):
        end = start + size - 1
        total_counts = {key: 0.0 for key in ("q_species", "q_t_alt", "q_other", "q_unresolved")}
        for row in group:
            total = number(row["n_resolved_quartets"]) + number(row["n_unresolved_quartets"])
            for key in total_counts:
                total_counts[key] += number(row[key]) * total
        denominator = sum(total_counts.values())
        distance_keys = ("A_distance", "A_relative", "mean_pair_distance_no_spretus")
        values = {key: np.array([number(row[key]) for row in group], dtype=float) for key in distance_keys}
        record = {
            "resolution_bp": size,
            "bin_start_bp": start,
            "bin_end_bp": end,
            "midpoint_bp": (start + end) / 2,
            "n_observed_windows": len(group),
            "expected_5kb_windows": size // 5_000,
            "observed_bp": observed_bp(group),
            "bin_bp": size,
        }
        record["coverage_fraction"] = record["observed_bp"] / size
        record["observed_window_fraction"] = len(group) / record["expected_5kb_windows"]
        for key, value in total_counts.items():
            record[key] = value / denominator if denominator else math.nan
        record["D"] = record["q_t_alt"] - record["q_species"]
        for key, value in values.items():
            valid = value[np.isfinite(value)]
            record[key] = float(np.mean(valid)) if len(valid) else math.nan
        output.append(record)
    return output


def segments(x):
    if len(x) < 2:
        return [np.arange(len(x))]
    gaps = np.diff(x)
    cutoff = max(float(np.median(gaps)) * 2.5, 0.01)
    return [part for part in np.split(np.arange(len(x)), np.flatnonzero(gaps > cutoff) + 1) if len(part)]


def shade(ax):
    for index, (start, end, _) in enumerate(interval_list()):
        color = INT_COLORS[index % len(INT_COLORS)]
        ax.axvspan(start / 1e6, end / 1e6, color=color, alpha=0.12)
        ax.axvline(start / 1e6, color=color, lw=0.55, ls="--")
        ax.axvline(end / 1e6, color=color, lw=0.55, ls="--")


def line(ax, x, y, color, label):
    for part in segments(x):
        ax.plot(x[part], y[part], color=color, lw=1.0, marker=".", ms=1.5,
                label=label if part[0] == 0 else None)


def save(fig, stem):
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def figure(rows, size, stem):
    x = np.array([number(row["midpoint_bp"]) / 1e6 for row in rows])
    fig, ax = plt.subplots(4, 1, figsize=(13, 10), sharex=True,
                           gridspec_kw={"height_ratios": [1, 1, 1, 1.35]},
                           constrained_layout=True)
    for axis in ax:
        shade(axis)
        axis.grid(axis="y", alpha=0.18)
        axis.set_xlim(x.min(), x.max())
    line(ax[0], x, np.array([number(row["D"]) for row in rows]), "#333333", "D = q_T - q_S")
    ax[0].axhline(0, color="#777", lw=0.6)
    ax[0].set_ylabel("D(w)")
    ax[0].set_title("(a) ALL_TIPS topology contrast", loc="left", weight="bold")
    line(ax[1], x, np.array([number(row["A_distance"]) for row in rows]), "#c44e52", "A_distance")
    ax[1].axhline(0, color="#777", lw=0.6)
    ax[1].set_ylabel("substitution units")
    ax[1].set_title("(b) t - standard patristic distance", loc="left", weight="bold")
    line(ax[2], x, np.array([number(row["A_relative"]) for row in rows]), "#55a868", "A_relative")
    ax[2].axhline(0, color="#777", lw=0.6)
    ax[2].set_ylabel("relative contrast")
    ax[2].set_title("(c) normalized distance contrast", loc="left", weight="bold")
    base = np.zeros(len(rows))
    for key, color, label in (("q_species", COLORS["q_species"], "q_S species history"),
                              ("q_t_alt", COLORS["q_t_alt"], "q_T t-associated alternative"),
                              ("q_other", COLORS["q_other"], "q_3 third")):
        value = np.array([number(row[key]) for row in rows])
        for part in segments(x):
            ax[3].fill_between(x[part], base[part], base[part] + value[part],
                               color=color, alpha=0.82, label=label if part[0] == 0 else None)
        base += value
    ax[3].set_ylim(0, 1)
    ax[3].set_ylabel("support")
    ax[3].set_title("(d) pooled ALL_TIPS three-topology support", loc="left", weight="bold")
    ax[3].set_xlabel("chromosome 17 position (Mb)")
    ax[3].legend(loc="upper right", ncol=3, fontsize=8, frameon=False)
    fig.suptitle(f"House-mouse t-complex: {size // 1000}-kb aggregation", fontsize=15, weight="bold")
    save(fig, stem)


def comparison_fig(by_resolution):
    fig, axes = plt.subplots(5, 3, figsize=(14, 13), sharex=True, constrained_layout=True)
    specs = (("D", "#333333", "D = q_T - q_S"),
             ("A_distance", "#c44e52", "A_distance"),
             ("A_relative", "#55a868", "A_relative"))
    for row_index, size in enumerate((5_000, *RESOLUTIONS)):
        records = by_resolution[size]
        x = np.array([number(row["midpoint_bp"]) / 1e6 for row in records])
        for col, (key, color, title) in enumerate(specs):
            axis = axes[row_index, col]
            shade(axis)
            axis.grid(axis="y", alpha=0.16)
            line(axis, x, np.array([number(row[key]) for row in records]), color, title)
            axis.axhline(0, color="#777", lw=0.5)
            if row_index == 0:
                axis.set_title(title, weight="bold")
            if col == 0:
                axis.set_ylabel("5 kb" if size == 5_000 else f"{size // 1000} kb")
            if row_index == 4:
                axis.set_xlabel("chromosome 17 position (Mb)")
    fig.suptitle("House-mouse t-complex: multi-resolution topology and distance comparison",
                 fontsize=15, weight="bold")
    save(fig, "house_mouse_t_complex_multiresolution_comparison")


def main():
    q = {row["locus_id"]: row for row in read(Q) if row["treatment"] == "ALL_TIPS"}
    distance = {row["locus_id"]: row for row in read(DIST)}
    metadata = {row["locus_id"]: row for row in read(META)}
    if set(q) != set(distance) or set(q) != set(metadata):
        raise RuntimeError("ALL_TIPS, distance, and metadata locus sets differ")
    rows = []
    for locus in sorted(q, key=lambda key: int(metadata[key]["start_bp"])):
        source = q[locus]
        dist = distance[locus]
        rows.append({
            "locus_id": locus, "start_bp": source["start_bp"], "end_bp": source["end_bp"],
            "midpoint_bp": source["midpoint_bp"], "n_resolved_quartets": source["n_resolved_quartets"],
            "n_unresolved_quartets": source["n_unresolved_quartets"],
            "q_species": source["q_species"], "q_t_alt": source["q_t_alt"],
            "q_other": source["q_other"], "q_unresolved": source["q_unresolved"],
            "D": number(source["q_t_alt"]) - number(source["q_species"]),
            "A_distance": dist["A_distance"], "A_relative": dist["A_relative"],
            "mean_pair_distance_no_spretus": dist["mean_pair_distance_no_spretus"],
        })
    by_resolution = {5_000: rows}
    summaries = []
    all_bins = []
    for size in RESOLUTIONS:
        binned = aggregate(rows, size)
        by_resolution[size] = binned
        fields = list(binned[0])
        output_name = (f"binned_{size // 1000}kb_tracks.tsv" if size != 250_000
                       else "multiresolution_binned_250kb_tracks.tsv")
        write(OUT / output_name, binned, fields)
        all_bins.extend(binned)
        values_by_key = {key: np.array([number(row[key]) for row in binned])
                         for key in ("D", "A_distance", "A_relative")}
        values = values_by_key["D"]
        summaries.append({
            "resolution_bp": size, "n_bins": len(binned),
            "n_observed_windows": sum(int(row["n_observed_windows"]) for row in binned),
            "mean_observed_window_fraction": float(np.mean([number(row["observed_window_fraction"]) for row in binned])),
            "min_observed_window_fraction": float(np.min([number(row["observed_window_fraction"]) for row in binned])),
            "n_incomplete_bins": sum(number(row["coverage_fraction"]) < 0.999999 for row in binned),
            "median_abs_step_D": float(np.median(np.abs(np.diff(values)))) if len(values) > 1 else math.nan,
            "median_abs_step_A_distance": float(np.median(np.abs(np.diff(values_by_key["A_distance"])))) if len(values) > 1 else math.nan,
            "median_abs_step_A_relative": float(np.median(np.abs(np.diff(values_by_key["A_relative"])))) if len(values) > 1 else math.nan,
            "D_range": float(np.nanmax(values) - np.nanmin(values)),
            "A_distance_range": float(np.nanmax(values_by_key["A_distance"]) - np.nanmin(values_by_key["A_distance"])),
            "A_relative_range": float(np.nanmax(values_by_key["A_relative"]) - np.nanmin(values_by_key["A_relative"])),
            "n_unique_D": len(np.unique(np.round(values, 12))),
            "description": "Pooled quartet counts; means of valid distance observations",
        })
        if size != 250_000:
            figure(binned, size, f"house_mouse_t_complex_topology_distance_{size // 1000}kb")
    comparison_fig(by_resolution)
    write(OUT / "multiresolution_summary.tsv", summaries, list(summaries[0]))
    (OUT / "multiresolution_methods.md").write_text("""# Multi-resolution methods

The native 5-kb rows are the frozen ALL_TIPS quartet and Stage 4 distance observations. For each 25-, 50-, 100-, and 250-kb bin, all source windows whose starts fall inside the common bin boundaries anchored at 5 Mb were included. Quartet support was pooled by reconstructing category counts from each window's support fraction times its total quartet count (`n_resolved_quartets + n_unresolved_quartets`), then dividing pooled counts by the pooled denominator. This preserves the source unresolved-quartet denominator. Distance features were arithmetic means of valid local-tree observations.

Each bin records observed window count, expected 5-kb window count, observed base-pair union, and coverage fractions. Missing sequence is retained as missing; no values are interpolated. Inversion boundaries are plotted at their frozen coordinates and are not snapped to bin edges. This is a figure-readability sensitivity analysis only; boundary detection was not rerun.
""")
    (OUT / "multiresolution_validation.md").write_text(f"""# Multi-resolution validation

- Native windows: {len(rows):,}.
- Resolutions generated: 25, 50, 100, and 250 kb.
- Common bin anchor: 5,000,000 bp.
- Source locus sets match across ALL_TIPS, Stage 4 distance, and metadata: {set(q) == set(distance) == set(metadata)}.
- Pooled quartet counts use the source total denominator including unresolved quartets.
- Distance aggregation uses means of valid observations.
- Unobserved genomic gaps remain explicit through bin-level coverage fields and line breaks in plots.

The original Stage 5 native 5-kb and 250-kb outputs remain preserved. The new 50-kb figure is the initial main-paper candidate; the summary table supports comparison with 25 and 100 kb based on observed-bin coverage and track roughness, without using inversion coordinates for selection.
""")
    (OUT / "multiresolution_report.md").write_text(f"""# Multi-resolution report

The native 5-kb track is retained for local detail. The 25-kb track preserves
more short-range variation but has the largest median step size in D(w)
({summaries[0]['median_abs_step_D']:.4f}) and is visibly noisy. The 100-kb
track is smoother ({summaries[2]['median_abs_step_D']:.4f}) but has fewer bins
and lower mean observed-window coverage ({summaries[2]['mean_observed_window_fraction']:.3f}).
The 250-kb track is the smoothest but attenuates the local spatial patterns.
The 50-kb track is therefore the initial main-paper candidate: its D(w) step
size is {summaries[1]['median_abs_step_D']:.4f}, with mean observed-window
coverage {summaries[1]['mean_observed_window_fraction']:.3f}. This choice is a
figure-readability judgment and does not use inversion coordinates or rerun
boundary detection.

All resolutions use pooled ALL_TIPS quartet counts, means of valid distance
observations, common bin boundaries, and explicit observed coverage. The
original native 5-kb and Stage 5 250-kb outputs are preserved separately.
""")
    manifest = {"stage": "5_multiresolution", "seed": SEED,
                "command": "python empirical/house_mouse_t_complex/scripts/06_multiresolution_unified_tracks.py",
                "resolutions_bp": list(RESOLUTIONS), "bin_anchor_bp": 5_000_000,
                "inputs": {str(path.relative_to(ROOT)): sha(path) for path in (Q, DIST, META, INT)},
                "preserved_outputs": True, "boundary_detection_rerun": False}
    (OUT / "multiresolution_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"n_windows": len(rows), "resolutions_bp": RESOLUTIONS,
                      "summary": str(OUT / "multiresolution_summary.tsv")}, indent=2))


if __name__ == "__main__":
    main()
