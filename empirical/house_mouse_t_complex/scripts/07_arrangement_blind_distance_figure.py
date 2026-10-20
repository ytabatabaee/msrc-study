#!/usr/bin/env python3
"""Add an arrangement-blind patristic-distance panel to house-mouse tracks."""
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
WINDOW = OUT / "window_unified_tracks.tsv"
Q = RES / "stage2_fixed_quartet_scan.tsv"
DIST = RES / "stage4_distance/stage4_window_distance.tsv"
META = ROOT / "data/house_mouse_t_complex/processed/house_mouse_t_complex_ml_5kb_metadata.tsv"
INT = RES / "t_complex_inversion_intervals.tsv"
SEED = 20261020
COLORS = {"D": "#333333", "A_distance": "#c44e52", "A_relative": "#55a868", "B": "#8172b2"}
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


def value(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def sha(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def intervals():
    return [(int(row["start_bp"]), int(row["end_bp"])) for row in read(INT)]


def segments(x):
    if len(x) < 2:
        return [np.arange(len(x))]
    gap = max(float(np.median(np.diff(x))) * 2.5, 0.01)
    return [part for part in np.split(np.arange(len(x)), np.flatnonzero(np.diff(x) > gap) + 1) if len(part)]


def shade(axis):
    for index, (start, end) in enumerate(intervals()):
        color = INT_COLORS[index % len(INT_COLORS)]
        axis.axvspan(start / 1e6, end / 1e6, color=color, alpha=0.12)
        axis.axvline(start / 1e6, color=color, ls="--", lw=0.55)
        axis.axvline(end / 1e6, color=color, ls="--", lw=0.55)


def plot_track(axis, x, y, color, label):
    for part in segments(x):
        axis.plot(x[part], y[part], color=color, lw=0.85, marker=".", ms=1.4,
                  label=label if part[0] == 0 else None)


def save(fig, stem):
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def figure(rows, resolution, stem):
    x = np.array([value(row["midpoint_bp"]) / 1e6 for row in rows])
    fig, axes = plt.subplots(4, 1, figsize=(13, 9.2), sharex=True, constrained_layout=True)
    for axis in axes:
        shade(axis)
        axis.grid(axis="y", alpha=0.18)
        axis.set_xlim(x.min(), x.max())
    specs = (("D", "D(w) = q_T - q_S", "ALL_TIPS topology contrast"),
             ("A_distance", "substitution units", "arrangement-informed A_distance"),
             ("A_relative", "relative contrast", "arrangement-informed A_relative"),
             ("B", "substitution units", "arrangement-blind mean pairwise distance"))
    for axis, (key, ylabel, title) in zip(axes, specs):
        plot_track(axis, x, np.array([value(row[key]) for row in rows]), COLORS[key], title)
        if key != "B":
            axis.axhline(0, color="#777", lw=0.6)
        axis.set_ylabel(ylabel)
        axis.set_title(title, loc="left", weight="bold")
    axes[-1].set_xlabel("chromosome 17 position (Mb)")
    axes[0].legend(frameon=False, loc="upper right", fontsize=8)
    fig.suptitle(f"House-mouse t-complex: topology and arrangement-blind distance ({resolution})",
                 fontsize=15, weight="bold")
    save(fig, stem)


def observed_bp(group):
    spans = sorted((int(row["start_bp"]), int(row["end_bp"]) + 1) for row in group)
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
    if current:
        total += current[1] - current[0]
    return total


def aggregate(rows, counts, size=250_000):
    grouped = defaultdict(list)
    for row in rows:
        start = int(row["start_bp"])
        grouped[5_000_000 + ((start - 5_000_000) // size) * size].append(row)
    output = []
    for start, group in sorted(grouped.items()):
        count_keys = ("q_species", "q_t_alt", "q_other", "q_unresolved")
        pooled = {key: 0.0 for key in count_keys}
        for row in group:
            total = counts[row["locus_id"]]
            for key in count_keys:
                pooled[key] += value(row[key]) * total
        denominator = sum(pooled.values())
        result = {"bin_start_bp": start, "bin_end_bp": start + size - 1,
                  "midpoint_bp": start + size / 2, "n_observed_windows": len(group),
                  "observed_bp": observed_bp(group), "bin_bp": size}
        result["coverage_fraction"] = result["observed_bp"] / size
        for key in count_keys:
            result[key] = pooled[key] / denominator
        result["D"] = result["q_t_alt"] - result["q_species"]
        for key in ("A_distance", "A_relative", "B"):
            vals = np.array([value(row[key]) for row in group])
            vals = vals[np.isfinite(vals)]
            result[key] = float(np.mean(vals)) if len(vals) else math.nan
        output.append(result)
    return output


def main():
    rows = read(WINDOW)
    source_distance = {row["locus_id"]: row for row in read(DIST)}
    source_q = {row["locus_id"]: row for row in read(Q) if row["treatment"] == "ALL_TIPS"}
    metadata = {row["locus_id"]: row for row in read(META)}
    if set(source_distance) != set(source_q) or {row["locus_id"] for row in rows} != set(source_q):
        raise RuntimeError("window, Stage 2, and Stage 4 locus sets differ")
    rows = sorted(rows, key=lambda row: int(row["start_bp"]))
    for row in rows:
        source = source_distance[row["locus_id"]]
        row["B"] = row["mean_pair_distance_no_spretus"]
        if row["start_bp"] != source["start_bp"] or row["end_bp"] != source["end_bp"]:
            raise RuntimeError(f"coordinate mismatch at {row['locus_id']}")
        if abs(value(row["A_distance"]) - value(source["A_distance"])) > 1e-15:
            raise RuntimeError(f"A_distance mismatch at {row['locus_id']}")
        if abs(value(row["A_relative"]) - value(source["A_relative"])) > 1e-15:
            raise RuntimeError(f"A_relative mismatch at {row['locus_id']}")
        if abs(value(row["B"]) - value(source["mean_pair_distance_no_spretus"])) > 1e-15:
            raise RuntimeError(f"arrangement-blind distance mismatch at {row['locus_id']}")
    counts = {locus: value(row["n_resolved_quartets"]) + value(row["n_unresolved_quartets"])
              for locus, row in source_q.items()}
    binned = aggregate(rows, counts, 250_000)
    binned_100 = aggregate(rows, counts, 100_000)
    # Separate filenames preserve the previously frozen Stage 5 aggregate table.
    write(OUT / "binned_250kb_tracks_arrangement_blind.tsv", binned, list(binned[0]))
    write(OUT / "binned_100kb_tracks_arrangement_blind.tsv", binned_100, list(binned_100[0]))
    figure(rows, "native 5-kb", "house_mouse_t_complex_topology_distance_blind_5kb")
    figure(binned_100, "100-kb aggregation", "house_mouse_t_complex_topology_distance_blind_100kb")
    figure(binned, "250-kb aggregation", "house_mouse_t_complex_topology_distance_blind_250kb")
    summary = []
    for key, label in (("D", "ALL_TIPS topology contrast"), ("A_distance", "A_distance"),
                       ("A_relative", "A_relative"), ("B", "mean_pair_distance_no_spretus")):
        values = np.array([value(row[key]) for row in rows])
        summary.append({"track": label, "n_windows": len(values), "mean": float(np.mean(values)),
                        "median": float(np.median(values)), "minimum": float(np.min(values)),
                        "maximum": float(np.max(values)), "sd": float(np.std(values)),
                        "units": "substitution units" if key in ("A_distance", "B") else "dimensionless"})
    write(OUT / "track_summary.tsv", summary, list(summary[0]))
    (OUT / "arrangement_blind_validation_report.md").write_text(f"""# Arrangement-blind track validation

- Unified window rows: {len(rows):,}; source Stage 2 ALL_TIPS rows: {len(source_q):,}; source Stage 4 rows: {len(source_distance):,}.
- Coordinate sets and ordering match: {all(rows[i]['locus_id'] == sorted(source_q, key=lambda k: int(metadata[k]['start_bp']))[i] for i in range(len(rows)))}.
- `A_distance`, `A_relative`, and `B=mean_pair_distance_no_spretus` match the frozen Stage 4 values to the 1e-15 check threshold.
- The arrangement-blind track uses all non-`M. spretus` pairwise distances and does not use t-versus-standard labels.
- Native 5-kb, separate 100-kb, and separate 250-kb figures were generated. The prior Stage 5 figures and `binned_250kb_tracks.tsv` were preserved.
- 100-kb and 250-kb aggregation use pooled ALL_TIPS quartet counts and arithmetic means of valid distance observations; bin coverage is recorded in the new arrangement-blind tables.
""")
    (OUT / "arrangement_blind_figure_caption.md").write_text("""**Figure. House-mouse t-complex topology and arrangement-blind divergence-distance tracks.** (a) The ALL_TIPS quartet contrast D(w)=q_T-q_S. (b) A_distance is the t-versus-standard patristic-distance contrast. (c) A_relative is its normalized form. (d) B(w) is the mean pairwise patristic distance among non-M. spretus samples and is arrangement-blind because it does not use t-versus-standard labels. All distance values are substitution-unit summaries from local ML trees, not calibrated times or coalescence times. Shading and dashed lines show the four frozen inversion intervals. Gaps in the observed windows are retained.
""")
    (OUT / "arrangement_blind_methods.md").write_text("""# Methods

The primary topology track is the frozen ALL_TIPS D(w)=q_T(w)-q_S(w). A_distance and A_relative are copied from the frozen Stage 4 local ML-tree table. The new arrangement-blind statistic B(w)=mean_pair_distance_no_spretus is the mean pairwise patristic distance among all non-M. spretus samples in each local tree; it does not use t-versus-standard labels. All patristic quantities remain in the original ML substitution units.

The native figure uses the observed 5-kb windows directly. The separate 250-kb table and figure use common bins anchored at 5 Mb, pooled quartet counts reconstructed from the frozen support fractions and total quartet counts, and arithmetic means of valid distance observations. Unobserved sequence is not interpolated.
""")
    manifest = {"stage": "5_arrangement_blind_distance_figure", "seed": SEED,
                "command": "python empirical/house_mouse_t_complex/scripts/07_arrangement_blind_distance_figure.py",
                "inputs": {str(path.relative_to(ROOT)): sha(path) for path in (WINDOW, Q, DIST, META, INT)},
                "arrangement_blind_feature": "mean_pair_distance_no_spretus",
                "units": "local ML patristic substitution units",
                "arrangement_blind_resolutions_bp": [100000, 250000],
                "existing_outputs_preserved": True}
    (OUT / "arrangement_blind_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"n_windows": len(rows), "n_bins_250kb": len(binned)}, indent=2))


if __name__ == "__main__":
    main()
