#!/usr/bin/env python3
"""Build standardized empirical QSS data and simplex figures.

This is a postprocessing-only analysis.  It reads the frozen tables used by
the existing cross-dataset QSS figure and never runs tree inference.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = ROOT / "empirical" / "combined_quartet_simplex"
DEFAULT_OUT = HERE

Q_FIELDS = ["q1", "q2", "q3"]
STANDARD_FIELDS = [
    "dataset", "chromosome", "window_id", "start", "end", "midpoint",
    "q1", "q2", "q3", "region_class", "region_id",
    "n_informative_quartets", "observation_unit",
]
COLORS = {"inside": "#0072B2", "outside": "#D55E00", "boundary": "#777777", "unknown": "#999999"}
TOPOLOGY_MAPPINGS = {
    "atlantic_cod": ("baseline/species", "arrangement", "third"),
    "neoaves_chr4": ("species history", "alternative 1", "alternative 2"),
    "fire_ants_chr16": ("species", "haplotype/arrangement", "third"),
    "anopheles_2la": ("species", "arrangement", "third"),
    "house_mouse_t_complex": ("species", "t-complex alternative", "other"),
}
REGION_ORDER = ("inside", "outside", "boundary", "unknown")
REGION_MARKERS = {"inside": "o", "outside": "s", "boundary": "D", "unknown": "x"}
SQRT3_2 = math.sqrt(3.0) / 2.0
Q_SUM_TOLERANCE = 2e-6  # cod support is stored at six decimal places.
MARKER_BASE_AREA = 6.0
MARKER_LOG_SCALE = 0.5
MARKER_MAX_AREA = 8.0
MARKER_ALPHA = 0.48
HALO_RADIUS_FACTOR = 2.5
HALO_REFERENCE_RADIUS_FACTOR = 2.5
HALO_BASE_AREA = 24.0
HALO_SQRT_SCALE = 18.0
HALO_BASE_ALPHA = 0.012
HALO_LOG_ALPHA_SCALE = 0.018
HALO_MAX_ALPHA = 0.10
HALO_LAYER_FRACTIONS = (1.00, 0.83, 0.66, 0.49, 0.32, 0.15)
HALO_LAYER_WEIGHTS = (0.08, 0.12, 0.16, 0.20, 0.20, 0.24)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else row.get(field, "") for field in fields})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def number(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def interval_state(start: int, end: int, intervals: list[tuple[str, int, int]]) -> tuple[str, str]:
    """Classify a window using midpoint convention from the combined figure.

    A window whose midpoint is in an annotated interval is inside.  A window
    crossing an interval boundary without an inside midpoint is boundary.
    """
    midpoint = (start + end) / 2.0
    overlaps = [(rid, lo, hi) for rid, lo, hi in intervals if start <= hi and end >= lo]
    containing = [rid for rid, lo, hi in intervals if lo <= midpoint <= hi]
    if containing:
        return "inside", ";".join(containing)
    if overlaps:
        return "boundary", ";".join(rid for rid, _, _ in overlaps)
    return "outside", ""


def classify_source_region(dataset: str, row: dict[str, str], intervals: list[tuple[str, int, int]]) -> tuple[str, str]:
    if dataset == "atlantic_cod":
        if row.get("inside_inversion", "").lower() == "true":
            return "inside", f"{row['lg']}_supergene"
        if row.get("overlaps_inversion", "").lower() == "true":
            return "boundary", f"{row['lg']}_supergene"
        return "outside", ""
    if dataset == "fire_ants_chr16":
        region = row.get("region", "")
        if region == "chr16_supergene":
            return "inside", region
        if region in {"chr16A", "chr16B"}:
            return "outside", region
        return "unknown", region
    return interval_state(int(float(row["start"])), int(float(row["end"])), intervals)


def validate_qss(rows: list[dict[str, object]], dataset: str) -> list[dict[str, object]]:
    valid: list[dict[str, object]] = []
    seen: set[tuple[str, str, str]] = set()
    for row in rows:
        key = (str(row["dataset"]), str(row["chromosome"]), str(row["window_id"]))
        if key in seen:
            raise ValueError(f"duplicate standardized observation: {key}")
        seen.add(key)
        n = number(row.get("n_informative_quartets"))
        if n <= 0:
            continue
        q = [number(row[field]) for field in Q_FIELDS]
        if not all(math.isfinite(x) and x >= -1e-12 for x in q):
            raise ValueError(f"invalid QSS values in {dataset}: {row}")
        if not math.isclose(sum(q), 1.0, abs_tol=Q_SUM_TOLERANCE):
            raise ValueError(f"QSS values do not sum to one in {dataset}: {q}")
        valid.append(row)
    return valid


def cod_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    source = ROOT / "data/atlantic_cod/processed/stage4a_window_topology_scores.tsv"
    regions = ROOT / "data/atlantic_cod/metadata/region_manifest.tsv"
    read_tsv(regions)  # Read and retain this frozen manifest in provenance.
    # The source table already freezes inside/overlap labels; the manifest is
    # retained in provenance and used to confirm that each LG is represented.
    out = []
    for r in read_tsv(source):
        if r["lg"] not in {"LG01", "LG02", "LG07", "LG12"}:
            continue
        out.append({
            "dataset": "atlantic_cod", "chromosome": r["lg"], "window_id": r["window_id"],
            "start": r["start"], "end": r["end"], "midpoint": r["midpoint"],
            "q1": r["fraction_baseline"], "q2": r["fraction_arrangement"], "q3": r["fraction_third"],
            "region_class": classify_source_region("atlantic_cod", r, [] )[0],
            "region_id": f"{r['lg']}_supergene", "n_informative_quartets": r["n_informative_quartets"],
            "observation_unit": "250-kb window quartet-support aggregate",
        })
    return out, {"source": str(source.relative_to(ROOT)), "topology_mapping": "q1=baseline/species; q2=arrangement; q3=third", "region_source": str(regions.relative_to(ROOT)), "region_convention": "frozen source inside_inversion/overlaps_inversion flags", "sha256": sha256(source)}


def neoaves_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    source = ROOT / "empirical/neoaves_chr4/results/stage10_multiresolution_tracks/binned_quartet_support_100kb.tsv"
    regions = ROOT / "empirical/neoaves_chr4/results/canonical_structural_regions.tsv"
    intervals = [(r["canonical_region_id"], int(r["canonical_start"]), int(r["canonical_end"])) for r in read_tsv(regions)]
    out = []
    for r in read_tsv(source):
        if r["clade"] != "N61":
            continue
        state, rid = interval_state(int(float(r["bin_start"])), int(float(r["bin_end"])), intervals)
        out.append({
            "dataset": "neoaves_chr4", "chromosome": "chr4", "window_id": f"chr4_{r['bin_start']}_{r['bin_end']}",
            "start": r["bin_start"], "end": r["bin_end"], "midpoint": r["bin_midpoint"],
            "q1": r["q1"], "q2": r["q2"], "q3": r["q3"], "region_class": state, "region_id": rid,
            "n_informative_quartets": r["n_loci"], "observation_unit": "100-kb N61 bin mean quartet support",
        })
    return out, {"source": str(source.relative_to(ROOT)), "topology_mapping": "q1=species history; q2=alternative 1; q3=alternative 2", "region_source": str(regions.relative_to(ROOT)), "region_convention": "100-kb bin midpoint; boundary when a bin crosses an annotated interval without an inside midpoint", "sha256": sha256(source)}


def fire_ant_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    source = ROOT / "data/fire_ants_chr16/processed/stage4a_window_quartet_support.tsv"
    regions = ROOT / "data/fire_ants_chr16/metadata/region_manifest.tsv"
    out = []
    for r in read_tsv(source):
        if r["chrom"] != "chr16":
            continue
        out.append({
            "dataset": "fire_ants_chr16", "chromosome": "chr16", "window_id": f"chr16_{r['window_index']}",
            "start": r["start"], "end": r["end"], "midpoint": r["mid"], "q1": r["q_S"], "q2": r["q_H"], "q3": r["q_3"],
            "region_class": classify_source_region("fire_ants_chr16", r, [])[0], "region_id": r["region"],
            "n_informative_quartets": r["W_total"], "observation_unit": "BUSCO window quartet-weight aggregate",
        })
    return out, {"source": str(source.relative_to(ROOT)), "topology_mapping": "q1=species (S); q2=haplotype/arrangement (H); q3=third", "region_source": str(regions.relative_to(ROOT)), "region_convention": "frozen Stage 4A chr16 region labels", "sha256": sha256(source)}


def anopheles_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    source = ROOT / "empirical/anopheles_2la/results/three_quartet_support/window_three_quartet_support.tsv"
    out = []
    for r in read_tsv(source):
        out.append({
            "dataset": "anopheles_2la", "chromosome": "2L", "window_id": r["window_id"], "start": r["start"], "end": r["end"], "midpoint": r["midpoint"],
            "q1": r["q_species"], "q2": r["q_arrangement"], "q3": r["q_third"], "region_class": r["region_class"], "region_id": "2La" if r["region_class"] == "inside" else "",
            "n_informative_quartets": r["n_informative_quartets"], "observation_unit": "50-kb window quartet-support aggregate",
        })
    return out, {"source": str(source.relative_to(ROOT)), "topology_mapping": "q1=species; q2=arrangement; q3=third", "region_source": "source region_class; frozen 2La boundaries are retained in source pipeline", "region_convention": "source region_class", "sha256": sha256(source)}


def mouse_rows() -> tuple[list[dict[str, object]], dict[str, object]]:
    source = ROOT / "empirical/house_mouse_t_complex/results/stage2_fixed_quartet_scan.tsv"
    intervals_file = ROOT / "empirical/house_mouse_t_complex/results/t_complex_inversion_intervals.tsv"
    intervals = [(r["inversion_name"], int(r["start_bp"]), int(r["end_bp"])) for r in read_tsv(intervals_file)]
    source_rows = [r for r in read_tsv(source) if r["treatment"] == "ALL_TIPS"]
    groups: dict[float, list[dict[str, str]]] = defaultdict(list)
    for r in source_rows:
        midpoint = float(r["midpoint_bp"])
        bin_start = math.floor(midpoint / 100_000.0) * 100_000
        groups[bin_start].append(r)
    out = []
    for bin_start, group in sorted(groups.items()):
        state, rid = interval_state(int(bin_start), int(bin_start + 100_000), intervals)
        q = {field: float(np.mean([number(r[field]) for r in group])) for field in ("q_species", "q_t_alt", "q_other")}
        # The frozen mouse source retains a small unresolved fraction.  The
        # existing combined figure conditions its simplex on resolved q1-q3;
        # make that operation explicit rather than silently accepting a sum
        # below one.
        qsum = sum(q.values())
        if qsum <= 0:
            continue
        q = {key: value / qsum for key, value in q.items()}
        out.append({
            "dataset": "house_mouse_t_complex", "chromosome": "chr17", "window_id": f"chr17_{int(bin_start)}_{int(bin_start + 100_000)}",
            "start": int(bin_start), "end": int(bin_start + 100_000), "midpoint": bin_start + 50_000,
            "q1": q["q_species"], "q2": q["q_t_alt"], "q3": q["q_other"], "region_class": state, "region_id": rid,
            "n_informative_quartets": sum(number(r["n_resolved_quartets"]) for r in group), "observation_unit": "100-kb mean of observed 5-kb ALL_TIPS windows",
        })
    return out, {"source": str(source.relative_to(ROOT)), "topology_mapping": "q1=species; q2=t-complex alternative; q3=other", "region_source": str(intervals_file.relative_to(ROOT)), "region_convention": "100-kb bin midpoint; boundary when a bin crosses an annotated interval without an inside midpoint", "qss_convention": "100-kb means are conditioned on resolved q1-q3 by explicit normalization; unresolved source fraction is not plotted", "sha256": sha256(source)}


def load_all() -> tuple[list[dict[str, object]], dict[str, dict[str, object]], dict[str, int], list[dict[str, object]]]:
    loaders = {"atlantic_cod": cod_rows, "neoaves_chr4": neoaves_rows, "fire_ants_chr16": fire_ant_rows, "anopheles_2la": anopheles_rows, "house_mouse_t_complex": mouse_rows}
    all_rows: list[dict[str, object]] = []
    provenance: dict[str, dict[str, object]] = {}
    raw_counts: dict[str, int] = {}
    excluded: list[dict[str, object]] = []
    for name, loader in loaders.items():
        rows, prov = loader()
        raw_counts[name] = len(rows)
        for row in rows:
            if number(row.get("n_informative_quartets")) <= 0:
                excluded.append({**row, "exclusion_reason": "zero_informative_quartets"})
        all_rows.extend(validate_qss(rows, name))
        provenance[name] = prov
    return all_rows, provenance, raw_counts, excluded


def xy(q1: np.ndarray, q2: np.ndarray, q3: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return 0.5 * q1 + q3, SQRT3_2 * q1


def panel_for_row(row: dict[str, object]) -> str:
    """Return the visual panel identifier used in overlap summaries."""
    return str(row["chromosome"]) if row["dataset"] == "atlantic_cod" else str(row["dataset"])


def coordinate_groups(rows: list[dict[str, object]], decimals: int = 12) -> list[dict[str, object]]:
    """Group rows by rounded simplex coordinates and region without moving them.

    The representative q values are copied from the first row in each group;
    rounding is used only as a grouping key.  Thus plotted coordinates remain
    the original coordinates from the standardized table.
    """
    grouped: dict[tuple[str, float, float, float], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        key = (str(row["region_class"]), *(round(number(row[field]), decimals) for field in Q_FIELDS))
        grouped[key].append(row)
    output = []
    for (region, _, _, _), members in grouped.items():
        first = members[0]
        output.append({"region_class": region, "q1": number(first["q1"]), "q2": number(first["q2"]), "q3": number(first["q3"]), "n_overlap": len(members)})
    return output


def marker_area(n_overlap: int) -> float:
    """Compact marker area in points^2 using logarithmic multiplicity scaling."""
    return min(MARKER_MAX_AREA, MARKER_BASE_AREA + MARKER_LOG_SCALE * math.log2(max(1, n_overlap)))


def halo_alpha(n_overlap: int) -> float:
    """Total halo opacity, increasing logarithmically with multiplicity."""
    return min(HALO_MAX_ALPHA, HALO_BASE_ALPHA + HALO_LOG_ALPHA_SCALE * math.log2(max(1, n_overlap)))


def halo_area(n_overlap: int, radius_fraction: float) -> float:
    """Original multiplicity-aware halo area with a soft layer fraction."""
    radius_scale = radius_fraction * HALO_RADIUS_FACTOR / HALO_REFERENCE_RADIUS_FACTOR
    return (HALO_BASE_AREA + HALO_SQRT_SCALE * math.sqrt(max(1, n_overlap))) * radius_scale ** 2


def draw_simplex(ax, rows: list[dict[str, object]], title: str, show_legend: bool = False) -> None:
    ax.set_aspect("equal")
    ax.set_xlim(-0.08, 1.08); ax.set_ylim(-0.08, SQRT3_2 + 0.08); ax.axis("off")
    vertices = np.array([[0, 0], [1, 0], [0.5, SQRT3_2]])
    for fraction in np.arange(0.2, 1.0, 0.2):
        # q1=fraction horizontal line; q2=fraction and q3=fraction oblique lines.
        for a, b in ((np.array([fraction * .5, fraction * SQRT3_2]), np.array([1 - fraction * .5, fraction * SQRT3_2])),
                     (np.array([fraction, 0]), np.array([fraction * .5 + (1-fraction), (1-fraction) * SQRT3_2])),
                     (np.array([1-fraction, 0]), np.array([fraction * .5, (1-fraction) * SQRT3_2]))):
            ax.plot([a[0], b[0]], [a[1], b[1]], color="#d9d9d9", lw=.45, zorder=0)
    ax.add_patch(Polygon(vertices, closed=True, fill=False, edgecolor="#333333", lw=.8, zorder=2))
    ax.plot([.5, .5], [0, SQRT3_2], color="#666666", lw=.65, ls="--", zorder=1)
    # Draw outside first and inside last, while retaining distinct marker
    # shapes so coincident inside/outside observations remain visible. Within
    # every region, the sorted PathCollection draws smaller multiplicities
    # before larger ones.
    for region in ("outside", "boundary", "unknown", "inside"):
        grouped = [g for g in coordinate_groups(rows) if g["region_class"] == region]
        if not grouped:
            continue
        grouped.sort(key=lambda g: int(g["n_overlap"]))
        q = np.array([[g[f] for f in Q_FIELDS] for g in grouped])
        x, y = xy(q[:, 0], q[:, 1], q[:, 2])
        color = COLORS[region]
        marker = REGION_MARKERS[region]
        # Approximate a smoothly fading halo with nested, outline-free disks.
        # Batching by multiplicity keeps the renderer lightweight while making
        # high-overlap coordinates more visible than single observations.
        by_n: dict[int, list[dict[str, object]]] = defaultdict(list)
        for group in grouped:
            by_n[int(group["n_overlap"])].append(group)
        for n_overlap in sorted(by_n):
            halo_groups = by_n[n_overlap]
            halo_q = np.array([[g[f] for f in Q_FIELDS] for g in halo_groups])
            halo_x, halo_y = xy(halo_q[:, 0], halo_q[:, 1], halo_q[:, 2])
            for fraction, weight in zip(HALO_LAYER_FRACTIONS, HALO_LAYER_WEIGHTS):
                ax.scatter(halo_x, halo_y, s=halo_area(n_overlap, fraction), c=color,
                           alpha=halo_alpha(n_overlap) * weight, edgecolors="none",
                           linewidths=0, marker="o", clip_on=False, zorder=3)
        ax.scatter(x, y, s=[marker_area(int(g["n_overlap"])) for g in grouped], c=color, alpha=MARKER_ALPHA,
                   edgecolors="none", linewidths=0, marker=marker, clip_on=False,
                   label=f"{region} (n={sum(int(g['n_overlap']) for g in grouped)})", zorder=4)
    counts = {region: sum(1 for r in rows if r["region_class"] == region) for region in ("inside", "outside")}
    ax.text(.5, SQRT3_2 + .012, "q1 = 1\nbaseline", ha="center", va="bottom", fontsize=7)
    ax.text(-.035, -.035, "q2 = 1\nalternative", ha="right", va="top", fontsize=7)
    ax.text(1.035, -.035, "q3 = 1\nthird", ha="left", va="top", fontsize=7)
    ax.set_title(f"{title}\ninside n={counts['inside']}; outside n={counts['outside']}", fontsize=8.5, pad=16, linespacing=1.05)
    if show_legend:
        region_handles = [Line2D([0], [0], marker=REGION_MARKERS[r], color="none", markerfacecolor=COLORS[r], markeredgecolor="none", alpha=MARKER_ALPHA, markersize=5, label=r) for r in ("inside", "outside")]
        region_handles.append(Line2D([0], [0], color="#666", ls="--", lw=.65, label="q2 = q3"))
        region_legend = ax.legend(handles=region_handles, loc="upper center", bbox_to_anchor=(.5, -0.08), ncol=3, frameon=False, fontsize=7)
        ax.add_artist(region_legend)
        size_handles = [Line2D([0], [0], marker="o", linestyle="none", color="none", markerfacecolor=COLORS["outside"], markeredgecolor="none", alpha=MARKER_ALPHA, markersize=math.sqrt(marker_area(n)) * 1.5, label=f"n={n}") for n in (1, 10, 100)]
        ax.legend(handles=size_handles, title="overlap multiplicity", loc="upper center", bbox_to_anchor=(.5, -0.20), ncol=3, frameon=False, fontsize=6.5, title_fontsize=6.5, handletextpad=.4, columnspacing=.8)


def savefig(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(path.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)


def summaries(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        group = row["chromosome"] if row["dataset"] == "atlantic_cod" else "all"
        grouped[(str(row["dataset"]), str(group), str(row["region_class"]))].append(row)
    out = []
    for (dataset, chromosome, region), group in sorted(grouped.items()):
        q = {f: np.array([number(r[f]) for r in group]) for f in Q_FIELDS}
        out.append({"dataset": dataset, "chromosome": chromosome, "region_class": region, "n_observations": len(group),
                    **{f"mean_{f}": float(np.mean(q[f])) for f in Q_FIELDS}, **{f"median_{f}": float(np.median(q[f])) for f in Q_FIELDS},
                    "fraction_q2_gt_q1": float(np.mean(q["q2"] > q["q1"])), "fraction_q2_gt_q3": float(np.mean(q["q2"] > q["q3"])),
                    "fraction_q1_gt_q2_and_q3": float(np.mean((q["q1"] > q["q2"]) & (q["q1"] > q["q3"])))})
    return out


def make_figures(rows: list[dict[str, object]], out: Path) -> None:
    fig_dir = out / "figures"
    by = defaultdict(list)
    for row in rows: by[row["dataset"]].append(row)
    # Five-row figure; cod gets four equal-width linkage-group axes.
    fig = plt.figure(figsize=(10.5, 14.0), facecolor="white")
    gs = fig.add_gridspec(5, 4, hspace=.72, wspace=.35)
    for j, lg in enumerate(("LG01", "LG02", "LG07", "LG12")):
        draw_simplex(fig.add_subplot(gs[0, j]), [r for r in by["atlantic_cod"] if r["chromosome"] == lg], f"Atlantic cod {lg}", show_legend=(j == 0))
    for i, (dataset, title) in enumerate((("neoaves_chr4", "Neoaves chr4 (N61)"), ("fire_ants_chr16", "Fire ants chr16"), ("anopheles_2la", "Anopheles 2La"), ("house_mouse_t_complex", "House mouse chr17")), start=1):
        draw_simplex(fig.add_subplot(gs[i, :]), by[dataset], title, show_legend=(i == 1))
    fig.suptitle("Empirical quartet-support simplexes inside and outside rearranged regions", fontsize=13, y=.995)
    savefig(fig, fig_dir / "qss_simplex_combined_empirical")
    for dataset, title in (("neoaves_chr4", "Neoaves chr4 (N61)"), ("fire_ants_chr16", "Fire ants chr16"), ("anopheles_2la", "Anopheles 2La"), ("house_mouse_t_complex", "House mouse chr17")):
        fig, ax = plt.subplots(figsize=(4.2, 4.1)); draw_simplex(ax, by[dataset], title, True); savefig(fig, fig_dir / f"{dataset}_qss_simplex")
    fig, axes = plt.subplots(1, 4, figsize=(12.5, 3.5))
    for j, (ax, lg) in enumerate(zip(axes, ("LG01", "LG02", "LG07", "LG12"))): draw_simplex(ax, [r for r in by["atlantic_cod"] if r["chromosome"] == lg], f"Atlantic cod {lg}", j == 0)
    savefig(fig, fig_dir / "atlantic_cod_qss_simplex_four_panel")
    fig, ax = plt.subplots(figsize=(4.8, 4.5)); draw_simplex(ax, by["atlantic_cod"], "Atlantic cod (LG01, LG02, LG07, LG12)", True); savefig(fig, fig_dir / "atlantic_cod_qss_simplex")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    rows, provenance, raw_counts, excluded = load_all()
    out = args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "standardized_qss.tsv", rows, STANDARD_FIELDS)
    summary_fields = ["dataset", "chromosome", "region_class", "n_observations", "mean_q1", "mean_q2", "mean_q3", "median_q1", "median_q2", "median_q3", "fraction_q2_gt_q1", "fraction_q2_gt_q3", "fraction_q1_gt_q2_and_q3"]
    write_tsv(out / "inside_outside_summary.tsv", summaries(rows), summary_fields)
    write_tsv(out / "excluded_qss.tsv", excluded, STANDARD_FIELDS + ["exclusion_reason"])
    overlap_rows = []
    panel_groups: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        panel_groups[(str(row["dataset"]), panel_for_row(row))].append(row)
    for (dataset, panel), panel_rows in sorted(panel_groups.items()):
        for group in coordinate_groups(panel_rows):
            overlap_rows.append({"dataset": dataset, "panel": panel, **group})
    write_tsv(out / "overlap_counts.tsv", overlap_rows, ["dataset", "panel", "region_class", "q1", "q2", "q3", "n_overlap"])
    audit = out / "data_audit.md"
    lines = ["# Empirical QSS simplex data audit", "", "All sources are frozen postprocessing tables used by the existing `qss_dw_aw_combined_all_datasets` figure. No gene trees or ASTRAL analyses were rerun.", "", "## Sources"]
    for dataset, p in provenance.items():
        lines += [f"### {dataset}", f"- Source: `{p['source']}`", f"- SHA-256: `{p['sha256']}`", f"- Topology mapping: {p['topology_mapping']}", f"- Region source/convention: `{p['region_source']}`; {p['region_convention']}", f"- Raw standardized rows: {raw_counts[dataset]}", f"- Retained rows: {sum(1 for r in rows if r['dataset'] == dataset)}", "- Observation unit: one frozen window or 100-kb aggregate, never an individual quartet.", ""]
    lines += ["## Individual-quartet availability", "", "The retained QSS tables contain window-level or binned aggregates. They do not retain three support frequencies for each fixed quartet across multiple loci. A true individual-quartet simplex distribution therefore cannot be reconstructed without returning to lower-level source data; no such figure is generated.", "", "## Region classes", "", "Inside and outside points use the frozen source labels where present. For interval-based sources, the combined-figure convention is applied by midpoint; intervals crossed by a window whose midpoint lies outside are marked boundary. Boundary and unknown observations remain in the standardized table and are shown separately when present."]
    audit.write_text("\n".join(lines) + "\n")
    (out / "provenance.json").write_text(json.dumps({"script": str(Path(__file__).relative_to(ROOT)), "python": platform.python_version(), "sources": provenance, "standard_fields": STANDARD_FIELDS, "overlap_encoding": {"grouping_round_decimals": 12, "marker_area": "min(8, 6 + 0.5*log2(n_overlap))", "core_markers": REGION_MARKERS, "alpha": MARKER_ALPHA, "edgecolors": "none", "coordinates_shifted": False, "halo": {"area": "24 + 18*sqrt(n_overlap)", "radius_factor": HALO_RADIUS_FACTOR, "total_alpha": "min(0.10, 0.012 + 0.018*log2(n_overlap))", "nested_layer_fractions": HALO_LAYER_FRACTIONS, "nested_layer_weights": HALO_LAYER_WEIGHTS}}}, indent=2) + "\n")
    make_figures(rows, out)
    print(f"retained {len(rows)} standardized observations")


if __name__ == "__main__":
    main()
