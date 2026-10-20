#!/usr/bin/env python3
"""Stage 9 fire-ant chromosome-16 patristic distance track.

This analysis uses the published RAxML local trees as trees with branch
lengths in substitution units.  It calculates a group-balanced contrast
between opposite and same supergene haplotypes and compares that track with
the frozen Stage-4A quartet D(w) track.  It deliberately does not interpret
patristic distances as coalescence times or assume a molecular clock.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
import platform
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/msrc-paper-mplconfig")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "fire_ants_chr16"
EMP = ROOT / "empirical" / "fire_ants_chr16"
PROCESSED = DATA / "processed"
RESULTS = EMP / "results" / "stage9_distance"
FIGURES = EMP / "figures"

TREE_TABLE = PROCESSED / "stage1_window_trees.tsv"
WINDOW_INDEX = PROCESSED / "stage1_window_index.tsv"
MAPPING = PROCESSED / "stage6c_individual_to_group.tsv"
GROUP_MAP = PROCESSED / "stage6c_group_label_map.tsv"
SUPPORT = PROCESSED / "stage4a_window_quartet_support.tsv"
REGION_MANIFEST = DATA / "metadata" / "region_manifest.tsv"
SPATIAL_SCRIPT = EMP / "scripts" / "04b_spatial_null_tests.py"

FOCAL = ["inv_mac_SB", "inv_mac_Sb", "richteri_SB", "richteri_Sb"]
OPPOSITE = [
    ("inv_mac_SB", "inv_mac_Sb"),
    ("inv_mac_SB", "richteri_Sb"),
    ("inv_mac_Sb", "richteri_SB"),
    ("richteri_SB", "richteri_Sb"),
]
SAME = [("inv_mac_SB", "richteri_SB"), ("inv_mac_Sb", "richteri_Sb")]
ALL_PAIRS = OPPOSITE + SAME
GROUP_FIELDS = {
    ("inv_mac_SB", "inv_mac_Sb"): "pair_inv_mac_SB__inv_mac_Sb",
    ("inv_mac_SB", "richteri_Sb"): "pair_inv_mac_SB__richteri_Sb",
    ("inv_mac_Sb", "richteri_SB"): "pair_inv_mac_Sb__richteri_SB",
    ("richteri_SB", "richteri_Sb"): "pair_richteri_SB__richteri_Sb",
    ("inv_mac_SB", "richteri_SB"): "pair_inv_mac_SB__richteri_SB",
    ("inv_mac_Sb", "richteri_Sb"): "pair_inv_mac_Sb__richteri_Sb",
}
SEED = 20261019
TOL = 1e-12


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({field: "" if row.get(field) is None else row.get(field) for field in fields})


def finite_mean(values: Iterable[float]) -> float | None:
    vals = [v for v in values if math.isfinite(v)]
    return statistics.fmean(vals) if vals else None


def parse_tree(text: str):
    return Phylo.read(io.StringIO(text.strip()), "newick")


def ancestor_depths(tree) -> dict[str, dict[object, float]]:
    """Return each tip's ancestor-to-root path lengths for fast LCA distances."""
    out: dict[str, dict[object, float]] = {}

    def visit(clade, depth: float, ancestors: dict[object, float]) -> None:
        current = dict(ancestors)
        current[clade] = depth
        if clade.is_terminal():
            if not clade.name:
                raise ValueError("tree contains an unnamed terminal")
            out[clade.name] = current
        else:
            for child in clade.clades:
                visit(child, depth + float(child.branch_length or 0.0), current)

    visit(tree.root, 0.0, {})
    return out


# The path dictionaries above are stored root-to-tip, so the first common
# ancestor is the root.  This helper computes the actual LCA without relying
# on global tree state and is used by distance_for_pair below.
def distance_for_pair(paths: dict[str, dict[object, float]], a: str, b: str) -> float:
    pa, pb = paths[a], paths[b]
    common = [node for node in pa if node in pb]
    if not common:
        raise ValueError(f"no common ancestor for {a}, {b}")
    lca = common[-1]
    return max(pa.values()) + max(pb.values()) - 2.0 * pa[lca]


def pair_summary(paths: dict[str, dict[object, float]], a_names: list[str], b_names: list[str]) -> tuple[float | None, int]:
    vals = [distance_for_pair(paths, a, b) for a in a_names for b in b_names if a != b]
    return finite_mean(vals), len(vals)


def load_spatial_functions():
    spec = importlib.util.spec_from_file_location("stage4b_spatial_null_tests", SPATIAL_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("could not load the frozen Stage-4B spatial-null module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def generic_circular_null(rows: list[dict[str, object]], value_field: str) -> dict[str, object]:
    rows = sorted(rows, key=lambda r: (r["mid"], r["window_index"]))
    mask = [r["region"] == "chr16_supergene" for r in rows]
    values = [float(r[value_field]) for r in rows]
    observed = statistics.fmean(v for v, inside in zip(values, mask) if inside) - statistics.fmean(v for v, inside in zip(values, mask) if not inside)
    null = []
    for shift in range(len(rows)):
        shifted = [values[(i - shift) % len(rows)] for i in range(len(rows))]
        delta = statistics.fmean(v for v, inside in zip(shifted, mask) if inside) - statistics.fmean(v for v, inside in zip(shifted, mask) if not inside)
        null.append(delta)
    n_ge = sum(v >= observed - TOL for v in null)
    return {"test": "exact_physical_order_circular_shift", "signal": value_field, "observed_delta": observed, "n_alignments": len(null), "rank": 1 + sum(v > observed + TOL for v in null), "n_ge_observed": n_ge, "p_one_sided": n_ge / len(null)}


def _delta_for_membership(values: dict[int, float], inside: tuple[int, ...], all_ids: list[int]) -> float:
    inside_set = set(inside)
    ivals = [values[i] for i in inside]
    ovals = [values[i] for i in all_ids if i not in inside_set]
    return statistics.fmean(ivals) - statistics.fmean(ovals)


def generic_coordinate_null(rows: list[dict[str, object]], value_field: str) -> tuple[dict[str, object], dict[str, object]]:
    """Generalize the Stage-4B coordinate-event null to any window signal."""
    rows = sorted(rows, key=lambda r: (r["mid"], r["window_index"]))
    span_rows = read_tsv(REGION_MANIFEST)
    span = next(r for r in span_rows if r["region_id"] == "chr16_supergene")
    observed_start, observed_end = float(span["coordinate_start"]), float(span["coordinate_end"])
    width = observed_end - observed_start
    domain_start, domain_end = min(float(r["start"]) for r in rows), max(float(r["end"]) for r in rows)
    min_start, max_start = domain_start, domain_end - width
    raw_events = sorted(set([float(r["mid"]) for r in rows] + [float(r["mid"]) - width for r in rows]))
    events = [x for x in raw_events if min_start <= x <= max_start]
    breaks = sorted(set([min_start, max_start, *events]))
    starts = {min_start, max_start, observed_start, *events}
    starts.update((a + b) / 2 for a, b in zip(breaks, breaks[1:]) if a < b)
    ids = [int(r["window_index"]) for r in rows]
    mids = {int(r["window_index"]): float(r["mid"]) for r in rows}
    values = {int(r["window_index"]): float(r[value_field]) for r in rows}
    keys: dict[tuple[int, ...], float] = {}
    for start in sorted(starts):
        key = tuple(i for i in ids if start <= mids[i] <= start + width)
        if key and len(key) < len(ids):
            keys.setdefault(key, start)
    observed_key = tuple(i for i in ids if observed_start <= mids[i] <= observed_end)
    if observed_key not in keys:
        keys[observed_key] = observed_start
    deltas = [_delta_for_membership(values, key, ids) for key in keys]
    observed = _delta_for_membership(values, observed_key, ids)
    n_ge = sum(v >= observed - TOL for v in deltas)
    summary = {"test": "unique_membership_state_coordinate_sensitivity", "signal": value_field, "observed_delta": observed, "n_alignments": len(deltas), "rank": 1 + sum(v > observed + TOL for v in deltas), "n_ge_observed": n_ge, "p_one_sided": n_ge / len(deltas), "observed_interval_start": observed_start, "observed_interval_end": observed_end, "interval_width": width, "coordinate_domain_start": domain_start, "coordinate_domain_end": domain_end}
    segment_rows = []
    extreme_length = 0.0
    total_length = max_start - min_start
    for left, right in zip(breaks, breaks[1:]):
        if right <= left:
            continue
        start = (left + right) / 2
        key = tuple(i for i in ids if start <= mids[i] <= start + width)
        if not key or len(key) == len(ids):
            continue
        delta = _delta_for_membership(values, key, ids)
        length = right - left
        if delta >= observed - TOL:
            extreme_length += length
        segment_rows.append({"start_min": left, "start_max": right, "segment_length": length, "delta": delta})
    weighted = {"test": "continuous_coordinate_length_weighted_sensitivity", "signal": value_field, "observed_delta": observed, "n_alignments": len(segment_rows), "rank": None, "n_ge_observed": None, "p_one_sided": extreme_length / total_length if total_length > 0 else None, "coordinate_domain_start": domain_start, "coordinate_domain_end": domain_end}
    return summary, weighted


def plot_series(ax, xs: list[float], ys: list[float], color: str) -> None:
    """Plot observed windows without drawing across large unobserved gaps."""
    if len(xs) < 2:
        ax.plot(xs, ys, marker="o", color=color)
        return
    gaps = [b - a for a, b in zip(xs, xs[1:]) if b > a]
    threshold = 3.0 * statistics.median(gaps) if gaps else math.inf
    start = 0
    for i, gap in enumerate([b - a for a, b in zip(xs, xs[1:])], start=0):
        if gap > threshold:
            ax.plot(xs[start:i + 1], ys[start:i + 1], color=color, linewidth=0.8, alpha=0.55)
            start = i + 1
    ax.plot(xs[start:], ys[start:], color=color, linewidth=0.8, alpha=0.55)


def load_data() -> tuple[list[dict[str, object]], dict[str, str], dict[str, str]]:
    trees = read_tsv(TREE_TABLE)
    index = {int(r["window_index"]): r for r in read_tsv(WINDOW_INDEX)}
    support = {int(r["window_index"]): r for r in read_tsv(SUPPORT)}
    mapping_rows = read_tsv(MAPPING)
    mapping = {r["sample"]: r["astral_group_safe"] for r in mapping_rows}
    if len(trees) != 213 or len(index) != 213 or len(support) != 213:
        raise ValueError("frozen fire-ant window tables must each contain 213 windows")
    rows: list[dict[str, object]] = []
    for tr in trees:
        wi = int(tr["window_index"])
        if wi not in index or wi not in support:
            raise ValueError(f"window {wi} missing from frozen index/support")
        row = {**index[wi], "tree_newick": tr["tree_newick"], "D": support[wi]["D"], "q_S": support[wi]["q_S"], "q_H": support[wi]["q_H"], "q_3": support[wi]["q_3"]}
        row["window_index"] = wi
        row["start"] = float(row["start"])
        row["end"] = float(row["end"])
        row["mid"] = (row["start"] + row["end"]) / 2.0
        rows.append(row)
    rows.sort(key=lambda r: int(r["window_index"]))
    return rows, mapping, {r["published_group"]: r["astral_safe_label"] for r in read_tsv(GROUP_MAP)}


def calculate(rows: list[dict[str, object]], mapping: dict[str, str]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    output, pair_audit = [], []
    nonfocal = [sample for sample, group in mapping.items() if group not in FOCAL]
    for row in rows:
        tree = parse_tree(str(row["tree_newick"]))
        paths = ancestor_depths(tree)
        groups = {group: [s for s, g in mapping.items() if g == group and s in paths] for group in FOCAL}
        valid = all(groups[g] for g in FOCAL)
        result = dict(row)
        result["status"] = "ok" if valid else "missing_focal_group"
        pair_means: dict[tuple[str, str], float | None] = {}
        pair_ns: dict[tuple[str, str], int] = {}
        for a, b in ALL_PAIRS:
            mean, n = pair_summary(paths, groups[a], groups[b]) if valid else (None, 0)
            pair_means[(a, b)], pair_ns[(a, b)] = mean, n
            result[GROUP_FIELDS[(a, b)] + "_mean"] = mean
            result[GROUP_FIELDS[(a, b)] + "_n"] = n
            pair_audit.append({"window_index": row["window_index"], "group_pair": f"{a}__{b}", "class": "opposite" if (a, b) in OPPOSITE else "same", "n_individual_pairs": n, "mean_distance": mean})
        opp = [pair_means[p] for p in OPPOSITE if pair_means[p] is not None]
        same = [pair_means[p] for p in SAME if pair_means[p] is not None]
        all_means = [pair_means[p] for p in ALL_PAIRS if pair_means[p] is not None]
        opp_weighted_vals, same_weighted_vals = [], []
        for pair in OPPOSITE:
            a, b = pair
            opp_weighted_vals.extend(distance_for_pair(paths, x, y) for x in groups[a] for y in groups[b])
        for pair in SAME:
            a, b = pair
            same_weighted_vals.extend(distance_for_pair(paths, x, y) for x in groups[a] for y in groups[b])
        mean_opp = finite_mean(opp)
        mean_same = finite_mean(same)
        mean_all = finite_mean(all_means)
        result["mean_opposite"] = mean_opp
        result["mean_same"] = mean_same
        result["mean_all_focal_pair_distances"] = mean_all
        result["A_distance"] = mean_opp - mean_same if mean_opp is not None and mean_same is not None else None
        result["A_relative"] = result["A_distance"] / mean_all if result["A_distance"] is not None and mean_all and mean_all != 0 else None
        weighted_opp, weighted_same = finite_mean(opp_weighted_vals), finite_mean(same_weighted_vals)
        weighted_all_vals = []
        for a, b in ALL_PAIRS:
            weighted_all_vals.extend(distance_for_pair(paths, x, y) for x in groups[a] for y in groups[b])
        weighted_all = finite_mean(weighted_all_vals)
        result["individual_weighted_opposite"] = weighted_opp
        result["individual_weighted_same"] = weighted_same
        result["A_distance_individual_weighted"] = weighted_opp - weighted_same if weighted_opp is not None and weighted_same is not None else None
        result["A_relative_individual_weighted"] = result["A_distance_individual_weighted"] / weighted_all if result["A_distance_individual_weighted"] is not None and weighted_all else None
        out_vals = [distance_for_pair(paths, x, y) for g in FOCAL for x in groups[g] for y in nonfocal if y in paths]
        outgroup_distance = finite_mean(out_vals)
        result["mean_focal_to_nonfocal_distance"] = outgroup_distance
        result["A_distance_outgroup_normalized"] = result["A_distance"] / outgroup_distance if result["A_distance"] is not None and outgroup_distance else None
        result["n_focal_individuals"] = sum(len(groups[g]) for g in FOCAL)
        result["n_nonfocal_individuals"] = len([s for s in nonfocal if s in paths])
        output.append(result)
    return output, pair_audit


def region_rows(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    return [r for r in rows if r["chrom"] == "chr16" and r["region"] in {"chr16A", "chr16_supergene", "chr16B"} and r.get("A_distance") is not None]


def plot_track(rows: list[dict[str, object]], out_pdf: Path, out_png: Path) -> None:
    r = sorted(region_rows(rows), key=lambda x: (x["mid"], x["window_index"]))
    fig, axes = plt.subplots(3, 1, figsize=(12, 8.5), sharex=True, constrained_layout=True)
    start, end = 11680438, 27917498
    xs = [float(x["mid"]) / 1e6 for x in r]
    colors = ["#b2182b" if x["region"] == "chr16_supergene" else "#2166ac" for x in r]
    for ax, key, label, color in [(axes[0], "A_distance", "A_distance: opposite − same (substitution units)", "#762a83"), (axes[1], "A_relative", "A_relative", "#1b7837"), (axes[2], "D", "D(w) = q_H − q_S", "#333333")]:
        ys = [float(x[key]) if x.get(key) is not None else math.nan for x in r]
        ax.scatter(xs, ys, c=colors, s=18, alpha=0.85, linewidths=0)
        plot_series(ax, xs, ys, color)
        ax.axvspan(start / 1e6, end / 1e6, color="#fdae61", alpha=0.15, zorder=0)
        ax.set_ylabel(label)
        ax.grid(alpha=0.2)
    axes[-1].set_xlabel("chromosome 16 coordinate (Mb); points are observed BUSCO-window midpoints")
    fig.suptitle("Fire-ant chromosome-16 patristic distance and quartet tracks")
    fig.savefig(out_pdf)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def write_text_outputs(rows: list[dict[str, object]], spatial: list[dict[str, object]], manifest: dict[str, object]) -> None:
    ok = [r for r in rows if r.get("A_distance") is not None]
    supergene = [r for r in ok if r["region"] == "chr16_supergene"]
    outside = [r for r in ok if r["region"] in {"chr16A", "chr16B"}]
    def summary(field: str, a: list[dict[str, object]], b: list[dict[str, object]]) -> tuple[float | None, float | None, float | None]:
        av, bv = finite_mean(float(r[field]) for r in a if r.get(field) is not None), finite_mean(float(r[field]) for r in b if r.get(field) is not None)
        return av, bv, av - bv if av is not None and bv is not None else None
    vals = {f: summary(f, supergene, outside) for f in ["A_distance", "A_relative", "A_distance_individual_weighted", "A_relative_individual_weighted", "A_distance_outgroup_normalized"]}
    report = [
        "# Fire-ant chromosome-16 patristic distance track (Stage 9)", "",
        "This analysis uses the 213 published RAxML local trees and the verified Stage-6C individual-to-group mapping. Branch-length distances are patristic distances in substitution units. They are not coalescence times and no molecular-clock assumption is made.", "",
        f"The primary contrast is group-balanced: each of the four opposite-haplotype group pairs is averaged within pair and then equally across pairs; each of the two same-haplotype cross-species pairs is treated analogously. The supergene/outside comparison uses {len(supergene)} and {len(outside)} chromosome-16 windows with complete focal groups.", "",
        "## Results", "",
    ]
    for f, (inside, out, delta) in vals.items():
        report.append(f"- `{f}`: supergene mean={inside}, outside mean={out}, difference={delta}.")
    report += ["", "The individual-pair-weighted values are a sensitivity analysis and can be influenced by unequal group sample sizes. The outgroup-normalized sensitivity divides the group-balanced distance contrast by the mean distance from focal individuals to nonfocal individuals observed in the same local tree.", "", "## Spatial null", ""]
    for row in spatial:
        report.append(f"- `{row['signal']}` / `{row['test']}`: observed supergene-minus-outside contrast={row['observed_delta']:.12g}; p={row['p_one_sided']:.12g} over {row['n_alignments']} alignments or coordinate states.")
    report += ["", "## Limitations", "", "The source trees are substitution-unit RAxML trees with irregular four-BUSCO window spacing. The track is descriptive and spatially tested against the frozen author-labeled supergene window mask; the interval is an observed BUSCO-window span rather than a nucleotide-resolved breakpoint annotation.", ""]
    (RESULTS / "stage9_report.md").write_text("\n".join(report))
    methods = """# Methods\n\nFor each published RAxML local tree, patristic distances were calculated between all individual tips belonging to the four focal Stage-6C groups: invicta/macdonaghi_SB, invicta/macdonaghi_Sb, richteri_SB, and richteri_Sb. Six unordered group pairs were retained. Opposite-haplotype pairs were the four SB/Sb combinations; same-haplotype pairs were the two cross-species SB/SB and Sb/Sb combinations. Distances were averaged within each group pair and then averaged equally across pairs, so unequal sample sizes did not determine the primary track. An individual-pair-weighted contrast was computed as a sensitivity analysis.\n\n`A_distance(w)` is the balanced opposite-haplotype mean minus the balanced same-haplotype mean. `A_relative(w)` divides this value by the equally weighted mean of the six focal group-pair distances, with zero denominators recorded as missing. An outgroup-normalized sensitivity divides `A_distance` by the mean focal-to-nonfocal patristic distance in that tree.\n\nThe existing Stage-4A `D(w)=q_H-q_S` values were joined by frozen window index. Chromosome-16 windows were plotted at their observed midpoints; no sequence was interpolated across absent coordinate intervals. The supergene/outside comparison used the existing Stage-4B physical-order exact circular-shift framework: the observed 52-window supergene mask was held fixed while the complete 96-window chromosome-16 signal was rotated through all alignments.\n"""
    (RESULTS / "stage9_methods_text.md").write_text(methods)
    (RESULTS / "stage9_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> None:
    global RESULTS
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=RESULTS)
    args = parser.parse_args()
    RESULTS = args.results_dir.resolve()
    RESULTS.mkdir(parents=True, exist_ok=True)
    rows, mapping, _ = load_data()
    calculated, pair_audit = calculate(rows, mapping)
    distance_fields = ["window_index", "chrom", "start", "end", "mid", "region", "D", "q_S", "q_H", "q_3", "status", "n_focal_individuals", "n_nonfocal_individuals", "mean_opposite", "mean_same", "mean_all_focal_pair_distances", "A_distance", "A_relative", "individual_weighted_opposite", "individual_weighted_same", "A_distance_individual_weighted", "A_relative_individual_weighted", "mean_focal_to_nonfocal_distance", "A_distance_outgroup_normalized"]
    for a, b in ALL_PAIRS:
        p = GROUP_FIELDS[(a, b)]
        distance_fields += [p + "_mean", p + "_n"]
    write_tsv(RESULTS / "stage9_window_patristic_distance.tsv", calculated, distance_fields)
    write_tsv(RESULTS / "stage9_pairwise_sample_sizes.tsv", pair_audit, ["window_index", "group_pair", "class", "n_individual_pairs", "mean_distance"])
    chr16 = region_rows(calculated)
    spatial = []
    for field in ["A_distance", "A_relative", "A_distance_individual_weighted", "A_relative_individual_weighted", "A_distance_outgroup_normalized"]:
        spatial.append(generic_circular_null(chr16, field))
        coordinate, weighted = generic_coordinate_null(chr16, field)
        spatial.extend([coordinate, weighted])
    write_tsv(RESULTS / "stage9_spatial_null.tsv", spatial, ["test", "signal", "observed_delta", "n_alignments", "rank", "n_ge_observed", "p_one_sided", "observed_interval_start", "observed_interval_end", "interval_width", "coordinate_domain_start", "coordinate_domain_end"])
    plot_track(calculated, FIGURES / "fire_ants_chr16_A_distance_track.pdf", FIGURES / "fire_ants_chr16_A_distance_track.png")
    outputs = [RESULTS / "stage9_window_patristic_distance.tsv", RESULTS / "stage9_pairwise_sample_sizes.tsv", RESULTS / "stage9_spatial_null.tsv", RESULTS / "stage9_report.md", RESULTS / "stage9_methods_text.md", FIGURES / "fire_ants_chr16_A_distance_track.pdf", FIGURES / "fire_ants_chr16_A_distance_track.png"]
    manifest = {"stage": "9_patristic_distance_track", "seed": SEED, "n_windows": len(rows), "n_chr16_windows": len(chr16), "inputs": {str(p.relative_to(ROOT)): sha256(p) for p in [TREE_TABLE, WINDOW_INDEX, MAPPING, GROUP_MAP, SUPPORT, REGION_MANIFEST]}, "outputs": {}, "python": sys.version, "platform": platform.platform(), "primary_units": "patristic substitution units", "clock_assumption": False, "spatial_null": "exact physical-order circular shifts of the 96-window chromosome-16 track with frozen 52-window supergene mask"}
    write_text_outputs(calculated, spatial, manifest)
    # The report and methods were written just before the manifest; refresh the
    # output hashes after writing all text artifacts.
    manifest["outputs"] = {str(p.relative_to(ROOT)): sha256(p) for p in outputs + [RESULTS / "stage9_manifest.json"] if p.exists() and p.name != "stage9_manifest.json"}
    (RESULTS / "stage9_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
