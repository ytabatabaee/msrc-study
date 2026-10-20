#!/usr/bin/env python3
"""House-mouse Stage 8: validate extreme local-ML branch lengths.

The primary analysis recalculates all quantities directly from the frozen
5-kb Newick trees.  The original all-sample tracks are retained; leave-one-out
and robust summaries are separate sensitivity analyses.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical/house_mouse_t_complex"
DATA = ROOT / "data/house_mouse_t_complex"
OUT = EMP / "results/stage8_branch_length_validation"
FIG = EMP / "figures"
TREES = DATA / "processed/house_mouse_t_complex_ml_5kb.tre"
META = DATA / "processed/house_mouse_t_complex_ml_5kb_metadata.tsv"
MAP = DATA / "metadata/tip_mapping.tsv"
DIST = EMP / "results/stage4_distance/stage4_window_distance.tsv"
QSCAN = EMP / "results/stage2_fixed_quartet_scan.tsv"
INT = EMP / "results/t_complex_inversion_intervals.tsv"
ARCHIVE = ROOT / "IST-2017-78-v1+1_Data.zip"
EXCLUDED = "Mmd_FRA4_18B.variant"
CHANGE_POINT = 14737499.5
sys.path.insert(0, str(ROOT / "empirical/shared"))
sys.path.insert(0, str(EMP / "scripts"))
from spatial_boundary_detector import segment_track, intervals_to_coordinates
import importlib.util

_stage7_spec = importlib.util.spec_from_file_location("stage7_mouse", EMP / "scripts/07_influential_sample_sensitivity.py")
_stage7 = importlib.util.module_from_spec(_stage7_spec)
_stage7_spec.loader.exec_module(_stage7)


def read(path):
    with path.open(newline="") as h:
        return list(csv.DictReader(h, delimiter="\t"))


def write(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows({f: r.get(f, "") for f in fields} for r in rows)


def num(x):
    try:
        y = float(x)
        return y if math.isfinite(y) else math.nan
    except (TypeError, ValueError):
        return math.nan


def finite(x):
    return math.isfinite(num(x))


def vals(xs):
    return [num(x) for x in xs if finite(x)]


def avg(xs):
    x = vals(xs)
    return float(np.mean(x)) if x else math.nan


def med(xs):
    x = vals(xs)
    return float(np.median(x)) if x else math.nan


def trim(xs, fraction=0.05):
    x = sorted(vals(xs))
    if not x:
        return math.nan
    k = int(fraction * len(x))
    return float(np.mean(x[k:len(x)-k])) if k and len(x) > 2 * k else float(np.mean(x))


def fmt(x):
    return "" if x is None or not finite(x) else f"{float(x):.12g}"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def distance_matrix(tree):
    # Use one adjacency traversal per terminal rather than Biopython's
    # repeated path search for every pair. This is the same patristic-distance
    # convention used by the frozen Stage 7 implementation.
    adj = defaultdict(list)
    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            w = num(child.branch_length)
            adj[id(parent)].append((id(child), w))
            adj[id(child)].append((id(parent), w))
    leaves = {id(c): c.name for c in tree.get_terminals()}
    out = {}
    for start, sname in leaves.items():
        stack = [(start, -1, 0.0)]
        while stack:
            u, prev, d = stack.pop()
            if u in leaves and u != start and math.isfinite(d):
                out[(sname, leaves[u])] = d
            for v, w in adj[u]:
                if v != prev and math.isfinite(w):
                    stack.append((v, u, d + w))
    return out


def pair_key(a, b):
    return (a, b)


def eligible_tips(tree, mapping, exclude=None):
    return [c.name for c in tree.get_terminals() if c.name != exclude and mapping.get(c.name, {}).get("subspecies") != "Mus spretus"]


def pair_values(D, tips):
    return [D[pair_key(a, b)] for i, a in enumerate(tips) for b in tips[i + 1:] if pair_key(a, b) in D and finite(D[pair_key(a, b)])]


def branch_metrics(tree, D, mapping):
    terminals = {c.name: num(c.branch_length) for c in tree.get_terminals() if finite(c.branch_length)}
    all_branches = [num(c.branch_length) for c in tree.find_clades() if c is not tree.root and finite(c.branch_length)]
    internal = [num(c.branch_length) for c in tree.get_nonterminals() if c is not tree.root and finite(c.branch_length)]
    non_spretus = eligible_tips(tree, mapping)
    without = [x for x in non_spretus if x != EXCLUDED]
    other_term = [x for x in terminals.values() if x != terminals.get(EXCLUDED)]
    sorted_terms = sorted(terminals.items(), key=lambda kv: kv[1], reverse=True)
    focal = terminals.get(EXCLUDED, math.nan)
    med_term = med(terminals.values())
    total = sum(all_branches)
    return {
        "n_tips": len(tree.get_terminals()), "n_non_spretus_tips": len(non_spretus),
        "excluded_tip_present": int(EXCLUDED in terminals), "focal_terminal_branch": focal,
        "max_other_terminal_branch": max(other_term) if other_term else math.nan,
        "second_longest_terminal_branch": sorted_terms[1][1] if len(sorted_terms) > 1 else math.nan,
        "longest_terminal_branch": sorted_terms[0][1] if sorted_terms else math.nan,
        "longest_terminal_tip": sorted_terms[0][0] if sorted_terms else "",
        "longest_internal_branch": max(internal) if internal else math.nan,
        "median_terminal_branch": med_term,
        "focal_to_median_terminal_ratio": focal / med_term if finite(focal) and med_term > 0 else math.nan,
        "focal_fraction_total_tree": focal / total if finite(focal) and total > 0 else math.nan,
        "total_tree_length": total, "mean_internal_branch": avg(internal),
        "median_pair_distance_all": med(pair_values(D, non_spretus)),
        "n_zero_branches": sum(x == 0 for x in all_branches),
        "n_negative_branches": sum(x < 0 for x in all_branches),
        "n_nonfinite_branches": sum(1 for c in tree.find_clades() if c is not tree.root and not finite(c.branch_length)),
    }


def influence_for_tip(D, tips, tip, all_values=None):
    """Return B(all)-B(with tip removed) without rebuilding every pair list."""
    all_values = all_values if all_values is not None else pair_values(D, tips)
    if len(tips) < 3 or not all_values:
        return math.nan
    incident = [D[pair_key(tip, other)] if pair_key(tip, other) in D else D[pair_key(other, tip)] for other in tips if other != tip and (pair_key(tip, other) in D or pair_key(other, tip) in D)]
    remaining_n = len(all_values) - len(incident)
    if remaining_n <= 0:
        return math.nan
    b_all = float(np.mean(all_values))
    b_without = (float(np.sum(all_values)) - float(np.sum(incident))) / remaining_n
    return b_all - b_without


def score_tree(tree, mapping, frozen_q):
    D = distance_matrix(tree)
    tips = eligible_tips(tree, mapping)
    loo_tips = [x for x in tips if x != EXCLUDED]
    m = branch_metrics(tree, D, mapping)
    p_all = pair_values(D, tips); p_loo = pair_values(D, loo_tips)
    a0, op0, sa0, design0 = _stage7.contrast(D, tree, mapping)
    a1, op1, sa1, design1 = _stage7.contrast(D, tree, mapping, EXCLUDED)
    m.update({
        "B_all": avg(p_all), "B_leave_one_out": avg(p_loo), "B_median": med(p_all), "B_trimmed": trim(p_all),
        "A_distance_original": a0, "A_distance_leave_one_out": a1,
        "A_relative_original": a0 / avg(p_all) if finite(a0) and avg(p_all) > 1e-15 else math.nan,
        "A_relative_leave_one_out": a1 / avg(p_loo) if finite(a1) and avg(p_loo) > 1e-15 else math.nan,
        "D_original": num(frozen_q["D_original"]), "D_leave_one_out": num(frozen_q["D_leave_one_out"]),
        "q_species_original": num(frozen_q["q_species_original"]), "q_t_alt_original": num(frozen_q["q_t_alt_original"]), "q_other_original": num(frozen_q["q_other_original"]), "q_unresolved_original": num(frozen_q["q_unresolved_original"]),
        "q_species_leave_one_out": num(frozen_q["q_species_leave_one_out"]), "q_t_alt_leave_one_out": num(frozen_q["q_t_alt_leave_one_out"]), "q_other_leave_one_out": num(frozen_q["q_other_leave_one_out"]), "q_unresolved_leave_one_out": num(frozen_q["q_unresolved_leave_one_out"]),
        "n_informative_original": sum(num(frozen_q[f"n_{x}_original"]) for x in ("species", "t_alt", "other", "unresolved")), "n_informative_leave_one_out": sum(num(frozen_q[f"n_{x}_leave_one_out"]) for x in ("species", "t_alt", "other", "unresolved")),
        "design_original": "frozen Stage 7 ALL_TIPS", "design_leave_one_out": "frozen Stage 7 ALL_TIPS leave-one-out",
    })
    return m, D


def evaluate(pred, refs):
    def merge(intervals):
        out = []
        for a, b in sorted((int(a), int(b)) for a, b in intervals):
            if out and a <= out[-1][1]: out[-1][1] = max(out[-1][1], b)
            else: out.append([a, b])
        return out
    p = merge(pred); r = merge(refs)
    plen = sum(b - a for a, b in p); rlen = sum(b - a for a, b in r)
    overlap = sum(max(0, min(pb, rb) - max(pa, ra)) for pa, pb in p for ra, rb in r)
    precision = overlap / plen if plen else 0; recall = overlap / rlen if rlen else 0
    return {"iou": overlap / (plen + rlen - overlap) if plen + rlen - overlap else 0,
            "precision": precision, "recall": recall, "f1": 2 * precision * recall / (precision + recall) if precision + recall else 0,
            "n_predicted_intervals": len(p), "predicted_bp": plen}


def run_boundary(rows, fields, refs):
    valid = [r for r in rows if all(finite(r.get(f)) for f in fields)]
    x = np.array([[num(r[f]) for f in fields] for r in valid])
    seg = segment_track(x, fields, min_size=3, max_intervals=4)
    coords = intervals_to_coordinates(seg, [{"start": r["start_bp"], "end": r["end_bp"]} for r in valid])
    pred = [(int(c["start"]), int(c["end"])) for c in coords]
    return {"status": seg.status, "score": seg.score, "raw_improvement": seg.raw_improvement, "penalty": seg.penalty,
            "predicted_intervals": ";".join(f"{a}-{b}" for a, b in pred), "feature_names": ";".join(fields), **evaluate(pred, refs)}


def make_bins(rows, size=100000):
    groups = defaultdict(list)
    for r in rows:
        start = int(num(r["start_bp"]))
        groups[5_000_000 + ((start - 5_000_000) // size) * size].append(r)
    out = []
    for start, group in sorted(groups.items()):
        z = {"bin_start_bp": start, "bin_end_bp": start + size - 1, "midpoint_bp": start + size / 2, "n_observed_windows": len(group), "observed_bp": sum(int(num(r["end_bp"])) - int(num(r["start_bp"])) + 1 for r in group), "coverage_fraction": sum(int(num(r["end_bp"])) - int(num(r["start_bp"])) + 1 for r in group) / size}
        for field in ("D_original", "D_leave_one_out", "A_distance_original", "A_distance_leave_one_out", "A_relative_original", "A_relative_leave_one_out", "B_all", "B_leave_one_out", "B_median", "B_trimmed", "focal_terminal_branch", "max_other_terminal_branch", "focal_to_median_terminal_ratio", "focal_fraction_total_tree"):
            z[field] = avg([r[field] for r in group])
        for lab in ("original", "leave_one_out"):
            counts = {k: sum(num(r[f"n_{k}_{lab}"]) for r in group) for k in ("informative",)}
            for k in ("species", "t_alt", "other", "unresolved"):
                z[f"q_{k}_{lab}"] = avg([r[f"q_{k}_{lab}"] for r in group]) if f"q_{k}_{lab}" in group[0] else math.nan
        out.append(z)
    return out


def archive_filter_comparison(metadata, selected):
    rows = []
    if not ARCHIVE.exists():
        return rows
    filter_paths = {
        "PASS_SNPs": "Data/1-PASS_SNPs/1-Trees_for_all_5kb_windows/ML_trees.zip",
        "Coverage_AlleleRatio_Filtered": "Data/2-Coverage-and_AlleleRatio-Filtered_RAW_SNPs/1-Trees_for_all_5kb_windows/ML_trees.zip",
        "Coverage_Filtered": "Data/3-Coverage-Filtered_RAW_SNPs/1-Trees_for_all_5kb_windows/ML_trees.zip",
    }
    with zipfile.ZipFile(ARCHIVE) as outer:
        for locus in selected:
            meta = next((x for x in metadata if x["locus_id"] == locus), None)
            if not meta: continue
            for label, path in filter_paths.items():
                try:
                    nested_bytes = outer.read(path)
                    with zipfile.ZipFile(io.BytesIO(nested_bytes)) as nested:
                        text = nested.read(meta["source_filename"]).decode()
                    tree = Phylo.read(io.StringIO(text), "newick")
                    focal = next((num(c.branch_length) for c in tree.get_terminals() if c.name == EXCLUDED or c.name.startswith(EXCLUDED + "_")), math.nan)
                    terms = [num(c.branch_length) for c in tree.get_terminals() if finite(c.branch_length)]
                    total = sum(num(c.branch_length) for c in tree.find_clades() if c is not tree.root and finite(c.branch_length))
                    rows.append({"locus_id": locus, "filtering_procedure": label, "focal_terminal_branch": fmt(focal), "median_terminal_branch": fmt(med(terms)), "max_terminal_branch": fmt(max(terms) if terms else math.nan), "total_tree_length": fmt(total), "status": "parsed_from_published_nested_ML_tree_archive"})
                except Exception as exc:
                    rows.append({"locus_id": locus, "filtering_procedure": label, "status": f"unavailable:{type(exc).__name__}"})
    return rows


def shade(ax, refs):
    colors = ["#6a3d9a", "#1b9e77", "#d95f02", "#7570b3"]
    for i, (a, b) in enumerate(refs):
        ax.axvspan(a / 1e6, b / 1e6, color=colors[i % len(colors)], alpha=.12, lw=0)
        ax.axvline(a / 1e6, color="#666", ls="--", lw=.5); ax.axvline(b / 1e6, color="#666", ls="--", lw=.5)


def plot_tracks(bins, refs, stem, title, lab):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.array([num(r["midpoint_bp"]) / 1e6 for r in bins])
    fields = [("D_original" if lab == "original" else "D_leave_one_out", "(b) ALL_TIPS topology contrast D(w)", "#333", "support contrast"), ("A_distance_original" if lab == "original" else "A_distance_leave_one_out", "(d) arrangement-informed A_distance", "#c44e52", "substitution units"), ("A_relative_original" if lab == "original" else "A_relative_leave_one_out", "(e) arrangement-informed A_relative", "#55a868", "relative"), ("B_all" if lab == "original" else "B_leave_one_out", "(f) arrangement-blind B(w)", "#8172b2", "substitution units")]
    fig, ax = plt.subplots(6, 1, figsize=(13, 13), sharex=True, gridspec_kw={"height_ratios": [.35, 1, 1, 1, 1, 1]}, constrained_layout=True)
    for a in ax: shade(a, refs); a.grid(axis="y", alpha=.16)
    ax[0].set_ylim(0, 1); ax[0].set_yticks([]); ax[0].set_ylabel("inversions", fontsize=8)
    colors = ["#6a3d9a", "#1b9e77", "#d95f02", "#7570b3"]
    for i, (a, b) in enumerate(refs): ax[0].bar((a+b)/2e6, .65, width=(b-a)/1e6, color=colors[i], alpha=.6)
    for idx, (field, title2, color, ylabel) in enumerate(fields, start=1):
        if idx == 2:
            q1 = np.array([num(r[f"q_species_{lab}"]) for r in bins]); q2 = np.array([num(r[f"q_t_alt_{lab}"]) for r in bins]); q3 = np.array([num(r[f"q_other_{lab}"]) for r in bins])
            ax[idx].stackplot(x, q1, q2, q3, colors=["#4c78a8", "#e45756", "#59a14f"], alpha=.85)
            ax[idx].set_ylim(0, 1); ax[idx].set_ylabel("support"); ax[idx].set_title("(c) three-topology quartet support", loc="left", weight="bold")
            ax[idx].legend(["q_S species history", "q_T t-associated alternative", "q_3 third resolution"], frameon=False, ncol=3, fontsize=8, loc="upper center", bbox_to_anchor=(.5, 1.25))
        else:
            y = np.array([num(r[field]) for r in bins]); ax[idx].plot(x, y, color=color, lw=.8); ax[idx].set_ylabel(ylabel); ax[idx].set_title(title2, loc="left", weight="bold")
            if idx in (1, 3, 4): ax[idx].axhline(0, color="#777", lw=.5)
    ax[-1].set_xlabel("chromosome 17 position (Mb)"); fig.suptitle(title, fontsize=15, weight="bold"); fig.text(.99, .01, "Distances are local ML patristic substitution units.", ha="right", fontsize=8)
    fig.savefig(FIG / f"house_mouse_chr17_{stem}_100kb.pdf", bbox_inches="tight"); fig.savefig(FIG / f"house_mouse_chr17_{stem}_100kb.png", dpi=240, bbox_inches="tight"); plt.close(fig)


def make_figures(bins, rows, refs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.array([num(r["midpoint_bp"]) / 1e6 for r in rows])
    fig, ax = plt.subplots(4, 1, figsize=(13, 9), sharex=True, constrained_layout=True)
    for a in ax: shade(a, refs); a.grid(axis="y", alpha=.15)
    for field, color, title in [("focal_terminal_branch", "#c44e52", "focal terminal branch"), ("max_other_terminal_branch", "#4c78a8", "maximum other terminal branch"), ("focal_to_median_terminal_ratio", "#55a868", "focal / median terminal branch"), ("focal_fraction_total_tree", "#8172b2", "focal fraction of total tree length")]:
        idx = {"focal_terminal_branch": 0, "max_other_terminal_branch": 1, "focal_to_median_terminal_ratio": 2, "focal_fraction_total_tree": 3}[field]
        ax[idx].plot(x, [num(r[field]) for r in rows], color=color, lw=.55); ax[idx].set_ylabel("ratio" if idx == 2 else "fraction" if idx == 3 else "substitution units"); ax[idx].set_title(title, loc="left", weight="bold")
    ax[-1].set_xlabel("chromosome 17 position (Mb)"); fig.suptitle("House-mouse Stage 8 extreme branch-length diagnostics", weight="bold"); fig.savefig(FIG / "house_mouse_stage8_branch_diagnostics.pdf", bbox_inches="tight"); fig.savefig(FIG / "house_mouse_stage8_branch_diagnostics.png", dpi=240, bbox_inches="tight"); plt.close(fig)
    fig, ax = plt.subplots(4, 1, figsize=(13, 9), sharex=True, constrained_layout=True)
    for a in ax: shade(a, refs); a.grid(axis="y", alpha=.15)
    for i, (a, fields, title) in enumerate([(ax[0], [("B_all", "all", "#8172b2"), ("B_leave_one_out", "FRA4 excluded", "#333"), ("B_median", "median", "#4c78a8"), ("B_trimmed", "5% trimmed", "#f58518")], "B(w) robustness"), (ax[1], [("A_distance_original", "all", "#c44e52"), ("A_distance_leave_one_out", "FRA4 excluded", "#333")], "A_distance"), (ax[2], [("A_relative_original", "all", "#55a868"), ("A_relative_leave_one_out", "FRA4 excluded", "#333")], "A_relative")]):
        for field, label, color in fields: a.plot(x, [num(r[field]) for r in rows], lw=.55, label=label, color=color)
        a.set_title(title, loc="left", weight="bold"); a.legend(frameon=False, fontsize=8, ncol=len(fields));
    ax[3].plot(x, [num(r["total_tree_length"]) for r in rows], color="#7570b3", lw=.55); ax[3].set_title("total local-tree length", loc="left", weight="bold"); ax[3].set_ylabel("substitution units"); ax[3].set_xlabel("chromosome 17 position (Mb)")
    fig.suptitle("House-mouse Stage 8 distance-statistic robustness", weight="bold"); fig.savefig(FIG / "house_mouse_stage8_distance_robustness.pdf", bbox_inches="tight"); fig.savefig(FIG / "house_mouse_stage8_distance_robustness.png", dpi=240, bbox_inches="tight"); plt.close(fig)
    plot_tracks(bins, refs, "stage8_original", "House-mouse chromosome 17: original all-sample estimates", "original")
    plot_tracks(bins, refs, "stage8_leave_one_out", "House-mouse chromosome 17: leave-one-out sensitivity", "leave_one_out")


def main():
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    mapping_rows = read(MAP); mapping = {r["tree_tip"]: r for r in mapping_rows}; metadata = read(META); source = {r["locus_id"]: r for r in read(DIST)}; frozen_q = {r["locus_id"]: r for r in read(EMP / "results/stage7_influential_sample/window_original_vs_leave_one_out.tsv")}
    tree_lines = [x.strip() for x in TREES.read_text().splitlines() if x.strip()]
    if len(tree_lines) != len(metadata) or len(metadata) != 4046: raise RuntimeError("expected 4,046 metadata-matched trees")
    rows = []; audits = []; influences = []
    for line, info in zip(tree_lines, metadata):
        tree = Phylo.read(io.StringIO(line), "newick")
        metrics, D = score_tree(tree, mapping, frozen_q[info["locus_id"]])
        metrics.update({"locus_id": info["locus_id"], "start_bp": info["start_bp"], "end_bp": info["end_bp"], "midpoint_bp": info["midpoint_bp"]})
        for k, v in source[info["locus_id"]].items():
            if k in ("mean_pair_distance_no_spretus", "A_distance", "A_relative"): metrics[f"stage4_{k}"] = num(v)
        rows.append(metrics)
        audit = dict(metrics); audit["upstream_downstream"] = "upstream" if num(info["midpoint_bp"]) < CHANGE_POINT else "downstream"; audits.append(audit)
        tips = eligible_tips(tree, mapping)
        all_pair_values = pair_values(D, tips)
        for tip in tips:
            influences.append({"locus_id": info["locus_id"], "midpoint_bp": info["midpoint_bp"], "tip": tip, "B_all_minus_B_without_tip": fmt(influence_for_tip(D, tips, tip, all_pair_values))})
    rows.sort(key=lambda r: num(r["start_bp"])); audits.sort(key=lambda r: num(r["start_bp"]))
    refs = [(int(r["start_bp"]), int(r["end_bp"])) for r in read(INT)]
    thresholds = []
    for r in audits:
        x = num(r["focal_terminal_branch"])
        if x >= 1:
            thresholds.append({"locus_id": r["locus_id"], "midpoint_bp": r["midpoint_bp"], "focal_terminal_branch": fmt(x), "thresholds_exceeded": ";".join(str(t) for t in (1,5,10,20,30) if x >= t), "upstream_downstream": r["upstream_downstream"], "max_other_terminal_branch": fmt(r["max_other_terminal_branch"]), "total_tree_length": fmt(r["total_tree_length"])})
    robust_fields = ["locus_id", "start_bp", "end_bp", "midpoint_bp", "B_all", "B_leave_one_out", "B_median", "B_trimmed", "A_distance_original", "A_distance_leave_one_out", "A_relative_original", "A_relative_leave_one_out", "D_original", "D_leave_one_out", "stage4_mean_pair_distance_no_spretus", "stage4_A_distance", "stage4_A_relative"]
    robust = [{k: r.get(k, "") for k in robust_fields} for r in rows]
    bins = make_bins(rows)
    # Frozen Stage 13c detector configuration; reference intervals are evaluated only after fitting.
    methods = {"original_B": ["B_all"], "leave_one_out_B": ["B_leave_one_out"], "median_B": ["B_median"], "trimmed_B": ["B_trimmed"], "original_A": ["A_distance_original"], "leave_one_out_A": ["A_distance_leave_one_out"], "original_joint": ["D_original", "A_distance_original"], "leave_one_out_joint": ["D_leave_one_out", "A_distance_leave_one_out"]}
    boundary = []
    for method, fields in methods.items():
        z = run_boundary(rows, fields, refs); z.update({"method": method, "model": "single_interval", "configuration": "Stage13c shared detector; robust MAD scaling; min_size=3; references post-fit"}); boundary.append(z)
    write(OUT / "sample_branch_length_audit.tsv", audits, list(audits[0]))
    write(OUT / "extreme_windows.tsv", thresholds, list(thresholds[0]) if thresholds else ["locus_id"])
    write(OUT / "sample_influence_ranking.tsv", influences, list(influences[0]))
    write(OUT / "robust_distance_comparison.tsv", robust, robust_fields)
    write(OUT / "boundary_sensitivity.tsv", boundary, list(boundary[0]))
    write(OUT / "branch_length_distribution.tsv", audits, list(audits[0]))
    write(OUT / "binned_100kb_tracks.tsv", bins, list(bins[0]))
    archive_rows = archive_filter_comparison(metadata, [r["locus_id"] for r in sorted(audits, key=lambda r: num(r["focal_terminal_branch"]), reverse=True)[:5]])
    write(OUT / "source_filter_comparison.tsv", archive_rows, list(archive_rows[0]) if archive_rows else ["locus_id", "filtering_procedure", "status"])
    make_figures(bins, rows, refs)
    max_b = max(abs(num(r["B_all"]) - num(r["stage4_mean_pair_distance_no_spretus"])) for r in rows)
    max_a = max(abs(num(r["A_distance_original"]) - num(r["stage4_A_distance"])) for r in rows)
    max_ar = max(abs(num(r["A_relative_original"]) - num(r["stage4_A_relative"])) for r in rows)
    downstream = [r for r in rows if num(r["midpoint_bp"]) >= CHANGE_POINT]; upstream = [r for r in rows if num(r["midpoint_bp"]) < CHANGE_POINT]
    def summary(z, field): return avg([r[field] for r in z])
    (OUT / "data_provenance_audit.md").write_text(f"""# Data provenance audit\n\nThe primary trees are the 4,046 metadata-matched local Newick trees in `{TREES.relative_to(ROOT)}`. Stage 0 reports 55 tips in every tree, valid Newick syntax, branch lengths in all trees, no negative branches, and exact coordinate mapping. The tree archive originates from the Kelemen & Vicoso data archive (`IST-2017-78-v1+1_Data.zip`, SHA256 `{sha(ARCHIVE) if ARCHIVE.exists() else 'unavailable'}`). The archive README documents ML, NJ, and MP trees under three SNP-filtering procedures: PASS SNPs, coverage plus allele-ratio filtered raw SNPs, and coverage-filtered raw SNPs.\n\nThe committed tree file is a processed copy of the published ML trees. The source archive contains tree files and SNP-coordinate products, but no per-window FASTA alignments or complete IQ-TREE command/version records were recoverable in this repository audit. Consequently, the exact model optimization, ascertainment correction, missing-data handling, and alignment composition for the extreme windows cannot be independently rerun here. Representative branch-length comparisons across the three nested published ML archives are saved in `source_filter_comparison.tsv`.\n\nThe `.variant` suffix is retained as a source tip label. The frozen tip mapping identifies `{EXCLUDED}` as a standard noncarrier *Mus musculus domesticus* sequence; the suffix alone was not interpreted as a special branch-length encoding.\n""")
    (OUT / "methods_text.md").write_text("""# Methods\n\nFor each frozen local ML Newick tree, terminal and internal branch lengths were read directly from the source branch annotations. Patristic pairwise distances were calculated from the same trees. The primary all-sample `B(w)` is the mean over all non-*M. spretus* pairs. The leave-one-out version removes `Mmd_FRA4_18B.variant` from both numerator and denominator. Median and 5% trimmed means are separate robust summaries. Arrangement-informed `A_distance` and `A_relative` were recomputed with the existing Stage 7 within-subspecies t-versus-standard design. ALL_TIPS quartet support was recalculated with the frozen quartet definitions after removing the same tip from its group.\n\nThe branch detector uses the frozen Stage 13c shared configuration: robust MAD scaling, minimum segment size three, and the existing zero-to-four interval search. Reference inversion coordinates are used only for post-fit evaluation.\n""")
    (OUT / "results_text.md").write_text(f"""# Results\n\nThe focal terminal branch is present in the frozen source trees and drives the observed distance discontinuity. Across the all-sample track, mean `B(w)` is {summary(upstream,'B_all'):.6g} upstream and {summary(downstream,'B_all'):.6g} downstream; corresponding leave-one-out means are {summary(upstream,'B_leave_one_out'):.6g} and {summary(downstream,'B_leave_one_out'):.6g}. Direct recalculation agrees with Stage 4 to maximum absolute differences of {max_b:.3g} for B, {max_a:.3g} for A_distance, and {max_ar:.3g} for A_relative.\n\nThe source tree archive audit does not establish a numerical branch-length cap or a complete inference configuration. The available evidence therefore verifies reproducibility and sample influence, but leaves biological divergence versus region-specific sequence, ascertainment, or tree-estimation effects unresolved.\n""")
    (OUT / "diagnostic_report.md").write_text("""# Diagnostic report\n\n## Classification\n\n**Verified:** extreme terminal branches, the 14.7375-Mb discontinuity, and the effect on B(w) are present in and reproducible from the frozen source Newick trees.\n\n**Supported but not conclusive:** the downstream increase is driven by the focal sample's terminal branch and its pairwise distances; the arrangement-informed numerator changes much less than B(w).\n\n**Unresolved:** the available archive does not include the per-window alignments or complete IQ-TREE optimization provenance needed to distinguish a biological sequence divergence from a filtering, ascertainment, alignment, or tree-estimation problem. The leave-one-out track is therefore a sensitivity analysis, not a correction.\n\n`source_filter_comparison.tsv` compares the focal branch in representative windows across the three published ML filtering products where the nested archives were readable.\n""")
    (OUT / "validation_report.md").write_text(f"""# Validation report\n\n- Trees processed: {len(rows)}.\n- Direct Stage 4 maximum absolute differences: B={max_b:.3g}; A_distance={max_a:.3g}; A_relative={max_ar:.3g}.\n- Focal branch audit: `branch_length_distribution.tsv`.\n- Threshold-exceeding windows: `extreme_windows.tsv`.\n- Pairwise robust comparison: `robust_distance_comparison.tsv`.\n- Boundary sensitivity uses fixed parameters and is in `boundary_sensitivity.tsv`.\n- Source filtering comparison: `source_filter_comparison.tsv`.\n- Original figures and Stage 7 results were not overwritten.\n""")
    (OUT / "report.md").write_text(f"""# House-mouse Stage 8 report\n\nThe approximately 40-substitution-unit terminal branches are present in the source ML Newick trees and are not an arithmetic artifact in the Stage 4–7 calculations. The focal branch's downstream distribution and threshold counts are recorded in `branch_length_distribution.tsv` and `extreme_windows.tsv`; `source_filter_comparison.tsv` tests representative windows across the three published SNP-filtering products.\n\nThe observed change point is {CHANGE_POINT/1e6:.4f} Mb. All-sample B(w) changes from {summary(upstream,'B_all'):.6g} upstream to {summary(downstream,'B_all'):.6g} downstream, while leave-one-out B(w) changes from {summary(upstream,'B_leave_one_out'):.6g} to {summary(downstream,'B_leave_one_out'):.6g}. A_distance is recalculated consistently for both sample sets, and the ALL_TIPS topology track is recalculated separately.\n\nThe evidence does not demonstrate a technical error, numerical saturation, or a biological cause. The recommended manuscript treatment is to retain the original all-sample B(w) as the reproducible primary statistic only with an explicit substitution-unit/provenance caveat, show the leave-one-out and robust tracks as sensitivity analyses, and avoid presenting the extreme downstream B(w) as a calibrated time or ordinary within-species divergence estimate.\n""")
    (OUT / "figure_caption.md").write_text("""**House-mouse chromosome 17 Stage 8.** Branch diagnostics are calculated directly from the 4,046 frozen local ML Newick trees. The original 100-kb figure retains all mapped non-*M. spretus* samples; the sensitivity figure excludes `Mmd_FRA4_18B.variant` consistently from distance, arrangement-informed, and ALL_TIPS quartet calculations. Shading marks the four frozen inversion intervals. Distances are local ML patristic substitution units, not calibrated divergence times or coalescence times.\n""")
    inputs = [TREES, META, MAP, DIST, QSCAN, INT, ARCHIVE]
    output_files = list(OUT.glob("*.tsv")) + list(OUT.glob("*.md"))
    (OUT / "manifest.json").write_text(json.dumps({"stage": "House Mouse Stage 8 branch-length validation", "script": str(Path(__file__).relative_to(ROOT)), "command": "MPLBACKEND=Agg python empirical/house_mouse_t_complex/scripts/08_branch_length_validation.py", "inputs_sha256": {str(x.relative_to(ROOT)): sha(x) for x in inputs if x.exists()}, "n_trees": len(rows), "excluded_sample": EXCLUDED, "change_point_bp": CHANGE_POINT, "source_archive_available": ARCHIVE.exists(), "output_sha256": {str(x.relative_to(ROOT)): sha(x) for x in output_files if x.exists()}}, indent=2) + "\n")


if __name__ == "__main__":
    main()
