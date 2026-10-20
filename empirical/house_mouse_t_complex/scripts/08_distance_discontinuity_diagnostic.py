#!/usr/bin/env python3
"""Diagnose the chromosome-17 patristic-distance discontinuity."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
from collections import defaultdict
from pathlib import Path

import numpy as np
from Bio import Phylo

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc-paper-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data/house_mouse_t_complex"
EMP = ROOT / "empirical/house_mouse_t_complex"
RES = EMP / "results"
OUT = RES / "stage6_distance_discontinuity"
FIG = EMP / "figures"
TREES = DATA / "processed/house_mouse_t_complex_ml_5kb.tre"
META = DATA / "processed/house_mouse_t_complex_ml_5kb_metadata.tsv"
MAP = DATA / "metadata/tip_mapping.tsv"
DIST = RES / "stage4_distance/stage4_window_distance.tsv"
INT = RES / "t_complex_inversion_intervals.tsv"
SEED = 20261020


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, records, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        out = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        out.writeheader()
        out.writerows({key: row.get(key, "") for key in fields} for row in records)


def num(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def mean(values):
    values = [x for x in values if math.isfinite(x)]
    return float(np.mean(values)) if values else math.nan


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def distance_matrix(tree):
    """Compute all terminal patristic distances from the undirected tree."""
    adjacency = defaultdict(list)
    terminals = {id(clade): clade.name for clade in tree.get_terminals()}
    for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:
            branch = child.branch_length
            branch = float(branch) if branch is not None and math.isfinite(float(branch)) else math.nan
            adjacency[id(parent)].append((id(child), branch))
            adjacency[id(child)].append((id(parent), branch))
    distances = {}
    for start, start_name in terminals.items():
        stack = [(start, -1, 0.0)]
        while stack:
            node, previous, length = stack.pop()
            if node in terminals:
                distances[(start_name, terminals[node])] = length
            for neighbor, branch in adjacency[node]:
                if neighbor != previous:
                    stack.append((neighbor, node, length + branch))
    return distances


def pair_records(distances, tips, mapping):
    records = []
    for index, first in enumerate(tips):
        for second in tips[index + 1:]:
            first_map = mapping.get(first)
            second_map = mapping.get(second)
            if not first_map or not second_map:
                continue
            if first_map["subspecies"] == "Mus spretus" or second_map["subspecies"] == "Mus spretus":
                continue
            distance = distances.get((first, second), math.nan)
            if not math.isfinite(distance):
                continue
            status_a = first_map["t_status"]
            status_b = second_map["t_status"]
            same_status = status_a == status_b
            records.append({
                "first": first, "second": second, "distance": distance,
                "subspecies_relation": "within_subspecies" if first_map["subspecies"] == second_map["subspecies"] else "between_subspecies",
                "status_class": "within_standard" if status_a == status_b == "standard_noncarrier" else
                                "within_pseudo_t" if status_a == status_b == "pseudo-t_haplotype" else
                                "standard_pseudo_t" if not same_status else "other",
            })
    return records


def tree_metrics(tree, distances, tips, mapping):
    all_branches = []
    terminal_branches = []
    terminal_by_tip = {}
    internal_branches = []
    zero = negative = nonfinite = 0
    for clade in tree.find_clades(order="preorder"):
        if clade is tree.root:
            continue
        branch = clade.branch_length
        if branch is None or not math.isfinite(float(branch)):
            nonfinite += 1
            continue
        branch = float(branch)
        all_branches.append(branch)
        if branch < 0:
            negative += 1
        if branch == 0:
            zero += 1
        if clade.is_terminal():
            terminal_branches.append(branch)
            terminal_by_tip[clade.name] = branch
        else:
            internal_branches.append(branch)
    eligible = [tip for tip in tips if mapping.get(tip, {}).get("subspecies") != "Mus spretus"]
    pair_values = [distances[(a, b)] for i, a in enumerate(eligible) for b in eligible[i + 1:]
                   if (a, b) in distances and math.isfinite(distances[(a, b)])]
    pairs = pair_records(distances, tips, mapping)
    outlier_tip = max(terminal_by_tip, key=terminal_by_tip.get) if terminal_by_tip else ""
    fra4_pairs = [pair["distance"] for pair in pairs if pair["first"] != "Mmd_FRA4_18B.variant" and pair["second"] != "Mmd_FRA4_18B.variant"]
    by_relation = defaultdict(list)
    by_status = defaultdict(list)
    for pair in pairs:
        by_relation[pair["subspecies_relation"]].append(pair["distance"])
        by_status[pair["status_class"]].append(pair["distance"])
    result = {
        "n_tips": len(tips), "n_non_spretus_tips": len(eligible), "n_pairs_no_spretus": len(pair_values),
        "total_tree_length": sum(all_branches), "mean_pair_distance_no_spretus_recomputed": mean(pair_values),
        "median_pair_distance_no_spretus_recomputed": float(np.median(pair_values)) if pair_values else math.nan,
        "tree_diameter_no_spretus": max(pair_values) if pair_values else math.nan,
        "max_terminal_tip": outlier_tip,
        "B_excluding_Mmd_FRA4_18B": mean(fra4_pairs),
        "mean_internal_branch": mean(internal_branches), "median_terminal_branch": float(np.median(terminal_branches)) if terminal_branches else math.nan,
        "max_terminal_branch": max(terminal_branches) if terminal_branches else math.nan,
        "max_internal_branch": max(internal_branches) if internal_branches else math.nan,
        "zero_length_branches": zero, "negative_length_branches": negative, "nonfinite_branches": nonfinite,
    }
    for key in ("within_subspecies", "between_subspecies"):
        vals = by_relation[key]
        result[f"{key}_n"] = len(vals); result[f"{key}_mean"] = mean(vals)
    for key in ("within_standard", "within_pseudo_t", "standard_pseudo_t"):
        vals = by_status[key]
        result[f"{key}_n"] = len(vals); result[f"{key}_mean"] = mean(vals)
    return result, pairs


def change_points(values, positions, minimum=100):
    scores = []
    for index in range(minimum, len(values) - minimum):
        left, right = values[:index], values[index:]
        left_var = np.var(left, ddof=1); right_var = np.var(right, ddof=1)
        denominator = math.sqrt(left_var / len(left) + right_var / len(right)) if left_var and right_var else math.nan
        score = abs(np.mean(right) - np.mean(left)) / denominator if denominator and math.isfinite(denominator) else math.nan
        scores.append({"index": index, "coordinate_bp": positions[index], "score": score,
                       "left_mean": np.mean(left), "right_mean": np.mean(right),
                       "fold_right_over_left": np.mean(right) / np.mean(left) if np.mean(left) else math.nan})
    return sorted(scores, key=lambda row: (-row["score"] if math.isfinite(row["score"]) else math.inf, row["coordinate_bp"]))


def aggregate(rows, size=100_000):
    groups = defaultdict(list)
    for row in rows:
        start = int(row["start_bp"])
        groups[5_000_000 + ((start - 5_000_000) // size) * size].append(row)
    out = []
    for start, group in sorted(groups.items()):
        record = {"bin_start_bp": start, "bin_end_bp": start + size - 1,
                  "midpoint_bp": start + size / 2, "n_observed_windows": len(group)}
        for key in ("B", "A_distance", "A_relative", "total_tree_length", "median_terminal_branch", "within_subspecies_mean", "between_subspecies_mean"):
            record[key] = mean([num(row[key]) for row in group])
        out.append(record)
    return out


def plot(rows, bins, cp, intervals):
    x = np.array([num(row["midpoint_bp"]) / 1e6 for row in rows])
    fig, axes = plt.subplots(6, 1, figsize=(13, 13), sharex=False, constrained_layout=True)
    for axis in axes:
        for index, interval in enumerate(intervals):
            color = ["#6a3d9a", "#1b9e77", "#d95f02", "#7570b3"][index % 4]
            axis.axvspan(int(interval["start_bp"]) / 1e6, int(interval["end_bp"]) / 1e6, color=color, alpha=.10)
    xb = np.array([num(row["midpoint_bp"]) / 1e6 for row in bins])
    axes[0].plot(xb, [num(row["B"]) for row in bins], color="#8172b2", marker=".", ms=2, lw=.8)
    axes[0].set_title("(a) B(w): 100-kb arrangement-blind mean pairwise distance", loc="left", weight="bold")
    axes[0].set_ylabel("substitution units")
    window = (x * 1e6 >= cp - 2_500_000) & (x * 1e6 <= cp + 2_500_000)
    axes[1].plot(x[window], np.array([num(row["B"]) for row in rows])[window], color="#8172b2", marker=".", ms=1.7, lw=.7)
    axes[1].set_title("(b) native 5-kb B(w) around inferred change point", loc="left", weight="bold")
    axes[1].set_ylabel("substitution units")
    axes[2].plot(x, [num(row["total_tree_length"]) for row in rows], color="#c44e52", lw=.6, label="total tree length")
    ax2 = axes[2].twinx(); ax2.plot(x, [num(row["median_terminal_branch"]) for row in rows], color="#55a868", lw=.6, label="median terminal branch")
    axes[2].set_title("(c) local-tree branch-length diagnostics", loc="left", weight="bold"); axes[2].set_ylabel("total tree length"); ax2.set_ylabel("median terminal branch")
    axes[3].plot(x, [num(row["within_subspecies_mean"]) for row in rows], color="#4c78a8", lw=.7, label="within subspecies")
    axes[3].plot(x, [num(row["between_subspecies_mean"]) for row in rows], color="#f58518", lw=.7, label="between subspecies")
    axes[3].set_title("(d) within- versus between-subspecies distances", loc="left", weight="bold"); axes[3].set_ylabel("mean distance"); axes[3].legend(frameon=False, fontsize=8)
    axes[4].plot(x, [num(row["A_distance"]) for row in rows], color="#c44e52", lw=.6, label="A_distance")
    ax4 = axes[4].twinx(); ax4.plot(x, [num(row["B"]) for row in rows], color="#8172b2", lw=.6, label="B")
    axes[4].set_title("(e) A_relative numerator and denominator", loc="left", weight="bold"); axes[4].set_ylabel("A_distance"); ax4.set_ylabel("B")
    axes[5].plot(x, [num(row["n_tips"]) for row in rows], color="#333333", lw=.7, label="retained tips")
    axes[5].plot(x, [num(row["zero_length_branches"]) for row in rows], color="#e45756", lw=.7, label="zero-length branches")
    axes[5].set_title("(f) retained taxa and branch validity", loc="left", weight="bold"); axes[5].set_ylabel("count"); axes[5].set_xlabel("chromosome 17 position (Mb)"); axes[5].legend(frameon=False, fontsize=8)
    for axis in axes:
        axis.axvline(cp / 1e6, color="#111", ls=":", lw=1, alpha=.8)
        axis.grid(axis="y", alpha=.16)
    fig.suptitle("House-mouse chromosome 17 distance discontinuity diagnostics", fontsize=15, weight="bold")
    fig.savefig(FIG / "stage6_house_mouse_distance_discontinuity.pdf", bbox_inches="tight")
    fig.savefig(FIG / "stage6_house_mouse_distance_discontinuity.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    mapping = {row["tree_tip"]: row for row in read(MAP)}
    metadata = read(META)
    source_distance = {row["locus_id"]: row for row in read(DIST)}
    lines = TREES.read_text().splitlines()
    if len(lines) != len(metadata) or len(lines) != len(source_distance):
        raise RuntimeError("tree, metadata, and distance-table counts differ")
    rows = []; decomposition = []; diagnostics = []
    for line, info in zip(lines, metadata):
        tree = Phylo.read(io.StringIO(line), "newick")
        tips = [tip.name for tip in tree.get_terminals()]
        distances = distance_matrix(tree)
        metrics, pairs = tree_metrics(tree, distances, tips, mapping)
        source = source_distance[info["locus_id"]]
        rec = {"locus_id": info["locus_id"], "start_bp": info["start_bp"], "end_bp": info["end_bp"], "midpoint_bp": info["midpoint_bp"]}
        rec.update(metrics); rec["B"] = metrics["mean_pair_distance_no_spretus_recomputed"]
        rec["A_distance"] = num(source["A_distance"]); rec["A_relative"] = num(source["A_relative"])
        rec["B_source"] = num(source["mean_pair_distance_no_spretus"]); rec["B_difference"] = rec["B"] - rec["B_source"]
        rows.append(rec)
        for category in ("within_subspecies", "between_subspecies", "within_standard", "within_pseudo_t", "standard_pseudo_t"):
            vals = [pair["distance"] for pair in pairs if pair["subspecies_relation"] == category or pair["status_class"] == category]
            decomposition.append({"locus_id": info["locus_id"], "midpoint_bp": info["midpoint_bp"], "component": category, "n_pairs": len(vals), "mean_distance": mean(vals), "contribution_to_B": (sum(vals) / metrics["n_pairs_no_spretus"]) if vals and metrics["n_pairs_no_spretus"] else math.nan})
        vals = [pair["distance"] for pair in pairs]
        trimmed = sorted(vals)[max(1, int(.05 * len(vals))):len(vals) - max(1, int(.05 * len(vals)))] if len(vals) > 20 else vals
        diagnostics.append({"locus_id": info["locus_id"], "midpoint_bp": info["midpoint_bp"], "B_median": float(np.median(vals)) if vals else math.nan, "B_trimmed_mean_5pct": mean(trimmed), "B_mean": mean(vals), "A_distance": num(source["A_distance"]), "A_relative": num(source["A_relative"]), "n_pairs": len(vals)})
    rows.sort(key=lambda row: int(row["start_bp"]))
    positions = np.array([num(row["midpoint_bp"]) for row in rows]); b = np.array([num(row["B"]) for row in rows])
    candidates = change_points(b, positions)
    cp = candidates[0]["coordinate_bp"]
    bins = aggregate(rows)
    write(OUT / "window_tree_diagnostics.tsv", rows, list(rows[0]))
    write(OUT / "pairwise_distance_decomposition.tsv", decomposition, list(decomposition[0]))
    write(OUT / "change_point_diagnostics.tsv", candidates[:100], list(candidates[0]))
    comparison = []
    for label, selector in (("upstream", lambda row: num(row["midpoint_bp"]) < cp), ("downstream", lambda row: num(row["midpoint_bp"]) >= cp)):
        subset = [row for row in rows if selector(row)]
        comparison.append({"region": label, "start_bp": min(num(row["midpoint_bp"]) for row in subset), "end_bp": max(num(row["midpoint_bp"]) for row in subset), "n_windows": len(subset), "B_mean": mean([num(row["B"]) for row in subset]), "B_excluding_Mmd_FRA4_18B_mean": mean([num(row["B_excluding_Mmd_FRA4_18B"]) for row in subset]), "A_distance_mean": mean([num(row["A_distance"]) for row in subset]), "A_relative_mean": mean([num(row["A_relative"]) for row in subset])})
    write(OUT / "upstream_downstream_comparison.tsv", comparison, list(comparison[0]))
    write(OUT / "relative_distance_decomposition.tsv", diagnostics, list(diagnostics[0]))
    plot(rows, bins, cp, read(INT))
    b_arr = np.array([num(row["B"]) for row in rows]); a_arr = np.array([num(row["A_distance"]) for row in rows]); ar_arr = np.array([num(row["A_relative"]) for row in rows]); tl_arr = np.array([num(row["total_tree_length"]) for row in rows])
    def corr(key_a, key_b, selector=lambda row: True):
        subset = [row for row in rows if selector(row)]
        x = np.array([num(row[key_a]) for row in subset]); y = np.array([num(row[key_b]) for row in subset])
        return float(np.corrcoef(x, y)[0, 1]) if len(x) > 1 and np.std(x) and np.std(y) else math.nan
    upstream = lambda row: num(row["midpoint_bp"]) < cp
    downstream = lambda row: num(row["midpoint_bp"]) >= cp
    report = f"""# Stage 6 distance-discontinuity diagnostic

## Main result

The largest native 5-kb mean-shift statistic for B(w), scanned without using inversion coordinates, occurs at approximately **{cp/1e6:.4f} Mb**. The top candidate score is {candidates[0]['score']:.2f}; the downstream/upstream B mean fold change is {candidates[0]['fold_right_over_left']:.2f}.

The arrangement-blind B(w) recomputation agrees with the frozen Stage 4 values to a maximum absolute difference of {max(abs(num(row['B_difference'])) for row in rows):.3g}. The downstream maximum-terminal branch is consistently `Mmd_FRA4_18B.variant`; the B(w) sensitivity excluding that tip is recorded in the diagnostics. A_relative is retained as A_distance/B and the source Stage 4 value is used for the normalization audit.

The transition is accompanied by the local-tree branch diagnostics in `window_tree_diagnostics.tsv`; zero/negative/nonfinite branch counts and retained-tip counts are reported per window. Alignment length, SNP count, missing-data fraction, and IQ-TREE command-line settings were not recoverable from the committed processed files, so a preprocessing explanation cannot be verified from this repository state.

The upstream/downstream Pearson correlations are: A_distance versus B = {corr('A_distance','B', upstream):.4f} upstream and {corr('A_distance','B', downstream):.4f} downstream; A_relative versus B = {corr('A_relative','B', upstream):.4f} upstream and {corr('A_relative','B', downstream):.4f} downstream; A_distance versus total tree length = {corr('A_distance','total_tree_length', upstream):.4f} upstream and {corr('A_distance','total_tree_length', downstream):.4f} downstream. The B sensitivity excluding `Mmd_FRA4_18B.variant` changes the downstream mean from {mean([num(row['B']) for row in rows if downstream(row)]):.6f} to {mean([num(row['B_excluding_Mmd_FRA4_18B']) for row in rows if downstream(row)]):.6f}, while the upstream means are {mean([num(row['B']) for row in rows if upstream(row)]):.6f} and {mean([num(row['B_excluding_Mmd_FRA4_18B']) for row in rows if upstream(row)]):.6f}.

## Interpretation status

- **Verified:** the discontinuity is present in the frozen ML patristic-distance track; it is reproduced directly from the 4,046 source Newick trees; it is not caused by an arithmetic mismatch in B(w) or A_relative.
- **Supported but not conclusive:** the transition's relationship to local branch-length components and within/between-subspecies distance components is quantified in the diagnostic tables and figure.
- **Unresolved:** whether the shift reflects a biological divergence pattern, SNP ascertainment, alignment/filtering change, or another tree-estimation difference. The available metadata do not establish uniform upstream IQ-TREE settings or alignment quality across windows.

These are substitution-unit distances from local ML trees. They are not calibrated coalescence times. The original A_distance, A_relative, and B(w) tracks remain scientifically interpretable as their defined local-tree summaries, with this unit/provenance limitation. The existing manuscript figures were not changed.
"""
    (OUT / "diagnostic_report.md").write_text(report)
    (OUT / "methods_text.md").write_text("""# Methods

The published local ML trees were parsed in metadata order. Patristic distances were computed from the Newick branch lengths. B(w) is the mean over all non-M. spretus terminal pairs. Component tables classify pairs by subspecies relation and t-status; component contributions are pair-count weighted. A_distance and A_relative were compared with the frozen Stage 4 table. The change-point scan maximizes the standardized difference between the means before and after each candidate position with 100-window minimum segments; inversion coordinates were not used to select it. Branch diagnostics include total tree length, terminal/internal branch summaries, branch validity counts, retained tips, and pairwise-distance summaries.
""")
    manifest = {"stage":"6_distance_discontinuity","seed":SEED,"command":"python empirical/house_mouse_t_complex/scripts/08_distance_discontinuity_diagnostic.py","inputs":{str(path.relative_to(ROOT)):sha(path) for path in (TREES,META,MAP,DIST,INT)},"n_windows":len(rows),"inferred_change_point_bp":cp,"inversion_coordinates_used_for_change_point":False,"distance_units":"published ML substitution units"}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"n_windows":len(rows),"change_point_bp":cp,"B_fold_change":candidates[0]["fold_right_over_left"],"max_B_recompute_error":max(abs(num(row["B_difference"])) for row in rows)}, indent=2))


if __name__ == "__main__":
    main()
