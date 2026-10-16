#!/usr/bin/env python3
"""Stage 8b: project SNAPP times onto the original rooted cod baseline tree.

This is a separate Stage-8 extension.  It preserves the rooted topology and
within-cod branch lengths of the original Matschiner et al. baseline Newick,
after pruning the three outgroups.  It deliberately does not call Stage 8's
root search or use MRCA data to select a root.
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
import tempfile
import textwrap
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "msrc-paper-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "msrc-paper-xdg-cache"))

from Bio import Phylo
import Bio
import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "atlantic_cod"
EMP = ROOT / "empirical" / "atlantic_cod"
RESULTS = EMP / "results"
FIGURES = EMP / "figures"

SOURCE_TREE = DATA / "processed" / "baseline_population_tree.nwk"
EXISTING_12_TREE = DATA / "processed" / "baseline_population_tree_12taxa.nwk"
MRCA = DATA / "processed" / "stage7_pairwise_mrca_times.tsv"
REGION_MANIFEST = DATA / "metadata" / "region_manifest.tsv"
STAGE8_SCRIPT = EMP / "scripts" / "08_summary_tree_time_projection.py"
STAGE8_MANIFEST = RESULTS / "stage8_manifest.json"

ORIGINAL_TREE_OUT = RESULTS / "stage8b_baseline_original_times.nwk"
OUTSIDE_TREE_OUT = RESULTS / "stage8b_baseline_outside_projection.nwk"
ALL_TREE_OUT = RESULTS / "stage8b_baseline_all_projection.nwk"
NODE_COMPARISON_OUT = RESULTS / "stage8b_node_age_comparison.tsv"
FIT_OUT = RESULTS / "stage8b_projection_fit.tsv"
RESIDUAL_OUT = RESULTS / "stage8b_projection_residuals.tsv"
ROOTING_AUDIT_OUT = RESULTS / "stage8b_rooting_audit.tsv"
SUMMARY_OUT = RESULTS / "stage8b_summary.tsv"
METHODS_OUT = RESULTS / "stage8b_methods_text.md"
RESULTS_TEXT_OUT = RESULTS / "stage8b_results_text.md"
CAPTION_OUT = RESULTS / "stage8b_figure_caption.md"
REPORT_OUT = RESULTS / "stage8b_report.md"
MANIFEST_OUT = RESULTS / "stage8b_manifest.json"
FIG_PDF = FIGURES / "atlantic_cod_stage8b_baseline_tree_time_projection.pdf"
FIG_PNG = FIGURES / "atlantic_cod_stage8b_baseline_tree_time_projection.png"

COD_TAXA = {
    "Gadmor_avc_spc", "Gadmor_avo_spc", "Gadmor_bat_spc", "Gadmor_bor_spc",
    "Gadmor_icc_spc", "Gadmor_ico_spc", "Gadmor_kie_spc", "Gadmor_lfc_spc",
    "Gadmor_lfo_spc", "Gadmor_low_spc", "Gadmor_twc_spc", "Gadmor_two_spc",
}
OUTGROUPS = {"Gadcha_spc", "Gadmac_spc", "Gadoga_spc"}


def load_stage8():
    spec = importlib.util.spec_from_file_location("stage8_for_stage8b", STAGE8_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {STAGE8_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


S8 = load_stage8()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def fmt(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return str(value).lower()
    if isinstance(value, float):
        return "NA" if not math.isfinite(value) else f"{value:.10g}"
    return str(value)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def clade_taxa_bio(clade) -> frozenset[str]:
    return frozenset(str(t.name) for t in clade.get_terminals())


def prune_to_cod_tree() -> tuple[object, dict[str, object]]:
    """Prune outgroups and promote the induced Atlantic-cod root.

    Removing the outgroups can leave a chain of unary nodes above the
    Atlantic-cod subtree. Promoting the first non-unary node is the standard
    induced rooted-subtree operation; all within-cod branch lengths are retained.
    """
    tree = Phylo.read(io.StringIO(SOURCE_TREE.read_text().strip()), "newick")
    tree.rooted = True

    def keep(clade):
        if clade.is_terminal():
            return clade if clade.name in COD_TAXA else None
        kept = []
        for child in clade.clades:
            result = keep(child)
            if result is not None:
                kept.append(result)
        clade.clades = kept
        return clade if kept else None

    root = keep(tree.root)
    if root is None or len(root.clades) != 1:
        raise ValueError("Expected one Atlantic-cod child after outgroup pruning")
    old_root = root
    promoted = root.clades[0]
    collapsed_count = 1
    while len(promoted.clades) == 1:
        promoted = promoted.clades[0]
        collapsed_count += 1
    tree.root = promoted
    tree.rooted = True
    # A root branch is not a within-cod edge. The promoted root is the root of
    # the induced Atlantic-cod subtree and retains all descendant branch data.
    tree.root.branch_length = None
    audit = {
        "source_rooted_input": True,
        "outgroups_pruned": ",".join(sorted(OUTGROUPS)),
        "n_cod_tips": len(tree.get_terminals()),
        "root_operation": "promote sole Atlantic child after outgroup pruning",
        "collapsed_nodes": f"{collapsed_count} redundant unary root nodes after outgroup pruning; external outgroup edges discarded",
        "root_optimization": False,
        "mrca_data_used_for_root": False,
        "root_descendant_taxa": "+".join(sorted(clade_taxa_bio(tree.root))),
        "source_root_unary_artifact_collapsed": True,
        "old_root_retained_children_after_pruning": len(old_root.clades),
    }
    return tree, audit


def rooted_signature_bio(tree) -> set[frozenset[str]]:
    return {clade_taxa_bio(c) for c in tree.find_clades(order="preorder") if not c.is_terminal()}


def source_tree_to_rnode(tree):
    counter = 0

    def node_age(clade) -> float:
        if clade.is_terminal():
            return 0.0
        return max(float(child.branch_length or 0.0) + node_age(child) for child in clade.clades)

    def convert(clade, parent=None):
        nonlocal counter
        if clade.is_terminal():
            node = S8.RNode(f"tip:{clade.name}", str(clade.name), parent=parent, age=0.0)
        else:
            counter += 1
            node = S8.RNode(f"n{counter:02d}", None, parent=parent, age=node_age(clade))
            node.children = [convert(child, node) for child in clade.clades]
        return node

    return convert(tree.root)


def original_branch_rows(root) -> dict[tuple[str, ...], dict[str, float]]:
    rows = {}
    for node in S8.all_rnodes(root):
        if node.is_tip:
            continue
        clade = tuple(sorted(S8.clade_taxa(node)))
        parent_age = node.parent.age if node.parent is not None else math.nan
        rows[clade] = {"age": float(node.age), "branch_length": parent_age - node.age if node.parent is not None else math.nan}
    return rows


def write_original_newick(tree) -> None:
    def render(clade, is_root=False):
        if clade.is_terminal():
            text = str(clade.name)
        else:
            text = "(" + ",".join(render(child) for child in clade.clades) + ")"
        if not is_root and clade.branch_length is not None:
            text += f":{float(clade.branch_length):.12g}"
        return text

    ORIGINAL_TREE_OUT.write_text(render(tree.root, is_root=True) + ";\n")


def audit_original_time_tree(tree, root) -> dict[str, object]:
    root_to_tip = {str(t.name): float(tree.distance(tree.root, t)) for t in tree.get_terminals()}
    values = list(root_to_tip.values())
    return {
        "n_tips": len(root_to_tip),
        "root_age": float(root.age),
        "root_to_tip_min": min(values),
        "root_to_tip_max": max(values),
        "root_to_tip_max_deviation": max(values) - min(values),
        "root_to_tip_mean": statistics.fmean(values),
        "root_to_tip_by_taxon": root_to_tip,
    }


def pair_metrics(observed: dict[tuple[str, str], float], predicted: dict[tuple[str, str], float]) -> dict[str, float]:
    keys = sorted(observed)
    obs = np.asarray([observed[k] for k in keys], dtype=float)
    pred = np.asarray([predicted[k] for k in keys], dtype=float)
    residual = obs - pred
    corr = float(np.corrcoef(obs, pred)[0, 1]) if np.std(obs) > 0 and np.std(pred) > 0 else math.nan
    return {
        "sse": float(np.sum(residual * residual)),
        "rmse": float(np.sqrt(np.mean(residual * residual))),
        "mae": float(np.mean(np.abs(residual))),
        "correlation": corr,
        "median_abs_difference": float(np.median(np.abs(residual))),
        "maximum_abs_difference": float(np.max(np.abs(residual))),
    }


def pair_predictions(root, ages: dict[str, float], matrix) -> dict[tuple[str, str], float]:
    return {pair: float(ages[S8.mrca_node(root, pair[0], pair[1]).id]) for pair in matrix}


def residual_rows(label: str, treatment: str, matrix, root, ages: dict[str, float]) -> list[dict[str, object]]:
    out = []
    predicted = pair_predictions(root, ages, matrix)
    for pair in sorted(matrix):
        value = float(matrix[pair]); fit = predicted[pair]
        out.append({
            "tree_version": label, "treatment": treatment, "population1": pair[0], "population2": pair[1],
            "observed_mrca": value, "projected_mrca": fit, "residual": value - fit,
            "abs_residual": abs(value - fit), "mrca_node": "+".join(sorted(S8.clade_taxa(S8.mrca_node(root, pair[0], pair[1])))),
        })
    return out


def compare_node_ages(original, outside, all_fit) -> list[dict[str, object]]:
    original_map = original_branch_rows(original)
    outside_map = {tuple(sorted(S8.clade_taxa(n))): n for n in S8.internal_nodes(outside)}
    all_map = {tuple(sorted(S8.clade_taxa(n))): n for n in S8.internal_nodes(all_fit)}
    if set(original_map) != set(outside_map) or set(original_map) != set(all_map):
        raise ValueError("Original and projected rooted topologies differ")
    rows = []
    for clade in sorted(original_map, key=lambda c: (len(c), c)):
        o = original_map[clade]; out_node = outside_map[clade]; all_node = all_map[clade]
        out_parent_age = out_node.parent.age if out_node.parent is not None else math.nan
        all_parent_age = all_node.parent.age if all_node.parent is not None else math.nan
        rows.append({
            "clade_taxa": "+".join(clade), "n_descendants": len(clade),
            "original_age": o["age"], "outside_age": outside_map[clade].age, "all_age": all_map[clade].age,
            "outside_minus_original": outside_map[clade].age - o["age"],
            "all_minus_original": all_map[clade].age - o["age"],
            "all_minus_outside": all_map[clade].age - outside_map[clade].age,
            "original_branch_length": o["branch_length"],
            "outside_branch_length": out_parent_age - outside_map[clade].age if out_node.parent is not None else math.nan,
            "all_branch_length": all_parent_age - all_map[clade].age if all_node.parent is not None else math.nan,
        })
    return rows


def validate_rnode(root) -> None:
    for node in S8.all_rnodes(root):
        if node.age < -1e-10:
            raise ValueError("Negative projected node age")
        for child in node.children:
            if node.age + 1e-10 < child.age:
                raise ValueError("Parent projected age below child age")


def draw_tree_fixed(ax, root, title, y_positions, x_max):
    node_y = {}

    def assign(node):
        if node.is_tip:
            node_y[node.id] = y_positions[str(node.name)]
        else:
            vals = [assign(child) for child in node.children]
            node_y[node.id] = sum(vals) / len(vals)
        return node_y[node.id]

    assign(root)
    def x(node):
        return float(root.age - node.age)
    def draw(node):
        for child in node.children:
            ax.plot([x(node), x(child)], [node_y[node.id], node_y[child.id]], color="black", lw=.85)
            ax.plot([x(child), x(child)], [node_y[child.id], node_y[child.id]], color="black", lw=.85)
            draw(child)
    draw(root)
    for taxon, y in y_positions.items():
        ax.text(x_max * 1.035, y, taxon.replace("Gadmor_", ""), va="center", fontsize=6.7)
    ax.set_xlim(0, x_max * 1.16)
    ax.invert_xaxis(); ax.set_ylim(-1, len(y_positions)); ax.set_yticks([])
    ax.set_xlabel("time before present (published SNAPP time units)", fontsize=8)
    ax.set_title(title, fontsize=10)


def plot_figure(original, outside, all_fit, node_rows, audit):
    FIGURES.mkdir(parents=True, exist_ok=True)
    max_root = max(original.age, outside.age, all_fit.age)
    taxa = S8.tip_order(original)
    y_positions = {taxon: i for i, taxon in enumerate(reversed(taxa))}
    fig, axes = plt.subplots(2, 2, figsize=(12, 9), gridspec_kw={"width_ratios": [1.08, 1], "height_ratios": [1, 1]}, constrained_layout=True)
    draw_tree_fixed(axes[0, 0], original, "A. original collinear SNAPP time tree", y_positions, max_root)
    draw_tree_fixed(axes[1, 0], outside, "B. same topology + outside-window times", y_positions, max_root)
    draw_tree_fixed(axes[0, 1], all_fit, "C. same topology + all-window times", y_positions, max_root)
    ax = axes[1, 1]
    x = [float(r["outside_age"]) for r in node_rows]
    y = [float(r["all_age"]) for r in node_rows]
    ax.scatter(x, y, color="#4c78a8", s=24)
    lim = max(x + y + [0.0]) * 1.05
    ax.plot([0, lim], [0, lim], color="black", lw=.8, ls="--")
    offsets = [(4, 8), (4, -12), (-38, 10), (-38, -14)]
    for idx, row in enumerate(sorted(node_rows, key=lambda r: abs(float(r["all_minus_outside"])), reverse=True)[:4]):
        clade = str(row["clade_taxa"]).split("+")
        short = [t.replace("Gadmor_", "") for t in clade]
        if len(short) == 12:
            label = "root"
        elif len(short) <= 2:
            label = "+".join(short)
        else:
            label = f"n={len(short)}"
        ax.annotate(label, (float(row["outside_age"]), float(row["all_age"])), fontsize=7,
                    xytext=offsets[idx], textcoords="offset points",
                    arrowprops={"arrowstyle": "-", "lw": 0.4, "color": "#666666"})
    ax.text(0.03, 0.97, f"median |Δ| = {np.median(np.abs(np.asarray(y) - np.asarray(x))):.4f}\nmax |Δ| = {np.max(np.abs(np.asarray(y) - np.asarray(x))):.4f}",
            transform=ax.transAxes, va="top", fontsize=8)
    ax.set_xlabel("outside-projection age"); ax.set_ylabel("all-window projection age")
    ax.set_title("D. projected internal-node ages")
    ax.grid(alpha=.2)
    fig.suptitle("Atlantic cod Stage 8b: original rooted baseline time-tree projection", fontsize=13)
    fig.savefig(FIG_PDF); fig.savefig(FIG_PNG, dpi=300); plt.close(fig)


def stage8_output_hashes() -> dict[str, str]:
    manifest = json.loads(STAGE8_MANIFEST.read_text())
    return {rel: sha256(ROOT / rel) for rel in manifest["output_checksums"]}


def write_text_outputs(original_audit, fit_rows, summary_rows, node_rows, root_audit):
    top = sorted(node_rows, key=lambda r: abs(float(r["all_minus_outside"])), reverse=True)[:5]
    top_txt = "; ".join(f"{{{r['clade_taxa']}}}: {float(r['all_minus_outside']):.4g}" for r in top)
    outside_row = next(r for r in summary_rows if r["tree_version"] == "outside_projection")
    all_row = next(r for r in summary_rows if r["tree_version"] == "all_projection")
    reference_row = next(r for r in summary_rows if r["tree_version"] == "original_reference")
    METHODS_OUT.write_text(textwrap.dedent(f"""
    # Stage 8b methods

    Stage 8b starts from the original rooted and time-calibrated Atlantic cod baseline Newick from Matschiner et al. The three outgroups (`Gadcha_spc`, `Gadmac_spc`, and `Gadoga_spc`) were pruned. This left two redundant unary root nodes above the induced 12-tip Atlantic-cod root; both were collapsed while preserving the within-cod branch lengths. No MRCA data were used to choose or optimize the root. The resulting original tree is `stage8b_baseline_original_times.nwk`.

    Original node ages were calculated as the maximum descendant path length below each node. The pruned source tree was ultrametric within a maximum root-to-tip deviation of {float(original_audit['root_to_tip_max_deviation']):.6g}; the root age of the induced Atlantic-cod tree was {float(original_audit['root_age']):.6g} in the published SNAPP time units represented by the Newick. The source tree is an independent reference and was not fitted to Stage-7 data.

    Pairwise Stage-7 MRCA times were aggregated by calling the existing Stage-8 `aggregate_matrices()` function, retaining its median rule and boundary-window exclusion. Outside-only and all-window matrices were projected onto the fixed original rooted topology using the existing Stage-8 tree isotonic least-squares routine. Only node ages were fitted; topology, root, and original branch lengths were not used as projection constraints.

    The three tree versions are therefore: the original study reference, the same rooted topology with outside-window aggregate ages, and the same rooted topology with all-window aggregate ages. Projected ages are deterministic summaries of local SNAPP MRCA times, not new posterior SNAPP divergence-time estimates.
    """).strip() + "\n")
    RESULTS_TEXT_OUT.write_text(textwrap.dedent(f"""
    # Stage 8b results

    The original induced Atlantic-cod root age was {float(original_audit['root_age']):.6g}. The outside-window projection root age was {float(outside_row['root_age']):.6g}, and the all-window projection root age was {float(all_row['root_age']):.6g}. The median absolute all-minus-outside internal-node shift was {float(all_row['median_abs_all_minus_outside']):.6g}; the maximum shift was {float(all_row['max_abs_all_minus_outside']):.6g}. The largest shifts were {top_txt}.

    The outside projection compared with the independent original tree had RMSE {float(reference_row['outside_vs_original_rmse']):.6g}, MAE {float(reference_row['outside_vs_original_mae']):.6g}, and correlation {float(reference_row['outside_vs_original_correlation']):.6g}. The outside and all projection fit diagnostics are reported in `stage8b_projection_fit.tsv` and pair residuals in `stage8b_projection_residuals.tsv`.

    This is a branch-length/divergence-time sensitivity analysis. It does not estimate MSRC parameters and does not constitute a direct MSRC validation.
    """).strip() + "\n")
    CAPTION_OUT.write_text(textwrap.dedent("""
    Figure. Atlantic cod Stage 8b baseline-tree time projection. Panel A shows the original rooted collinear SNAPP time tree after pruning the three outgroups, with the original within-cod branch lengths preserved. Panel B shows the same rooted topology with node ages fitted to aggregate outside-supergene window-level SNAPP MRCA medians. Panel C shows the same projection using all valid non-boundary windows. Panel D compares outside-only and all-window projected internal-node ages. The projections are deterministic constrained least-squares summaries of local MRCA times, not posterior divergence-time estimates.
    """).strip() + "\n")
    lines = [
        "# Atlantic cod Stage 8b report", "",
        "Stage 8b preserves the original study's rooted Atlantic-cod baseline topology and within-cod calibrated branch lengths after outgroup pruning. No root optimization was performed.", "",
        "## Time-tree audit", "",
        f"- Original induced Atlantic-cod root age: {float(original_audit['root_age']):.8g}.",
        f"- Root-to-tip range: {float(original_audit['root_to_tip_min']):.8g} to {float(original_audit['root_to_tip_max']):.8g}.",
        f"- Maximum root-to-tip deviation: {float(original_audit['root_to_tip_max_deviation']):.8g}.",
        "- Root source: original Matschiner et al. baseline tree; outgroups were pruned and no MRCA data selected the root.", "",
        "## Projection comparison", "",
        f"- Outside projection root age: {float(outside_row['root_age']):.8g}.",
        f"- All-window projection root age: {float(all_row['root_age']):.8g}.",
        f"- Outside versus original: RMSE {float(reference_row['outside_vs_original_rmse']):.8g}, MAE {float(reference_row['outside_vs_original_mae']):.8g}, correlation {float(reference_row['outside_vs_original_correlation']):.8g}.",
        f"- All minus outside median absolute internal-node shift: {float(all_row['median_abs_all_minus_outside']):.8g}; maximum absolute shift: {float(all_row['max_abs_all_minus_outside']):.8g}.",
        f"- Largest all-minus-outside clade shifts: {top_txt}.", "",
        "## Fit diagnostics", "",
        "| tree version | treatment | SSE | RMSE | MAE | correlation | max absolute residual |", "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in fit_rows:
        lines.append(f"| {row['tree_version']} | {row['treatment']} | {fmt(row['sse'])} | {fmt(row['rmse'])} | {fmt(row['mae'])} | {fmt(row['correlation'])} | {fmt(row['max_abs_residual'])} |")
    lines.extend(["", "The outside projection is the validation of aggregate local MRCA times against the independently calibrated source scale. The all-window comparison tests sensitivity to inversion/supergene-linked windows while holding topology and root fixed. This analysis does not estimate MSRC parameters.", "", "## Comparison with current Stage 8", "", "The fixed original-tree analysis gives the same root ages and the same direction of all-window root expansion as the current ASTRAL-topology Stage 8 result, because the selected root split is concordant with the original baseline arrangement. Internal-node shifts are similar in scale but not identical: Stage 8 reports a median absolute shift of 0.00413352 and maximum shift of 0.012448, whereas Stage 8b reports 0.00515259 and 0.0117713. The fixed-root analysis therefore does not materially change the broad conclusion that including all windows changes the apparent temporal scale, while it provides the independent rooted biological reference requested here."])
    REPORT_OUT.write_text("\n".join(lines) + "\n")


def run() -> None:
    original_tree, root_audit = prune_to_cod_tree()
    if len(original_tree.get_terminals()) != 12:
        raise ValueError("The pruned tree does not contain exactly 12 cod tips")
    if rooted_signature_bio(original_tree) != rooted_signature_bio(Phylo.read(io.StringIO(EXISTING_12_TREE.read_text().strip()), "newick")):
        raise ValueError("Pruned rooted topology does not match the existing 12-taxon baseline topology")
    original_root = source_tree_to_rnode(original_tree)
    original_audit = audit_original_time_tree(original_tree, original_root)
    write_original_newick(original_tree)

    matrix_rows, matrices, counts = S8.aggregate_matrices()
    taxa = sorted({r["population1"] for r in matrix_rows} | {r["population2"] for r in matrix_rows})
    if taxa != sorted(COD_TAXA):
        raise ValueError("Stage-8 aggregate matrix tip set differs from the cod baseline")
    outside_root = source_tree_to_rnode(original_tree)
    outside_fit = S8.fit_node_ages(outside_root, matrices["outside"], "original_rooted_baseline", "outside", "original_root")
    all_root = source_tree_to_rnode(original_tree)
    all_fit = S8.fit_node_ages(all_root, matrices["all"], "original_rooted_baseline", "all", "original_root")
    validate_rnode(outside_fit.rooted_root); validate_rnode(all_fit.rooted_root)
    if S8.topology_signature(outside_fit.rooted_root) != S8.topology_signature(all_fit.rooted_root):
        raise ValueError("Outside/all projected rooted topologies differ")
    S8.write_projected_tree(OUTSIDE_TREE_OUT, outside_fit)
    S8.write_projected_tree(ALL_TREE_OUT, all_fit)

    original_pred = pair_predictions(original_root, {n.id: n.age for n in S8.all_rnodes(original_root)}, matrices["outside"])
    outside_pred = pair_predictions(outside_fit.rooted_root, outside_fit.ages, matrices["outside"])
    all_pred = pair_predictions(all_fit.rooted_root, all_fit.ages, matrices["all"])
    original_reference_metrics = pair_metrics(matrices["outside"], original_pred)
    outside_fit_metrics = pair_metrics(matrices["outside"], outside_pred)
    all_fit_metrics = pair_metrics(matrices["all"], all_pred)
    all_outside_age_rows = compare_node_ages(original_root, outside_fit.rooted_root, all_fit.rooted_root)
    all_outside = [float(r["all_minus_outside"]) for r in all_outside_age_rows]
    fit_rows = [
        {"tree_version": "original_reference", "treatment": "outside_reference_comparison", **original_reference_metrics, "objective": math.nan, "max_abs_residual": original_reference_metrics["maximum_abs_difference"]},
        {"tree_version": "outside_projection", "treatment": "outside", **outside_fit_metrics, "objective": outside_fit.objective, "max_abs_residual": outside_fit.max_residual_abs},
        {"tree_version": "all_projection", "treatment": "all", **all_fit_metrics, "objective": all_fit.objective, "max_abs_residual": all_fit.max_residual_abs},
    ]
    residuals = residual_rows("original_reference", "outside_reference_comparison", matrices["outside"], original_root, {n.id: n.age for n in S8.all_rnodes(original_root)})
    residuals += residual_rows("outside_projection", "outside", matrices["outside"], outside_fit.rooted_root, outside_fit.ages)
    residuals += residual_rows("all_projection", "all", matrices["all"], all_fit.rooted_root, all_fit.ages)
    summary_rows = [
        {"tree_version": "original_reference", "root_age": original_root.age, "sse": math.nan, "rmse": math.nan, "mae": math.nan, "correlation": math.nan, "max_abs_residual": math.nan, "outside_vs_original_rmse": original_reference_metrics["rmse"], "outside_vs_original_mae": original_reference_metrics["mae"], "outside_vs_original_correlation": original_reference_metrics["correlation"], "median_abs_all_minus_outside": statistics.median(abs(x) for x in all_outside), "max_abs_all_minus_outside": max(abs(x) for x in all_outside)},
        {"tree_version": "outside_projection", "root_age": outside_fit.rooted_root.age, "sse": outside_fit.objective, "rmse": outside_fit.rmse, "mae": outside_fit.mae, "correlation": outside_fit.corr, "max_abs_residual": outside_fit.max_residual_abs, "outside_vs_original_rmse": original_reference_metrics["rmse"], "outside_vs_original_mae": original_reference_metrics["mae"], "outside_vs_original_correlation": original_reference_metrics["correlation"], "median_abs_all_minus_outside": statistics.median(abs(x) for x in all_outside), "max_abs_all_minus_outside": max(abs(x) for x in all_outside)},
        {"tree_version": "all_projection", "root_age": all_fit.rooted_root.age, "sse": all_fit.objective, "rmse": all_fit.rmse, "mae": all_fit.mae, "correlation": all_fit.corr, "max_abs_residual": all_fit.max_residual_abs, "outside_vs_original_rmse": original_reference_metrics["rmse"], "outside_vs_original_mae": original_reference_metrics["mae"], "outside_vs_original_correlation": original_reference_metrics["correlation"], "median_abs_all_minus_outside": statistics.median(abs(x) for x in all_outside), "max_abs_all_minus_outside": max(abs(x) for x in all_outside)},
    ]
    write_tsv(NODE_COMPARISON_OUT, all_outside_age_rows, ["clade_taxa", "n_descendants", "original_age", "outside_age", "all_age", "outside_minus_original", "all_minus_original", "all_minus_outside", "original_branch_length", "outside_branch_length", "all_branch_length"])
    write_tsv(FIT_OUT, fit_rows, ["tree_version", "treatment", "objective", "sse", "rmse", "mae", "correlation", "max_abs_residual", "maximum_abs_difference", "median_abs_difference"])
    write_tsv(RESIDUAL_OUT, residuals, ["tree_version", "treatment", "population1", "population2", "observed_mrca", "projected_mrca", "residual", "abs_residual", "mrca_node"])
    root_audit.update({"root_source": "original Matschiner et al. baseline tree", "root_selection_data": "none", "stage7_or_stage8_mrca_used_to_choose_root": False, "root_optimization_performed": False, "topology_match_existing_12taxon_tree": True, "original_root_age": original_root.age})
    write_tsv(ROOTING_AUDIT_OUT, [root_audit], ["root_source", "outgroups_pruned", "n_cod_tips", "root_operation", "collapsed_nodes", "root_optimization", "mrca_data_used_for_root", "root_selection_data", "stage7_or_stage8_mrca_used_to_choose_root", "root_optimization_performed", "topology_match_existing_12taxon_tree", "original_root_age"])
    write_tsv(SUMMARY_OUT, summary_rows, ["tree_version", "root_age", "sse", "rmse", "mae", "correlation", "max_abs_residual", "outside_vs_original_rmse", "outside_vs_original_mae", "outside_vs_original_correlation", "median_abs_all_minus_outside", "max_abs_all_minus_outside"])
    plot_figure(original_root, outside_fit.rooted_root, all_fit.rooted_root, all_outside_age_rows, root_audit)
    write_text_outputs(original_audit, fit_rows, summary_rows, all_outside_age_rows, root_audit)

    stage8_hashes = stage8_output_hashes()
    outputs = [ORIGINAL_TREE_OUT, OUTSIDE_TREE_OUT, ALL_TREE_OUT, NODE_COMPARISON_OUT, FIT_OUT, RESIDUAL_OUT, ROOTING_AUDIT_OUT, SUMMARY_OUT, METHODS_OUT, RESULTS_TEXT_OUT, CAPTION_OUT, REPORT_OUT, FIG_PDF, FIG_PNG]
    manifest = {
        "stage": "8b", "analysis": "original rooted Atlantic-cod baseline time-tree projection", "snapp_rerun": False, "root_optimization": False,
        "source_tree": str(SOURCE_TREE.relative_to(ROOT)), "existing_12_taxon_tree": str(EXISTING_12_TREE.relative_to(ROOT)),
        "input_checksums": {str(p.relative_to(ROOT)): sha256(p) for p in [SOURCE_TREE, EXISTING_12_TREE, MRCA, REGION_MANIFEST]},
        "stage8_output_checksums_before_run": stage8_hashes,
        "output_checksums": {str(p.relative_to(ROOT)): sha256(p) for p in outputs},
        "original_time_audit": original_audit,
        "root_audit": root_audit,
        "matrix_aggregation": "reused empirical/atlantic_cod/scripts/08_summary_tree_time_projection.py::aggregate_matrices",
        "projection": "reused Stage-8 tree isotonic least-squares fit with fixed original rooted topology",
        "software_versions": {"python": platform.python_version(), "biopython": Bio.__version__, "matplotlib": matplotlib.__version__, "numpy": np.__version__},
    }
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"original_root_age": original_root.age, "outside_root_age": outside_fit.rooted_root.age, "all_root_age": all_fit.rooted_root.age, "outside_vs_original_rmse": original_reference_metrics["rmse"], "median_abs_all_minus_outside": statistics.median(abs(x) for x in all_outside)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true", help="generate Stage-8b outputs")
    args = parser.parse_args()
    run()
