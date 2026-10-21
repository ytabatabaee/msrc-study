#!/usr/bin/env python3
"""Stage 8c: run the published TCMM implementation on Atlantic-cod trees.

This stage is deliberately independent of Stage 8b.  It uses the original
TCMM command-line implementation from the pinned external repository and
preserves its raw output.  Published 250-kb window MCC-like trees are the
reference trees; per-window posterior samples were not available locally.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean, median

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "atlantic_cod"
EMP = ROOT / "empirical" / "atlantic_cod"
RESULTS = EMP / "results" / "stage8c_tcmm"
RAW = RESULTS / "raw"
FIGURES = EMP / "figures"

BASELINE = DATA / "processed" / "baseline_population_tree_12taxa.nwk"
QUERY_TREE = RAW / "tcmm_query_baseline_12taxa_root_stem_removed.nwk"
WINDOWS = DATA / "processed" / "cod_window_trees.tsv"
PAIRWISE = DATA / "processed" / "stage8_pairwise_mrca_matrices.tsv"
REGIONS = DATA / "metadata" / "region_manifest.tsv"
STAGE8B_NODE = EMP / "results" / "stage8b_node_age_comparison.tsv"
STAGE8B_FIT = EMP / "results" / "stage8b_projection_fit.tsv"

TAXA = {
    "Gadmor_avc_spc", "Gadmor_avo_spc", "Gadmor_bat_spc", "Gadmor_bor_spc",
    "Gadmor_icc_spc", "Gadmor_ico_spc", "Gadmor_kie_spc", "Gadmor_lfc_spc",
    "Gadmor_lfo_spc", "Gadmor_low_spc", "Gadmor_twc_spc", "Gadmor_two_spc",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def parse_tree(newick: str) -> Tree:
    return Phylo.read(__import__("io").StringIO(newick.strip()), "newick")


def terminals(tree: Tree) -> set[str]:
    return {str(c.name) for c in tree.get_terminals()}


def descendants(clade: Clade) -> tuple[str, ...]:
    return tuple(sorted(str(c.name) for c in clade.get_terminals()))


def topology_keys(tree: Tree) -> set[tuple[str, ...]]:
    return {descendants(c) for c in tree.get_nonterminals()}


def depths(tree: Tree) -> dict[Clade, float]:
    return {c: float(v or 0.0) for c, v in tree.depths().items()}


def audit_tree(tree: Tree) -> dict[str, object]:
    d = depths(tree)
    tips = tree.get_terminals()
    root_height = max((d[c] for c in tips), default=0.0)
    tip_values = [d[c] for c in tips]
    edge_values = [float(c.branch_length) for c in tree.find_clades() if c is not tree.root and c.branch_length is not None]
    negatives = [x for x in edge_values if x < 0]
    return {
        "n_taxa": len(tips),
        "taxa": terminals(tree),
        "n_edges": len(edge_values),
        "n_negative_edges": len(negatives),
        "minimum_edge": min(edge_values) if edge_values else 0.0,
        "finite_edges": all(math.isfinite(x) for x in edge_values),
        "root_height_max": root_height,
        "root_to_tip_min": min(tip_values) if tip_values else 0.0,
        "root_to_tip_max": max(tip_values) if tip_values else 0.0,
        "root_to_tip_deviation": (max(tip_values) - min(tip_values)) if tip_values else 0.0,
        "valid_nonnegative": all(x >= -1e-10 for x in edge_values),
    }


def node_rows(tree: Tree, treatment: str, reference_id: str, lam: float) -> list[dict[str, object]]:
    d = depths(tree)
    root_height = max((d[c] for c in tree.get_terminals()), default=0.0)
    rows = []
    for c in tree.get_nonterminals():
        key = "+".join(descendants(c))
        age = root_height - d[c]
        rows.append({
            "treatment": treatment,
            "lambda": lam,
            "reference_id": reference_id,
            "clade_taxa": key,
            "n_descendants": len(c.get_terminals()),
            "age_metric_from_root_height": age,
            "branch_length": float(c.branch_length or 0.0),
        })
    return rows


def pairwise_values(tree: Tree) -> dict[tuple[str, str], float]:
    d = depths(tree)
    root_height = max((d[c] for c in tree.get_terminals()), default=0.0)
    result = {}
    names = sorted(terminals(tree))
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            m = tree.common_ancestor(a, b)
            result[(a, b)] = root_height - d[m]
    return result


def patristic_values(tree: Tree) -> dict[tuple[str, str], float]:
    names = sorted(terminals(tree))
    return {(a, b): float(tree.distance(a, b)) for i, a in enumerate(names) for b in names[i + 1:]}


def load_baseline() -> Tree:
    return Phylo.read(str(BASELINE), "newick")


def prepare_query_tree() -> Path:
    """Preserve the induced rooted topology while removing the source root stem.

    The 12-taxon source Newick retains a unary external stem inherited from the
    outgroup-containing tree. Stage 8b's validated induced baseline removes
    that stem. This is a topology-preserving input normalization, not rerooting
    and not an MRCA-driven operation.
    """
    tree = load_baseline()
    tree.root.branch_length = 0.0
    QUERY_TREE.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(QUERY_TREE), "newick")
    return QUERY_TREE


def extract_references() -> tuple[list[dict[str, str]], dict[str, Path]]:
    rows = read_tsv(WINDOWS)
    outdir = RAW / "reference_trees"
    outdir.mkdir(parents=True, exist_ok=True)
    refs: list[dict[str, str]] = []
    paths: dict[str, Path] = {}
    for treatment in ("outside", "all"):
        selected = [r for r in rows if treatment == "all" or r["inside_inversion"].lower() == "false"]
        path = outdir / f"reference_{treatment}.trees"
        with path.open("w") as f:
            for r in selected:
                tree = parse_tree(r["tree_newick"])
                audit = audit_tree(tree)
                valid = audit["taxa"] == TAXA and bool(audit["finite_edges"])
                reason = "valid" if valid else "taxon_or_finite_length_mismatch"
                if valid:
                    f.write(r["tree_newick"].strip() + "\n")
                refs.append({
                    "treatment": treatment,
                    "window_id": r["window_id"],
                    "lg": r["lg"],
                    "start": r["start"],
                    "end": r["end"],
                    "inside_inversion": r["inside_inversion"],
                    "source_file": r["source_file"],
                    "tree_type": r["tree_type"],
                    "valid_for_tcmm": str(valid).lower(),
                    "validation_reason": reason,
                })
        paths[treatment] = path
    write_tsv(RESULTS / "tcmm_reference_tree_manifest.tsv", refs, list(refs[0]))
    return refs, paths


def run_tcmm(tcmm_repo: Path, python: str, query: Path, references: Path, treatment: str, lam: float) -> dict[str, Path | float]:
    tag = str(lam).replace(".", "p")
    outbase = RAW / "tcmm_runs" / f"{treatment}_lambda_{tag}"
    outbase.parent.mkdir(parents=True, exist_ok=True)
    existing = {"treatment": treatment, "lambda": lam, "base": outbase, "trees": Path(str(outbase) + ".trees"), "branches": Path(str(outbase) + ".branches"), "rates": Path(str(outbase) + ".rates"), "log": Path(str(outbase) + ".log")}
    if all(Path(existing[key]).exists() for key in ("trees", "branches", "rates", "log")):
        return existing
    cmd = [python, str(tcmm_repo / "multiple_tree_matching.py"), "-i", str(query), "-r", str(references), "-o", str(outbase), "-l", str(lam), "-g", str(outbase) + ".log"]
    completed = subprocess.run(cmd, capture_output=True, text=True, check=True)
    (outbase.with_suffix(".command.txt")).write_text(" ".join(cmd) + "\n\nSTDOUT\n" + completed.stdout + "\nSTDERR\n" + completed.stderr)
    return {"treatment": treatment, "lambda": lam, "base": outbase, "trees": Path(str(outbase) + ".trees"), "branches": Path(str(outbase) + ".branches"), "rates": Path(str(outbase) + ".rates"), "log": Path(str(outbase) + ".log")}


def load_tree_lines(path: Path) -> list[Tree]:
    return [parse_tree(line) for line in path.read_text().splitlines() if line.strip()]


def summarize_tcmm(run: dict[str, object], baseline_keys: set[tuple[str, ...]]) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    trees = load_tree_lines(run["trees"])
    node = []
    audits = []
    topology_ok = True
    for i, tree in enumerate(trees, 1):
        audit = audit_tree(tree)
        topology_ok &= terminals(tree) == TAXA and topology_keys(tree) == baseline_keys
        audits.append({"treatment": run["treatment"], "lambda": run["lambda"], "reference_index": i, **{k: v for k, v in audit.items() if k != "taxa"}, "topology_preserved": str(topology_keys(tree) == baseline_keys).lower(), "taxa_preserved": str(terminals(tree) == TAXA).lower()})
        node.extend(node_rows(tree, f"tcmm_{run['treatment']}", str(i), float(run["lambda"])))
    return node, audits, {"n_trees": len(trees), "topology_preserved": topology_ok}


def observed_matrix(treatment: str) -> dict[tuple[str, str], float]:
    rows = read_tsv(PAIRWISE)
    key = "median_outside" if treatment == "outside" else "median_all"
    return {(min(r["population1"], r["population2"]), max(r["population1"], r["population2"])): float(r[key]) for r in rows}


def fit_diagnostic(trees: list[Tree], observed: dict[tuple[str, str], float], treatment: str, reference_kind: str) -> dict[str, object]:
    pred_rows = []
    predictions = []
    for i, tree in enumerate(trees, 1):
        pred = pairwise_values(tree)
        predictions.append(pred)
        residual = [pred[k] - observed[k] for k in observed]
        pred_rows.append({"treatment": treatment, "reference_kind": reference_kind, "reference_index": i, "sse": sum(x*x for x in residual), "rmse": math.sqrt(mean(x*x for x in residual)), "mae": mean(abs(x) for x in residual), "max_abs_residual": max(abs(x) for x in residual)})
    return {
        "treatment": treatment,
        "reference_kind": reference_kind,
        "n_reference_trees": len(pred_rows),
        "mean_sse": mean(r["sse"] for r in pred_rows),
        "mean_rmse": mean(r["rmse"] for r in pred_rows),
        "mean_mae": mean(r["mae"] for r in pred_rows),
        "max_abs_residual": max(r["max_abs_residual"] for r in pred_rows),
        "predicted_pairwise_mean_rmse": math.sqrt(mean((mean(pred[k] for pred in predictions) - observed[k]) ** 2 for k in observed)),
    }, pred_rows


def write_figures(node_rows_all: list[dict[str, object]], audits: list[dict[str, object]], fit_rows: list[dict[str, object]], runs: dict[str, list[Tree]]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/msrc-paper-matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGURES.mkdir(parents=True, exist_ok=True)
    stage8b = {r["clade_taxa"]: r for r in read_tsv(STAGE8B_NODE)}
    clades = sorted(stage8b, key=lambda x: (x.count("+") + 1, x))

    # Figure A: keep calibrated time scales and raw TCMM metric scales in
    # separate panels because the latter are not generally chronograms.
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    root_row = next((stage8b[c] for c in stage8b if int(c.count("+") + 1) == 12), None)
    axes[0].bar(["original", "8b outside", "8b all"], [float(root_row["original_age"]), float(root_row["outside_age"]), float(root_row["all_age"])])
    axes[0].set_ylabel("calibrated root age")
    for treatment, color in (("outside", "tab:blue"), ("all", "tab:orange")):
        vals = [float(a["root_height_max"]) for a in audits if a["treatment"] == treatment and float(a["lambda"]) == 0.0]
        axes[1].boxplot(vals, positions=[1 if treatment == "outside" else 2], widths=0.5, patch_artist=True, boxprops={"facecolor": color})
    axes[1].set_xticks([1, 2], ["TCMM outside", "TCMM all"]); axes[1].set_ylabel("raw TCMM root metric"); axes[1].set_title("not calibrated time")
    fig.tight_layout()
    for ext in ("pdf", "png"): fig.savefig(FIGURES / f"atlantic_cod_stage8c_tree_comparison.{ext}", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharex=True)
    ax = axes[0]
    series = [("original", "original_age"), ("stage8b outside", "outside_age"), ("stage8b all", "all_age")]
    for label, col in series:
        ax.plot(range(len(clades)), [float(stage8b[c][col]) for c in clades], marker="o", ms=3, lw=1, label=label)
    ax.set_xlabel("internal node (descendant-set order)"); ax.set_ylabel("calibrated age"); ax.legend(fontsize=8); ax.tick_params(axis="x", labelbottom=False); fig.tight_layout()
    ax = axes[1]
    for label, treatment in (("TCMM outside", "outside"), ("TCMM all", "all")):
        vals = defaultdict(list)
        for r in node_rows_all:
            if r["treatment"] == f"tcmm_{treatment}" and float(r["lambda"]) == 0.0:
                vals[r["clade_taxa"]].append(float(r["age_metric_from_root_height"]))
        ax.plot(range(len(clades)), [mean(vals[c]) if vals[c] else math.nan for c in clades], marker=".", ms=3, lw=1, label=label)
    ax.set_xlabel("internal node (descendant-set order)"); ax.set_ylabel("TCMM metric height"); ax.legend(fontsize=8); ax.tick_params(axis="x", labelbottom=False)
    fig.tight_layout()
    for ext in ("pdf", "png"): fig.savefig(FIGURES / f"atlantic_cod_stage8c_node_age_comparison.{ext}", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5))
    outside = {r["clade_taxa"]: r for r in node_rows_all if r["treatment"] == "tcmm_outside" and float(r["lambda"]) == 0.0}
    allv = {r["clade_taxa"]: r for r in node_rows_all if r["treatment"] == "tcmm_all" and float(r["lambda"]) == 0.0}
    xs, ys, stage8b_shift, tcmm_shift = [], [], [], []
    for c in clades:
        if c in outside and c in allv:
            xs.append(mean(float(r["age_metric_from_root_height"]) for r in [outside[c]])); ys.append(mean(float(r["age_metric_from_root_height"]) for r in [allv[c]]))
    for c in clades:
        if c in stage8b:
            stage8b_shift.append(float(stage8b[c]["all_age"]) - float(stage8b[c]["outside_age"]))
        if c in outside and c in allv:
            tcmm_shift.append(float(allv[c]["age_metric_from_root_height"]) - float(outside[c]["age_metric_from_root_height"]))
    ax.plot(range(len(stage8b_shift)), stage8b_shift, "o-", ms=3, lw=1, label="Stage 8b")
    ax.plot(range(len(tcmm_shift)), tcmm_shift, ".-", ms=4, lw=1, label="TCMM")
    ax.axhline(0, color="k", lw=0.8); ax.set_xlabel("internal node (descendant-set order)"); ax.set_ylabel("all minus outside shift"); ax.legend(); ax.tick_params(axis="x", labelbottom=False); fig.tight_layout()
    for ext in ("pdf", "png"): fig.savefig(FIGURES / f"atlantic_cod_stage8c_all_minus_outside.{ext}", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    for treatment, color in (("outside", "tab:blue"), ("all", "tab:orange")):
        vals = [r["rmse"] for r in fit_rows if r["treatment"] == treatment and r["reference_kind"].startswith("tcmm") and "rmse" in r]
        if vals: ax.hist(vals, bins=12, alpha=0.6, label=f"TCMM {treatment}", color=color)
    ax.set_xlabel("pairwise MRCA-metric RMSE"); ax.set_ylabel("reference-tree count"); ax.legend(); fig.tight_layout()
    for ext in ("pdf", "png"): fig.savefig(FIGURES / f"atlantic_cod_stage8c_fit_diagnostics.{ext}", dpi=220)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tcmm-repo", type=Path, default=Path("/tmp/tcmm-stage8c"))
    p.add_argument("--tcmm-python", default=sys.executable)
    p.add_argument("--lambdas", nargs="+", type=float, default=[0.0, 0.001])
    args = p.parse_args()
    RESULTS.mkdir(parents=True, exist_ok=True)
    if not (args.tcmm_repo / "multiple_tree_matching.py").exists(): raise SystemExit("Pinned TCMM repository not found; pass --tcmm-repo")
    baseline = load_baseline()
    baseline_audit = audit_tree(baseline)
    query_tree = prepare_query_tree()
    query_audit = audit_tree(Phylo.read(str(query_tree), "newick"))
    refs, ref_paths = extract_references()
    runs = []
    for treatment in ("outside", "all"):
        for lam in args.lambdas:
            runs.append(run_tcmm(args.tcmm_repo, args.tcmm_python, query_tree, ref_paths[treatment], treatment, lam))
    baseline_keys = topology_keys(baseline)
    node_rows_all, audits_all, run_audits = [], [], []
    fit_summary, fit_individual = [], []
    run_trees = {}
    for run in runs:
        nodes, audits, run_audit = summarize_tcmm(run, baseline_keys)
        node_rows_all.extend(nodes); audits_all.extend(audits); run_audits.append({**run_audit, "treatment": run["treatment"], "lambda": run["lambda"], "trees": str(run["trees"].relative_to(ROOT))})
        trees = load_tree_lines(run["trees"]); run_trees[f"{run['treatment']}_{run['lambda']}"] = trees
        summary, individual = fit_diagnostic(trees, observed_matrix(run["treatment"]), run["treatment"], f"tcmm_lambda_{run['lambda']}")
        fit_summary.append(summary); fit_individual.extend(individual)

    stage8b_rows = read_tsv(STAGE8B_NODE)
    comparison = []
    for r in stage8b_rows:
        clade = r["clade_taxa"]
        for treatment in ("outside", "all"):
            vals = [float(x["age_metric_from_root_height"]) for x in node_rows_all if x["treatment"] == f"tcmm_{treatment}" and float(x["lambda"]) == 0.0 and x["clade_taxa"] == clade]
            if vals:
                tcmm_age = mean(vals)
                comparison.append({"clade_taxa": clade, "n_descendants": r["n_descendants"], "original_age": r["original_age"], "stage8b_outside_age": r["outside_age"], "stage8b_all_age": r["all_age"], "tcmm_outside_age": tcmm_age if treatment == "outside" else "", "tcmm_all_age": tcmm_age if treatment == "all" else "", "comparison_treatment": treatment})
    # Combine one row per clade with both primary TCMM estimates.
    joined = {}
    for row in comparison:
        target = joined.setdefault(row["clade_taxa"], {})
        for key, value in row.items():
            if value != "":
                target[key] = value
    final_nodes = []
    for clade, row in joined.items():
        out = float(row.get("tcmm_outside_age", 0.0) or 0.0); allv = float(row.get("tcmm_all_age", 0.0) or 0.0)
        final_nodes.append({"clade_taxa": clade, "n_descendants": row["n_descendants"], "original_age": row["original_age"], "stage8b_outside_age": row["stage8b_outside_age"], "stage8b_all_age": row["stage8b_all_age"], "tcmm_outside_age": out, "tcmm_all_age": allv, "stage8b_all_minus_outside": float(row["stage8b_all_age"]) - float(row["stage8b_outside_age"]), "tcmm_all_minus_outside": allv - out, "tcmm_minus_stage8b_outside": out - float(row["stage8b_outside_age"]), "tcmm_minus_stage8b_all": allv - float(row["stage8b_all_age"]), "percentage_tcmm_shift": (100.0 * (allv - out) / out) if out else math.nan})
    fields = list(final_nodes[0]) if final_nodes else []
    write_tsv(RESULTS / "tcmm_node_age_comparison.tsv", final_nodes, fields)
    branch_rows = []
    for clade in sorted(joined):
        def branch_for(treatment: str) -> float:
            vals = [float(r["branch_length"]) for r in node_rows_all if r["treatment"] == f"tcmm_{treatment}" and r["clade_taxa"] == clade]
            return mean(vals) if vals else math.nan
        branch_rows.append({"clade_taxa": clade, "tcmm_outside_branch_length_mean": branch_for("outside"), "tcmm_all_branch_length_mean": branch_for("all"), "tcmm_all_minus_outside_branch_length": branch_for("all") - branch_for("outside")})
    write_tsv(RESULTS / "tcmm_branch_length_comparison.tsv", branch_rows, list(branch_rows[0]))
    write_tsv(RESULTS / "tcmm_ultrametric_audit.tsv", audits_all, list(audits_all[0]))
    fit_fields = sorted({key for row in fit_summary + fit_individual for key in row})
    write_tsv(RESULTS / "tcmm_fit_diagnostics.tsv", fit_summary + fit_individual, fit_fields)
    sensitivity = []
    for a in run_audits:
        audit_rows = [x for x in audits_all if x["treatment"] == a["treatment"] and x["lambda"] == a["lambda"]]
        sensitivity.append({"treatment": a["treatment"], "lambda": a["lambda"], "n_reference_trees": a["n_trees"], "topology_preserved": a["topology_preserved"], "raw_negative_edge_trees": sum(int(x["n_negative_edges"]) > 0 for x in audit_rows), "minimum_raw_edge": min(float(x["minimum_edge"]) for x in audit_rows), "maximum_root_to_tip_deviation": max(float(x["root_to_tip_deviation"]) for x in audit_rows), "mean_root_height": mean(float(x["root_height_max"]) for x in audit_rows)})
    write_tsv(RESULTS / "tcmm_sensitivity.tsv", sensitivity, list(sensitivity[0]))
    write_figures(node_rows_all, audits_all, fit_individual, run_trees)

    tcmm_commit = subprocess.check_output(["git", "-C", str(args.tcmm_repo), "rev-parse", "HEAD"], text=True).strip()
    try:
        from importlib.metadata import version
        package_versions = {name: version(name) for name in ("biopython", "matplotlib", "numpy", "cvxpy", "treeswift")}
    except Exception:
        package_versions = {}
    output_checksums = {}
    for output in RESULTS.rglob("*"):
        if output.is_file() and output.name not in {"manifest.json", "tcmm_environment.json"}:
            output_checksums[str(output.relative_to(ROOT))] = sha256(output)
    manifest = {
        "stage": "8c",
        "analysis": "published TCMM projection onto fixed Atlantic-cod baseline topology",
        "tcmm_repository": "https://github.com/shayesteh99/TCMM",
        "tcmm_commit": tcmm_commit,
        "tcmm_entrypoint": "multiple_tree_matching.py",
        "tcmm_objective": "original implementation; nonnegative edge lengths; metric matching with optional query-length regularization",
        "lambda_values": args.lambdas,
        "query_tree": str(BASELINE.relative_to(ROOT)),
        "tcmm_query_tree": str(QUERY_TREE.relative_to(ROOT)),
        "query_normalization": "set only the induced-tree external root stem to zero; retain rooted topology, labels, and all within-tree branch lengths; no rerooting or MRCA-based root selection",
        "reference_source": str(WINDOWS.relative_to(ROOT)),
        "reference_type": "published 250-kb Source Data Fig. 4 MCC-like window trees; not per-window posterior samples",
        "reference_counts": {t: sum(1 for r in refs if r["treatment"] == t and r["valid_for_tcmm"] == "true") for t in ("outside", "all")},
        "original_inputs": {str(p.relative_to(ROOT)): sha256(p) for p in (BASELINE, WINDOWS, PAIRWISE, REGIONS)},
        "output_checksums_before_manifest": output_checksums,
        "baseline_audit": {k: v for k, v in baseline_audit.items() if k != "taxa"},
        "query_audit": {k: v for k, v in query_audit.items() if k != "taxa"},
        "stage8b_preserved": True,
        "command_python": args.tcmm_python,
        "platform": platform.platform(),
        "python": sys.version,
        "package_versions": package_versions,
        "limitations": ["No per-window SNAPP posterior tree samples were available locally.", "TCMM raw outputs are preserved; tiny negative numerical edges and non-ultrametricity are reported, not silently repaired.", "Metric-tree output is not automatically interpreted as calibrated divergence time."],
    }
    (RESULTS / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    (RESULTS / "tcmm_environment.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
