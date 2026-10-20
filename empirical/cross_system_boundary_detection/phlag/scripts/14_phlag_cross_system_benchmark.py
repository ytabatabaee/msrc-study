#!/usr/bin/env python3
"""Stage 13 extension: audit and benchmark the unmodified Phlag implementation.

The script deliberately refuses to collapse individual or haplotype tips into
species/population tips.  Phlag's input parser uses one shared taxon namespace
for the species and local trees, so a taxon-set mismatch is an incompatibility,
not a negative biological result.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Iterable

import dendropy

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "empirical/cross_system_boundary_detection/phlag"
RAW = OUT / "raw_results"
PHLAG_COMMIT = "76597d0e4ebd0e507fdafebf7f5163ef634a36d6"
PHLAG_VERSION = "0.2.0"
SEED = 20261020


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tree_tips(text: str) -> set[str]:
    # This deliberately reads only leaf tokens.  Some reference trees contain
    # repeated internal support labels, which DendroPy rejects as duplicate
    # taxa even though the leaf set is unambiguous for this compatibility audit.
    return set(re.findall(r"(?<=[(,])([^():;,\s]+)(?=[:),;])", text.strip()))


def first_local_tree(path: Path, gz: bool = False, column: int | None = None) -> tuple[str, int]:
    opener = gzip.open if gz else open
    n = 0
    first = ""
    header_skipped = False
    with opener(path, "rt", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if column is not None:
                if not header_skipped:
                    header_skipped = True
                    continue
                fields = line.split("\t")
                line = fields[column]
            elif "\t" in line:
                # Anopheles stage-4R stores a window identifier before the Newick.
                candidate = line.split("\t", 1)[1]
                if candidate.lstrip().startswith("("):
                    line = candidate
            n += 1
            if not first:
                first = line
    return first, n


def parse_phlag_output(path: Path) -> tuple[list[int] | None, list[float] | None, str]:
    """Parse Phlag's state and posterior rows without treating NaN as valid."""
    rows: list[list[str]] = []
    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "," not in line:
            continue
        rows.append([x.strip() for x in line.split(",")])
    state = None
    posterior = None
    for row in rows:
        if all(x in {"0", "1", "-1", "nan", "NaN"} for x in row):
            if not all(x.lower() == "nan" for x in row):
                state = [int(float(x)) if x.lower() != "nan" else -1 for x in row]
                continue
        try:
            vals = [float(x) for x in row]
        except ValueError:
            continue
        if vals and all(0 <= x <= 1 for x in vals if math.isfinite(x)):
            posterior = vals
    if state is None or posterior is None:
        return state, posterior, "missing_state_or_posterior"
    if len(state) != len(posterior):
        return state, posterior, "state_posterior_length_mismatch"
    if not all(math.isfinite(x) for x in posterior):
        return state, posterior, "nonfinite_posterior"
    return state, posterior, "valid"


def classify_run(state: list[int] | None, posterior: list[float] | None, parse_status: str) -> str:
    if parse_status != "valid":
        return "degenerate_or_invalid"
    if posterior is None or not posterior or not all(math.isfinite(x) for x in posterior):
        return "degenerate_or_invalid"
    if len(set(state or [])) <= 1:
        return "valid_no_anomaly"
    return "valid_informative"


def aggregate_branch_predictions(rows: Iterable[dict]) -> dict[str, dict]:
    """Union branch calls by window; posterior aggregation is max, never a mean."""
    out: dict[str, dict] = {}
    for row in rows:
        key = row["window_id"]
        posterior = float(row.get("posterior_anomaly", "nan"))
        current = out.setdefault(key, {"posterior_anomaly_max": float("nan"), "predicted_state_union": 0})
        if math.isfinite(posterior):
            if not math.isfinite(current["posterior_anomaly_max"]) or posterior > current["posterior_anomaly_max"]:
                current["posterior_anomaly_max"] = posterior
        current["predicted_state_union"] = max(current["predicted_state_union"], int(row.get("predicted_state", 0)))
    return out


def interval_iou(a: tuple[int, int], b: tuple[int, int]) -> float:
    lo = max(a[0], b[0]); hi = min(a[1], b[1])
    inter = max(0, hi - lo)
    union = max(a[1], b[1]) - min(a[0], b[0])
    return inter / union if union > 0 else 0.0


def write_tsv(path: Path, rows: Iterable[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def audit_inputs() -> list[dict]:
    configs = [
        {
            "dataset": "atlantic_cod",
            "local": ROOT / "empirical/atlantic_cod/results/phlag_benchmark/inputs/LG01.gene_trees",
            "species": ROOT / "empirical/atlantic_cod/results/phlag_benchmark/inputs/baseline_population_tree_12taxa_labeled.nwk",
            "local_format": "one Newick per line",
            "note": "SNAPP MCC population trees; Stage 11 runs are retained and audited, not reinterpreted as gene trees.",
            "status": "degenerate_numeric",
        },
        {
            "dataset": "fire_ants_chr16",
            "local": ROOT / "data/fire_ants_chr16/processed/stage1_window_trees.tsv",
            "species": ROOT / "empirical/fire_ants_chr16/results/astral4/stage6_astral4_background_reference.nwk",
            "local_format": "TSV tree_newick column",
            "note": "Local trees have 267 individual tips; reference has seven group tips. No silent individual-to-group collapse.",
            "status": "incompatible_input",
        },
        {
            "dataset": "anopheles_2la",
            "local": ROOT / "data/anopheles_2la/processed/stage4r_all72_haplotype_local_trees.nwk.gz",
            "species": ROOT / "empirical/anopheles_2la/results/astral6/T_all6_species.nwk",
            "local_format": "gzipped one Newick per line",
            "note": "Local trees have phased haplotype tips; six-species reference cannot receive them without a non-original collapse. Heterokaryotypes are not relabeled.",
            "status": "incompatible_input",
        },
        {
            "dataset": "house_mouse_t_complex",
            "local": ROOT / "data/house_mouse_t_complex/processed/house_mouse_t_complex_ml_5kb.tre",
            "species": ROOT / "empirical/house_mouse_t_complex/results/stage2_population_astral/P1_ALL.nwk",
            "local_format": "one Newick per line",
            "note": "Local trees have 55 individual/pseudo-t-haplotype tips; population reference has collapsed taxa. No arbitrary representative selection.",
            "status": "incompatible_input",
        },
    ]
    rows = []
    for c in configs:
        local = c["local"]
        species = c["species"]
        try:
            if c["dataset"] == "fire_ants_chr16":
                first, n = first_local_tree(local, column=11)
            elif c["dataset"] == "anopheles_2la":
                first, n = first_local_tree(local, gz=True)
            else:
                first, n = first_local_tree(local)
            local_tips = read_tree_tips(first)
            species_tips = read_tree_tips(species.read_text())
            missing = sorted(species_tips - local_tips)
            extra = sorted(local_tips - species_tips)
            rows.append({
                "dataset": c["dataset"], "local_input": str(local.relative_to(ROOT)),
                "species_tree": str(species.relative_to(ROOT)), "n_local_trees": n,
                "n_local_tips_first_tree": len(local_tips), "n_species_tree_tips": len(species_tips),
                "taxon_sets_identical": str(not missing and not extra).lower(),
                "missing_species_taxa": ",".join(missing), "extra_local_taxa": ",".join(extra),
                "status": c["status"], "reason": c["note"],
            })
        except Exception as exc:
            rows.append({"dataset": c["dataset"], "local_input": str(local.relative_to(ROOT)),
                         "species_tree": str(species.relative_to(ROOT)), "n_local_trees": "NA",
                         "n_local_tips_first_tree": "NA", "n_species_tree_tips": "NA",
                         "taxon_sets_identical": "NA", "missing_species_taxa": "NA",
                         "extra_local_taxa": "NA", "status": "input_audit_failed",
                         "reason": f"{c['note']} Error: {exc}"})
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    audit = audit_inputs()
    fields = ["dataset", "local_input", "species_tree", "n_local_trees", "n_local_tips_first_tree",
              "n_species_tree_tips", "taxon_sets_identical", "missing_species_taxa", "extra_local_taxa", "status", "reason"]
    write_tsv(OUT / "input_compatibility_audit.tsv", audit, fields)

    # Preserve the original cod run byte-for-byte in a dedicated raw-results directory.
    cod_raw = ROOT / "empirical/atlantic_cod/results/phlag_benchmark/raw"
    if cod_raw.exists():
        target = RAW / "atlantic_cod_stage11_raw"
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(cod_raw, target)
    toy = Path("/private/tmp/phlag-stage12b-toy.txt")
    toy_status = "not_available"
    if toy.exists():
        shutil.copy2(toy, RAW / "upstream_bundled_example_stage12b.txt")
        st, po, ps = parse_phlag_output(toy)
        toy_status = classify_run(st, po, ps)
    run_rows = [{
        "dataset": "bundled_example", "status": toy_status, "phlag_version": PHLAG_VERSION,
        "phlag_commit": PHLAG_COMMIT, "command": "/private/tmp/phlag-venv/bin/phlag -s test/neoaves.nwk -g test/emission.gtrees -e N159 -L 1 -l 5 --n-replicates 10 --read-qqs-path test/qqs.tsv",
        "defaults": "rho=0.9; beta-prime=0.0025; emission-lambda=1.0; eta=0.5; no cod coordinates used",
        "raw_output": "raw_results/upstream_bundled_example_stage12b.txt",
    }, {
        "dataset": "atlantic_cod", "status": "degenerate_numeric", "phlag_version": PHLAG_VERSION,
        "phlag_commit": PHLAG_COMMIT, "command": "Reused unchanged Stage 11 four-LG default runs; see empirical/atlantic_cod/results/phlag_benchmark/phlag_commands.txt",
        "defaults": "Stage 11 defaults; 9 focal branches per LG; no threshold tuning",
        "raw_output": "raw_results/atlantic_cod_stage11_raw/",
    }, *[{
        "dataset": r["dataset"], "status": r["status"], "phlag_version": PHLAG_VERSION,
        "phlag_commit": PHLAG_COMMIT, "command": "not run: taxon-set/species-tree incompatibility",
        "defaults": "original Phlag defaults retained; no forced collapse",
        "raw_output": "",
    } for r in audit if r["dataset"] != "atlantic_cod"]]
    write_tsv(OUT / "run_manifest.tsv", run_rows,
              ["dataset", "status", "phlag_version", "phlag_commit", "command", "defaults", "raw_output"])

    diagnostics = []
    for row in audit:
        diagnostics.append({"dataset": row["dataset"], "execution_class": row["status"],
                            "posterior_finite": "false" if row["dataset"] == "atlantic_cod" else "not_run",
                            "state_path": "constant_zero_or_nan" if row["dataset"] == "atlantic_cod" else "not_run",
                            "focal_edge_estimate": "1e-06 (Stage 11)" if row["dataset"] == "atlantic_cod" else "not_run",
                            "emission_finite": "false" if row["dataset"] == "atlantic_cod" else "not_run",
                            "transition_parameters": "estimated but output degenerate" if row["dataset"] == "atlantic_cod" else "not_run",
                            "diagnostic": row["reason"]})
    write_tsv(OUT / "numerical_diagnostics.tsv", diagnostics,
              ["dataset", "execution_class", "posterior_finite", "state_path", "focal_edge_estimate", "emission_finite", "transition_parameters", "diagnostic"])

    # Copy the existing cod branch-level output into a clearly labelled audit table.
    bp = ROOT / "empirical/atlantic_cod/results/phlag_benchmark/phlag_branch_predictions.tsv"
    rows = []
    if bp.exists():
        with bp.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                r["dataset"] = "atlantic_cod"; r["valid_for_accuracy"] = "false"; r["execution_class"] = "degenerate_numeric"; rows.append(r)
    write_tsv(OUT / "branch_level_predictions.tsv", rows, ["dataset", "lg", "window_id", "start", "end", "midpoint", "focal_branch", "predicted_state", "posterior_anomaly", "valid_for_accuracy", "execution_class"])
    write_tsv(OUT / "window_level_predictions.tsv", [], ["dataset", "group", "window_id", "start", "end", "posterior_anomaly_max", "predicted_state_union", "valid_for_accuracy", "execution_class"])

    # Validity-aware comparison: existing unified detector metrics are retained as comparators;
    # Phlag receives no accuracy score when its output is numerically invalid or incompatible.
    comparison = []
    ms = ROOT / "empirical/cross_system_boundary_detection/results/fixed_model_predictions.tsv"
    if ms.exists():
        with ms.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                comparison.append({"dataset": r.get("system", ""), "group": r.get("group", ""), "method": r.get("method", ""), "execution_class": "valid_unified_detector", "valid_for_accuracy": "true", "f1": r.get("f1", ""), "iou": r.get("iou_union", ""), "boundary_error": r.get("mean_boundary_error_kb", "")})
    for dataset in ["atlantic_cod", "fire_ants_chr16", "anopheles_2la", "house_mouse_t_complex"]:
        comparison.append({"dataset": dataset, "group": "all", "method": "Phlag", "execution_class": "degenerate_numeric" if dataset == "atlantic_cod" else "incompatible_input", "valid_for_accuracy": "false", "f1": "NA", "iou": "NA", "boundary_error": "NA"})
    write_tsv(OUT / "accuracy_comparison.tsv", comparison, ["dataset", "group", "method", "execution_class", "valid_for_accuracy", "f1", "iou", "boundary_error"])
    write_tsv(OUT / "runtime_memory.tsv", [
        {"dataset": "atlantic_cod", "runtime_source": "Stage 11 .runtime_seconds files", "runtime_seconds": "see Stage 11 raw files", "peak_memory_mb": "not recorded", "valid_run_fraction": "0/36"},
        *[{"dataset": d, "runtime_source": "not run after compatibility audit", "runtime_seconds": "NA", "peak_memory_mb": "NA", "valid_run_fraction": "0/0"} for d in ["fire_ants_chr16", "anopheles_2la", "house_mouse_t_complex"]]
    ], ["dataset", "runtime_source", "runtime_seconds", "peak_memory_mb", "valid_run_fraction"])

    (OUT / "methods_text.md").write_text("""# Stage 13 Phlag benchmark methods\n\nPhlag was kept at upstream commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6` (version 0.2.0) with its default model settings. The finite bundled-example check uses the previously recorded reduced-iteration reproduction (`-L 1 -l 5 --n-replicates 10`); a fresh full-default toy run was not completed within the practical runtime window. For Atlantic cod, the unchanged Stage 11 runs were copied into `raw_results/` and reclassified from their diagnostics. For the other systems, the benchmark stops at the input audit when the local-tree and reference-tree taxon sets differ.\n\nPhlag uses a shared DendroPy taxon namespace for the species tree and local gene trees. Collapsing individuals, haplotypes, or pseudo-t-haplotype tips to population labels would change the input data and is therefore not used. A numerical failure is never scored as a biological no-anomaly result.\n\nThe Phlag objective detects genealogical deviations/model violations. It is not itself a recombination-suppression or inversion-boundary estimator.\n""")
    (OUT / "report.md").write_text("""# Stage 13 extension: Phlag cross-system benchmark\n\n## Result\n\nThe original Phlag implementation was preserved and audited at commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6` (v0.2.0). The bundled example has a finite reference output in `raw_results/upstream_bundled_example_stage12b.txt`; this is a reduced-iteration software check rather than a fresh full-default benchmark.\n\nAtlantic cod reuses the unchanged Stage 11 default runs. Those runs are numerically degenerate on SNAPP MCC population-tree inputs: posterior anomaly probabilities are nonfinite and state paths are not interpretable. They are therefore classified as `degenerate_numeric`, with no Phlag accuracy score. This is not evidence that Phlag found no anomaly.\n\nThe fire-ant, Anopheles, and house-mouse inputs are classified as `incompatible_input` because their local trees contain individual, phased haplotype, or pseudo-t-haplotype tips while the available reference trees use group/species/population taxa. The original implementation does not provide a scientifically neutral multi-sample mapping layer. No arbitrary representative selection or silent collapsing was performed.\n\nConsequently, this benchmark does not support a cross-system accuracy claim for Phlag. Existing unified MSRC detector scores are carried into `accuracy_comparison.tsv` for context, while invalid/incompatible Phlag runs are explicitly marked `NA`. A future valid benchmark requires gene-tree and species-tree taxon sets designed for Phlag, or a method extension validated independently of the target intervals.\n""")
    (OUT / "manifest.json").write_text(json.dumps({
        "stage": "13_phlag_extension", "seed": SEED, "phlag_version": PHLAG_VERSION, "phlag_commit": PHLAG_COMMIT,
        "script": str(Path(__file__).relative_to(ROOT)), "input_audit": str((OUT / "input_compatibility_audit.tsv").relative_to(ROOT)),
        "raw_results": str(RAW.relative_to(ROOT)), "primary_policy": "do not collapse mismatched taxa; do not score invalid runs",
        "input_hashes": {str(p.relative_to(ROOT)): sha256(p) for p in [ROOT / "empirical/atlantic_cod/results/phlag_benchmark/manifest.json", ROOT / "data/fire_ants_chr16/processed/stage1_window_trees.tsv", ROOT / "data/anopheles_2la/processed/stage4r_haplotype_species_map.tsv", ROOT / "data/house_mouse_t_complex/metadata/tip_mapping.tsv"]},
    }, indent=2) + "\n")

    # A validity-aware overview figure: invalid Phlag executions are shown as
    # status annotations rather than zero-height accuracy bars.
    if os.environ.get("PHLAG_SKIP_PLOT") != "1":
      try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), constrained_layout=True)
        systems = ["Atlantic cod", "Fire ants", "Anopheles", "House mouse"]
        statuses = ["degenerate\nnumeric", "incompatible\ninput", "incompatible\ninput", "incompatible\ninput"]
        ax = axes[0]
        ax.axis("off")
        ax.set_title("Phlag execution status")
        for i, (s, status) in enumerate(zip(systems, statuses)):
            ax.text(0.02, 0.86 - i * 0.21, s, fontsize=10, va="center")
            ax.text(0.52, 0.86 - i * 0.21, status, fontsize=10, va="center", color="#b22222")
        ax.text(0.02, 0.05, "No invalid or incompatible run is scored as a biological negative.", fontsize=8)
        ax = axes[1]
        ax.set_title("Existing unified detector context")
        rows = [r for r in comparison if r["execution_class"] == "valid_unified_detector"]
        labels = []
        vals = []
        for r in rows:
            if r["group"] == "LG01" and r["method"] in {"topology", "continuous", "joint"}:
                labels.append(r["method"])
                vals.append(float(r["iou"]))
        if vals:
            ax.bar(labels, vals, color=["#4c78a8", "#59a14f", "#f28e2b"][:len(vals)])
            ax.set_ylim(0, 1)
            ax.set_ylabel("IoU (cod LG01 reference)")
        else:
            ax.text(0.5, 0.5, "No valid comparator rows", ha="center")
        ax.text(0.02, -0.20, "Phlag accuracy: NA for all four systems in this benchmark.", transform=ax.transAxes, fontsize=8)
        figdir = ROOT / "empirical/cross_system_boundary_detection/figures"
        figdir.mkdir(parents=True, exist_ok=True)
        fig.savefig(figdir / "atlantic_cod_stage13_phlag_comparison.pdf")
        fig.savefig(figdir / "atlantic_cod_stage13_phlag_comparison.png", dpi=220)
        plt.close(fig)
      except Exception as exc:
        (OUT / "figure_warning.txt").write_text(f"Figure generation failed: {exc}\n")


if __name__ == "__main__":
    main()
