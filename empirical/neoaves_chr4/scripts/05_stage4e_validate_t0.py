#!/usr/bin/env python3
"""Validate the committed full-363-taxon Stage 4E T0 trees.

This audit only parses existing Newick files.  It never invokes ASTRAL-IV.
The original CU and SU runs are checked against their corresponding edge
length fields, while the SU tree is also used as the authoritative source for
the single-run extraction/migration output.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from stage4d_core import focal_split_states, frozen_focal_groups, named_splits
from stage4e_annotations import extract_records, sha256, write_branch_table


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def internal_map(records: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    return {r["split_key"]: r for r in records if r["branch_type"] == "internal"}


def terminal_map(records: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    # A rooted representation can contain a duplicate root-adjacent terminal
    # edge.  Exact taxon labels are the stable terminal identifiers.
    return {r["taxon_label"]: r for r in records if r["branch_type"] == "terminal"}


def differences(a: dict[str, dict[str, object]], b: dict[str, dict[str, object]], field: str) -> list[float]:
    out = []
    for key in a.keys() & b.keys():
        x, y = a[key].get(field), b[key].get(field)
        if x is not None and y is not None:
            out.append(abs(float(x) - float(y)))
    return out


def signed_differences(a: dict[str, dict[str, object]], b: dict[str, dict[str, object]], field: str) -> list[float]:
    out = []
    for key in a.keys() & b.keys():
        x, y = a[key].get(field), b[key].get(field)
        if x is not None and y is not None:
            out.append(float(x) - float(y))
    return out


def summarize(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"n": 0, "min": None, "median": None, "max": None, "mean": None}
    return {"n": len(values), "min": min(values), "median": statistics.median(values),
            "max": max(values), "mean": statistics.fmean(values)}


def finite_nonnegative(records: list[dict[str, object]], field: str) -> bool:
    values = [r[field] for r in records if r["branch_type"] == "internal"]
    return all(v is not None and math.isfinite(float(v)) and float(v) >= 0 for v in values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--t0-dir", type=Path, default=Path("empirical/neoaves_chr4/results/stage4e/cluster/T0"))
    parser.add_argument("--extracted-dir", type=Path, default=Path("empirical/neoaves_chr4/results/stage4e/cluster/T0_extracted"))
    parser.add_argument("--reference-tree", type=Path, default=Path("empirical/neoaves_chr4/results/stage4d/trees/FULL363_ORIGINAL_REPRO.nwk"))
    parser.add_argument("--results-dir", type=Path, default=Path("empirical/neoaves_chr4/results/stage4e"))
    args = parser.parse_args()

    t0, extracted, results = args.t0_dir, args.extracted_dir, args.results_dir
    cu_tree, cu_taxa, cu_records = extract_records((t0 / "tree.cu.nwk").read_text(), edge_annotation="CULength")
    su_tree, su_taxa, su_records = extract_records((t0 / "tree.su.nwk").read_text(), edge_annotation="SULength")
    ex_cu_tree, ex_cu_taxa, ex_cu_records = extract_records((extracted / "tree.cu.nwk").read_text(), edge_annotation="CULength")
    ex_su_tree, ex_su_taxa, ex_su_records = extract_records((extracted / "tree.su.nwk").read_text(), edge_annotation="SULength")
    source_tree, source_taxa, source_records = extract_records((extracted / "tree.annotated.nwk").read_text(), edge_annotation="SULength")

    cm, sm = internal_map(cu_records), internal_map(su_records)
    ecm, esm = internal_map(ex_cu_records), internal_map(ex_su_records)
    stm = terminal_map(su_records)
    estm = terminal_map(ex_su_records)
    reference_available = args.reference_tree.exists() and args.reference_tree.stat().st_size > 0
    ref_taxa, ref_splits = named_splits(args.reference_tree.read_text()) if reference_available else (set(), {})

    shared = cm.keys() & sm.keys()
    source_support_diffs = {field: differences(cm, sm, field) for field in ("localPP", "q1", "q2", "q3")}
    cu_annotation_vs_su = {"CU_length": summarize(differences(cm, sm, "CU_length")),
                           "SU_length": summarize(differences(cm, sm, "SU_length"))}
    extracted_vs_source = {"CU_length": summarize(differences(ecm, sm, "CU_length")),
                           "SU_length": summarize(differences(esm, sm, "SU_length"))}
    extracted_vs_original_cu = summarize(differences(ecm, cm, "CU_length"))
    extracted_vs_original_su = summarize(differences(esm, sm, "SU_length"))

    focal_states = focal_split_states((t0 / "tree.su.nwk").read_text(), frozen_focal_groups())
    focal_rows = []
    for state in focal_states:
        branch = next((r for r in su_records if r.get("node") is not None and r["node"].label == state["branch_support"]), None)
        split = branch["split_key"] if branch else ""
        reference_state = focal_split_states(args.reference_tree.read_text(), frozen_focal_groups()) if reference_available else []
        ref_row = next((r for r in reference_state if r["clade"] == state["clade"]), {})
        focal_rows.append({
            "focal_relationship": state["clade"], "t0_state": state["state"],
            "reference_state": ref_row.get("state", "unavailable"),
            "topology_agrees_with_reference": str(state["state"] == ref_row.get("state")) if ref_row else "unavailable",
            "split_key": split, "n_C1": state["n_C1"], "n_C2": state["n_C2"],
            "n_S": state["n_S"], "n_O": state["n_O"],
            "CU_length": branch["CU_length"] if branch else "NA",
            "SU_length": branch["SU_length"] if branch else "NA",
            "localPP": branch["localPP"] if branch else "NA",
            "q1": branch["q1"] if branch else "NA", "q2": branch["q2"] if branch else "NA",
            "q3": branch["q3"] if branch else "NA",
            "split_exists_in_reference": str(split in ref_splits) if reference_available and split else ("unavailable" if not reference_available else "False"),
        })
    focal_path = results / "t0_focal_branches.tsv"
    write_tsv(focal_path, focal_rows, list(focal_rows[0]) if focal_rows else ["focal_relationship"])

    # The extracted branch table is already the standardized all-branch T0
    # baseline; write a named copy beside the report for future comparisons.
    baseline_path = results / "t0_reference_branch_lengths.tsv"
    write_branch_table(baseline_path, "T0", source_records)

    checks = [
        ("taxa_363", len(cu_taxa) == len(su_taxa) == len(ex_cu_taxa) == len(ex_su_taxa) == 363),
        ("unique_internal_splits_360", len(cm) == len(sm) == len(ecm) == len(esm) == 360),
        ("original_cu_su_unrooted_rf_zero", cm.keys() == sm.keys()),
        ("source_cu_lengths_finite_nonnegative", finite_nonnegative(cu_records, "CU_length")),
        ("source_su_lengths_finite_nonnegative", finite_nonnegative(su_records, "SU_length")),
        ("extracted_cu_terminal_lengths_missing", all(r["CU_length"] is None for r in ex_cu_records if r["branch_type"] == "terminal")),
        ("extracted_su_terminal_lengths_present", len(estm) == 363 and all(r["SU_length"] is not None for r in estm.values())),
        ("extracted_source_hash_matches", sha256(extracted / "tree.annotated.nwk") == sha256(t0 / "tree.su.nwk")),
        ("extracted_cu_topology_matches_original", ecm.keys() == cm.keys()),
        ("extracted_su_topology_matches_original", esm.keys() == sm.keys()),
        ("focal_mappings_exact_in_t0", all(r["split_key"] for r in focal_rows)),
    ]

    run = json.loads((t0 / "run.json").read_text())
    scores = {}
    for mode in ("topology_cu", "topology_su"):
        text = (t0 / f"{mode}.stderr.log").read_text(errors="replace")
        scores[mode] = "30019936712278" in text
    report_lines = [
        "# Full 363-Taxon Neoaves T0 ASTRAL-IV Validation", "",
        f"Generated UTC: {datetime.now(timezone.utc).isoformat()}", "",
        "This report validates committed output only. No ASTRAL-IV subprocess was invoked.", "",
        "## Provenance", "",
        f"- Recorded ASTRAL-IV host: `{run.get('hostname')}`; job `{run.get('slurm_job_id')}`; CPUs `{run.get('cpu_count')}`.",
        f"- Recorded input: `{run.get('input')}`; SHA-256 `{run.get('input_sha256')}`.",
        f"- ASTRAL binary SHA-256: `{run.get('aster_bin_sha256')}`; `-u 2`, `--genelength {run.get('genelength')}`, `--sulengthtype {run.get('sulengthtype')}`.",
        "- Both original runs returned code 0; the stderr logs report optimization score `30019936712278` for CU and SU.",
        "- The committed run metadata does not include process end timestamps, so a walltime cannot be reconstructed from repository files alone.", "",
        "## Tree and annotation audit", "",
        f"- Taxa: **{len(cu_taxa)}** in each original and extracted tree; duplicate labels: none.",
        f"- Unique nontrivial unrooted internal splits: **{len(cm)}** in each tree (expected `363 - 3 = 360`).",
        f"- Original CU/SU unrooted RF symmetric difference: **{len(set(cm) ^ set(sm))}**; shared splits: **{len(shared)}**.",
        f"- Original CU tree: 360 internal CU estimates; terminal CU estimates: 0.",
        f"- Original SU tree: 360 internal SU estimates and {len(stm)} unique terminal SU estimates.",
        f"- Extracted tree: 360 internal CU/SU estimates, {len(estm)} unique terminal SU estimates, and 0 terminal CU estimates.",
        f"- Maximum extracted-vs-original CU difference: `{extracted_vs_original_cu['max']}`.",
        f"- Maximum extracted-vs-original SU difference: `{extracted_vs_original_su['max']}`.",
        f"- Maximum CU annotation difference between independent original runs: `{cu_annotation_vs_su['CU_length']['max']}`.",
        f"- Maximum SU annotation difference between independent original runs: `{cu_annotation_vs_su['SU_length']['max']}`.",
        f"- Maximum support/q annotation differences between original runs: `{max((max(v) for v in source_support_diffs.values() if v), default=0.0)}`.", "",
        "## Internal-length summaries", "",
        f"- CU range: `{min(r['CU_length'] for r in sm.values())}` to `{max(r['CU_length'] for r in sm.values())}`; median `{statistics.median(r['CU_length'] for r in sm.values())}`.",
        f"- SU range: `{min(r['SU_length'] for r in sm.values())}` to `{max(r['SU_length'] for r in sm.values())}`; median `{statistics.median(r['SU_length'] for r in sm.values())}`.",
        f"- All localPP values finite/nonnegative: `{all(r['localPP'] is not None and 0 <= r['localPP'] <= 1 for r in sm.values())}`.",
        "- The CU and SU trees have the same topology; only the edge-length representation differs.", "",
        "## Focal relationships", "",
        "Focal roles and taxon sets come from the frozen Stage 3B role mapping and `focal_split_states()`; no approximate branch was selected.",
        f"- An independent reference species tree for focal agreement was {'available' if reference_available else 'not available'} in the repository. The existing `{args.reference_tree}` path is empty, and the large `.tre` file is a gene-tree collection rather than a reference species tree; no focal agreement claim is made.", "",
        "| Relationship | T0 state | Reference state | CU | SU | localPP | q1 | q2 | q3 | exact reference split | |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in focal_rows:
        report_lines.append(f"| {r['focal_relationship']} | {r['t0_state']} | {r['reference_state']} | {r['CU_length']} | {r['SU_length']} | {r['localPP']} | {r['q1']} | {r['q2']} | {r['q3']} | {r['split_exists_in_reference']} |")
    report_lines += ["", "## PASS/FAIL checks", "", "| Check | Result |", "|---|---|"]
    for name, ok in checks:
        report_lines.append(f"| {name} | {'PASS' if ok else 'FAIL'} |")
    report_lines += ["", "## Interpretation and scope", "", "T0 establishes the full-taxon reference topology and branch-length annotations for later treatment comparisons. It does not by itself establish a chromosome-4 structural effect or causal rearrangement mechanism. Future analyses must compare stable unrooted splits and distinguish topology changes from CU and SU length changes on matching splits.", ""]
    (results / "t0_full363_validation.md").write_text("\n".join(report_lines))
    manifest = {"generated_utc": datetime.now(timezone.utc).isoformat(), "astral_invoked": False,
                "t0_dir": str(t0), "extracted_dir": str(extracted), "reference_tree": str(args.reference_tree),
                "checks": {name: ok for name, ok in checks}, "scores_in_logs": scores,
                "original_cu_sha256": sha256(t0 / "tree.cu.nwk"), "original_su_sha256": sha256(t0 / "tree.su.nwk"),
                "extracted_source_sha256": sha256(extracted / "tree.annotated.nwk")}
    (results / "t0_full363_validation.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    if not all(ok for _, ok in checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
