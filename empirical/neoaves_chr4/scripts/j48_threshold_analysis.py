#!/usr/bin/env python3
"""Bounded validation and threshold preparation for completed Neoaves J48 runs.

This module deliberately does not read the large gene-tree collections.  It
validates the five completed, 48-taxon ASTRAL result trees, uses the existing
small QQS summary for an explicitly labelled quartet-level analysis, and
prepares a deterministic Bridges-2 mixture manifest.  New ASTRAL runs are
handled by the separate SLURM wrapper.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import re
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "empirical" / "neoaves_chr4"
RESULTS = BASE / "results" / "stage4e"
CLUSTER = RESULTS / "cluster"
DATA = BASE / "data" / "stage4e"
MANIFEST = BASE / "results" / "stage4d_locus_manifest.tsv"
QQS = BASE / "results" / "stage6_genomic_tracks" / "window_quartet_support.tsv"
TAXA = BASE / "data" / "stage4d" / "jarvis48_taxa.tsv"
TREATMENTS = ["J48_T0", "J48_T_PNAS", "J48_T_STRUCT_PRIMARY", "J48_T_BLOCK_500K", "J48_NONCHR4"]
EPSILONS = (0.00, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.75, 1.00)
FOCALS = ("Columbea", "N61", "N62")
FOCAL_ROLES = {
    "Columbea": (("Columba_livia", "Pterocles_gutturalis"), ("Phoenicopterus_ruber",), ("Tauraco_erythrolophus", "Cuculus_canorus", "Ciconia_maguari"), ("Gallus_gallus",)),
    "N61": (("Tauraco_erythrolophus", "Cuculus_canorus"), ("Columba_livia", "Pterocles_gutturalis"), ("Ciconia_maguari",), ("Phoenicopterus_ruber", "Gallus_gallus")),
    "N62": (("Columba_livia",), ("Pterocles_gutturalis",), ("Tauraco_erythrolophus", "Cuculus_canorus"), ("Phoenicopterus_ruber", "Ciconia_maguari", "Gallus_gallus")),
}
PNAS_INTERVALS = ((25030000, 32670000), (33510000, 34470000), (44130000, 56810000))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def boolish(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


@dataclass
class Node:
    name: str = ""
    children: tuple["Node", ...] = ()


def parse_newick(text: str) -> Node:
    """Parse the simple ASTRAL Newick form without depending on dendropy."""
    text = text.strip().rstrip(";")
    pos = 0

    def skip_annotation() -> None:
        nonlocal pos
        while pos < len(text) and text[pos] == "[":
            end = text.find("]", pos + 1)
            if end < 0:
                raise ValueError("unterminated Newick annotation")
            pos = end + 1

    def label() -> str:
        nonlocal pos
        start = pos
        while pos < len(text) and text[pos] not in ",()":
            if text[pos] == ":":
                break
            pos += 1
        value = text[start:pos].strip()
        if pos < len(text) and text[pos] == ":":
            pos += 1
            while pos < len(text) and text[pos] not in ",()":
                pos += 1
        skip_annotation()
        return value

    def subtree() -> Node:
        nonlocal pos
        if pos >= len(text):
            raise ValueError("unexpected end of Newick")
        if text[pos] == "(":
            pos += 1
            children = [subtree()]
            while pos < len(text) and text[pos] == ",":
                pos += 1
                children.append(subtree())
            if pos >= len(text) or text[pos] != ")":
                raise ValueError("unclosed Newick subtree")
            pos += 1
            name = label()
            return Node(name=name, children=tuple(children))
        return Node(name=label())

    root = subtree()
    if pos != len(text):
        raise ValueError(f"unparsed Newick suffix at {pos}")
    return root


def leaf_names(node: Node) -> set[str]:
    if not node.children:
        return {node.name}
    result: set[str] = set()
    for child in node.children:
        result.update(leaf_names(child))
    return result


def split_keys(node: Node, all_taxa: set[str]) -> set[str]:
    keys: set[str] = set()

    def visit(current: Node) -> set[str]:
        if not current.children:
            return {current.name}
        left: set[str] = set()
        for child in current.children:
            left.update(visit(child))
        if 1 < len(left) < len(all_taxa) - 1:
            other = all_taxa - left
            chosen = min(tuple(sorted(left)), tuple(sorted(other)))
            keys.add("|".join(chosen))
        return left

    visit(node)
    return keys


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_taxa() -> set[str]:
    rows = read_tsv(TAXA)
    names = set()
    for row in rows:
        names.add(row.get("stiller2024_taxon") or row.get("taxon") or row.get("jarvis_taxon") or "")
    return {name for name in names if name}


def branch_rows(treatment: str) -> dict[str, dict[str, str]]:
    rows = read_tsv(CLUSTER / treatment / "branch_lengths.tsv")
    return {row["split_key"]: row for row in rows if row.get("branch_type") == "internal"}


def validate_results() -> tuple[dict[str, dict[str, object]], dict[tuple[str, str], int]]:
    expected = read_taxa()
    summaries: dict[str, dict[str, object]] = {}
    splits: dict[str, set[str]] = {}
    for treatment in TREATMENTS:
        directory = CLUSTER / treatment
        exit_code = (directory / "exit_code.txt").read_text().strip()
        if exit_code != "0":
            raise ValueError(f"{treatment}: exit_code.txt is {exit_code!r}")
        tree_path = directory / "tree.annotated.nwk"
        tree = parse_newick(tree_path.read_text())
        taxa = leaf_names(tree)
        if taxa != expected:
            raise ValueError(f"{treatment}: taxon set differs from jarvis48_taxa.tsv")
        split_set = split_keys(tree, taxa)
        rows = branch_rows(treatment)
        if split_set != set(rows):
            raise ValueError(f"{treatment}: annotated tree and branch_lengths.tsv split sets differ")
        for row in rows.values():
            for field in ("CU_length", "SU_length", "localPP", "q1", "q2", "q3"):
                if row[field] != "NA":
                    float(row[field])
            q = [float(row[f]) for f in ("q1", "q2", "q3") if row[f] != "NA"]
            if q and abs(sum(q) - 1.0) > 2e-5:
                raise ValueError(f"{treatment}: quartet frequencies do not sum to one")
        provenance = json.loads((directory / "provenance.json").read_text())
        run = json.loads((directory / "run.json").read_text())
        summaries[treatment] = {
            "taxa": len(taxa), "internal_splits": len(split_set),
            "input": provenance.get("input", ""),
            "input_sha256": provenance.get("input_sha256", ""),
            "annotated_sha256": sha256(tree_path),
            "runtime_seconds": run.get("inference_wall_seconds", run.get("elapsed_seconds", run.get("runtime_seconds", "NA"))),
            "branch_rows": len(rows),
        }
        splits[treatment] = split_set
    pairwise: dict[tuple[str, str], int] = {}
    for i, left in enumerate(TREATMENTS):
        for right in TREATMENTS[i + 1:]:
            pairwise[(left, right)] = len(splits[left] ^ splits[right])
    return summaries, pairwise


def write_topology_reports(summaries: dict[str, dict[str, object]], pairwise: dict[tuple[str, str], int]) -> None:
    out = RESULTS / "j48_topology_comparison.tsv"
    fields = ["left", "right", "rf_symmetric_difference", "left_internal_splits", "right_internal_splits"]
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for (left, right), rf in pairwise.items():
            writer.writerow({"left": left, "right": right, "rf_symmetric_difference": rf,
                             "left_internal_splits": summaries[left]["internal_splits"],
                             "right_internal_splits": summaries[right]["internal_splits"]})
    t0 = split_keys(parse_newick((CLUSTER / "J48_T0" / "tree.annotated.nwk").read_text()), read_taxa())
    groups = []
    for treatment in TREATMENTS[1:]:
        current = split_keys(parse_newick((CLUSTER / treatment / "tree.annotated.nwk").read_text()), read_taxa())
        groups.append((treatment, sorted(t0 - current), sorted(current - t0)))
    lines = ["# Completed J48 topology comparison", "",
             "The RF value below is the symmetric difference of canonical unrooted internal splits.",
             "The arbitrary Newick root is not used as a biological root.", "",
             "| Treatment | Taxa | Internal splits | Input SHA-256 | Runtime |", "|---|---:|---:|---|---:|"]
    for treatment in TREATMENTS:
        s = summaries[treatment]
        lines.append(f"| {treatment} | {s['taxa']} | {s['internal_splits']} | `{s['input_sha256']}` | {s['runtime_seconds']} |")
    lines += ["", "## Topology groups", "", "`J48_T0` and `J48_T_STRUCT_PRIMARY` share one unrooted topology; `J48_T_PNAS`, `J48_T_BLOCK_500K`, and `J48_NONCHR4` share the other. Their RF distance is 6, i.e. three internal splits in each direction.", "", "## Changed splits relative to J48_T0", ""]
    for treatment, removed, added in groups:
        if removed or added:
            lines.append(f"### {treatment}")
            lines.append(f"- Removed from T0: {len(removed)} split(s).")
            lines.append(f"- Added: {len(added)} split(s).")
            for split in removed:
                lines.append(f"  - removed `{split}`")
            for split in added:
                lines.append(f"  - added `{split}`")
    (RESULTS / "j48_topology_comparison.md").write_text("\n".join(lines) + "\n")


def focal_rows() -> list[dict[str, str]]:
    reference = {row["focal_relationship"]: row for row in read_tsv(RESULTS / "t0_focal_branches.tsv")}
    all_taxa = read_taxa()
    result: list[dict[str, str]] = []
    for treatment in TREATMENTS:
        rows = branch_rows(treatment)
        for focal in FOCALS:
            ref = reference[focal]
            projected = set(ref["split_key"].split("|")) & all_taxa
            other = all_taxa - projected
            projected_key = "|".join(min(tuple(sorted(projected)), tuple(sorted(other))))
            row = rows.get(projected_key)
            mapping_status = "exact_projected_reference_split" if row else "reference_split_absent_or_363_taxon_scope"
            result.append({"treatment": treatment, "focal": focal, "reference_state": ref["t0_state"],
                           "reference_split_key": ref["split_key"], "observed_split_key": projected_key if row else "NA",
                           "candidate_edge_count": "1" if row else "0", "mapping_status": mapping_status,
                           "branch_present": str(bool(row)).lower(), "state": ref["t0_state"] if row else "absent",
                           "CU_length": row["CU_length"] if row else "NA",
                           "SU_length": row["SU_length"] if row else "NA",
                           "localPP": row["localPP"] if row else "NA",
                           "q1": row["q1"] if row else "NA", "q2": row["q2"] if row else "NA",
                           "q3": row["q3"] if row else "NA"})
    return result


def write_focal_reports() -> None:
    rows = focal_rows()
    fields = list(rows[0])
    with (RESULTS / "j48_focal_branch_comparison.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    lines = ["# J48 focal branch comparison", "", "Focal branches are matched by the exact T0 canonical unrooted split. Missing rows are reported as absent; no replacement edge is substituted.", "", "| Focal | Treatment | Present | CU | SU | localPP | q1 | q2 | q3 |", "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for row in rows:
        lines.append("| {focal} | {treatment} | {branch_present} | {CU_length} | {SU_length} | {localPP} | {q1} | {q2} | {q3} |".format(**row))
    (RESULTS / "j48_focal_branch_comparison.md").write_text("\n".join(lines) + "\n")


def aggregate_qqs() -> list[dict[str, object]]:
    manifest = {row["original_gene_id"]: row for row in read_tsv(MANIFEST)}
    values: dict[tuple[str, str], list[list[float]]] = {}
    for row in read_tsv(QQS):
        locus = manifest.get(row["Gene"])
        if not locus:
            continue
        partition = "affected_chr4" if boolish(locus["is_chr4"]) else "background_nonchr4"
        key = (row["clade"], partition)
        values.setdefault(key, []).append([float(row["q1_species_history"]), float(row["q2_alternative_1"]), float(row["q3_alternative_2"])])
    result = []
    for clade in FOCALS:
        for partition in ("background_nonchr4", "affected_chr4"):
            observations = values.get((clade, partition), [])
            mean = [sum(v[i] for v in observations) / len(observations) for i in range(3)] if observations else [math.nan] * 3
            result.append({"clade": clade, "partition": partition, "n_loci": len(observations), "q1": mean[0], "q2": mean[1], "q3": mean[2], "estimand": "normalized_QQS_locus_mean"})
    return result


def crossover(qb: list[float], qr: list[float], i: int, j: int) -> dict[str, object]:
    if any(math.isnan(x) for x in qb + qr):
        return {"status": "missing_endpoint", "epsilon": "NA", "pair": f"q{i+1}=q{j+1}"}
    numerator = qb[i] - qb[j]
    denominator = numerator + qr[j] - qr[i]
    if abs(denominator) < 1e-12:
        return {"status": "denominator_zero", "epsilon": "NA", "pair": f"q{i+1}=q{j+1}"}
    eps = numerator / denominator
    status = "valid_interior" if 0.0 < eps < 1.0 else "outside_[0,1]"
    if abs(eps) < 1e-12 or abs(eps - 1) < 1e-12:
        status = "endpoint"
    return {"status": status, "epsilon": eps, "pair": f"q{i+1}=q{j+1}"}


def write_mixture_reports() -> None:
    rows = aggregate_qqs()
    by = {(r["clade"], r["partition"]): r for r in rows}
    output = RESULTS / "j48_quartet_mixture.tsv"
    fields = ["clade", "background_n_loci", "affected_n_loci", "q1_B", "q2_B", "q3_B", "q1_R", "q2_R", "q3_R", "sample_weight_R", "estimand", "pair", "epsilon_star", "status", "third_topology_at_epsilon"]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader()
        for clade in FOCALS:
            b = by[(clade, "background_nonchr4")]; r = by[(clade, "affected_chr4")]
            qb = [b[f"q{i}"] for i in (1, 2, 3)]; qr = [r[f"q{i}"] for i in (1, 2, 3)]
            sample_weight = r["n_loci"] / (b["n_loci"] + r["n_loci"]) if b["n_loci"] + r["n_loci"] else math.nan
            for i, j in ((0, 1), (0, 2), (1, 2)):
                c = crossover(qb, qr, i, j); third = "NA"
                if isinstance(c["epsilon"], float) and 0 < c["epsilon"] < 1:
                    third_i = ({0, 1, 2} - {i, j}).pop()
                    third = (1 - c["epsilon"]) * qb[third_i] + c["epsilon"] * qr[third_i]
                writer.writerow({"clade": clade, "background_n_loci": b["n_loci"], "affected_n_loci": r["n_loci"],
                                 "q1_B": qb[0], "q2_B": qb[1], "q3_B": qb[2], "q1_R": qr[0], "q2_R": qr[1], "q3_R": qr[2],
                                 "sample_weight_R": sample_weight, "estimand": b["estimand"], "pair": c["pair"],
                                 "epsilon_star": c["epsilon"], "status": c["status"], "third_topology_at_epsilon": third})
    manifest_rows = list(read_tsv(MANIFEST))[1:] if False else read_tsv(MANIFEST)
    chr4 = sum(boolish(r["is_chr4"]) for r in manifest_rows)
    total = len(manifest_rows)
    physical_span = max((int(r["end"]) for r in manifest_rows if r["chromosome"] == "chr4"), default=0)
    summary = RESULTS / "j48_threshold_summary.tsv"
    with summary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["quantity", "value", "definition"], delimiter="\t"); writer.writeheader()
        writer.writerow({"quantity": "chr4_locus_fraction", "value": chr4 / total, "definition": "6262/63430 locus-manifest fraction; distinct from epsilon"})
        writer.writerow({"quantity": "chr4_locus_count", "value": chr4, "definition": "stage4d locus manifest"})
        writer.writerow({"quantity": "total_locus_count", "value": total, "definition": "stage4d locus manifest"})
        writer.writerow({"quantity": "chr4_coordinate_extent_end", "value": physical_span, "definition": "maximum chr4 end in locus manifest; not a genome-wide affected fraction"})
        writer.writerow({"quantity": "astral_transition", "value": "pending_bridges2", "definition": "no new inference was run locally"})
    lines = ["# J48 quartet mixture and threshold analysis", "", "The quartet summary uses the existing `normalized_QQS` window summary, averaged across available loci. It is not ASTRAL branch-quartet frequency and does not establish an MSRC mechanism.", "", f"The locus manifest contains {total} loci, including {chr4} chromosome-4 loci ({chr4 / total:.6f}). This sampled-locus fraction is reported separately from the mixture weight epsilon.", "", "For each focal quartet, `Q(epsilon)=(1-epsilon)Q_B+epsilon Q_R`; pairwise crossover values are in `j48_quartet_mixture.tsv`. A valid interior pairwise equality is only a candidate quartet crossover: the third topology must be checked, and it is not a species-tree transition.", "", "The current local deliverables do not claim an empirical full-tree transition. The coarse mixture grid and matched controls are prepared for Bridges-2."]
    (RESULTS / "j48_threshold_analysis.md").write_text("\n".join(lines) + "\n")


def write_experiment_manifest() -> None:
    path = RESULTS / "j48_threshold_experiment_manifest.tsv"
    fields = ["mixture_id", "kind", "epsilon_requested", "seed", "source_background", "source_affected", "target_records", "status", "command"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader()
        writer.writerow({"mixture_id": "J48_BASELINE_T0", "kind": "exact_completed_baseline", "epsilon_requested": "observed_manifest_fraction", "seed": "NA",
                         "source_background": "J48_T0", "source_affected": "J48_T0", "target_records": 62945,
                         "status": "reference_only", "command": "use existing empirical/neoaves_chr4/data/stage4e/J48_T0.tre"})
        for index, epsilon in enumerate(EPSILONS):
            mixture = f"J48_MIX_E{epsilon:.2f}".replace(".", "p")
            source = "background_stratum" if epsilon == 0.0 else ("affected_chr4_stratum" if epsilon == 1.0 else "streaming_stratified_resample")
            writer.writerow({"mixture_id": mixture, "kind": "coarse_grid", "epsilon_requested": f"{epsilon:.2f}", "seed": 74000 + index,
                             "source_background": "J48_NONCHR4", "source_affected": "J48_T0_chr4_stratum", "target_records": 62945,
                             "status": "prepared_not_run", "command": f"python scripts/05_stage4e_prepare_threshold_mixtures.py --mixture-id {mixture} --epsilon {epsilon:.2f} --seed {74000 + index}"})
        for replicate in range(1, 6):
            writer.writerow({"mixture_id": f"J48_CTRL_RANDOM_{replicate:02d}", "kind": "random_genome_deletion_control", "epsilon_requested": "NA", "seed": 75000 + replicate,
                             "source_background": "J48_T0", "source_affected": "NA", "target_records": "matched", "status": "prepared_not_run",
                             "command": f"python scripts/05_stage4e_prepare_threshold_mixtures.py --control random_genome_deletion --replicate {replicate} --seed {75000 + replicate}"})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare-only", action="store_true", help="write reports without reading large gene-tree collections")
    args = parser.parse_args()
    _ = args
    summaries, pairwise = validate_results()
    write_topology_reports(summaries, pairwise)
    write_focal_reports()
    write_mixture_reports()
    write_experiment_manifest()
    print("validated five completed J48 trees; wrote bounded reports and Bridges-2 manifest")


if __name__ == "__main__":
    main()
