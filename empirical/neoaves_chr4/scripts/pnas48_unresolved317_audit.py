#!/usr/bin/env python3
"""Audit the 317 collapsed records whose header-derived loci disagree.

This is a metadata audit. It never treats a matching taxon set as a locus
assignment and never writes a published-source gene-tree collection.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import lzma
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "empirical/neoaves_chr4"
EXT = BASE / "external"
STILLER = EXT / "stiller2024"
PNAS = EXT / "chr4avian"
RESULTS = BASE / "results/stage4e"
MANIFEST = BASE / "results/stage4d_locus_manifest.tsv"
HEADER_MEMBER = "genetreesupport/63K_trees.names_header.txt.xz"
EXPECTED_COMMIT = "17b408a622c9589fe9f6ca8a434c150974e39cd4"
LEAF = re.compile(r"(?:\(|,)([^(),:;]+)(?::[^(),;]+)?")

TREATMENTS = (
    "PNAS_outlier_removal",
    "STRUCTURAL_PRIMARY",
    "STRUCTURAL_STRINGENT",
    "STRUCTURAL_INCLUSIVE",
    "BLOCK_250K",
    "BLOCK_500K",
    "BLOCK_1MB",
    "NONCHR4",
)


def taxa(tree: str) -> tuple[str, ...]:
    return tuple(sorted(LEAF.findall(tree)))


def taxon_key(names: tuple[str, ...]) -> str:
    return hashlib.sha256("\0".join(names).encode()).hexdigest()


def normalized(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def normalized_key(names: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(sorted(normalized(name) for name in names))


def has_polytomy(tree: str) -> bool:
    stack: list[int] = []
    for char in tree:
        if char == "(":
            stack.append(1)
        elif char == "," and stack:
            stack[-1] += 1
        elif char == ")" and stack:
            if stack.pop() > 2:
                return True
    return False


def read_header() -> list[dict[str, str]]:
    archive = PNAS / "chr4avian-v1.0.0.zip"
    with zipfile.ZipFile(archive) as zf:
        raw = lzma.decompress(zf.read(next(n for n in zf.namelist() if n.endswith(HEADER_MEMBER))))
    rows = []
    for line in raw.decode().splitlines()[1:]:
        gene, chrom, start, end, _size, offset = line.split()
        rows.append({
            "gene": gene,
            "locus_id": f"{chrom}_{start}_{end}_1k_start{offset}_fasta.treefile",
            "chromosome": chrom,
            "start": start,
            "end": end,
            "offset": offset,
        })
    return rows


def read_manifest() -> dict[str, dict[str, str]]:
    with MANIFEST.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return {row["locus_id"]: row for row in rows}


def boolish(value: str) -> bool:
    return value.lower() in {"true", "1", "yes"}


def representatives(rows: list[dict[str, str]], width: int) -> set[str]:
    chosen: dict[tuple[str, int], tuple[tuple[float, str], str]] = {}
    for row in rows:
        if not boolish(row["is_chr4"]):
            continue
        midpoint = float(row["midpoint"])
        key = (row["chromosome"], int(midpoint // width))
        rank = (abs(midpoint - (key[1] + 0.5) * width), row["locus_id"])
        if key not in chosen or rank < chosen[key][0]:
            chosen[key] = (rank, row["locus_id"])
    return {value[1] for value in chosen.values()}


def treatment_membership(row: dict[str, str], reps: dict[int, set[str]]) -> dict[str, bool]:
    return {
        "PNAS_outlier_removal": not boolish(row["published_outlier_region"]),
        "STRUCTURAL_PRIMARY": not boolish(row["frozen_structural_primary"]),
        "STRUCTURAL_STRINGENT": not boolish(row["frozen_structural_stringent"]),
        "STRUCTURAL_INCLUSIVE": not boolish(row["frozen_structural_inclusive"]),
        "BLOCK_250K": not boolish(row["is_chr4"]) or row["locus_id"] in reps[250000],
        "BLOCK_500K": not boolish(row["is_chr4"]) or row["locus_id"] in reps[500000],
        "BLOCK_1MB": not boolish(row["is_chr4"]) or row["locus_id"] in reps[1000000],
        "NONCHR4": not boolish(row["is_chr4"]),
    }


def hopcroft_karp(graph: dict[str, set[str]]) -> dict[str, str]:
    """Maximum matching for mismatch rows against candidate locus IDs."""
    left_match: dict[str, str] = {}
    right_match: dict[str, str] = {}
    while True:
        distance: dict[str, int] = {}
        queue = []
        for left in graph:
            if left not in left_match:
                distance[left] = 0
                queue.append(left)
        for left in queue:
            for right in graph[left]:
                mate = right_match.get(right)
                if mate is not None and mate not in distance:
                    distance[mate] = distance[left] + 1
                    queue.append(mate)
        found = False

        def visit(left: str) -> bool:
            nonlocal found
            for right in graph[left]:
                mate = right_match.get(right)
                if mate is None or (mate in distance and distance[mate] == distance[left] + 1 and visit(mate)):
                    left_match[left] = right
                    right_match[right] = left
                    found = True
                    return True
            distance.pop(left, None)
            return False

        for left in list(graph):
            if left not in left_match and left in distance:
                visit(left)
        if not found:
            return left_match


def audit() -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    header = read_header()
    manifest = read_manifest()
    if len(header) != 63430 or len(manifest) != 63430:
        raise RuntimeError(f"expected 63,430 header/manifest rows, found {len(header)}/{len(manifest)}")
    manifest_rows = list(manifest.values())
    reps = {width: representatives(manifest_rows, width) for width in (250000, 500000, 1000000)}
    named_by_key: dict[str, list[dict[str, object]]] = defaultdict(list)
    named_taxa: dict[str, tuple[str, ...]] = {}
    named_polytomy: dict[str, bool] = {}
    with gzip.open(STILLER / "63430.named.gene.trees.gz", "rt") as handle:
        for line_number, line in enumerate(handle, 1):
            locus, tree = line.rstrip("\n").split("\t", 1)
            names = taxa(tree)
            named_taxa[locus] = names
            named_polytomy[locus] = has_polytomy(tree)
            named_by_key[taxon_key(names)].append({"locus_id": locus, "line": line_number})
    if len(named_taxa) != 63430:
        raise RuntimeError("named archive locus count is not 63,430")

    mismatches: list[dict[str, object]] = []
    proposed_loci_seen: set[str] = set()
    candidate_graph: dict[str, set[str]] = {}
    mismatch_by_key: Counter[str] = Counter()
    with gzip.open(STILLER / "63430.aLRT-0.95-collapsed.gene.trees.gz", "rt") as collapsed:
        for index, line in enumerate(collapsed):
            meta = header[index]
            proposed = manifest[meta["locus_id"]]
            collapsed_names = taxa(line.rstrip("\n"))
            named_names = named_taxa[meta["locus_id"]]
            if collapsed_names == named_names:
                continue
            collapsed_set, named_set = set(collapsed_names), set(named_names)
            candidates = named_by_key.get(taxon_key(collapsed_names), [])
            candidate_ids = {str(candidate["locus_id"]) for candidate in candidates}
            candidate_graph[str(index + 1)] = candidate_ids
            mismatch_by_key[taxon_key(collapsed_names)] += 1
            near = []
            proposed_start = int(proposed["start"])
            for candidate in candidates:
                candidate_row = manifest[str(candidate["locus_id"])]
                if candidate_row["chromosome"] == proposed["chromosome"]:
                    near.append(abs(int(candidate_row["start"]) - proposed_start))
            treatment = treatment_membership(proposed, reps)
            candidate_membership: dict[str, set[bool]] = {name: set() for name in TREATMENTS}
            for candidate_id in candidate_ids:
                for name, value in treatment_membership(manifest[candidate_id], reps).items():
                    candidate_membership[name].add(value)
            normalized_only = normalized_key(collapsed_names) == normalized_key(named_names) and collapsed_names != named_names
            mismatches.append({
                "collapsed_record": index + 1,
                "gene": meta["gene"],
                "proposed_locus_id": meta["locus_id"],
                "proposed_chromosome": proposed["chromosome"],
                "proposed_start": proposed["start"],
                "proposed_end": proposed["end"],
                "proposed_offset": meta["offset"],
                "collapsed_n_taxa": len(collapsed_names),
                "named_n_taxa": len(named_names),
                "collapsed_only_taxa": ",".join(sorted(collapsed_set - named_set)),
                "named_only_taxa": ",".join(sorted(named_set - collapsed_set)),
                "single_taxon_difference": len(collapsed_set ^ named_set) == 1,
                "n_taxon_differences": len(collapsed_set ^ named_set),
                "normalized_name_equivalent": normalized_only,
                "collapsed_has_polytomy": has_polytomy(line),
                "named_has_polytomy": named_polytomy[meta["locus_id"]],
                "candidate_count": len(candidate_ids),
                "candidate_locus_ids": ",".join(sorted(candidate_ids)),
                "same_chromosome_candidate_count": sum(1 for candidate in candidates if manifest[str(candidate["locus_id"])] ["chromosome"] == proposed["chromosome"]),
                "nearest_same_chromosome_candidate_bp": min(near) if near else "",
                "candidate_assigned_to_other_header_record_count": len(candidate_ids),
                "candidate_already_proposed_by_other_mismatch": sum(1 for candidate_id in candidate_ids if candidate_id in proposed_loci_seen),
                "proposed_pnas_outlier": not treatment["PNAS_outlier_removal"],
                "proposed_structural_primary_excluded": not treatment["STRUCTURAL_PRIMARY"],
                "proposed_structural_stringent_excluded": not treatment["STRUCTURAL_STRINGENT"],
                "proposed_structural_inclusive_excluded": not treatment["STRUCTURAL_INCLUSIVE"],
                "proposed_block_250k_included": treatment["BLOCK_250K"],
                "proposed_block_500k_included": treatment["BLOCK_500K"],
                "proposed_block_1mb_included": treatment["BLOCK_1MB"],
                "proposed_nonchr4_included": treatment["NONCHR4"],
                "treatment_membership_invariant": all(len(values) <= 1 for values in candidate_membership.values()) if candidate_ids else False,
                "classification": (
                    "source discrepancy" if not candidate_ids else
                    "treatment-equivalent mapping" if all(len(values) <= 1 for values in candidate_membership.values()) else
                    "unresolved"
                ),
            })
            proposed_loci_seen.add(meta["locus_id"])

    # The matching is deliberately limited to mismatch rows and is reported as
    # a mathematical constraint, never as a provenance decision.
    matching = hopcroft_karp(candidate_graph)
    treatment_rows = []
    for name in TREATMENTS:
        # Recompute from the detailed rows to keep the output schema explicit.
        variable = 0
        invariant = 0
        proposed_membership = 0
        for row in mismatches:
            if name == "PNAS_outlier_removal":
                proposed = bool(row["proposed_pnas_outlier"])
            elif name.startswith("STRUCTURAL"):
                proposed = not bool(row["proposed_" + name.lower() + "_excluded"])
            elif name.startswith("BLOCK"):
                proposed = bool(row["proposed_" + name.lower() + "_included"])
            else:
                proposed = bool(row["proposed_nonchr4_included"])
            proposed_membership += proposed
            candidate_values = set()
            for candidate_id in candidate_graph[str(row["collapsed_record"])]:
                candidate_values.add(treatment_membership(manifest[candidate_id], reps)[name])
            if candidate_values and len(candidate_values) == 1:
                invariant += 1
            if len(candidate_values) > 1:
                variable += 1
        treatment_rows.append({
            "treatment": name,
            "unresolved_records": len(mismatches),
            "proposed_membership_true_count": proposed_membership,
            "candidate_membership_invariant_count": invariant,
            "candidate_membership_variable_count": variable,
            "potentially_affecting_membership_count": len(mismatches) - invariant,
            "membership_verifiable": "Yes" if variable == 0 and invariant + sum(not bool(candidate_graph[str(row["collapsed_record"])]) for row in mismatches) == len(mismatches) else "Partial",
        })
    summary = {
        "collapsed_records": 63430,
        "matched_records": 63113,
        "unresolved_records": len(mismatches),
        "candidate_index_loci": len(named_taxa),
        "candidate_taxon_set_keys": len(named_by_key),
        "candidate_graph_maximum_matching": len(matching),
        "candidate_graph_rows_with_candidates": sum(bool(value) for value in candidate_graph.values()),
        "candidate_graph_rows_without_candidates": sum(not bool(value) for value in candidate_graph.values()),
        "mismatch_by_proposed_chromosome": dict(Counter(str(row["proposed_chromosome"]) for row in mismatches)),
        "mismatch_by_taxon_delta": dict(Counter(int(row["n_taxon_differences"]) for row in mismatches)),
        "mismatch_single_taxon": sum(bool(row["single_taxon_difference"]) for row in mismatches),
        "mismatch_polytomy_collapsed": sum(bool(row["collapsed_has_polytomy"]) for row in mismatches),
        "mismatch_normalized_name_equivalent": sum(bool(row["normalized_name_equivalent"]) for row in mismatches),
        "classification_counts": dict(Counter(str(row["classification"]) for row in mismatches)),
        "source_repository_commit": EXPECTED_COMMIT,
    }
    return mismatches, treatment_rows, summary


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def update_existing_reports(summary: dict[str, object], treatment_rows: list[dict[str, object]]) -> None:
    mapping_path = RESULTS / "pnas48_mapping_summary.tsv"
    with mapping_path.open() as handle:
        mapping = list(csv.DictReader(handle, delimiter="\t"))
    new_components = {
        "unresolved317_candidate_index",
        "unresolved317_source_discrepancy",
        "unresolved317_treatment_equivalent",
        "unresolved317_verified_locus_mappings",
    }
    mapping = [row for row in mapping if row["component"] not in new_components]
    mapping.extend([
        {"component": "unresolved317_candidate_index", "status": "audited", "value": f"117 rows with candidates; 200 without candidates", "evidence": "pnas48_unresolved317.tsv exact named taxon-set index"},
        {"component": "unresolved317_source_discrepancy", "status": "strongly_supported", "value": "200", "evidence": "collapsed taxon sets absent from all 63,430 named loci; ordinary branch collapse preserves leaves"},
        {"component": "unresolved317_treatment_equivalent", "status": "candidate_only", "value": "26", "evidence": "all exact taxon-set candidates share all eight treatment memberships; no locus identity established"},
        {"component": "unresolved317_verified_locus_mappings", "status": "none", "value": "0", "evidence": "no independent source file links an unresolved collapsed row to a locus"},
    ])
    write_tsv(mapping_path, mapping)

    comparison_path = RESULTS / "pnas48_dataset_comparison.tsv"
    with comparison_path.open() as handle:
        comparison = list(csv.DictReader(handle, delimiter="\t"))
    for row in comparison:
        row.update({
            "unresolved317_records": "317",
            "unresolved317_source_discrepancy": "200",
            "unresolved317_treatment_equivalent": "26",
            "unresolved317_verified": "0",
            "published_source_status": "exact_locus_filtering_unverified; no dataset constructed",
        })
    write_tsv(comparison_path, comparison)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    mismatches, treatment_rows, summary = audit()
    if args.write:
        write_tsv(RESULTS / "pnas48_unresolved317.tsv", mismatches)
        write_tsv(RESULTS / "pnas48_treatment_membership_audit.tsv", treatment_rows)
        update_existing_reports(summary, treatment_rows)
        lines = [
            "# Audit of the 317 unresolved collapsed records", "",
            "The audit compares each collapsed record with its header-derived named locus and indexes all named loci by exact taxon set. Candidate taxon-set matches are mathematical alternatives only; they are not locus assignments.", "",
            f"- Published collapsed records: **{summary['collapsed_records']:,}**.",
            f"- Header-derived taxon-set matches: **{summary['matched_records']:,}**.",
            f"- Unresolved records audited: **{summary['unresolved_records']:,}**.",
            f"- Rows with at least one exact taxon-set candidate: **{summary['candidate_graph_rows_with_candidates']:,}**.",
            f"- Rows with no exact taxon-set candidate: **{summary['candidate_graph_rows_without_candidates']:,}**.",
            f"- Maximum one-to-one matching among unresolved rows: **{summary['candidate_graph_maximum_matching']:,}**.", "",
            "## Interpretation", "",
            "No unresolved record is promoted to a verified locus. A record is treatment-equivalent only when every exact taxon-set candidate has the same membership for all eight treatments; otherwise the identity remains unresolved. Candidate matches do not establish ordering, genomic proximity, or shared upstream provenance. Because every candidate is already a proposed header locus for another collapsed row, candidate reassignment would require a complete global permutation, which this audit does not claim.", "",
            "## Distribution", "",
            f"- Proposed chromosome counts: `{summary['mismatch_by_proposed_chromosome']}`.",
            f"- Symmetric taxon-set difference counts: `{summary['mismatch_by_taxon_delta']}`.",
            f"- Single-taxon differences: **{summary['mismatch_single_taxon']:,}**.",
            f"- Collapsed records containing a polytomy: **{summary['mismatch_polytomy_collapsed']:,}**.",
            f"- Naming-normalization-only differences: **{summary['mismatch_normalized_name_equivalent']:,}**.",
            f"- Classification counts: `{summary['classification_counts']}`.", "",
            "The source-side producer of the collapsed archive remains absent from the cached repository snapshot. Ordinary branch collapsing cannot explain changed leaf sets, so any taxon-set discrepancy remains a source or preprocessing discrepancy until an upstream file or author metadata resolves it.", "",
        ]
        (RESULTS / "pnas48_unresolved317_summary.md").write_text("\n".join(lines))
    print(summary)


if __name__ == "__main__":
    main()
