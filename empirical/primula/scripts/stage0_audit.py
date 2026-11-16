#!/usr/bin/env python3
"""Audit deposited Primula Newick trees without inventing genomic coordinates."""
from __future__ import annotations

import csv
import sys
import statistics
from collections import Counter
import re
from pathlib import Path
from typing import Iterable

try:
    from Bio import Phylo
    from io import StringIO
except ImportError as exc:  # pragma: no cover - environment diagnostic
    raise SystemExit("Biopython is required: python3 -m pip install biopython") from exc

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data/primula/raw/dryad_2022"
RESULTS = ROOT / "empirical/primula/results"
TREE_FILE = RAW / "PrimPrim_phylogeny_allGF.trees"
PARTITION_PATTERN = re.compile(r"DNA,\s*.+?_(?P<chrom>\d+):(?P<start>\d+)-(?P<end>\d+)\.min4\.fasta\s*=\s*(?P<site_start>\d+)-(?P<site_end>\d+)")


def split_newick_records(text: str) -> list[str]:
    """Split concatenated Newick trees at unquoted, top-level semicolons."""
    records: list[str] = []
    start = 0
    quote = False
    bracket_depth = 0
    for i, char in enumerate(text):
        if char == "'":
            quote = not quote
        elif not quote:
            if char == "[":
                bracket_depth += 1
            elif char == "]":
                bracket_depth = max(0, bracket_depth - 1)
            elif char == ";" and bracket_depth == 0:
                record = text[start:i + 1].strip()
                if record:
                    records.append(record)
                start = i + 1
    tail = text[start:].strip()
    if tail:
        records.append(tail)
    return records


def inspect_newick(record: str) -> dict[str, object]:
    tree = Phylo.read(StringIO(record), "newick")
    tips = [tip.name for tip in tree.get_terminals()]
    if any(name is None or name == "" for name in tips):
        raise ValueError("unnamed terminal")
    lengths = [clade.branch_length for clade in tree.find_clades(order="preorder") if clade != tree.root]
    present = sum(x is not None for x in lengths)
    return {
        "taxa": tuple(sorted(tips)),
        "n_taxa": len(tips),
        "branch_lengths_present": present > 0,
        "branch_lengths_complete": bool(lengths) and present == len(lengths),
        "polytomies": sum(len(clade.clades) > 2 for clade in tree.get_nonterminals() if clade != tree.root),
        "root_degree": len(tree.root.clades),
        "branch_length_count": len(lengths),
    }


def reconcile_partition_counts(tree_count: int, partition_count: int) -> str:
    if tree_count != partition_count:
        return f"unresolved: {partition_count} FASTA partition blocks versus {tree_count} unlabeled Newick records ({partition_count-tree_count} count difference); no association established"
    return "counts agree; tree-record ordering still requires source verification"


def unresolved_mapping_rows(tree_count: int) -> list[dict[str, object]]:
    """Represent every Newick ordinal without inventing a coordinate or source ID."""
    return [{"tree_index": i, "source_tree_identifier": "NA",
             "chromosome": "NA", "start_bp": "NA", "end_bp": "NA",
             "mapping_evidence": "Newick record has no embedded ID; partition files map FASTA site blocks only",
             "validation_status": "unresolved_no_tree_coordinate_crosswalk"}
            for i in range(1, tree_count + 1)]


def parse_partition_files() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Read physical windows and alignment ranges without assigning tree IDs."""
    rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    for chrom in range(1, 12):
        path = RAW / f"PrimPrim_phylogeny_chr{chrom:03d}.partition.txt"
        if not path.exists():
            continue
        chrom_rows = []
        for line_number, line in enumerate(path.read_text().splitlines(), 1):
            match = PARTITION_PATTERN.search(line)
            if not match:
                continue
            data = {k: int(v) for k, v in match.groupdict().items()}
            if data["chrom"] != chrom:
                raise ValueError(f"{path.name}:{line_number}: chromosome label mismatch")
            if data["end"] - data["start"] + 1 != 50_000:
                raise ValueError(f"{path.name}:{line_number}: not a 50-kb interval")
            row = {"chromosome": f"chr{chrom:03d}", "partition_row": len(chrom_rows) + 1,
                   "start_bp": data["start"], "end_bp": data["end"],
                   "alignment_site_start": data["site_start"], "alignment_site_end": data["site_end"],
                   "partition_file": path.name, "tree_index": "NA",
                   "mapping_status": "partition metadata only; tree-record association unresolved"}
            rows.append(row)
            chrom_rows.append(data)
        physical = sorted(chrom_rows, key=lambda x: x["start"])
        site_ranges = sorted(chrom_rows, key=lambda x: x["site_start"])
        gaps = [(b["start"] - a["start"]) for a, b in zip(physical, physical[1:])]
        summaries.append({"chromosome": f"chr{chrom:03d}", "partition_rows": len(chrom_rows),
                          "minimum_start_bp": min((r["start"] for r in chrom_rows), default="NA"),
                          "maximum_end_bp": max((r["end"] for r in chrom_rows), default="NA"),
                          "all_intervals_50kb": all(r["end"] - r["start"] + 1 == 50_000 for r in chrom_rows),
                          "alignment_ranges_contiguous": all(a["site_end"] + 1 == b["site_start"] for a, b in zip(site_ranges, site_ranges[1:])),
                          "physical_100kb_gaps": sum(g == 100_000 for g in gaps),
                          "physical_gaps_gt_100kb": sum(g > 100_000 for g in gaps),
                          "partition_line_order_is_physical_order": all(a["start"] < b["start"] for a, b in zip(chrom_rows, chrom_rows[1:]))})
    return rows, summaries


def validate_explicit_mapping(rows: Iterable[dict[str, str]], expected_trees: int) -> None:
    """Reject missing, duplicated, unordered, or invalid coordinates."""
    rows = list(rows)
    indices = [int(row["tree_index"]) for row in rows]
    if indices != list(range(1, expected_trees + 1)):
        raise ValueError("mapping must contain each 1-based tree index exactly once in order")
    for row in rows:
        if not row.get("chromosome") or row.get("mapping_status") != "verified":
            raise ValueError("each coordinate must have an explicit verified status")
        start, end = int(row["start_bp"]), int(row["end_bp"])
        if start < 1 or end < start:
            raise ValueError("invalid physical coordinates")


def quartets(tree, quartet: tuple[str, str, str, str]) -> tuple[tuple[str, str], tuple[str, str]] | None:
    """Return the unrooted split for a fully resolved quartet, else None."""
    q = set(quartet)
    present = {x.name for x in tree.get_terminals()}
    if not q <= present:
        return None
    # Find the edge whose removal separates two focal taxa on each side.
    splits = set()
    for clade in tree.find_clades(order="preorder"):
        if clade == tree.root:
            continue
        side = {x.name for x in clade.get_terminals()} & q
        if len(side) == 2:
            other = q - side
            splits.add(tuple(sorted((tuple(sorted(side)), tuple(sorted(other))))))
    return next(iter(splits)) if len(splits) == 1 else None


def patristic(tree, a: str, b: str) -> float:
    value = tree.distance(a, b)
    return float(value)


def main() -> int:
    if not TREE_FILE.exists():
        print(f"Stage 0 blocked: missing {TREE_FILE}", file=sys.stderr)
        print("Acquire the small Dryad tree and partition files; do not infer coordinates from tree order.", file=sys.stderr)
        return 2
    RESULTS.mkdir(parents=True, exist_ok=True)
    records = split_newick_records(TREE_FILE.read_text())
    counts: Counter[tuple[str, ...]] = Counter()
    inventory = []
    failures = []
    branch_lengths = []
    root_degrees: Counter[int] = Counter()
    parsed_taxa: dict[int, tuple[str, ...]] = {}
    for index, record in enumerate(records, 1):
        try:
            tree = Phylo.read(StringIO(record), "newick")
            tips = [tip.name for tip in tree.get_terminals()]
            if any(name is None or name == "" for name in tips):
                raise ValueError("unnamed terminal")
            taxa = tuple(sorted(tips))
            lengths = [clade.branch_length for clade in tree.find_clades(order="preorder") if clade != tree.root]
            present = [value for value in lengths if value is not None]
            branch_lengths.extend(present)
            root_degrees[len(tree.root.clades)] += 1
            polytomies = sum(len(clade.clades) > 2 for clade in tree.get_nonterminals() if clade != tree.root)
            counts[taxa] += 1
            parsed_taxa[index] = taxa
            inventory.append({"tree_index": index, "taxon_set_id": "", "n_taxa": len(taxa),
                              "branch_lengths_present": bool(present),
                              "branch_lengths_complete": bool(lengths) and len(present) == len(lengths),
                              "n_nonroot_edges": len(lengths), "root_degree": len(tree.root.clades),
                              "rooting": "Newick root placement; biological rooting convention not inferred",
                              "polytomies": polytomies, "parse_status": "usable", "notes": "none"})
        except Exception as exc:  # retain malformed records in audit
            failures.append((index, str(exc)))
            inventory.append({"tree_index": index, "taxon_set_id": "", "n_taxa": "NA",
                              "branch_lengths_present": "NA", "branch_lengths_complete": "NA",
                              "n_nonroot_edges": "NA", "root_degree": "NA", "rooting": "unknown",
                              "polytomies": "NA", "parse_status": "malformed", "notes": str(exc)})
    taxon_ids = {taxa: f"T{i:03d}" for i, taxa in enumerate(sorted(counts), 1)}
    for row in inventory:
        if row["parse_status"] == "usable":
            row["taxon_set_id"] = taxon_ids[parsed_taxa[int(row["tree_index"])]]
    out = RESULTS / "gene_tree_inventory.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(inventory[0]) if inventory else ["tree_index"], delimiter="\t")
        writer.writeheader()
        writer.writerows(inventory)
    (RESULTS / "tree_taxon_sets.tsv").write_text("taxon_set_id\tn_taxa\tn_trees\ttaxa\n" + "".join(
        f"{taxon_ids[taxa]}\t{len(taxa)}\t{n}\t{','.join(taxa)}\n" for taxa, n in sorted(counts.items())))
    complete = sum(row["branch_lengths_complete"] is True for row in inventory)
    usable = len(records) - len(failures)
    summary = [
        ("records", len(records), "Newick records split from the deposited tree file"),
        ("usable_trees", usable, "Trees parsed by Biopython"),
        ("malformed_trees", len(failures), "Records not parsed as Newick"),
        ("taxon_sets", len(counts), "Distinct terminal-label sets"),
        ("taxa_per_tree", min((r["n_taxa"] for r in inventory if r["parse_status"] == "usable"), default="NA"), "minimum across usable trees"),
        ("taxa_per_tree_max", max((r["n_taxa"] for r in inventory if r["parse_status"] == "usable"), default="NA"), "maximum across usable trees"),
        ("complete_branch_length_trees", complete, "All nonroot edges have numeric branch lengths"),
        ("branch_length_edges", len(branch_lengths), "Nonroot branch-length observations; substitution-distance scale"),
        ("branch_length_min", min(branch_lengths) if branch_lengths else "NA", "IQ-TREE fitted branch-length scale"),
        ("branch_length_median", statistics.median(branch_lengths) if branch_lengths else "NA", "IQ-TREE fitted branch-length scale"),
        ("branch_length_mean", statistics.mean(branch_lengths) if branch_lengths else "NA", "IQ-TREE fitted branch-length scale"),
        ("branch_length_max", max(branch_lengths) if branch_lengths else "NA", "IQ-TREE fitted branch-length scale"),
        ("zero_length_edges", sum(value == 0 for value in branch_lengths), "Zero-valued branch lengths"),
        ("negative_length_edges", sum(value < 0 for value in branch_lengths), "Negative-valued branch lengths"),
        ("unrooted_root_degree_3_trees", root_degrees.get(3, 0), "Newick representation has 3-way root node; treated as unrooted, not a polytomy"),
        ("nonroot_polytomy_trees", sum(bool(r["polytomies"]) for r in inventory if r["parse_status"] == "usable"), "Trees with a nonroot node of degree >2"),
    ]
    with (RESULTS / "tree_quality_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value", "interpretation"])
        writer.writerows(summary)
    partition_rows, partition_summaries = parse_partition_files()
    if partition_rows:
        fields = list(partition_rows[0])
        with (RESULTS / "partition_window_inventory.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(partition_rows)
        with (RESULTS / "partition_coordinate_summary.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(partition_summaries[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(partition_summaries)
    with (RESULTS / "mapping_reconciliation.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["tree_records", "partition_windows", "difference", "mapping_status"])
        writer.writerow([len(records), len(partition_rows), len(partition_rows) - len(records),
                         reconcile_partition_counts(len(records), len(partition_rows))])
    mapping = unresolved_mapping_rows(len(records))
    with (RESULTS / "tree_coordinate_mapping.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(mapping[0]) if mapping else ["tree_index"], delimiter="\t")
        writer.writeheader()
        writer.writerows(mapping)
    print(f"records={len(records)} usable={usable} malformed={len(failures)}")
    print(f"partition_windows={len(partition_rows)} tree_partition_difference={len(partition_rows)-len(records)}")
    print(f"taxon_sets={len(counts)} complete_branch_trees={complete} inventory={out}")
    if failures:
        print(f"malformed indices: {[i for i, _ in failures]}")
    print("Coordinate mapping is not inferred from tree order; validate explicit mapping separately.")
    for path in RESULTS.glob("*.tsv"):
        path.write_text(path.read_text().replace("\r\n", "\n"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
