#!/usr/bin/env python3
"""Audit deposited Primula Newick trees without inventing genomic coordinates."""
from __future__ import annotations

import csv
import sys
from collections import Counter
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
        "polytomies": sum(len(clade.clades) > 2 for clade in tree.get_nonterminals()),
        "root_degree": len(tree.root.clades),
        "branch_length_count": len(lengths),
    }


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
    text = TREE_FILE.read_text()
    records = split_newick_records(text)
    counts: Counter[tuple[str, ...]] = Counter()
    inventory = []
    failures = []
    for index, record in enumerate(records, 1):
        try:
            info = inspect_newick(record)
            counts[info["taxa"]] += 1
            inventory.append({"tree_index": index, "taxon_set_id": "", "n_taxa": info["n_taxa"],
                              "branch_lengths_present": info["branch_lengths_present"],
                              "branch_lengths_complete": info["branch_lengths_complete"],
                              "rooting": "Newick root placement; biological rooting convention not inferred",
                              "polytomies": info["polytomies"], "parse_status": "usable", "notes": ""})
        except Exception as exc:  # retain malformed records in audit
            failures.append((index, str(exc)))
            inventory.append({"tree_index": index, "taxon_set_id": "", "n_taxa": "NA",
                              "branch_lengths_present": "NA", "branch_lengths_complete": "NA",
                              "rooting": "unknown", "polytomies": "NA", "parse_status": "malformed", "notes": str(exc)})
    taxon_ids = {taxa: f"T{i:03d}" for i, taxa in enumerate(sorted(counts), 1)}
    for row in inventory:
        if row["parse_status"] == "usable":
            try:
                row["taxon_set_id"] = taxon_ids[inspect_newick(records[int(row["tree_index"]) - 1])["taxa"]]
            except Exception:
                pass
    out = RESULTS / "gene_tree_inventory.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(inventory[0]) if inventory else ["tree_index"], delimiter="\t")
        writer.writeheader()
        writer.writerows(inventory)
    (RESULTS / "tree_taxon_sets.tsv").write_text("taxon_set_id\tn_taxa\tn_trees\ttaxa\n" + "".join(
        f"{taxon_ids[taxa]}\t{len(taxa)}\t{n}\t{','.join(taxa)}\n" for taxa, n in sorted(counts.items())))
    print(f"records={len(records)} usable={len(records)-len(failures)} malformed={len(failures)}")
    print(f"taxon_sets={len(counts)} inventory={out}")
    if failures:
        print(f"malformed indices: {[i for i, _ in failures]}")
    print("Coordinate mapping is not inferred from tree order; validate explicit mapping separately.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
