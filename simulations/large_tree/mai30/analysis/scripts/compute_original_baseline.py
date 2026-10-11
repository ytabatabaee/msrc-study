#!/usr/bin/env python3
"""Validate and summarize selected original Mai D1 true gene trees.

The topology calculations are implemented directly on Newick clades so that
the normalized RF convention is explicit: unrooted RF divided by 2(n-3) for
fully bifurcating n-taxon trees.  No original gene-tree result is treated as
an MSRC result.
"""

from __future__ import annotations

import csv
import math
import random
import re
import sys
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "datasets" / "raw"
PROCESSED = ROOT / "datasets" / "processed"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from population_mapping import map_simphy_haploid_to_msrc  # noqa: E402


@dataclass
class ParsedTree:
    taxa: frozenset[str]
    splits: frozenset[frozenset[str]]
    root_height: float
    n_internal: int


def canonical_split(side: frozenset[str], taxa: frozenset[str]) -> frozenset[str] | None:
    other = taxa - side
    if len(side) < 2 or len(other) < 2:
        return None
    a = tuple(sorted(side))
    b = tuple(sorted(other))
    return frozenset(side) if a <= b else frozenset(other)


def parse_newick(text: str) -> ParsedTree:
    tokens = re.findall(r"[(),:;]|[^(),:;\s]+", text.strip())
    pos = 0

    def subtree() -> tuple[frozenset[str], list[frozenset[str]], float, int]:
        nonlocal pos
        child_data = []
        if tokens[pos] == "(":
            pos += 1
            child_data.append(subtree())
            while tokens[pos] == ",":
                pos += 1
                child_data.append(subtree())
            if tokens[pos] != ")":
                raise ValueError("Malformed Newick: expected ')'" )
            pos += 1
            if pos < len(tokens) and tokens[pos] not in ",);:":
                pos += 1  # internal label
            is_internal = True
        else:
            label = tokens[pos]
            if label in ",);:":
                raise ValueError("Malformed Newick: missing taxon label")
            pos += 1
            child_data.append((frozenset([label]), [], 0.0, 0))
            is_internal = False
        length = 0.0
        if pos < len(tokens) and tokens[pos] == ":":
            pos += 1
            length = float(tokens[pos])
            if not math.isfinite(length) or length < 0:
                raise ValueError("Negative or non-finite branch length")
            pos += 1
        taxa = frozenset().union(*(x[0] for x in child_data))
        splits = [split for x in child_data for split in x[1]]
        if is_internal:
            splits.append(taxa)
        heights = [x[2] for x in child_data]
        height = (max(heights) if heights else 0.0) + length
        n_internal = sum(x[3] for x in child_data) + int(is_internal)
        return taxa, splits, height, n_internal

    taxa, raw_splits, root_height, n_internal = subtree()
    if pos < len(tokens) and tokens[pos] == ";":
        pos += 1
    if pos != len(tokens):
        raise ValueError(f"Unparsed Newick tokens: {tokens[pos:]}")
    splits = frozenset(
        split for raw in raw_splits
        if (split := canonical_split(raw, taxa)) is not None
    )
    return ParsedTree(taxa, splits, root_height, n_internal)


def rf(tree: ParsedTree, species: ParsedTree) -> tuple[int, float]:
    distance = len(tree.splits ^ species.splits)
    maximum = 2 * (len(species.taxa) - 3)
    return distance, distance / maximum


def quartet_key(tree: ParsedTree, quartet: tuple[str, str, str, str]) -> frozenset[str]:
    q = frozenset(quartet)
    for split in tree.splits:
        a = split & q
        b = q - a
        if len(a) == 2 and len(b) == 2:
            left = tuple(sorted(a))
            right = tuple(sorted(b))
            return frozenset(a) if left <= right else frozenset(b)
    raise ValueError("Tree does not resolve the sampled quartet")


def condition_replicates() -> list[tuple[str, str]]:
    out = []
    for condition in sorted(p.name for p in RAW.glob("outgroup.0.species.*") if p.is_dir()):
        for rep in sorted(p.name for p in (RAW / condition).glob("[0-9][0-9][0-9]") if (p / "truegenetrees").exists()):
            out.append((condition, rep))
    return out


def summarize_replicate(condition: str, rep: str, quartet_count: int = 0) -> dict[str, str]:
    species_path = RAW / condition / rep / "s_tree.trees"
    genes_path = RAW / condition / rep / "truegenetrees"
    species = parse_newick(species_path.read_text())
    gene_rows = []
    rf_values = []
    concordance_values = []
    gene_splits: list[ParsedTree] = []
    with genes_path.open() as handle:
        for index, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            gene = parse_newick(line)
            if gene.taxa != species.taxa:
                raise ValueError(f"{genes_path}:{index}: taxon set differs from species tree")
            if gene.n_internal != len(species.taxa) - 1:
                raise ValueError(f"{genes_path}:{index}: gene tree is not fully bifurcating")
            distance, normalized = rf(gene, species)
            rf_values.append(normalized)
            concordance_values.append(sum(split in gene.splits for split in species.splits) / len(species.splits))
            gene_splits.append(gene)
    if len(gene_splits) != 500:
        raise ValueError(f"{genes_path}: expected 500 trees, found {len(gene_splits)}")
    quartets = []
    if quartet_count:
        rng = random.Random(20261011)
        taxa = sorted(species.taxa)
        for _ in range(quartet_count):
            quartet = tuple(rng.sample(taxa, 4))
            truth = quartet_key(species, quartet)
            concordant = sum(quartet_key(gene, quartet) == truth for gene in gene_splits)
            quartets.append(concordant / len(gene_splits))
    return {
        "condition": condition,
        "replicate": rep,
        "n_taxa": str(len(species.taxa)),
        "n_genes": str(len(gene_splits)),
        "mean_normalized_rf": f"{sum(rf_values) / len(rf_values):.10g}",
        "exact_species_topology_fraction": f"{sum(r == 0 for r in rf_values) / len(rf_values):.10g}",
        "mean_species_split_concordance": f"{sum(concordance_values) / len(concordance_values):.10g}",
        "mean_sampled_quartet_concordance": f"{sum(quartets) / len(quartets):.10g}" if quartets else "",
        "gene_tree_rf_min": f"{min(rf_values):.10g}",
        "gene_tree_rf_median": f"{sorted(rf_values)[len(rf_values)//2]:.10g}",
        "gene_tree_rf_max": f"{max(rf_values):.10g}",
        "gene_tree_branch_lengths": "present; units not used for coalescent comparison until SimPhy gene-tree units are independently confirmed",
    }


def main() -> None:
    all_rows = []
    for condition, rep in condition_replicates():
        all_rows.append(summarize_replicate(condition, rep))
    all_rows.sort(key=lambda row: float(row["mean_normalized_rf"]))
    # Quantiles of measured discordance, restricted to the 30-leaf no-outgroup data.
    indices = [0, len(all_rows) // 3, (2 * len(all_rows)) // 3, len(all_rows) - 1]
    selected = [all_rows[i] for i in indices]
    selected_keys = {(row["condition"], row["replicate"]) for row in selected}
    parameter_by_key = {}
    with (PROCESSED / "mai30_replicate_parameters.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            parameter_by_key[(row["condition"], row["replicate"])] = row
    selected_rows = []
    for condition, rep in sorted(selected_keys):
        row = summarize_replicate(condition, rep, quartet_count=1000)
        parameter = parameter_by_key[(condition, rep)]
        mapping = map_simphy_haploid_to_msrc(int(float(parameter["haploid_effective_population_size"])))
        row.update({
            "species_tree_path": parameter["species_tree"],
            "original_gene_trees_path": str(RAW / condition / rep / "truegenetrees"),
            "species_tree_sha256": parameter["species_tree_sha256"],
            "simphy_command": parameter["command"],
            "generations": parameter["generations"],
            "haploid_effective_population_size": parameter["haploid_effective_population_size"],
            "msrc_diploid_ne_exact": f"{mapping.msrc_diploid_ne_exact:.10g}",
            "msrc_diploid_ne_integer": str(mapping.msrc_diploid_ne_integer),
            "integer_mapping_relative_rate_error": f"{mapping.relative_rate_error:.10g}",
            "species_tree_root_height_generations": parameter["root_height_generations"],
            "gene_tree_units_note": "original gene-tree branch lengths are retained; topology statistics only until units are confirmed",
        })
        selected_rows.append(row)
    fields = list(selected_rows[0])
    out = PROCESSED / "mai30_original_gene_tree_baseline.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(selected_rows)
    manifest = PROCESSED / "mai30_original_gene_tree_selection.tsv"
    with manifest.open("w", newline="") as handle:
        fields = ["selection_rank", "selection_rule", *selected_rows[0].keys()]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for rank, row in enumerate(selected, start=1):
            enriched = next(x for x in selected_rows if x["condition"] == row["condition"] and x["replicate"] == row["replicate"])
            writer.writerow({"selection_rank": rank, "selection_rule": "discordance quantile among 30-leaf no-outgroup replicates", **enriched})
    print(f"audited {len(all_rows)} replicates")
    print(f"selected {len(selected_rows)} representatives")
    for row in selected_rows:
        print(row["condition"], row["replicate"], row["mean_normalized_rf"], row["mean_sampled_quartet_concordance"])


if __name__ == "__main__":
    main()
