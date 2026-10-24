#!/usr/bin/env python3
"""Select a deterministic quartet panel without using gene-tree outcomes."""
from __future__ import annotations

import csv
import itertools
import sys
from pathlib import Path

import yaml

MSRC_ROOT = Path("/Users/ytabatabaee/Desktop/msrc-sim")
ROOT = Path(__file__).resolve().parents[5]
MAI = ROOT / "simulations/large_tree/mai30"
OUT = MAI / "analysis/origin_placement"
sys.path.insert(0, str(MSRC_ROOT))
from msrcsim import SpeciesTree  # noqa: E402

TREES = [
    ("outgroup.0.species.1.5.genes.1.5", "007"),
    ("outgroup.0.species.1.5.genes.1.5", "042"),
    ("outgroup.0.species.1.5.genes.1.5", "074"),
    ("outgroup.0.species.5.genes.5", "005"),
]
SCENARIOS = ("A_existing_internal", "B_deep_internal")


def config_for(condition: str, replicate: str) -> dict:
    path = MAI / "datasets/configs/mai30_msrc_pilot" / f"{condition}__replicate_{replicate}__msrc_no_suppression__freq_0.5.yaml"
    return yaml.safe_load(path.read_text())


def mrca(tree: SpeciesTree, names: set[str]):
    nodes = {node.name: node for node in tree.nodes()}
    current = [nodes[n] for n in names]
    ancestors = []
    x = current[0]
    while x is not None:
        ancestors.append(x)
        x = x.parent
    for candidate in ancestors:
        if all((lambda y: candidate in y)([node for node in _ancestors(n)]) for n in current):
            return candidate
    raise ValueError("no MRCA")


def _ancestors(node):
    out = []
    while node is not None:
        out.append(node)
        node = node.parent
    return out


def quartet_split(tree: SpeciesTree, quartet: tuple[str, ...]) -> tuple[frozenset[str], frozenset[str]]:
    nodes = {node.name: node for node in tree.nodes()}
    for left in itertools.combinations(quartet, 2):
        right = tuple(n for n in quartet if n not in left)
        if mrca(tree, set(left)) is not mrca(tree, set(right)):
            a, b = frozenset(left), frozenset(right)
            return (a, b) if tuple(sorted(a)) <= tuple(sorted(b)) else (b, a)
    raise ValueError(f"unresolved quartet {quartet}")


def panel_for(tree: SpeciesTree, scenario: dict, per_bin: int = 2) -> list[dict]:
    descendants = set(tree.descendants(scenario["origin_branch"])) if scenario["origin_branch"] != tree.root.name else set(tree.taxa)
    bins: dict[int, list[dict]] = {2: [], 3: [], 4: []}
    for quartet in itertools.combinations(tree.taxa, 4):
        exposed = len(descendants & set(quartet))
        if exposed not in bins:
            continue
        left, right = quartet_split(tree, quartet)
        all_mrca = mrca(tree, set(quartet))
        pair_mrca = max(mrca(tree, set(left)).age, mrca(tree, set(right)).age)
        internal_length = float(all_mrca.age - pair_mrca)
        bins[exposed].append({"quartet": ",".join(quartet), "exposure_count": exposed,
                              "induced_internal_length": internal_length})
    rows = []
    for exposed in (2, 3, 4):
        candidates = sorted(bins[exposed], key=lambda r: (r["induced_internal_length"], r["quartet"]))
        chosen = candidates[:per_bin] + candidates[-per_bin:]
        seen = set()
        for rank, row in enumerate(chosen, 1):
            if row["quartet"] in seen:
                continue
            seen.add(row["quartet"])
            row = dict(row)
            row["selection_rule"] = f"geometry-only exposure={exposed}; {'short' if rank <= per_bin else 'long'} induced internal branch"
            rows.append(row)
    return rows


def main() -> None:
    rows = []
    scenario_file = OUT / "origin_scenarios.tsv"
    with scenario_file.open() as handle:
        scenarios = list(csv.DictReader(handle, delimiter="\t"))
    for condition, replicate in TREES:
        config = config_for(condition, replicate)
        for scenario_name in SCENARIOS:
            scenario = next(r for r in scenarios if r["condition"] == condition and r["replicate"] == replicate and r["scenario"] == scenario_name)
            tree = SpeciesTree(config["species_tree"]["newick"], config["species_tree"]["default_effective_population_size"], float(scenario["root_extension"]), ultrametric_tolerance=1e-6)
            for row in panel_for(tree, scenario):
                rows.append({"condition": condition, "replicate": replicate, "scenario": scenario_name, **row})
    out = OUT / "selected_quartets.tsv"
    with out.open("w", newline="") as handle:
        fields = list(rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    print(f"wrote {len(rows)} selected quartet rows to {out}")


if __name__ == "__main__":
    main()
