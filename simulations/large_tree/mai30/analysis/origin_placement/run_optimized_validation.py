#!/usr/bin/env python3
"""Validated follow-up for origin placement using the optimized simulator.

Forward histories are shared across rho values.  Genealogy-level results are
only emitted when explicitly requested; structural exposure and persistence
are never interpreted as gene-tree topology results.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import yaml

MSRC_ROOT = Path("/Users/ytabatabaee/Desktop/msrc-sim")
ROOT = Path(__file__).resolve().parents[5]
MAI = ROOT / "simulations/large_tree/mai30"
OUT = MAI / "analysis/origin_placement"
sys.path.insert(0, str(MSRC_ROOT))
sys.path.insert(0, str(MAI / "analysis/scripts"))

from compute_original_baseline import parse_newick, quartet_key  # noqa: E402
from msrcsim import SpeciesTree  # noqa: E402
from msrcsim.population_process import simulate_population_history  # noqa: E402
from msrcsim.rearrangement import Rearrangement  # noqa: E402
from msrcsim.structured_coalescent import simulate_genealogy  # noqa: E402

TREES = [
    ("outgroup.0.species.1.5.genes.1.5", "007"),
    ("outgroup.0.species.1.5.genes.1.5", "042"),
    ("outgroup.0.species.1.5.genes.1.5", "074"),
    ("outgroup.0.species.5.genes.5", "005"),
]
RHOS = [0.0, 1e-4, 1e-3, 0.002, 0.01, 0.02, 0.1, 1.0]
SCENARIOS = ["A_existing_internal", "B_deep_internal", "C_root_adjacent", "D_older_root_extended"]
SEED = 20261022


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def config_for(condition: str, replicate: str) -> dict:
    path = MAI / "datasets/configs/mai30_msrc_pilot" / f"{condition}__replicate_{replicate}__msrc_no_suppression__freq_0.5.yaml"
    return yaml.safe_load(path.read_text())


def sample_quartets(tree: SpeciesTree, scenario: dict, per_bin: int = 2) -> list[tuple[str, str, str, str]]:
    descendants = set(tree.descendants(scenario["origin_branch"]))
    bins: dict[int, list[tuple[str, str, str, str]]] = {2: [], 3: [], 4: []}
    for quartet in itertools.combinations(tree.taxa, 4):
        count = len(descendants & set(quartet))
        if count in bins and len(bins[count]) < per_bin:
            bins[count].append(quartet)
    return [quartet for count in (2, 3, 4) for quartet in bins[count]]


def scenario_rows(condition: str, replicate: str, config: dict, tree: SpeciesTree) -> list[dict]:
    source = config["rearrangement"]
    candidates = []
    for branch in tree.branches.values():
        descendants = tree.descendants(branch.child_node)
        if branch.branch_id != tree.root.name and 4 <= len(descendants) <= len(tree.taxa) - 2:
            duration = branch.older_age - branch.younger_age
            if duration >= 2:
                candidates.append((branch.older_age, len(descendants), branch.branch_id, duration))
    _, _, deep_branch, deep_duration = max(candidates)
    root = tree.root.name
    rows = [
        ("A_existing_internal", source["origin_branch"], int(source["origin_time_from_branch_start"]), tree.root_extension),
        ("B_deep_internal", deep_branch, max(1, int(deep_duration - min(50, max(1, deep_duration - 1)))), tree.root_extension),
        ("C_root_adjacent", root, max(1, int(tree.root_extension - 50)), tree.root_extension),
        ("D_older_root_extended", root, 2, 2 * tree.root_extension),
    ]
    result = []
    for scenario, branch, origin_time, root_extension in rows:
        active_tree = SpeciesTree(config["species_tree"]["newick"], config["species_tree"]["default_effective_population_size"], root_extension, ultrametric_tolerance=1e-6)
        result.append({"condition": condition, "replicate": replicate, "scenario": scenario,
                       "origin_branch": branch, "origin_time_from_branch_start": origin_time,
                       "root_extension": root_extension, "descendant_taxa": len(active_tree.descendants(branch)) if branch != root else len(active_tree.taxa)})
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--histories", type=int, default=1000)
    parser.add_argument("--genealogy-histories", type=int, default=0,
                        help="0 disables backward simulation; use a prespecified count for quartet validation")
    parser.add_argument("--scenarios", nargs="+", choices=SCENARIOS,
                        default=["A_existing_internal", "B_deep_internal", "C_root_adjacent"])
    parser.add_argument("--rhos", nargs="+", type=float, default=RHOS)
    parser.add_argument("--tree-indices", nargs="+", type=int, default=[0, 1, 2, 3])
    parser.add_argument("--output-prefix", type=Path, default=OUT / "optimized_validation")
    args = parser.parse_args()
    forward_rows, pattern_rows, quartet_rows = [], [], []
    scenario_rows_all = []
    benchmark_rows = []
    for tree_index, (condition, replicate) in enumerate(TREES):
        if tree_index not in args.tree_indices:
            continue
        config = config_for(condition, replicate)
        base_tree = SpeciesTree(config["species_tree"]["newick"], config["species_tree"]["default_effective_population_size"], config["species_tree"]["root_extension"], ultrametric_tolerance=1e-6)
        for scenario_index, scenario in enumerate(scenario_rows(condition, replicate, config, base_tree)):
            if scenario["scenario"] not in args.scenarios:
                continue
            scenario_rows_all.append(scenario)
            tree = SpeciesTree(config["species_tree"]["newick"], config["species_tree"]["default_effective_population_size"], scenario["root_extension"], ultrametric_tolerance=1e-6)
            rr = config["rearrangement"]
            rearrangement = Rearrangement("origin_validation", rr["type"], scenario["origin_branch"], scenario["origin_time_from_branch_start"], int(rr["initial_copy_count"]), 0.0)
            selected = sample_quartets(tree, scenario)
            counts = Counter()
            started = time.perf_counter()
            for history_index in range(args.histories):
                history = simulate_population_history(tree, rearrangement, np.random.default_rng(SEED + tree_index * 10_000_000 + scenario_index * 100_000 + history_index), config)
                origin = history.by_branch[scenario["origin_branch"]]
                polymorphic = [0 < record.copy_count_A1 < record.population_chromosomes for record in origin]
                terminal = {taxon: int(np.random.default_rng(SEED + 777 + history_index).random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
                forward_rows.append({"condition": condition, "replicate": replicate, "scenario": scenario["scenario"], "history": history_index,
                                     "origin_polymorphic_fraction": float(np.mean(polymorphic)),
                                     "origin_end_polymorphic": int(polymorphic[-1]),
                                     "ever_fixed": int(any(record.status == "fixed" for record in origin)),
                                     "ever_lost": int(any(record.status == "lost" for record in origin)),
                                     "terminal_polymorphic": int(any(0 < history.terminal_frequency(t) < 1 for t in tree.taxa))})
                for quartet in selected:
                    pattern = "".join(str(terminal[taxon]) for taxon in quartet)
                    pattern_rows.append({"condition": condition, "replicate": replicate, "scenario": scenario["scenario"], "history": history_index, "quartet": ",".join(quartet), "pattern": pattern})
                if history_index >= args.genealogy_histories:
                    continue
                for rho in args.rhos:
                    sampled = terminal
                    result = simulate_genealogy(history_index, tree, history, sampled, rho, 1.0,
                                                np.random.default_rng(SEED + 1_000_000_000 + tree_index * 10_000_000 + scenario_index * 100_000 + history_index * 100 + int(rho * 1000)), False)
                    gene = parse_newick(result.newick)
                    for quartet in selected:
                        actual = quartet_key(gene, quartet)
                        species = quartet_key(parse_newick(config["species_tree"]["newick"]), quartet)
                        label = "species_tree" if actual == species else "discordant"
                        counts[(rho, ",".join(quartet), "".join(str(terminal[taxon]) for taxon in quartet), label)] += 1
            benchmark_rows.append({"condition": condition, "replicate": replicate, "scenario": scenario["scenario"], "histories": args.histories, "genealogy_histories": min(args.genealogy_histories, args.histories), "seconds": time.perf_counter() - started})
            for (rho, quartet, pattern, label), count in sorted(counts.items()):
                quartet_rows.append({"condition": condition, "replicate": replicate, "scenario": scenario["scenario"], "rho": rho, "quartet": quartet, "pattern": pattern, "topology_class": label, "count": count})
    prefix = args.output_prefix
    write_tsv(prefix.with_name(prefix.name + "_forward.tsv"), forward_rows)
    write_tsv(prefix.with_name(prefix.name + "_patterns.tsv"), pattern_rows)
    write_tsv(prefix.with_name(prefix.name + "_quartets.tsv"), quartet_rows)
    write_tsv(prefix.with_name(prefix.name + "_benchmark.tsv"), benchmark_rows)
    write_tsv(prefix.with_name(prefix.name + "_scenarios.tsv"), scenario_rows_all)
    (prefix.with_name(prefix.name + "_manifest.json")).write_text(json.dumps({
        "status": "executed", "histories": args.histories, "genealogy_histories": args.genealogy_histories,
        "rho_grid": args.rhos, "seed": SEED, "simulator": str(MSRC_ROOT),
        "forward_histories_reused_across_rho": True,
        "quartets_are_induced_from_full_30_taxon_genealogies": True,
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
