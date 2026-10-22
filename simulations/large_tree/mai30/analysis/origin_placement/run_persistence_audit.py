#!/usr/bin/env python3
"""Forward-history audit for all origin scenarios.

This deliberately excludes backward genealogies.  It measures unconditional
Wright--Fisher persistence and terminal arrangement sorting with the same
authoritative simulator process used by the genealogy runner.
"""

from __future__ import annotations

import csv
import itertools
from pathlib import Path
import sys

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[5]
MAI = ROOT / "simulations/large_tree/mai30"
MSRC_ROOT = Path("/Users/ytabatabaee/Desktop/msrc-sim")
sys.path.insert(0, str(MSRC_ROOT))
sys.path.insert(0, str(MAI / "analysis/scripts"))
from msrcsim import SpeciesTree  # noqa: E402
from msrcsim.population_process import simulate_population_history  # noqa: E402
from msrcsim.rearrangement import Rearrangement  # noqa: E402
from run_origin_discovery import scenario_rows, load_base, TREES, SEED, N_HISTORIES  # noqa: E402
from compute_original_baseline import parse_newick, quartet_key  # noqa: E402


def write(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def main():
    out = HERE
    persistence = []
    patterns = {}
    for tree_index, (condition, replicate) in enumerate(TREES):
        base = load_base(condition, replicate)
        base_n = int(base["species_tree"]["default_effective_population_size"])
        base_tree = SpeciesTree(base["species_tree"]["newick"], base_n, base["species_tree"]["root_extension"], ultrametric_tolerance=1e-6)
        reference_tree = parse_newick(base["species_tree"]["newick"])
        scenarios = scenario_rows(condition, replicate, base, base_tree)
        for scenario_index, scenario in enumerate(scenarios):
            tree = SpeciesTree(base["species_tree"]["newick"], base_n, scenario["root_extension"], ultrametric_tolerance=1e-6)
            rr = Rearrangement("origin_audit", "inversion", scenario["origin_branch"], int(scenario["origin_time_from_branch_start"]), int(base["rearrangement"]["initial_copy_count"]), 0.0)
            for history_index in range(N_HISTORIES):
                seed = SEED + 10_000_000 + tree_index * 1_000_000 + scenario_index * 10_000 + history_index
                history = simulate_population_history(tree, rr, np.random.default_rng(seed), base)
                records = history.by_branch[scenario["origin_branch"]]
                persistence.append({
                    "condition": condition, "replicate": replicate, "scenario": scenario["scenario"], "history": history_index,
                    "origin_branch": scenario["origin_branch"], "origin_age": scenario["origin_age"],
                    "origin_branch_polymorphic_fraction": sum(0 < r.copy_count_A1 < r.population_chromosomes for r in records) / len(records),
                    "origin_branch_end_polymorphic": int(0 < records[-1].copy_count_A1 < records[-1].population_chromosomes),
                    "origin_branch_ever_fixed": int(any(r.status == "fixed" for r in records)),
                    "origin_branch_ever_lost": int(any(r.status == "lost" for r in records)),
                    "root_terminal_polymorphic": int(0 < history.frequency_at(tree.root.name, tree.root.age + tree.root_extension) < 1),
                    "status": "executed_forward_history",
                })
                sampled = {taxon: int(np.random.default_rng(seed + 77).random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
                # Use every quartet for terminal-pattern exposure; this is
                # cheap because it does not require gene-tree simulation.
                for quartet in itertools.combinations(tree.taxa, 4):
                    pattern = "".join(str(sampled[taxon]) for taxon in quartet)
                    n = pattern.count("1")
                    if n in (0, 4):
                        category = "uniform_4_0"
                    elif n in (1, 3):
                        category = "3_1"
                    else:
                        ones = {taxon for taxon, bit in zip(quartet, pattern) if bit == "1"}
                        side = set(quartet_key(reference_tree, quartet))
                        category = "concordant_2_2" if ones in (side, set(quartet) - side) else "discordant_2_2"
                    key = (condition, replicate, scenario["scenario"], pattern, category)
                    patterns[key] = patterns.get(key, 0) + 1
    write(out / "persistence_results.tsv", persistence)
    pattern_rows = []
    for key, count in sorted(patterns.items()):
        c, r, s, p, category = key
        total = sum(v for (cc, rr, ss, _p, _cat), v in patterns.items() if (cc, rr, ss) == (c, r, s))
        pattern_rows.append({"condition": c, "replicate": r, "scenario": s, "pattern": p, "category": category, "count": count, "probability": count / total, "support": "all_27405_quartets_per_history"})
    write(out / "arrangement_patterns.tsv", pattern_rows)
    print(f"wrote {len(persistence)} forward-history rows and {len(pattern_rows)} pattern rows")


if __name__ == "__main__":
    main()
