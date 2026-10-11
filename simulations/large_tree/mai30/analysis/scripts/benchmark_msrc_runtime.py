#!/usr/bin/env python3
"""Benchmark scaled forward histories and full-tree generation.

The benchmark never relabels a scaled run as an exact Mai simulation.  It
records the scale factor and the coalescent-unit preservation explicitly.
"""

from __future__ import annotations

import csv
import resource
import sys
from pathlib import Path
from time import perf_counter

import numpy as np

MSRC_ROOT = Path("/Users/ytabatabaee/Desktop/msrc-sim")
STUDY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(MSRC_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from msrcsim import SpeciesTree  # noqa: E402
from msrcsim.arbitrary import rearrangement_from_config  # noqa: E402
from msrcsim.population_process import simulate_population_history  # noqa: E402
from msrcsim.structured_coalescent import simulate_genealogy  # noqa: E402
from population_mapping import map_simphy_haploid_to_msrc  # noqa: E402


CONDITION = "outgroup.0.species.1.5.genes.1.5"
REPLICATE = "042"
SCALES = (0.0001, 0.0005, 0.001, 0.002)
GENE_TREES_FOR_BENCHMARK = 20


def scale_newick(newick: str, factor: float) -> str:
    import re
    return re.sub(r":([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)", lambda m: f":{float(m.group(1)) * factor:.12g}", newick)


def main() -> None:
    params = {}
    with (STUDY_ROOT / "datasets/processed/mai30_replicate_parameters.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            params[(row["condition"], row["replicate"])] = row
    row = params[(CONDITION, REPLICATE)]
    source_newick = (STUDY_ROOT / row["species_tree"]).read_text().strip()
    mapping = map_simphy_haploid_to_msrc(int(float(row["haploid_effective_population_size"])))
    records = []
    for index, scale in enumerate(SCALES):
        n_scaled = max(1, int(round(mapping.msrc_diploid_ne_integer * scale)))
        scaled_newick = scale_newick(source_newick, scale)
        root_extension_cu = 5.0
        root_extension_generations = root_extension_cu * 2.0 * n_scaled
        tree = SpeciesTree(scaled_newick, n_scaled, root_extension_generations, ultrametric_tolerance=1e-6)
        config = {"population_process": {"model": "wright_fisher"}}
        rearrangement = rearrangement_from_config({
            "rearrangement": {
                "id": "benchmark_inv",
                "type": "inversion",
                "origin_branch": tree.root.name,
                "origin_time_from_branch_start": max(1, int(root_extension_generations // 2)),
                "initial_copy_count": 1,
                "selection": {"coefficient": 0.0},
            }
        })
        rng = np.random.default_rng(20261011 + index)
        history_start = perf_counter()
        history = simulate_population_history(tree, rearrangement, rng, config)
        history_seconds = perf_counter() - history_start
        sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
        gene_start = perf_counter()
        results = [simulate_genealogy(i, tree, history, sampled, 0.0, 1.0, rng, record_events=False) for i in range(GENE_TREES_FOR_BENCHMARK)]
        gene_seconds = perf_counter() - gene_start
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        rss_mb = rss / (1024 * 1024) if sys.platform != "darwin" else rss / (1024 * 1024)
        records.append({
            "condition": CONDITION,
            "replicate": REPLICATE,
            "scale_factor_c": scale,
            "coalescent_units_preserved": "yes; branch generations and N scaled together",
            "source_haploid_ne": mapping.simphy_haploid_ne,
            "source_species_root_generations": row["root_height_generations"],
            "scaled_msrc_diploid_ne": n_scaled,
            "scaled_species_root_generations": tree.root.age,
            "root_extension_coalescent_units": root_extension_cu,
            "root_extension_generations": root_extension_generations,
            "n_population_records": len(history.records),
            "n_gene_trees": len(results),
            "history_runtime_seconds": history_seconds,
            "genealogy_runtime_seconds": gene_seconds,
            "seconds_per_gene_tree": gene_seconds / len(results),
            "estimated_500_gene_tree_seconds": gene_seconds / len(results) * 500.0,
            "peak_rss_mb_process": rss_mb,
            "population_process": "wright_fisher",
            "complete_genealogies": all(len(result.coalescences) == len(tree.taxa) - 1 for result in results),
        })
    out = STUDY_ROOT / "analysis/results/mai30_runtime_benchmark.tsv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    print(out)
    for record in records:
        print(record["scale_factor_c"], record["history_runtime_seconds"], record["estimated_500_gene_tree_seconds"])


if __name__ == "__main__":
    main()
