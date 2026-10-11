#!/usr/bin/env python3
"""Generate the coalescent-timescale-scaled MSC controls for the four Mai reps."""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path
from time import perf_counter

import numpy as np

MSRC_ROOT = Path("/Users/ytabatabaee/Desktop/msrc-sim")
STUDY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(MSRC_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from compute_original_baseline import parse_newick, quartet_key, rf  # noqa: E402
from msrcsim import SpeciesTree, simulate_msc_genealogy  # noqa: E402
from population_mapping import map_simphy_haploid_to_msrc  # noqa: E402


SELECTED = [
    ("outgroup.0.species.1.5.genes.1.5", "007"),
    ("outgroup.0.species.1.5.genes.1.5", "042"),
    ("outgroup.0.species.1.5.genes.1.5", "074"),
    ("outgroup.0.species.5.genes.5", "005"),
]
SCALE = 0.001
ROOT_CU = 20.0
N_TREES = 500


def scale_newick(newick: str) -> str:
    return re.sub(
        r":([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)",
        lambda m: f":{float(m.group(1)) * SCALE:.12g}",
        newick,
    )


def parameter_map() -> dict[tuple[str, str], dict[str, str]]:
    with (STUDY_ROOT / "datasets/processed/mai30_replicate_parameters.tsv").open() as handle:
        return {(r["condition"], r["replicate"]): r for r in csv.DictReader(handle, delimiter="\t")}


def summarize(generated, species, original):
    generated_rf = [rf(g, species)[1] for g in generated]
    original_rf = [rf(g, species)[1] for g in original]
    generated_split = [sum(s in g.splits for s in species.splits) / len(species.splits) for g in generated]
    original_split = [sum(s in g.splits for s in species.splits) / len(species.splits) for g in original]
    return {
        "scaled_msc_mean_normalized_rf": float(np.mean(generated_rf)),
        "original_mai_mean_normalized_rf": float(np.mean(original_rf)),
        "scaled_minus_original_mean_rf": float(np.mean(generated_rf) - np.mean(original_rf)),
        "scaled_msc_exact_topology_fraction": float(np.mean(np.asarray(generated_rf) == 0)),
        "original_mai_exact_topology_fraction": float(np.mean(np.asarray(original_rf) == 0)),
        "scaled_msc_mean_species_split_concordance": float(np.mean(generated_split)),
        "original_mai_mean_species_split_concordance": float(np.mean(original_split)),
        "scaled_msc_rf_q05": float(np.quantile(generated_rf, 0.05)),
        "scaled_msc_rf_q50": float(np.quantile(generated_rf, 0.50)),
        "scaled_msc_rf_q95": float(np.quantile(generated_rf, 0.95)),
        "original_mai_rf_q05": float(np.quantile(original_rf, 0.05)),
        "original_mai_rf_q50": float(np.quantile(original_rf, 0.50)),
        "original_mai_rf_q95": float(np.quantile(original_rf, 0.95)),
    }


def main() -> None:
    params = parameter_map()
    rows = []
    out_root = STUDY_ROOT / "datasets/processed/mai30_scaled_msc_controls"
    for selected_index, (condition, rep) in enumerate(SELECTED):
        parameter = params[(condition, rep)]
        source_species = parse_newick((STUDY_ROOT / parameter["species_tree"]).read_text())
        original = [parse_newick(line) for line in (
            STUDY_ROOT / "datasets/raw" / condition / rep / "truegenetrees"
        ).read_text().splitlines() if line.strip()]
        mapping = map_simphy_haploid_to_msrc(int(float(parameter["haploid_effective_population_size"])))
        scaled_exact = mapping.msrc_diploid_ne_exact * SCALE
        n_scaled = max(1, int(round(scaled_exact)))
        species_newick = scale_newick((STUDY_ROOT / parameter["species_tree"]).read_text().strip())
        tree = SpeciesTree(
            species_newick,
            n_scaled,
            ROOT_CU * 2.0 * n_scaled,
            ultrametric_tolerance=1e-6,
        )
        seed = 20300000 + selected_index
        rng = np.random.default_rng(seed)
        generated = []
        start = perf_counter()
        out = out_root / condition / rep
        out.mkdir(parents=True, exist_ok=True)
        with (out / "msc_gene_trees.nwk").open("w") as handle:
            for locus_id in range(N_TREES):
                result = simulate_msc_genealogy(locus_id, tree, rng, record_events=True)
                parsed = parse_newick(result.newick)
                if parsed.taxa != source_species.taxa:
                    raise ValueError(f"taxon set mismatch: {condition}/{rep}")
                if parsed.n_internal != len(source_species.taxa) - 1:
                    raise ValueError(f"non-bifurcating tree: {condition}/{rep}")
                if len(result.coalescences) != len(source_species.taxa) - 1:
                    raise ValueError(f"incomplete genealogy: {condition}/{rep}")
                generated.append(parsed)
                handle.write(result.newick + "\n")
        runtime = perf_counter() - start
        summary = summarize(generated, tree_from_source(source_species), original)
        summary.update({
            "condition": condition,
            "replicate": rep,
            "scale_factor": SCALE,
            "root_extension_coalescent_units": ROOT_CU,
            "root_extension_generations": ROOT_CU * 2.0 * n_scaled,
            "simphy_haploid_ne": mapping.simphy_haploid_ne,
            "msrc_diploid_ne_exact": mapping.msrc_diploid_ne_exact,
            "scaled_ne_exact": scaled_exact,
            "scaled_ne_integer": n_scaled,
            "scaled_ne_relative_rounding_error": (n_scaled - scaled_exact) / scaled_exact,
            "seed": seed,
            "n_taxa": len(source_species.taxa),
            "n_gene_trees": N_TREES,
            "runtime_seconds": runtime,
            "seconds_per_gene_tree": runtime / N_TREES,
            "output_contract": "PASS",
            "output_path": str(out / "msc_gene_trees.nwk"),
        })
        with (out / "validation_summary.json").open("w") as handle:
            json.dump(summary, handle, indent=2)
        rows.append(summary)
    out_tsv = STUDY_ROOT / "analysis/results/mai30_scaled_msc_controls.tsv"
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out_tsv}")


def tree_from_source(source):
    """Return the same topology object used by the source-summary machinery."""
    return source


if __name__ == "__main__":
    main()
