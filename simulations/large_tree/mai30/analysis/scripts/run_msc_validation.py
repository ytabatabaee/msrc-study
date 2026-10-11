#!/usr/bin/env python3
"""Generate matched ordinary MSC controls for selected Mai replicates."""

from __future__ import annotations

import csv
import json
import random
import sys
from pathlib import Path
from time import perf_counter

import numpy as np

MSRC_ROOT = Path("/Users/ytabatabaee/Desktop/msrc-sim")
STUDY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(MSRC_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from msrcsim import SpeciesTree, simulate_msc_genealogy  # noqa: E402
from compute_original_baseline import parse_newick, quartet_key, rf  # noqa: E402
from population_mapping import map_simphy_haploid_to_msrc  # noqa: E402


SELECTED = [
    ("outgroup.0.species.1.5.genes.1.5", "007"),
    ("outgroup.0.species.1.5.genes.1.5", "042"),
    ("outgroup.0.species.1.5.genes.1.5", "074"),
    ("outgroup.0.species.5.genes.5", "005"),
]
ROOT_CU = (5, 10, 20, 40)


def parameter_map() -> dict[tuple[str, str], dict[str, str]]:
    with (STUDY_ROOT / "datasets/processed/mai30_replicate_parameters.tsv").open() as handle:
        return {(r["condition"], r["replicate"]): r for r in csv.DictReader(handle, delimiter="\t")}


def original_data(condition: str, rep: str):
    raw = STUDY_ROOT / "datasets/raw" / condition / rep
    species = parse_newick((raw / "s_tree.trees").read_text())
    genes = [parse_newick(line) for line in (raw / "truegenetrees").read_text().splitlines() if line.strip()]
    return species, genes


def quartet_sample(taxa: list[str], count: int = 1000) -> list[tuple[str, str, str, str]]:
    rng = random.Random(20261011)
    return [tuple(rng.sample(taxa, 4)) for _ in range(count)]


def summarize(generated, species, original, quartets):
    generated_rf = [rf(g, species)[1] for g in generated]
    original_rf = [rf(g, species)[1] for g in original]
    gen_split = [sum(s in g.splits for s in species.splits) / len(species.splits) for g in generated]
    orig_split = [sum(s in g.splits for s in species.splits) / len(species.splits) for g in original]
    gen_quartet = []
    orig_quartet = []
    for quartet in quartets:
        truth = quartet_key(species, quartet)
        gen_quartet.append(sum(quartet_key(g, quartet) == truth for g in generated) / len(generated))
        orig_quartet.append(sum(quartet_key(g, quartet) == truth for g in original) / len(original))
    gen_sorted = sorted(generated_rf)
    orig_sorted = sorted(original_rf)
    ecdf_distance = sum(abs(a - b) for a, b in zip(gen_sorted, orig_sorted)) / len(gen_sorted)
    return {
        "generated_mean_normalized_rf": float(np.mean(generated_rf)),
        "original_mean_normalized_rf": float(np.mean(original_rf)),
        "mean_normalized_rf_difference": float(np.mean(generated_rf) - np.mean(original_rf)),
        "generated_exact_topology_fraction": float(np.mean(np.array(generated_rf) == 0)),
        "original_exact_topology_fraction": float(np.mean(np.array(original_rf) == 0)),
        "generated_mean_species_split_concordance": float(np.mean(gen_split)),
        "original_mean_species_split_concordance": float(np.mean(orig_split)),
        "generated_mean_sampled_quartet_concordance": float(np.mean(gen_quartet)),
        "original_mean_sampled_quartet_concordance": float(np.mean(orig_quartet)),
        "sorted_rf_ecdf_l1_distance": float(ecdf_distance),
        "generated_rf_q05": float(np.quantile(generated_rf, 0.05)),
        "generated_rf_q50": float(np.quantile(generated_rf, 0.50)),
        "generated_rf_q95": float(np.quantile(generated_rf, 0.95)),
        "original_rf_q05": float(np.quantile(original_rf, 0.05)),
        "original_rf_q50": float(np.quantile(original_rf, 0.50)),
        "original_rf_q95": float(np.quantile(original_rf, 0.95)),
    }


def main() -> None:
    params = parameter_map()
    rows = []
    out_root = STUDY_ROOT / "datasets/processed/msc_validation"
    for selected_index, (condition, rep) in enumerate(SELECTED):
        parameter = params[(condition, rep)]
        source_species, original = original_data(condition, rep)
        mapping = map_simphy_haploid_to_msrc(int(float(parameter["haploid_effective_population_size"])))
        species_newick = (STUDY_ROOT / parameter["species_tree"]).read_text().strip()
        quartets = quartet_sample(sorted(source_species.taxa))
        for root_cu in ROOT_CU:
            root_extension = root_cu * 2.0 * mapping.msrc_diploid_ne_integer
            # SimPhy writes branch lengths to 8 decimal places.  The Mai trees
            # are ultrametric within that serialization precision; retain their
            # lengths and relax only the parser's equality tolerance.
            tree = SpeciesTree(species_newick, mapping.msrc_diploid_ne_integer, root_extension, ultrametric_tolerance=1e-6)
            seed = 20261011 + selected_index * 100000 + root_cu
            rng = np.random.default_rng(seed)
            generated = []
            start = perf_counter()
            for locus_id in range(500):
                result = simulate_msc_genealogy(locus_id, tree, rng, record_events=True)
                parsed = parse_newick(result.newick)
                if parsed.taxa != source_species.taxa:
                    raise ValueError(f"Generated taxon set mismatch: {condition}/{rep}/R{root_cu}")
                if parsed.n_internal != len(source_species.taxa) - 1:
                    raise ValueError(f"Generated tree is not bifurcating: {condition}/{rep}/R{root_cu}")
                if len(result.coalescences) != len(source_species.taxa) - 1:
                    raise ValueError(f"Incomplete genealogy: {condition}/{rep}/R{root_cu}")
                generated.append(parsed)
            runtime = perf_counter() - start
            out = out_root / condition / rep / f"root_{root_cu}cu"
            out.mkdir(parents=True, exist_ok=True)
            # Recreate the validated stream once for export, avoiding event objects in memory.
            export_rng = np.random.default_rng(seed)
            with (out / "msc_gene_trees.nwk").open("w") as handle:
                for locus_id in range(500):
                    handle.write(simulate_msc_genealogy(locus_id, tree, export_rng, record_events=False).newick + "\n")
            summary = summarize(generated, source_species, original, quartets)
            summary.update({
                "condition": condition,
                "replicate": rep,
                "root_extension_coalescent_units": root_cu,
                "root_extension_generations": root_extension,
                "seed": seed,
                "n_taxa": len(source_species.taxa),
                "n_gene_trees": len(generated),
                "runtime_seconds": runtime,
                "seconds_per_gene_tree": runtime / len(generated),
                "msrc_diploid_ne_exact": mapping.msrc_diploid_ne_exact,
                "msrc_diploid_ne_integer": mapping.msrc_diploid_ne_integer,
                "integer_mapping_relative_rate_error": mapping.relative_rate_error,
                "output_contract": "PASS",
                "output_path": str(out / "msc_gene_trees.nwk"),
            })
            with (out / "validation_summary.json").open("w") as handle:
                json.dump(summary, handle, indent=2)
            rows.append(summary)
    out_tsv = STUDY_ROOT / "analysis/results/mai30_msc_validation.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out_tsv}")


if __name__ == "__main__":
    main()
