#!/usr/bin/env python3
"""Prepare, but do not execute, the computationally scaled MSRC pilot."""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import yaml

MSRC_ROOT = Path("/Users/ytabatabaee/Desktop/msrc-sim")
STUDY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(MSRC_ROOT))
from msrcsim import SpeciesTree  # noqa: E402
from population_mapping import map_simphy_haploid_to_msrc  # noqa: E402


SELECTED = [
    ("outgroup.0.species.1.5.genes.1.5", "007"),
    ("outgroup.0.species.1.5.genes.1.5", "042"),
    ("outgroup.0.species.1.5.genes.1.5", "074"),
    ("outgroup.0.species.5.genes.5", "005"),
]
TREATMENTS = {
    "msrc_no_suppression": 1.0,
    "msrc_moderate_suppression": 0.5,
    "msrc_strong_suppression": 0.1,
}
SCALE = 0.001
N_LOCI = 20
N_HISTORY_REPLICATES = 25
INITIAL_FREQUENCIES = {"freq_0.5": 0.5, "freq_0.2": 0.2}


def scale_newick(newick: str, factor: float) -> str:
    return re.sub(r":([0-9]+(?:\.[0-9]+)?(?:[eE][+-]?[0-9]+)?)", lambda m: f":{float(m.group(1)) * factor:.12g}", newick)


def choose_origin(tree: SpeciesTree) -> tuple[str, int]:
    candidates = []
    for branch in tree.branches.values():
        node = next(n for n in tree.nodes() if n.name == branch.child_node)
        if node.is_tip() or node is tree.root:
            continue
        descendants = tree.descendants(node.name)
        duration = branch.older_age - branch.younger_age
        if len(descendants) >= 3 and duration >= 4:
            candidates.append((duration, branch.branch_id))
    if not candidates:
        return tree.root.name, max(0, int(tree.root_extension // 2))
    duration, branch_id = min(candidates)
    return branch_id, max(1, min(int(duration // 2), int(duration)))


def main() -> None:
    params = {}
    with (STUDY_ROOT / "datasets/processed/mai30_replicate_parameters.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            params[(row["condition"], row["replicate"])] = row
    config_dir = STUDY_ROOT / "datasets/configs/mai30_msrc_pilot"
    config_dir.mkdir(parents=True, exist_ok=True)
    seed_rows = []
    for index, (condition, replicate) in enumerate(SELECTED):
        row = params[(condition, replicate)]
        mapping = map_simphy_haploid_to_msrc(int(float(row["haploid_effective_population_size"])))
        scaled_exact = mapping.msrc_diploid_ne_exact * SCALE
        n_scaled = max(1, int(round(scaled_exact)))
        newick = scale_newick((STUDY_ROOT / row["species_tree"]).read_text().strip(), SCALE)
        root_extension = 5 * 2 * n_scaled
        tree = SpeciesTree(newick, n_scaled, root_extension, ultrametric_tolerance=1e-6)
        origin_branch, origin_time = choose_origin(tree)
        for treatment_index, (treatment, fraction) in enumerate(TREATMENTS.items()):
          for frequency_index, (frequency_label, initial_frequency) in enumerate(INITIAL_FREQUENCIES.items()):
            name = f"{condition}__replicate_{replicate}__{treatment}__{frequency_label}"
            path = config_dir / f"{name}.yaml"
            initial_copy_count = int(round(2 * n_scaled * initial_frequency))
            config = {
                "mode": "arbitrary_tree",
                "seed": 20261011 + index * 100000 + treatment_index,
                "num_loci": N_LOCI,
                "species_tree": {
                    "newick": newick,
                    "root_extension": root_extension,
                    "default_effective_population_size": n_scaled,
                    "ultrametric_tolerance": 1e-6,
                },
                "rearrangement": {
                    "id": "mai_pilot_inv",
                    "type": "inversion",
                    "origin_branch": origin_branch,
                    "origin_time_from_branch_start": origin_time,
                    "initial_copy_count": initial_copy_count,
                    "selection": {"model": "genic", "coefficient": 0.0},
                },
                "recombination": {
                    "baseline_rate": 0.02,
                    "effective_cross_arrangement_fraction": fraction,
                },
                "genome": {
                    "length_bp": N_LOCI,
                    "rearrangement_interval": {"start": 0, "end": N_LOCI, "id": "mai_pilot_inv"},
                    "breakpoint_rate_per_bp": 0.0,
                },
                "population_process": {"model": "wright_fisher"},
                "ordinary_msc_outside": True,
                "output": {
                    "directory": str(STUDY_ROOT / "datasets/processed/mai30_msrc_pilot" / name),
                    "plots": False,
                },
                "provenance": {
                    "source": "Mai et al. 2017 D1 fixed backbone",
                    "condition": condition,
                    "replicate": replicate,
                    "computational_scale_factor": SCALE,
                    "scale_interpretation": "coalescent-timescale-preserving approximation; not original-scale MSRC",
                    "simphy_haploid_ne": mapping.simphy_haploid_ne,
                    "msrc_diploid_ne_exact": mapping.msrc_diploid_ne_exact,
                    "msrc_diploid_ne_integer": mapping.msrc_diploid_ne_integer,
                    "scaled_effective_population_size_exact": scaled_exact,
                    "scaled_effective_population_size_integer": n_scaled,
                    "scaled_population_rounding_relative_error": (n_scaled - scaled_exact) / scaled_exact,
                    "initial_frequency_requested": initial_frequency,
                    "initial_copy_count": initial_copy_count,
                    "cross_arrangement_fraction": fraction,
                    "independent_history_replicates_planned": N_HISTORY_REPLICATES,
                    "gene_trees_per_history": N_LOCI,
                },
            }
            with path.open("w") as handle:
                yaml.safe_dump(config, handle, sort_keys=False)
            for history_index in range(N_HISTORY_REPLICATES):
                seed_rows.append({
                    "config": str(path.relative_to(STUDY_ROOT)),
                    "treatment": treatment,
                    "history_index": history_index,
                    "seed": config["seed"] + history_index * 1000,
                    "initial_frequency": initial_frequency,
                    "execution_command": f"python -m msrcsim.arbitrary_cli --config {path.relative_to(STUDY_ROOT)}",
                })
    with (config_dir / "SEEDS.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(seed_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(seed_rows)
    (config_dir / "README.md").write_text(
        "These configs are prepared but intentionally not executed. They use "
        "c=0.001, 20 loci per history, and three independent histories per "
        "treatment. The MSC control is generated by the direct MSC wrapper, "
        "because the arbitrary-tree runner always constructs a forward history.\n"
    )
    print(f"wrote {len(SELECTED) * len(TREATMENTS) * len(INITIAL_FREQUENCIES)} configs to {config_dir}")


if __name__ == "__main__":
    main()
