#!/usr/bin/env python3
"""Create the planned 300-replicate, three-treatment Mai30 manifest.

This writes metadata only.  It does not launch simulations.  The manifest is
the stable handoff between the authenticated SimPhy parameters and the
resumable SLURM collection runner.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import yaml

STUDY_ROOT = Path(__file__).resolve().parents[2]
MSRC_ROOT = STUDY_ROOT.parents[3] / "msrc-sim"
if "MSRC_SIM_ROOT" in __import__("os").environ:
    MSRC_ROOT = Path(__import__("os").environ["MSRC_SIM_ROOT"])
sys.path.insert(0, str(MSRC_ROOT))
from msrcsim import SpeciesTree  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from population_mapping import map_simphy_haploid_to_msrc  # noqa: E402


CONDITIONS = (
    "outgroup.0.species.0.15.genes.0.15",
    "outgroup.0.species.1.5.genes.1.5",
    "outgroup.0.species.5.genes.5",
)
TREATMENTS = (
    ("msrc_no_suppression", 1.0),
    ("msrc_moderate_suppression", 0.5),
    ("msrc_strong_suppression", 0.1),
)
SCALE = 0.001
HISTORIES = 25
LOCI_PER_HISTORY = 20
FREQUENCY = 0.5


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


def relative_source(path: str) -> str:
    candidate = Path(path)
    try:
        return str(candidate.relative_to(STUDY_ROOT))
    except ValueError:
        return str(candidate)


def main() -> None:
    parameter_path = STUDY_ROOT / "datasets/processed/mai30_replicate_parameters.tsv"
    params = {}
    with parameter_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["condition"] in CONDITIONS:
                params[(row["condition"], row["replicate"])] = row

    rows = []
    task_id = 0
    for condition_index, condition in enumerate(CONDITIONS):
        for replicate_number in range(1, 101):
            replicate = f"{replicate_number:03d}"
            row = params[(condition, replicate)]
            mapping = map_simphy_haploid_to_msrc(int(float(row["haploid_effective_population_size"])))
            n_scaled = max(1, int(round(mapping.msrc_diploid_ne_exact * SCALE)))
            species_path = STUDY_ROOT / row["species_tree"]
            newick = species_path.read_text().strip()
            root_extension = 10 * n_scaled
            tree = SpeciesTree(newick, n_scaled, root_extension, ultrametric_tolerance=1e-6)
            origin_branch, origin_time = choose_origin(tree)
            for treatment_index, (treatment, fraction) in enumerate(TREATMENTS):
                collection = f"{condition}__replicate_{replicate}__{treatment}__freq_0.5"
                rows.append({
                    "task_id": task_id,
                    "condition": condition,
                    "replicate": replicate,
                    "treatment": treatment,
                    "cross_arrangement_fraction": fraction,
                    "initial_frequency": FREQUENCY,
                    "n_histories": HISTORIES,
                    "loci_per_history": LOCI_PER_HISTORY,
                    "scale_factor": SCALE,
                    "simphy_haploid_ne": mapping.simphy_haploid_ne,
                    "msrc_diploid_ne_exact": mapping.msrc_diploid_ne_exact,
                    "scaled_effective_population_size": n_scaled,
                    "scaled_root_extension_generations": root_extension,
                    "origin_branch": origin_branch,
                    "origin_time_from_branch_start": origin_time,
                    "seed_base": 20300000 + condition_index * 10000000 + replicate_number * 10000 + treatment_index,
                    "species_tree_path": relative_source(row["species_tree"]),
                    "collection_id": collection,
                    "output_root": "datasets/processed/mai30_production_runs",
                    "status": "planned",
                })
                task_id += 1

    if len(rows) != 900:
        raise ValueError(f"expected 900 collection tasks, found {len(rows)}")
    out = STUDY_ROOT / "datasets/processed/mai30_production_manifest.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    (out.with_suffix(".yaml")).write_text(yaml.safe_dump({
        "conditions": list(CONDITIONS),
        "replicates_per_condition": 100,
        "treatments": {name: fraction for name, fraction in TREATMENTS},
        "scale_factor": SCALE,
        "initial_frequency": FREQUENCY,
        "histories_per_collection": HISTORIES,
        "loci_per_history": LOCI_PER_HISTORY,
        "total_collection_tasks": len(rows),
        "total_gene_trees_if_complete": len(rows) * HISTORIES * LOCI_PER_HISTORY,
        "execution": "planned only; submit SLURM array after pilot/resource review",
    }, sort_keys=False))
    print(f"wrote {out} ({len(rows)} collection tasks)")


if __name__ == "__main__":
    main()
