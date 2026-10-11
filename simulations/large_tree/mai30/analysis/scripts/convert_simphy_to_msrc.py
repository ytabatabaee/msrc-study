#!/usr/bin/env python3
"""Convert an authenticated Mai/SimPhy record into an msrc-sim config.

This converter deliberately requires demographic values as explicit inputs. It
does not infer effective population sizes, generation time, root extension, or
units from condition names or from Newick branch lengths.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import yaml
from population_mapping import map_simphy_haploid_to_msrc


def read_branch_parameters(path: Path | None) -> dict[str, dict[str, int]]:
    if path is None:
        return {}
    rows: dict[str, dict[str, int]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            branch = row["branch_id"]
            if not row.get("effective_population_size"):
                raise ValueError(f"Missing effective_population_size for branch {branch}")
            rows[branch] = {"effective_population_size": int(row["effective_population_size"])}
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--species-tree", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--condition", required=True)
    parser.add_argument("--replicate", required=True)
    parser.add_argument("--simphy-haploid-ne", type=int, required=True)
    parser.add_argument("--root-extension-generations", type=float, required=True)
    parser.add_argument("--branch-parameters", type=Path)
    parser.add_argument("--origin-branch", required=True)
    parser.add_argument("--origin-time-from-branch-start", type=int, required=True)
    parser.add_argument("--initial-copy-count", type=int, required=True)
    parser.add_argument("--num-loci", type=int, default=500)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--baseline-rate", type=float, required=True)
    parser.add_argument("--effective-cross-arrangement-fraction", type=float, required=True)
    args = parser.parse_args()

    newick = args.species_tree.read_text().strip()
    if not newick.endswith(";"):
        raise ValueError("Species tree must be Newick and end with ';'")
    if args.num_loci <= 0:
        raise ValueError("--num-loci must be positive")
    if args.simphy_haploid_ne <= 0:
        raise ValueError("--simphy-haploid-ne must be positive")
    mapping = map_simphy_haploid_to_msrc(args.simphy_haploid_ne)

    config = {
        "mode": "arbitrary_tree",
        "seed": args.seed,
        "num_loci": args.num_loci,
        "species_tree": {
            "newick": newick,
            "root_extension": args.root_extension_generations,
            "default_effective_population_size": mapping.msrc_diploid_ne_integer,
            "branch_parameters": read_branch_parameters(args.branch_parameters),
        },
        "rearrangement": {
            "id": "mai_inv_1",
            "type": "inversion",
            "origin_branch": args.origin_branch,
            "origin_time_from_branch_start": args.origin_time_from_branch_start,
            "initial_copy_count": args.initial_copy_count,
            "selection": {"model": "genic", "coefficient": 0.0},
        },
        "recombination": {
            "baseline_rate": args.baseline_rate,
            "effective_cross_arrangement_fraction": args.effective_cross_arrangement_fraction,
        },
        "output": {"directory": f"datasets/processed/{args.condition}/replicate_{args.replicate}"},
        "provenance": {
            "source": "Mai et al. 2017 D1 SimPhy dataset",
            "condition": args.condition,
            "replicate": args.replicate,
            "branch_length_units": "generations; verified from SimPhy log before use",
            "population_size_convention": "SimPhy haploid N_H mapped to msrc-sim diploid N_H/2",
            "simphy_haploid_effective_population_size": mapping.simphy_haploid_ne,
            "msrc_diploid_effective_population_size_exact": mapping.msrc_diploid_ne_exact,
            "msrc_diploid_effective_population_size_integer": mapping.msrc_diploid_ne_integer,
            "integer_mapping_relative_rate_error": mapping.relative_rate_error,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)


if __name__ == "__main__":
    main()
