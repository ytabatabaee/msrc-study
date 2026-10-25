"""Small independent Monte Carlo validation for the pectinate recursion.

This script uses the authoritative ``msrc-sim`` checkout for simulation and
the operator implementation in this directory for expectations.  It is
intended for smoke tests and reproducible cluster jobs; it never reuses a
frequency history as if it were an independent locus.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

try:
    from .unbalanced_species_tree import unbalanced_species_tree_quartet_probabilities
except ImportError:  # pragma: no cover
    from unbalanced_species_tree import unbalanced_species_tree_quartet_probabilities


def _se(p: float, n: int) -> float:
    return float(np.sqrt(max(p * (1.0 - p), 0.0) / n))


def _tree(msrc_sim: Path, M: int, sister: float, intermediate: float, root: float, selection: float):
    sys.path.insert(0, str(msrc_sim))
    from msrcsim.species_tree import SpeciesTree

    ne = M // 2
    newick = f"(((1:0,2:0)N12:{sister},3:{sister})N123:{intermediate},4:{sister + intermediate})ROOT;"
    branch_parameters = {
        name: {"effective_population_size": ne, "selection_coefficient": selection}
        for name in ("ROOT", "N123", "N12", "1", "2", "3", "4")
    }
    return SpeciesTree(newick, ne, root, branch_parameters=branch_parameters)


def simulate_cell(msrc_sim: Path, M: int, k0: int, rho: float, sister: float,
                  intermediate: float, root: float, selection: float,
                  replicates: int, seed: int) -> np.ndarray:
    sys.path.insert(0, str(msrc_sim))
    from msrcsim.moran import simulate_frequency_history
    from msrcsim.rearrangement import Rearrangement
    from msrcsim.structured_coalescent import simulate_genealogy

    tree = _tree(msrc_sim, M, sister, intermediate, root, selection)
    rearrangement = Rearrangement("pectinate_validation", "inversion", "ROOT", 0, k0, selection)
    rng = np.random.default_rng(seed)
    counts = np.zeros(3, dtype=int)
    for locus in range(replicates):
        history = simulate_frequency_history(tree, rearrangement, rng)
        sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
        result = simulate_genealogy(locus, tree, history, sampled, rho, 1.0, rng, record_events=False)
        counts[int(result.topology_index)] += 1
    return counts


def run(args: argparse.Namespace) -> tuple[list[dict[str, object]], dict[str, object]]:
    if args.smoke:
        M_values = (10,)
        fractions = (0.5,)
        rhos = (0.0, 0.17)
        sisters = (0.3,)
        intermediates = (0.0, 0.2)
        roots = (0.0, 0.5)
        replicates = args.replicates
    else:
        M_values = tuple(args.M)
        fractions = tuple(args.fraction)
        rhos = tuple(args.rho)
        sisters = tuple(args.sister_length)
        intermediates = tuple(args.intermediate_length)
        roots = tuple(args.root_length)
        replicates = args.replicates

    rows: list[dict[str, object]] = []
    for M in M_values:
        for fraction in fractions:
            k0 = int(round(M * fraction))
            pi = np.zeros(M + 1)
            pi[k0] = 1.0
            for rho in rhos:
                for sister in sisters:
                    for intermediate in intermediates:
                        for root in roots:
                            analytic = unbalanced_species_tree_quartet_probabilities(
                                M_root=M, M_intermediate=M, M_sister=M,
                                M_terminal=(M, M, M, M), rho=rho,
                                root_length=root, intermediate_length=intermediate,
                                sister_length=sister, terminal_lengths=(0.0,) * 4,
                                pi_old=pi, selection=args.selection,
                            )
                            counts = simulate_cell(
                                args.msrc_sim, M, k0, rho, sister, intermediate, root,
                                args.selection, replicates,
                                args.seed + M + k0 + round(100 * rho) + round(100 * sister)
                                + round(100 * intermediate) + round(100 * root),
                            )
                            observed = counts / replicates
                            for r in range(3):
                                error = float(observed[r] - analytic.quartet_probabilities[r])
                                se = _se(float(observed[r]), replicates)
                                rows.append({
                                    "M": M, "fraction": fraction, "rho": rho,
                                    "sister_length": sister, "intermediate_length": intermediate,
                                    "root_length": root, "topology": f"q{r + 1}",
                                    "analytic": analytic.quartet_probabilities[r],
                                    "observed": observed[r], "absolute_error": abs(error),
                                    "mc_se": se, "standardized_error": error / se if se else np.nan,
                                    "replicates": replicates,
                                })
    residuals = np.asarray([float(row["absolute_error"]) for row in rows])
    standardized = np.asarray([
        abs(float(row["standardized_error"])) for row in rows
        if np.isfinite(row["standardized_error"])
    ])
    summary = {
        "cells": len(rows) // 3,
        "replicates_per_cell": replicates,
        "rmse": float(np.sqrt(np.mean(residuals ** 2))),
        "maximum_absolute_error": float(np.max(residuals)),
        "fraction_within_2_se": float(np.mean(standardized <= 2.0)),
        "fraction_within_3_se": float(np.mean(standardized <= 3.0)),
    }
    return rows, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--msrc-sim", type=Path, default=Path("/Users/ytabatabaee/Desktop/msrc-sim"))
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results" / "unbalanced_species_tree_validation.tsv")
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261025)
    parser.add_argument("--selection", type=float, default=0.0)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--M", type=int, nargs="+", default=[10, 20, 40])
    parser.add_argument("--fraction", type=float, nargs="+", default=[0.25, 0.5, 0.75])
    parser.add_argument("--rho", type=float, nargs="+", default=[0.0, 0.1, 0.5])
    parser.add_argument("--sister-length", type=float, nargs="+", default=[0.2, 0.5])
    parser.add_argument("--intermediate-length", type=float, nargs="+", default=[0.0, 0.3])
    parser.add_argument("--root-length", type=float, nargs="+", default=[0.0, 0.5])
    args = parser.parse_args()
    rows, summary = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    summary_path = args.output.with_name(args.output.stem + "_summary.tsv")
    summary_path.write_text("metric\tvalue\n" + "\n".join(f"{k}\t{v}" for k, v in summary.items()) + "\n")
    print(summary)


if __name__ == "__main__":
    main()
