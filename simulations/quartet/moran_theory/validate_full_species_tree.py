"""Validate the full finite-M species-tree recursion against msrc-sim.

The Monte Carlo unit is one independently simulated forward Moran history and
one genealogy conditional on that history. Histories are never reused as if
they were independent loci.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

try:
    from .full_species_tree import full_species_tree_quartet_probabilities
except ImportError:
    from full_species_tree import full_species_tree_quartet_probabilities


def _write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        if not rows:
            return
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _se(probability: float, n: int) -> float:
    return float(np.sqrt(max(probability * (1.0 - probability), 0.0) / n))


def _direct_cell(msrc_sim: Path, M: int, k0: int, rho: float, T_L: float, T_A: float, reps: int, seed: int):
    sys.path.insert(0, str(msrc_sim))
    from msrcsim.moran import simulate_frequency_history
    from msrcsim.rearrangement import Rearrangement
    from msrcsim.species_tree import SpeciesTree
    from msrcsim.structured_coalescent import simulate_genealogy

    tree = SpeciesTree(
        f"((1:0,2:0)AB:{T_L},(3:0,4:0)CD:{T_L})ROOT;",
        M // 2,
        T_A,
    )
    rr = Rearrangement("quartet", "inversion", "ROOT", 0, k0, 0.0)
    rng = np.random.default_rng(seed)
    topology_counts = np.zeros(3, dtype=int)
    sister_q1 = 0
    for locus in range(reps):
        history = simulate_frequency_history(tree, rr, rng)
        sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
        result = simulate_genealogy(locus, tree, history, sampled, rho, 1.0, rng, record_events=False)
        topology_counts[int(result.topology_index)] += 1
        if any(coal.population in {"AB", "CD"} for coal in result.coalescences):
            sister_q1 += 1
    return topology_counts, sister_q1


def _analytic(M: int, k0: int, rho: float, T_L: float, T_A: float):
    pi = np.zeros(M + 1)
    pi[k0] = 1.0
    result = full_species_tree_quartet_probabilities(
        M_root=M,
        M_left=M,
        M_right=M,
        M_terminal=(M, M, M, M),
        rho=rho,
        root_length=T_A,
        sister_lengths=(T_L, T_L),
        terminal_lengths=(0.0, 0.0, 0.0, 0.0),
        pi_old=pi,
    )
    return result


def run_validation(args: argparse.Namespace):
    rows = []
    for M in args.M:
        for fraction in args.initial_fraction:
            k0 = int(round(M * fraction))
            for rho in args.rho:
                for T_L in args.sister_length:
                    for T_A in args.root_length:
                        analytic = _analytic(M, k0, rho, T_L, T_A)
                        counts, sister_count = _direct_cell(
                            args.msrc_sim, M, k0, rho, T_L, T_A, args.replicates,
                            args.seed + M + k0 + round(1000 * rho) + round(100 * T_L) + round(10 * T_A),
                        )
                        observed = counts / args.replicates
                        for r in range(3):
                            se = _se(observed[r], args.replicates)
                            error = observed[r] - analytic.quartet_probabilities[r]
                            rows.append({
                                "M": M, "k0": k0, "initial_fraction": fraction,
                                "rho": rho, "T_L": T_L, "T_A": T_A,
                                "quantity": f"q{r + 1}", "analytic": analytic.quartet_probabilities[r],
                                "empirical": observed[r], "absolute_error": abs(error),
                                "mc_se": se, "standardized_residual": error / se if se else np.nan,
                                "replicates": args.replicates, "history_unit": "independent_forward_history",
                            })
                        observed_sister = sister_count / args.replicates
                        expected_sister = analytic.diagnostics["sister_q1"]
                        se = _se(observed_sister, args.replicates)
                        error = observed_sister - expected_sister
                        rows.append({
                            "M": M, "k0": k0, "initial_fraction": fraction,
                            "rho": rho, "T_L": T_L, "T_A": T_A,
                            "quantity": "sister_q1", "analytic": expected_sister,
                            "empirical": observed_sister, "absolute_error": abs(error),
                            "mc_se": se, "standardized_residual": error / se if se else np.nan,
                            "replicates": args.replicates, "history_unit": "independent_forward_history",
                        })
    finite = np.asarray([abs(float(r["standardized_residual"])) for r in rows if np.isfinite(r["standardized_residual"])])
    errors = np.asarray([float(r["absolute_error"]) for r in rows])
    summary = {
        "cells": len(args.M) * len(args.initial_fraction) * len(args.rho) * len(args.sister_length) * len(args.root_length),
        "replicates_per_cell": args.replicates,
        "rmse": float(np.sqrt(np.mean(errors**2))),
        "maximum_absolute_error": float(np.max(errors)),
        "fraction_within_2_se": float(np.mean(finite <= 2.0)),
        "fraction_within_3_se": float(np.mean(finite <= 3.0)),
    }
    return rows, summary


def sensitivity_rows(args: argparse.Namespace) -> list[dict]:
    M = args.M[0]
    k0 = int(round(M * args.initial_fraction[0]))
    pi = np.zeros(M + 1)
    pi[k0] = 1.0
    rho_grid = np.linspace(0.0, args.sensitivity_max_rho, args.sensitivity_points)
    qs = []
    for rho in rho_grid:
        qs.append(_analytic(M, k0, float(rho), args.sister_length[0], args.root_length[0]).quartet_probabilities)
    qs = np.asarray(qs)
    rows = []
    for i, rho in enumerate(rho_grid):
        derivative = np.gradient(qs[:, 0], rho_grid)[i] if len(rho_grid) > 1 else np.nan
        rows.append({"rho": rho, "q1": qs[i, 0], "q2": qs[i, 1], "q3": qs[i, 2], "dq1_drho": derivative})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--msrc-sim", type=Path, default=Path("/Users/ytabatabaee/Desktop/msrc-sim"))
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results")
    parser.add_argument("--replicates", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261021)
    parser.add_argument("--M", type=int, nargs="+", default=[20, 40])
    parser.add_argument("--initial-fraction", type=float, nargs="+", default=[0.25, 0.5, 0.75])
    parser.add_argument("--rho", type=float, nargs="+", default=[0.0, 0.1, 0.5])
    parser.add_argument("--sister-length", type=float, nargs="+", default=[0.25, 0.5])
    parser.add_argument("--root-length", type=float, nargs="+", default=[0.25, 0.5])
    parser.add_argument("--sensitivity-max-rho", type=float, default=1.0)
    parser.add_argument("--sensitivity-points", type=int, default=21)
    args = parser.parse_args()
    rows, summary = run_validation(args)
    _write_tsv(args.output / "full_species_tree_validation.tsv", rows)
    _write_tsv(args.output / "rho_sensitivity.tsv", sensitivity_rows(args))
    (args.output / "full_species_tree_validation_summary.tsv").write_text(
        "metric\tvalue\n" + "\n".join(f"{key}\t{value}" for key, value in summary.items()) + "\n"
    )
    print(summary)


if __name__ == "__main__":
    main()
