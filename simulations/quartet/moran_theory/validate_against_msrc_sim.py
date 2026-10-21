"""Compare the finite-M operator with the authoritative msrc-sim Moran path.

The simulator is imported from the separate ``msrc-sim`` checkout.  This
script is deliberately standalone so publication-scale replicate counts can be
run without making them part of the unit-test suite.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

from moran_quartet_operator import PATTERNS, integrated_quartet_probabilities


def _write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _se(x: float, n: int) -> float:
    return float(np.sqrt(max(x * (1.0 - x), 0.0) / n))


def _scenario(msrc_sim: Path, M: int, k0: int, rho: float, T: float, reps: int, seed: int):
    sys.path.insert(0, str(msrc_sim))
    from msrcsim.moran import simulate_frequency_history
    from msrcsim.rearrangement import Rearrangement
    from msrcsim.species_tree import SpeciesTree
    from msrcsim.structured_coalescent import simulate_genealogy

    ne = M // 2
    tree = SpeciesTree("(1:0,2:0,3:0,4:0)ROOT;", ne, T)
    rr = Rearrangement("quartet", "inversion", "ROOT", 0, k0, 0.0)
    rng = np.random.default_rng(seed)
    top_counts = np.zeros(3, dtype=int)
    pattern_counts = np.zeros(16, dtype=int)
    joint_counts = np.zeros((16, 3), dtype=int)
    for i in range(reps):
        history = simulate_frequency_history(tree, rr, rng)
        sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
        pattern = "".join(str(sampled[taxon]) for taxon in tree.taxa)
        zi = PATTERNS.index(pattern)
        result = simulate_genealogy(i, tree, history, sampled, rho, 1.0, rng, record_events=False)
        ri = int(result.topology_index)
        if not 0 <= ri < 3:
            raise RuntimeError(f"msrc-sim returned invalid topology index {ri}")
        top_counts[ri] += 1
        pattern_counts[zi] += 1
        joint_counts[zi, ri] += 1
    return top_counts, pattern_counts, joint_counts


def run(args: argparse.Namespace) -> tuple[list[dict], dict]:
    if args.smoke:
        Ms, fractions, rhos, durations = (20,), (0.5,), (0.0, 0.1), (0.5,)
    else:
        Ms = tuple(args.M)
        fractions = tuple(args.fraction)
        rhos = tuple(args.rho)
        durations = tuple(args.duration)
    rows = []
    all_errors = []
    for M in Ms:
        for fraction in fractions:
            k0 = int(round(M * fraction))
            if not 0 <= k0 <= M:
                raise ValueError("initial fraction must produce a count in [0,M]")
            pi = np.zeros(M + 1)
            pi[k0] = 1.0
            for rho in rhos:
                for multiplier in durations:
                    T = multiplier * M
                    expected = integrated_quartet_probabilities(M, rho, T, pi)
                    tops, patterns, joint = _scenario(args.msrc_sim, M, k0, rho, T, args.replicates, args.seed + M + k0 + round(100 * rho) + round(10 * multiplier))
                    n = args.replicates
                    observed_q = tops / n
                    observed_pattern = patterns / n
                    for r in range(3):
                        se = _se(observed_q[r], n)
                        err = abs(observed_q[r] - expected.quartet_probabilities[r])
                        all_errors.append(err)
                        rows.append({"M": M, "k0": k0, "initial_fraction": fraction, "rho": rho, "T": T,
                                     "quantity": f"q{r+1}", "pattern": "", "observed": observed_q[r],
                                     "analytic": expected.quartet_probabilities[r], "absolute_error": err,
                                     "mc_se": se, "standardized_error": err / se if se else np.nan,
                                     "replicates": n, "pattern_observations": n})
                    for zi, z in enumerate(PATTERNS):
                        se = _se(observed_pattern[zi], n)
                        err = abs(observed_pattern[zi] - expected.pattern_probabilities[zi])
                        all_errors.append(err)
                        rows.append({"M": M, "k0": k0, "initial_fraction": fraction, "rho": rho, "T": T,
                                     "quantity": "pattern", "pattern": z, "observed": observed_pattern[zi],
                                     "analytic": expected.pattern_probabilities[zi], "absolute_error": err,
                                     "mc_se": se, "standardized_error": err / se if se else np.nan,
                                     "replicates": n, "pattern_observations": n})
                        if patterns[zi] >= args.min_pattern_observations:
                            conditional_observed = joint[zi] / patterns[zi]
                            for r in range(3):
                                se = _se(conditional_observed[r], patterns[zi])
                                err = abs(conditional_observed[r] - expected.conditional_quartet_probabilities[zi, r])
                                all_errors.append(err)
                                rows.append({"M": M, "k0": k0, "initial_fraction": fraction, "rho": rho, "T": T,
                                             "quantity": f"P(Q{r+1}|Z)", "pattern": z, "observed": conditional_observed[r],
                                             "analytic": expected.conditional_quartet_probabilities[zi, r], "absolute_error": err,
                                             "mc_se": se, "standardized_error": err / se if se else np.nan,
                                             "replicates": n, "pattern_observations": int(patterns[zi])})
    errors = np.asarray(all_errors, dtype=float)
    standardized = np.asarray([abs(float(row["standardized_error"])) for row in rows if np.isfinite(row["standardized_error"])])
    summary = {
        "cells": len(Ms) * len(fractions) * len(rhos) * len(durations),
        "replicates_per_cell": args.replicates,
        "reported_rows": len(rows),
        "rmse": float(np.sqrt(np.mean(errors**2))),
        "maximum_absolute_error": float(np.max(errors)),
        "fraction_within_2_se": float(np.mean(standardized <= 2.0)),
        "fraction_within_3_se": float(np.mean(standardized <= 3.0)),
    }
    return rows, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--msrc-sim", type=Path, default=Path("/Users/ytabatabaee/Desktop/msrc-sim"))
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results")
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--min-pattern-observations", type=int, default=25)
    parser.add_argument("--seed", type=int, default=20261021)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--M", type=int, nargs="+", default=[20, 40, 80])
    parser.add_argument("--fraction", type=float, nargs="+", default=[0.25, 0.5, 0.75])
    parser.add_argument("--rho", type=float, nargs="+", default=[0.0, 0.02, 0.1, 0.5])
    parser.add_argument("--duration", type=float, nargs="+", default=[0.25, 0.5, 1.0])
    args = parser.parse_args()
    rows, summary = run(args)
    _write_tsv(args.output / "moran_quartet_validation.tsv", rows)
    (args.output / "moran_quartet_validation_summary.tsv").write_text(
        "metric\tvalue\n" + "\n".join(f"{k}\t{v}" for k, v in summary.items()) + "\n"
    )
    print(summary)


if __name__ == "__main__":
    main()

