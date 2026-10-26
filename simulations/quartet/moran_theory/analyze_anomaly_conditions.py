"""Exact operator calculations supporting anomaly-condition statements."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from .moran_quartet_operator import PATTERNS
from .unbalanced_species_tree import unbalanced_species_tree_conditional_quartet_probabilities
from .validate_unbalanced_conditional import balanced_result, unbalanced_result


def base_condition(**overrides):
    row = {"condition": "base", "p0": .5, "rho": .17, "sister_length": .3,
           "intermediate_length": .2, "root_length": .5, "M_root": 10,
           "M_intermediate": 10, "M_sister": 10, "M_terminal": (10, 10, 10, 10),
           "selection": 0.0, "k0": 5}
    row.update(overrides)
    row["k0"] = int(round(row["M_root"] * row["p0"]))
    return row


def conditional_q(condition, pattern="0101", balanced=False):
    result = balanced_result(condition) if balanced else unbalanced_result(condition)
    return result.conditional_quartet_probabilities[PATTERNS.index(pattern)]


def zero_switching_bounds() -> list[dict[str, object]]:
    """Rigorous conservative sufficient bounds for crossed patterns."""
    rows = []
    for p in (.1, .25, .5, .75, .9):
        for pattern, background in (("0101", 0), ("1010", 1), ("0110", 1), ("1001", 0)):
            # In the intermediate branch, the matching discordant pair has
            # rate c_b=1/(M p_b).  If it coalesces, its topology is fixed.
            pb = p if background else 1.0 - p
            threshold = 10.0 * pb * np.log(2.0)
            rows.append({"pattern": pattern, "p": p, "matching_background": background,
                         "M": 10, "pair_rate": 1.0 / (10.0 * pb),
                         "sufficient_intermediate_length": threshold,
                         "bound_argument": "1-exp(-c_b*T_I) > exp(-c_b*T_I)"})
    return rows


def root_thresholds() -> list[dict[str, object]]:
    rows = []
    for intermediate in (0.0, .1, .2, .4, .8):
        def f(rho):
            return conditional_q(base_condition(intermediate_length=intermediate, rho=rho))[1] - conditional_q(base_condition(intermediate_length=intermediate, rho=rho))[0]
        lo, hi = 0.0, 100.0
        if f(lo) * f(hi) < 0:
            root = brentq(f, lo, hi, xtol=1e-9)
            rows.append({"parameter": "rho", "intermediate_length": intermediate,
                         "root_q2_equals_q1": root, "existence": "numerically bracketed"})
        else:
            rows.append({"parameter": "rho", "intermediate_length": intermediate,
                         "root_q2_equals_q1": "NaN", "existence": "no sign change on [0,100]"})
    return rows


def region_grid(n: int = 7) -> list[dict[str, object]]:
    rows = []
    pattern = "0101"
    values_rho = np.linspace(0.0, 10.0, n)
    values_I = np.linspace(0.0, .8, n)
    for rho in values_rho:
        for intermediate in values_I:
            condition = base_condition(rho=float(rho), intermediate_length=float(intermediate))
            q = conditional_q(condition, pattern)
            qb = conditional_q(condition, pattern, balanced=True)
            rows.append({"rho": rho, "intermediate_length": intermediate,
                         "gamma_unbalanced": max(q[1], q[2]) - q[0],
                         "delta_unbalanced": q[1] - q[2],
                         "gamma_balanced": max(qb[1], qb[2]) - qb[0],
                         "delta_balanced": qb[1] - qb[2],
                         "anomalous_unbalanced": int(max(q[1], q[2]) > q[0]),
                         "anomalous_balanced": int(max(qb[1], qb[2]) > qb[0]),
                         "pattern_probability": float(unbalanced_result(condition).pattern_probabilities[PATTERNS.index(pattern)])})
    return rows


def simplex_rows() -> list[dict[str, object]]:
    rows = []
    for intermediate in np.linspace(0.0, .8, 9):
        q = conditional_q(base_condition(intermediate_length=float(intermediate)), "0101")
        rows.append({"intermediate_length": intermediate, "q1": q[0], "q2": q[1], "q3": q[2]})
    return rows


def write(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).parent / "results")
    parser.add_argument("--grid-size", type=int, default=7)
    args = parser.parse_args()
    write(args.output_dir / "anomaly_zero_switching_bounds.tsv", zero_switching_bounds())
    write(args.output_dir / "anomaly_switching_thresholds.tsv", root_thresholds())
    write(args.output_dir / "anomaly_region_grid.tsv", region_grid(args.grid_size))
    write(args.output_dir / "anomaly_simplex_path.tsv", simplex_rows())
    print("wrote anomaly bounds, thresholds, region grid, and simplex path")


if __name__ == "__main__":
    main()
