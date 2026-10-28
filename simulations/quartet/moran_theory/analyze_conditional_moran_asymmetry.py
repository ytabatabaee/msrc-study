"""Targeted exact audit of dominant-probability/asymmetry coordinates.

The calculation is intentionally a small certificate, not a global search.
It evaluates the exact finite-M conditional recursion at symmetric and
asymmetric pectinate parameter settings.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from .unbalanced_species_tree import unbalanced_species_tree_conditional_quartet_probabilities


def evaluate_pectinate(*, rho, root, intermediate, sister, terminal_lengths, M=4, K=2):
    result = unbalanced_species_tree_conditional_quartet_probabilities(
        M_root=M, M_intermediate=M, M_sister=M, M_terminal=(M,) * 4,
        rho=rho, root_length=root, intermediate_length=intermediate,
        sister_length=sister, terminal_lengths=terminal_lengths,
        pi_old=np.eye(M + 1)[K], selection=0.0,
    )
    q = np.asarray(result.conditional_quartet_probabilities[5], dtype=float)
    return q, float(result.pattern_probabilities[5])


CASES = (
    ("symmetric_zero_lower", dict(rho=0.0, root=1.0, intermediate=0.0,
                                   sister=0.0, terminal_lengths=(0.0,) * 4)),
    ("r_below_half", dict(rho=0.1, root=1.0, intermediate=5.0,
                          sister=0.0, terminal_lengths=(0.0, 0.3, 0.3, 0.0))),
    ("r_above_half", dict(rho=0.1, root=1.0, intermediate=5.0,
                           sister=0.0, terminal_lengths=(0.3, 0.0, 0.0, 0.3))),
    ("dominant_probability_long_root", dict(rho=0.1, root=5.0, intermediate=0.2,
                                             sister=0.0, terminal_lengths=(0.0, 0.1, 0.2, 0.3))),
)


def certificate_rows():
    rows = []
    for name, params in CASES:
        q, pz = evaluate_pectinate(**params)
        d = float(q[1])
        r = float(q[0] / (q[0] + q[2]))
        rows.append({
            "case": name, **params,
            "pattern": "0101", "q1": q[0], "q2": q[1], "q3": q[2],
            "d_q2": d, "r_q1_over_q1q3": r, "P_Z": pz,
            "normalization_error": float(q.sum() - 1.0),
        })
    return rows


def write_tsv(path: Path):
    rows = certificate_rows()
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ("case", "pattern", "rho", "root", "intermediate", "sister",
              "terminal_lengths", "q1", "q2", "q3", "d_q2",
              "r_q1_over_q1q3", "P_Z", "normalization_error")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row[key] for key in fields})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(
        "simulations/quartet/moran_theory/results/conditional_moran_asymmetry/"
        "asymmetry_certificates.tsv"
    ))
    args = parser.parse_args()
    write_tsv(args.output)
    print(args.output)


if __name__ == "__main__":  # pragma: no cover
    main()
