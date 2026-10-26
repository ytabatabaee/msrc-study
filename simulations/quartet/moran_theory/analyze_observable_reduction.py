"""Small exact diagnostics for unconditional observable reduction.

The principal diagnostic is the killed two-lineage Moran count process.  With
zero switching and a pair in background 0, unresolved survival is governed by
Q_M - diag(c_0(k)).  Its interior spectrum supplies independent exponentials
whose number grows with M.  This is a structural obstruction to a universal
fixed-dimensional linear closure, not a Monte Carlo result.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from .moran_quartet_operator import moran_generator
from .full_species_tree import full_species_tree_quartet_probabilities


def killed_pair_spectrum(M: int) -> np.ndarray:
    Q = moran_generator(M).toarray()
    c0 = np.array([1.0 / max(M - k, 1) for k in range(M + 1)])
    # Boundary counts are absorbing in the Moran process.  The killed
    # interior dynamics relevant to a segregating branch are the principal
    # block on k=1,...,M-1.
    return np.linalg.eigvals(Q[1:M, 1:M] - np.diag(c0[1:M]))


def run(output: Path) -> list[dict[str, object]]:
    rows = []
    for M in (4, 6, 8, 10, 12, 16, 20, 24, 32):
        eig = killed_pair_spectrum(M)
        eig = np.sort_complex(eig)
        min_gap = float(np.min(np.abs(np.diff(eig)))) if len(eig) > 1 else np.nan
        rows.append({"diagnostic": "killed_pair_interior_modes", "M": M,
                     "modes": len(eig), "minimum_eigenvalue_gap": min_gap})
    M = 10
    fixed = np.eye(M + 1)[5]
    mixture = np.zeros(M + 1); mixture[[2, 8]] = .5
    kwargs = dict(M_root=M, M_left=M, M_right=M, M_terminal=(M,) * 4,
                  rho=.17, root_length=.5, sister_lengths=(.3, .1),
                  terminal_lengths=(.2, .4, .1, .3), selection=0.0)
    q_fixed = full_species_tree_quartet_probabilities(pi_old=fixed, **kwargs).quartet_probabilities
    q_mix = full_species_tree_quartet_probabilities(pi_old=mixture, **kwargs).quartet_probabilities
    rows.append({"diagnostic": "same_mean_species_tree_difference", "M": M,
                 "modes": np.nan, "minimum_eigenvalue_gap": float(np.max(np.abs(q_fixed - q_mix)))})
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("diagnostic", "M", "modes", "minimum_eigenvalue_gap"), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results" / "observable_reduction_diagnostics.tsv")
    args = parser.parse_args(); run(args.output); print(args.output)


if __name__ == "__main__":
    main()
