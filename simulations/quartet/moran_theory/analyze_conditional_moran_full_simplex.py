"""Small exact certificate for the proved part of conditional coverage.

This script deliberately evaluates the vertex-closure family only.  It is not
a global optimizer and does not treat a collection of parameter settings as a
mixture model.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from .full_species_tree import full_species_tree_quartet_probabilities
from .unbalanced_species_tree import (
    unbalanced_species_tree_conditional_quartet_probabilities,
)
from .moran_quartet_operator import PATTERNS


PATTERN_TO_TOPOLOGY = {"0011": 0, "0101": 1, "0110": 2}


def vertex_closure_rows(
    *,
    M: int = 4,
    K: int = 2,
    root_lengths: tuple[float, ...] = (0.0, 1.0, 8.0),
) -> list[dict[str, float | int | str]]:
    """Evaluate the zero-lower-branch, zero-switching vertex family."""
    pi = np.zeros(M + 1)
    pi[K] = 1.0
    rows: list[dict[str, float | int | str]] = []
    for shape in ("balanced", "pectinate"):
        for pattern, selected in PATTERN_TO_TOPOLOGY.items():
            for root_length in root_lengths:
                common = dict(
                    rho=0.0,
                    root_length=root_length,
                    terminal_lengths=(0.0, 0.0, 0.0, 0.0),
                    pi_old=pi,
                    selection=0.0,
                )
                pattern_index = PATTERNS.index(pattern)
                if shape == "balanced":
                    result = full_species_tree_quartet_probabilities(
                        M_root=M,
                        M_left=M,
                        M_right=M,
                        M_terminal=(M, M, M, M),
                        sister_lengths=(0.0, 0.0),
                        **common,
                    )
                    conditional = result.conditional_quartet_probabilities
                    pattern_probability = result.pattern_probabilities[pattern_index]
                else:
                    result = unbalanced_species_tree_conditional_quartet_probabilities(
                            M_root=M, M_intermediate=M, M_sister=M,
                            M_terminal=(M, M, M, M),
                            intermediate_length=0.0, sister_length=0.0,
                            **common,
                        )
                    conditional = result.conditional_quartet_probabilities
                    pattern_probability = result.pattern_probabilities[pattern_index]
                q = np.asarray(conditional[pattern_index])
                rows.append({
                    "shape": shape,
                    "pattern": pattern,
                    "selected_topology": selected + 1,
                    "root_length": root_length,
                    "pattern_probability": float(pattern_probability),
                    "q1": float(q[0]), "q2": float(q[1]), "q3": float(q[2]),
                    "selected_q": float(q[selected]),
                    "normalization_error": float(q.sum() - 1.0),
                })
    return rows


def write_certificate(path: Path) -> None:
    rows = vertex_closure_rows()
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(
        "simulations/quartet/moran_theory/results/conditional_moran_full_simplex/"
        "vertex_closure_certificate.tsv"
    ))
    args = parser.parse_args()
    write_certificate(args.output)
    print(args.output)


if __name__ == "__main__":  # pragma: no cover
    main()
