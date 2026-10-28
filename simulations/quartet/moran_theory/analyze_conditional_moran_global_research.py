"""Exact certificate for the radial conditional Moran quartet families.

The calculation tests the symmetric zero-switching, zero-lower-branch
subfamily. It records the radial path from the uniform quartet law toward the
vertex favored by each paired arrangement pattern.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from .full_species_tree import full_species_tree_quartet_probabilities
from .moran_quartet_operator import PATTERNS
from .unbalanced_species_tree import unbalanced_species_tree_conditional_quartet_probabilities


PAIRED = {"0011": 0, "0101": 1, "0110": 2}


def radial_rows(shape: str, root_lengths=(0.0, 0.5, 1.0, 4.0, 8.0), M=4, K=2):
    pi = np.eye(M + 1)[K]
    rows = []
    for pattern, selected in PAIRED.items():
        for root in root_lengths:
            if shape == "balanced":
                result = full_species_tree_quartet_probabilities(
                    M_root=M, M_left=M, M_right=M, M_terminal=(M,) * 4,
                    rho=0.0, root_length=root, sister_lengths=(0.0, 0.0),
                    terminal_lengths=(0.0,) * 4, pi_old=pi,
                )
            elif shape == "pectinate":
                result = unbalanced_species_tree_conditional_quartet_probabilities(
                    M_root=M, M_intermediate=M, M_sister=M,
                    M_terminal=(M,) * 4, rho=0.0, root_length=root,
                    intermediate_length=0.0, sister_length=0.0,
                    terminal_lengths=(0.0,) * 4, pi_old=pi,
                )
            else:
                raise ValueError(shape)
            index = PATTERNS.index(pattern)
            q = np.asarray(result.conditional_quartet_probabilities[index])
            rows.append({
                "shape": shape, "pattern": pattern, "selected_topology": selected + 1,
                "root_length": root, "q1": q[0], "q2": q[1], "q3": q[2],
                "selected_q": q[selected],
                "nonselected_difference": float(q[(selected + 1) % 3] - q[(selected + 2) % 3]),
                "normalization_error": float(q.sum() - 1.0),
                "pattern_probability": float(result.pattern_probabilities[index]),
            })
    return rows


def write_certificate(path: Path):
    rows = radial_rows("balanced") + radial_rows("pectinate")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(
        "simulations/quartet/moran_theory/results/conditional_moran_global_research/"
        "radial_family_certificate.tsv"
    ))
    args = parser.parse_args()
    write_certificate(args.output)
    print(args.output)


if __name__ == "__main__":  # pragma: no cover
    main()
