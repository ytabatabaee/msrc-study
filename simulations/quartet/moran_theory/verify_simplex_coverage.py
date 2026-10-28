"""Certificates and adversarial checks for quartet-simplex coverage claims."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .two_regime_species_tree import RegimeBranch, TwoRegimeModel, balanced_two_regime_quartet_probabilities


def deterministic_pattern_kernel(bit: int) -> np.ndarray:
    """Force a daughter lineage into background ``bit`` from either regime."""
    return np.array([[1.0, 0.0], [1.0, 0.0]]) if bit == 0 else np.array([[0.0, 1.0], [0.0, 1.0]])


def crossed_pattern_quartet(pattern: str, root_duration: float) -> np.ndarray:
    """M=1, p_L=0, p_H=1, rho=0, zero lower branches, fixed pattern.

    For patterns 0011, 0101, and 0110, the matching pair of backgrounds has
    two possible coalescing pairs. Each pair has rate one. If neither pair
    coalesces during the finite root branch, the ancestral unresolved quartet
    contributes uniformly. Hence the favored topology has
    1 - (2/3) exp(-2T), and the other two have (1/3) exp(-2T).
    """
    if pattern not in ("0011", "0101", "0110") or root_duration < 0:
        raise ValueError("pattern must be one of the three paired patterns")
    b = lambda t: RegimeBranch(t, 0.0, 0.0, 1)
    ks = [deterministic_pattern_kernel(int(x)) for x in pattern]
    result = balanced_two_regime_quartet_probabilities(
        model=TwoRegimeModel(0.0, 1.0, 0.0), root_branch=b(root_duration),
        left_pair_branch=b(0), right_pair_branch=b(0),
        terminal_branches=(b(0), b(0), b(0), b(0)), initial_regime=np.array([1.0, 0.0]),
        H_root_left=np.eye(2), H_root_right=np.eye(2),
        H_left_1=ks[0], H_left_2=ks[1], H_right_3=ks[2], H_right_4=ks[3])
    return result.quartet_probabilities


def theoretical_edge_vector(favored: int, T: float) -> np.ndarray:
    q = np.full(3, np.exp(-2 * T) / 3)
    q[favored] = 1 - 2 * np.exp(-2 * T) / 3
    return q


def run(output: Path):
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for pattern, favored in (("0011", 0), ("0101", 1), ("0110", 2)):
        for T in (0.0, 1.0, 4.0, 8.0, 12.0):
            observed = crossed_pattern_quartet(pattern, T)
            expected = theoretical_edge_vector(favored, T)
            rows.append((pattern, favored + 1, T, *observed, float(np.max(np.abs(observed - expected)))))
    with (output / "edge_certificates.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t")
        w.writerow(("pattern", "favored_topology", "T_root", "q1", "q2", "q3", "max_formula_error")); w.writerows(rows)
    return rows


if __name__ == "__main__":
    run(Path(__file__).parent / "results" / "simplex_coverage")
