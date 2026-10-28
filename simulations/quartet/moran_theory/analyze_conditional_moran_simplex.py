"""Exact finite-M Moran certificates for arrangement-conditional quartets."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .full_species_tree import full_species_tree_quartet_probabilities
from .unbalanced_species_tree import unbalanced_species_tree_conditional_quartet_probabilities


def conditional_case(z: str, T: float, M: int = 10, K: int | None = None, rho: float = 0.0):
    if z not in ("0011", "0101", "0110"):
        raise ValueError("z must be one of the three paired patterns")
    K = M // 2 if K is None else K
    result = full_species_tree_quartet_probabilities(
        M_root=M, M_left=M, M_right=M, M_terminal=(M,) * 4, rho=rho,
        root_length=T, sister_lengths=(0.0, 0.0), terminal_lengths=(0.0,) * 4,
        pi_old=np.eye(M + 1)[K])
    index = int(z, 2)
    return result.conditional_quartet_probabilities[index], float(result.pattern_probabilities[index])


def run(output: Path):
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for z, favored in (("0011", 0), ("0101", 1), ("0110", 2)):
        for M, K in ((4, 1), (10, 5)):
            for T in (0.0, 1.0, 5.0, 20.0):
                q, pz = conditional_case(z, T, M, K)
                rows.append((z, M, K, T, pz, *q, int(np.argmax(q)) + 1, float(q[favored])))
    with (output / "vertex_certificates.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t")
        w.writerow(("z", "M", "K", "T_root", "P_Z", "q1", "q2", "q3", "largest_topology", "favored_q")); w.writerows(rows)
    pectinate_rows = []
    for z, favored in (("0011", 0), ("0101", 1), ("0110", 2)):
        result = unbalanced_species_tree_conditional_quartet_probabilities(
            M_root=4, M_intermediate=4, M_sister=4, M_terminal=(4,) * 4, rho=0,
            root_length=20.0, intermediate_length=0.0, sister_length=0.0,
            terminal_lengths=(0.0,) * 4, pi_old=np.eye(5)[1])
        q = result.conditional_quartet_probabilities[int(z, 2)]
        pectinate_rows.append((z, *q, result.pattern_probabilities[int(z, 2)], int(np.argmax(q)) + 1, favored + 1))
    with (output / "pectinate_vertex_certificates.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t"); w.writerow(("z", "q1", "q2", "q3", "P_Z", "largest_topology", "expected_topology")); w.writerows(pectinate_rows)
    # A local rank-2 diagnostic for one fixed observed pattern. This is only a
    # local reachability check, not a global coverage result.
    def q(x):
        result = full_species_tree_quartet_probabilities(
            M_root=10, M_left=10, M_right=10, M_terminal=(10,) * 4, rho=.1,
            root_length=float(x[0]), sister_lengths=(.2, .3),
            terminal_lengths=(float(x[1]), .1, .2, .3), pi_old=np.eye(11)[5])
        return result.conditional_quartet_probabilities[int("0101", 2), :2]
    x = np.array([1.0, .2]); h = 1e-4
    J = np.column_stack(((q(x + np.array([h, 0])) - q(x - np.array([h, 0]))) / (2 * h),
                         (q(x + np.array([0, h])) - q(x - np.array([0, h]))) / (2 * h)))
    with (output / "local_jacobian.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t"); w.writerow(("rank", "singular_values", "determinant")); w.writerow((np.linalg.matrix_rank(J), ",".join(map(str, np.linalg.svd(J, compute_uv=False))), np.linalg.det(J)))
    targets = (("center", (1 / 3, 1 / 3, 1 / 3)), ("interior", (.2, .5, .3)),
               ("near_Q2", (.05, .9, .05)), ("interior_alt", (.7, .1, .2)))
    with (output / "target_status.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t"); w.writerow(("target", "q1", "q2", "q3", "status"))
        for name, target in targets: w.writerow((name, *target, "unresolved_no_global_inverse"))
    return rows


if __name__ == "__main__":
    run(Path(__file__).parent / "results" / "conditional_moran_simplex")
