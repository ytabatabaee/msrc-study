"""Exact finite-M Moran-integrated quartet law for ``(((1,2),3),4)``.

The implementation follows the hierarchical species-tree model used by the
balanced recursion: daughter counts are conditionally independent given their
shared parental count, and genealogy states are propagated with the transpose
dual operators.  The three-lineage operator stops when the first pair among
1, 2, and 3 coalesces, because lineage 4 cannot alter the resulting unrooted
quartet split.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse.linalg import expm_multiply

try:
    from .full_species_tree import inheritance_kernel, _validate_pi
    from .moran_quartet_operator import (
        BACKGROUND_INDEX,
        BACKGROUND_STATES,
        G1,
        PAIR_BACKGROUND_STATES,
        PATTERNS,
        TRIPLE_BACKGROUND_INDEX,
        dual_lift,
        joint_operator_n,
        moran_generator,
        residual_vectors,
    )
except ImportError:  # pragma: no cover
    from full_species_tree import inheritance_kernel, _validate_pi
    from moran_quartet_operator import (
        BACKGROUND_INDEX,
        BACKGROUND_STATES,
        G1,
        PAIR_BACKGROUND_STATES,
        PATTERNS,
        TRIPLE_BACKGROUND_INDEX,
        dual_lift,
        joint_operator_n,
        moran_generator,
        residual_vectors,
    )


TRIPLE_TRANSIENT_COUNT = 8
TRIPLE_STATE_COUNT = 11
PAIR_SURVIVOR_COUNT = 4


@dataclass(frozen=True)
class UnbalancedSpeciesTreeResult:
    """Results and diagnostics for the pectinate quartet recursion."""

    quartet_probabilities: np.ndarray
    early_topology_contributions: np.ndarray
    root_topology_contributions: np.ndarray
    root_entry_unresolved: np.ndarray
    diagnostics: dict[str, float]


def _one_lineage_terminal_values(M: int) -> np.ndarray:
    out = np.zeros((M + 1) * 2, dtype=float)
    for k in range(M + 1):
        p = k / M
        out[2 * k:2 * k + 2] = (1.0 - p, p)
    return out


def _one_lineage_at_node(
    M_node: int,
    M_terminal: int,
    terminal_length: float,
    rho: float,
    selection: float,
) -> np.ndarray:
    terminal = dual_lift(
        M_terminal, rho, terminal_length,
        _one_lineage_terminal_values(M_terminal), 1, selection=selection,
    ).reshape(M_terminal + 1, 2)
    return inheritance_kernel(M_node, M_terminal) @ terminal


def _pair_at_node(
    M_node: int,
    M_sister: int,
    M_terminal: tuple[int, int],
    terminal_lengths: tuple[float, float],
    sister_length: float,
    rho: float,
    selection: float,
) -> np.ndarray:
    """Pair state at the younger end of the intermediate branch.

    The returned array is conditional on the count at that node and has the
    six G2 states ``00,01,10,11,C0,C1``.
    """
    u1 = dual_lift(
        M_terminal[0], rho, terminal_lengths[0],
        _one_lineage_terminal_values(M_terminal[0]), 1, selection=selection,
    ).reshape(M_terminal[0] + 1, 2)
    u2 = dual_lift(
        M_terminal[1], rho, terminal_lengths[1],
        _one_lineage_terminal_values(M_terminal[1]), 1, selection=selection,
    ).reshape(M_terminal[1] + 1, 2)
    H1 = inheritance_kernel(M_sister, M_terminal[0])
    H2 = inheritance_kernel(M_sister, M_terminal[1])
    young = np.zeros((M_sister + 1, 6), dtype=float)
    for ks in range(M_sister + 1):
        for b1, b2 in PAIR_BACKGROUND_STATES:
            young[ks, 2 * b1 + b2] = (H1[ks] @ u1[:, b1]) * (H2[ks] @ u2[:, b2])
    old = dual_lift(
        M_sister, rho, sister_length, young.reshape(-1), 2,
        selection=selection,
    ).reshape(M_sister + 1, 6)
    return inheritance_kernel(M_node, M_sister) @ old


def _three_lineage_at_root(
    M_root: int,
    M_intermediate: int,
    pair_at_node: np.ndarray,
    lineage3_at_node: np.ndarray,
    intermediate_length: float,
    rho: float,
    selection: float,
) -> np.ndarray:
    """Return G3 states conditional on the count at the root node."""
    young = np.zeros((M_intermediate + 1, TRIPLE_STATE_COUNT), dtype=float)
    for k in range(M_intermediate + 1):
        # A sister-pair coalescence already fixes 12|34.  It is seeded into
        # the G3 absorbing state so the intermediate count process is retained.
        young[k, 8] = pair_at_node[k, 4] + pair_at_node[k, 5]
        for pair_index, pair_bits in enumerate(PAIR_BACKGROUND_STATES[:PAIR_SURVIVOR_COUNT]):
            for b3 in (0, 1):
                state = pair_bits + (b3,)
                young[k, TRIPLE_BACKGROUND_INDEX[state]] = (
                    pair_at_node[k, pair_index] * lineage3_at_node[k, b3]
                )
    old = dual_lift(
        M_intermediate, rho, intermediate_length, young.reshape(-1), 3,
        selection=selection,
    ).reshape(M_intermediate + 1, TRIPLE_STATE_COUNT)
    return inheritance_kernel(M_root, M_intermediate) @ old


def unbalanced_species_tree_quartet_probabilities(
    *,
    M_root: int,
    M_intermediate: int,
    M_sister: int,
    M_terminal: tuple[int, int, int, int],
    rho: float,
    root_length: float,
    intermediate_length: float,
    sister_length: float,
    terminal_lengths: tuple[float, float, float, float],
    pi_old: np.ndarray,
    selection: float = 0.0,
) -> UnbalancedSpeciesTreeResult:
    """Evaluate the exact hierarchical finite-M law for ``(((1,2),3),4)``.

    ``root_length`` is the older root ancestral branch, ``intermediate_length``
    is the branch ancestral to 1,2,3, and ``sister_length`` is the branch
    ancestral to 1,2.  Durations use the simulator's generation scale, so the
    same-background pair rate is ``1/M`` at count boundaries.
    """
    sizes = (M_root, M_intermediate, M_sister, *M_terminal)
    if any(int(M) <= 0 for M in sizes):
        raise ValueError("all population chromosome counts must be positive")
    if rho < 0 or any(float(x) < 0 for x in (root_length, intermediate_length, sister_length, *terminal_lengths)):
        raise ValueError("rho and branch lengths must be nonnegative")
    if len(M_terminal) != 4 or len(terminal_lengths) != 4:
        raise ValueError("M_terminal and terminal_lengths must have length four")
    pi = _validate_pi(pi_old, M_root)

    pair_node = _pair_at_node(
        M_intermediate, M_sister, M_terminal[:2], terminal_lengths[:2],
        sister_length, rho, selection,
    )
    lineage3_node = _one_lineage_at_node(
        M_intermediate, M_terminal[2], terminal_lengths[2], rho, selection,
    )
    triple_root = _three_lineage_at_root(
        M_root, M_intermediate, pair_node, lineage3_node,
        intermediate_length, rho, selection,
    )
    lineage4_root = _one_lineage_at_node(
        M_root, M_terminal[3], terminal_lengths[3], rho, selection,
    )

    # Count distribution at the younger end of the root branch, used only for
    # events already resolved below the root branch.
    pi_root = np.asarray(expm_multiply(
        moran_generator(M_root, selection).transpose() * root_length, pi,
    ))
    early = pi_root @ triple_root[:, 8:]

    # Build the four-lineage unresolved entry vector at the root node.
    F = np.zeros((M_root + 1, 19), dtype=float)
    for k in range(M_root + 1):
        for z3, i3 in TRIPLE_BACKGROUND_INDEX.items():
            if i3 >= 8:
                continue
            for b4 in (0, 1):
                z4 = z3 + (b4,)
                F[k, BACKGROUND_INDEX[z4]] = triple_root[k, i3] * lineage4_root[k, b4]

    A4 = joint_operator_n(M_root, rho, 4, selection)
    propagated = expm_multiply(A4 * root_length, F.reshape(-1))
    initial = np.column_stack([np.kron(pi, h) for h in residual_vectors()])
    root = np.asarray(initial.T @ propagated)
    q = np.asarray(early + root)

    # Separate the part of q1 seeded by the sister branch for diagnostics.
    pair_coalesced_young = pair_node[:, 4] + pair_node[:, 5]
    pair_coalesced_old = expm_multiply(
        moran_generator(M_intermediate, selection).transpose() * intermediate_length,
        pair_coalesced_young,
    )
    pair_coalesced_root = inheritance_kernel(M_root, M_intermediate) @ pair_coalesced_old
    sister_q1 = float(pi_root @ pair_coalesced_root)
    unresolved_entry_mass = float(np.dot(pi_root, F[:, :16].sum(axis=1)))
    diagnostics = {
        "sister_q1": sister_q1,
        "intermediate_or_sister_early_mass": float(early.sum()),
        "root_mass": float(root.sum()),
        "root_entry_unresolved_mass": unresolved_entry_mass,
        "probability_sum": float(q.sum()),
        "minimum_probability": float(q.min()),
        "early_plus_root_error": float(q.sum() - 1.0),
    }
    return UnbalancedSpeciesTreeResult(q, early, root, F, diagnostics)


__all__ = ["UnbalancedSpeciesTreeResult", "unbalanced_species_tree_quartet_probabilities"]
