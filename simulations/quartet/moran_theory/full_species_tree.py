"""Exact finite-M Moran MSRC quartet law for ``S=((1,2),(3,4))``.

The recursion keeps the common parental Moran count at the root speciation
event until the left and right sister-clade contributions have been multiplied.
This is the finite-population dependence induced by msrc-sim's binomial
daughter initialization.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import comb

import numpy as np
from scipy.sparse.linalg import expm_multiply

try:
    from .moran_quartet_operator import (
        BACKGROUND_STATES,
        G1,
        G2,
        G4,
        PATTERNS,
        PAIR_BACKGROUND_STATES,
        TOPOLOGY_NAMES,
        dual_lift,
        joint_operator_n,
        moran_generator,
        residual_vectors,
    )
except ImportError:  # pragma: no cover - supports direct script execution
    from moran_quartet_operator import (
        BACKGROUND_STATES,
        G1,
        G2,
        G4,
        PATTERNS,
        PAIR_BACKGROUND_STATES,
        TOPOLOGY_NAMES,
        dual_lift,
        joint_operator_n,
        moran_generator,
        residual_vectors,
    )


PAIR_STATE_NAMES = ("00", "01", "10", "11", "C0", "C1")
PAIR_SURVIVOR_INDICES = (0, 1, 2, 3)


@dataclass(frozen=True)
class FullSpeciesTreeResult:
    quartet_probabilities: np.ndarray
    joint: np.ndarray
    pattern_probabilities: np.ndarray
    conditional_quartet_probabilities: np.ndarray
    left_root_states: np.ndarray
    right_root_states: np.ndarray
    root_survivor_entry: np.ndarray
    diagnostics: dict[str, float]


def inheritance_kernel(M_parent: int, M_daughter: int) -> np.ndarray:
    """Return the msrc-sim default daughter-count kernel.

    Rows are parent counts ``k=0,...,M_parent`` and columns are daughter
    counts ``a=0,...,M_daughter``.  This is conditional independence of each
    daughter given the parental count, not marginal independence.
    """
    if M_parent <= 0 or M_daughter <= 0:
        raise ValueError("population chromosome counts must be positive")
    H = np.zeros((M_parent + 1, M_daughter + 1), dtype=float)
    for k in range(M_parent + 1):
        p = k / M_parent
        for a in range(M_daughter + 1):
            H[k, a] = comb(M_daughter, a) * p**a * (1.0 - p) ** (M_daughter - a)
    return H


def _validate_pi(pi_old: np.ndarray, M: int) -> np.ndarray:
    pi = np.asarray(pi_old, dtype=float)
    if pi.shape != (M + 1,) or np.any(pi < -1e-14) or not np.isclose(pi.sum(), 1.0):
        raise ValueError("pi_old must be a probability vector of length M_root+1")
    return pi


def _one_lineage_terminal_values(M: int, sampled_bit: int | None = None) -> np.ndarray:
    """Young-boundary values, optionally conditional on one sampled bit."""
    out = np.zeros((M + 1) * 2, dtype=float)
    for k in range(M + 1):
        p = k / M
        if sampled_bit is None:
            out[2 * k:2 * k + 2] = (1.0 - p, p)
        else:
            out[2 * k + int(sampled_bit)] = p if sampled_bit else 1.0 - p
    return out


def _pair_terminal_values(
    M_sister: int,
    M_1: int,
    M_2: int,
    terminal_lengths: tuple[float, float],
    rho: float,
    selection: float,
    terminal_bits: tuple[int, int] | None = None,
) -> np.ndarray:
    """Join two terminal-lifted lineages at the younger sister boundary."""
    u1 = dual_lift(
        M_1, rho, terminal_lengths[0], _one_lineage_terminal_values(M_1, None if terminal_bits is None else terminal_bits[0]), 1,
        selection=selection,
    ).reshape(M_1 + 1, 2)
    u2 = dual_lift(
        M_2, rho, terminal_lengths[1], _one_lineage_terminal_values(M_2, None if terminal_bits is None else terminal_bits[1]), 1,
        selection=selection,
    ).reshape(M_2 + 1, 2)
    H1 = inheritance_kernel(M_sister, M_1)
    H2 = inheritance_kernel(M_sister, M_2)
    out = np.zeros((M_sister + 1, 6), dtype=float)
    for ks in range(M_sister + 1):
        for b1, b2 in PAIR_BACKGROUND_STATES:
            out[ks, 2 * b1 + b2] = (
                H1[ks] @ u1[:, b1]
            ) * (H2[ks] @ u2[:, b2])
    return out


def _sister_root_states(
    M_root: int,
    M_sister: int,
    M_1: int,
    M_2: int,
    terminal_lengths: tuple[float, float],
    sister_length: float,
    rho: float,
    selection: float,
    terminal_bits: tuple[int, int] | None = None,
) -> np.ndarray:
    """Return left/right pair states conditional on each root-node count."""
    young = _pair_terminal_values(M_sister, M_1, M_2, terminal_lengths, rho, selection, terminal_bits)
    old = dual_lift(
        M_sister, rho, sister_length, young.reshape(-1), 2, selection=selection
    ).reshape(M_sister + 1, 6)
    return inheritance_kernel(M_root, M_sister) @ old


def _fixed_one_lineage_terminal(M: int, T: float, rho: float, k: int) -> np.ndarray:
    return np.asarray(expm_multiply(G1(k, M, rho).transpose() * T, np.array([1 - k / M, k / M])))


def _fixed_sister_root_states(
    M: int,
    terminal_lengths: tuple[float, float],
    sister_length: float,
    rho: float,
    k: int,
) -> np.ndarray:
    u1 = _fixed_one_lineage_terminal(M, terminal_lengths[0], rho, k)
    u2 = _fixed_one_lineage_terminal(M, terminal_lengths[1], rho, k)
    young = np.zeros(6)
    for b1, b2 in PAIR_BACKGROUND_STATES:
        young[2 * b1 + b2] = u1[b1] * u2[b2]
    return np.asarray(expm_multiply(G2(k, M, rho).transpose() * sister_length, young))


def _root_quartet_from_survivors(
    M_root: int,
    root_length: float,
    rho: float,
    pi_old: np.ndarray,
    survivor_values: np.ndarray,
    selection: float,
) -> np.ndarray:
    """Apply the 19-state root operator to conditional four-lineage mass."""
    A = joint_operator_n(M_root, rho, 4, selection)
    propagated = expm_multiply(A * root_length, survivor_values.reshape(-1))
    initial = np.column_stack([np.kron(pi_old, h) for h in residual_vectors()])
    return initial.T @ propagated


def full_species_tree_quartet_probabilities(
    *,
    M_root: int,
    M_left: int,
    M_right: int,
    M_terminal: tuple[int, int, int, int],
    rho: float,
    root_length: float,
    sister_lengths: tuple[float, float],
    terminal_lengths: tuple[float, float, float, float],
    pi_old: np.ndarray,
    selection: float = 0.0,
) -> FullSpeciesTreeResult:
    """Return the exact quartet law for the four-taxon species tree.

    All branch durations run forward from older to younger time for Moran
    counts. Genealogies are evaluated backward through each branch with the
    transpose/dual operators. The root count at the root speciation event is
    retained while left and right clade contributions are multiplied.
    """
    if rho < 0 or root_length < 0 or any(x < 0 for x in sister_lengths + terminal_lengths):
        raise ValueError("rho and branch lengths must be nonnegative")
    if len(M_terminal) != 4 or len(terminal_lengths) != 4 or len(sister_lengths) != 2:
        raise ValueError("terminal and sister parameter tuples have the wrong length")
    pi = _validate_pi(pi_old, M_root)
    pi_root = np.asarray(expm_multiply(moran_generator(M_root, selection).transpose() * root_length, pi))
    pair_bits = tuple("".join(str(x) for x in bits) for bits in PAIR_BACKGROUND_STATES)
    left_by_bits = {
        bits: _sister_root_states(
            M_root, M_left, M_terminal[0], M_terminal[1], terminal_lengths[:2],
            sister_lengths[0], rho, selection, tuple(int(x) for x in bits),
        )
        for bits in pair_bits
    }
    right_by_bits = {
        bits: _sister_root_states(
            M_root, M_right, M_terminal[2], M_terminal[3], terminal_lengths[2:],
            sister_lengths[1], rho, selection, tuple(int(x) for x in bits),
        )
        for bits in pair_bits
    }
    # Summing over the four possible sampled pair patterns recovers the
    # unconditional sister outputs used by the original recursion.
    left = sum(left_by_bits.values())
    right = sum(right_by_bits.values())
    left_survivors = left[:, :4].sum(axis=1)
    right_survivors = right[:, :4].sum(axis=1)

    survivor_values_by_pattern = np.zeros((M_root + 1, 19, 16), dtype=float)
    sister_q1_by_pattern = np.zeros(16, dtype=float)
    for zi, z in enumerate(PATTERNS):
        left_z = left_by_bits[z[:2]]
        right_z = right_by_bits[z[2:]]
        left_survivor_z = left_z[:, :4].sum(axis=1)
        right_survivor_z = right_z[:, :4].sum(axis=1)
        for k in range(M_root + 1):
            for il, zl in enumerate(PAIR_BACKGROUND_STATES):
                for ir, zr in enumerate(PAIR_BACKGROUND_STATES):
                    background = BACKGROUND_STATES.index(zl + zr)
                    survivor_values_by_pattern[k, background, zi] = left_z[k, il] * right_z[k, ir]
        left_total_z = left_z.sum(axis=1)
        right_total_z = right_z.sum(axis=1)
        sister_q1_by_pattern[zi] = float(
            np.dot(pi_root, left_total_z * right_total_z - left_survivor_z * right_survivor_z)
        )

    A = joint_operator_n(M_root, rho, 4, selection)
    propagated = expm_multiply(A * root_length, survivor_values_by_pattern.reshape((M_root + 1) * 19, 16))
    initial = np.column_stack([np.kron(pi, h) for h in residual_vectors()])
    root_joint = (initial.T @ propagated).T
    joint = np.asarray(root_joint, dtype=float)
    joint[:, 0] += sister_q1_by_pattern
    pattern = joint.sum(axis=1)
    q = joint.sum(axis=0)
    conditional = np.divide(
        joint,
        pattern[:, None],
        out=np.full_like(joint, np.nan),
        where=pattern[:, None] > 0,
    )
    sister_q1 = float(sister_q1_by_pattern.sum())
    diagnostics = {
        "sister_q1": sister_q1,
        "neither_sister_coalesces": float(np.dot(pi_root, left_survivors * right_survivors)),
        "left_survivor_mean": float(np.dot(pi_root, left_survivors)),
        "right_survivor_mean": float(np.dot(pi_root, right_survivors)),
    }
    survivor_entry = np.zeros((M_root + 1, 19), dtype=float)
    for k in range(M_root + 1):
        for il, zl in enumerate(PAIR_BACKGROUND_STATES):
            for ir, zr in enumerate(PAIR_BACKGROUND_STATES):
                survivor_entry[k, BACKGROUND_STATES.index(zl + zr)] = left[k, il] * right[k, ir]
    return FullSpeciesTreeResult(q, joint, pattern, conditional, left, right, survivor_entry, diagnostics)


def fixed_p_half_species_tree_quartet_probabilities(
    *,
    M: int,
    rho: float,
    root_length: float,
    sister_lengths: tuple[float, float],
    terminal_lengths: tuple[float, float, float, float],
) -> np.ndarray:
    """Fixed-count ``p=1/2`` special case for regression against old theory."""
    if M % 2:
        raise ValueError("M must be even")
    k = M // 2
    left = _fixed_sister_root_states(M, terminal_lengths[:2], sister_lengths[0], rho, k)
    right = _fixed_sister_root_states(M, terminal_lengths[2:], sister_lengths[1], rho, k)
    survivor = np.zeros(19)
    for il, zl in enumerate(PAIR_BACKGROUND_STATES):
        for ir, zr in enumerate(PAIR_BACKGROUND_STATES):
            survivor[BACKGROUND_STATES.index(zl + zr)] = left[il] * right[ir]
    root = np.asarray(expm_multiply(G4(k, M, rho).transpose() * root_length, survivor))
    q = np.array([h @ root for h in residual_vectors()])
    q[0] += 1.0 - left[:4].sum() * right[:4].sum()
    return q


def fixed_p_half_species_tree_conditional(
    z: str,
    *,
    M: int,
    rho: float,
    root_length: float,
    sister_lengths: tuple[float, float],
    terminal_lengths: tuple[float, float, float, float],
) -> np.ndarray:
    """Fixed-count conditional law for a specified terminal pattern.

    This is the regression interface for the old ``H_m``/``A_m`` formulas.
    The unconditional fixed-``p=1/2`` quartet law averages this function over
    all 16 patterns and is generally uniform by symmetry.
    """
    if M % 2 or len(z) != 4 or set(z) - set("01"):
        raise ValueError("M must be even and z must be a four-bit pattern")
    k = M // 2

    def pair(bits: str, lengths: tuple[float, float], sister: float) -> np.ndarray:
        u1 = np.asarray(expm_multiply(G1(k, M, rho).transpose() * lengths[0], np.eye(2)[int(bits[0])]))
        u2 = np.asarray(expm_multiply(G1(k, M, rho).transpose() * lengths[1], np.eye(2)[int(bits[1])]))
        young = np.zeros(6)
        for b1, b2 in PAIR_BACKGROUND_STATES:
            young[2 * b1 + b2] = u1[b1] * u2[b2]
        return np.asarray(expm_multiply(G2(k, M, rho).transpose() * sister, young))

    left = pair(z[:2], terminal_lengths[:2], sister_lengths[0])
    right = pair(z[2:], terminal_lengths[2:], sister_lengths[1])
    survivor = np.zeros(19)
    for il, zl in enumerate(PAIR_BACKGROUND_STATES):
        for ir, zr in enumerate(PAIR_BACKGROUND_STATES):
            survivor[BACKGROUND_STATES.index(zl + zr)] = left[il] * right[ir]
    root = np.asarray(expm_multiply(G4(k, M, rho).transpose() * root_length, survivor))
    q = np.array([h @ root for h in residual_vectors()])
    q[0] += 1.0 - left[:4].sum() * right[:4].sum()
    return q
