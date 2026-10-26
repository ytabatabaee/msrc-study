"""Exact balanced and pectinate species-tree recursions for a two-regime model.

The regime process is a separate stochastic model from the finite-M Moran
process.  A regime has fixed background frequency p_L or p_H and switches with
a two-state CTMC on each branch.  Regime inheritance is represented by an
explicit parent-to-child stochastic matrix at each speciation.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
from scipy.linalg import expm
from scipy.sparse import block_diag, csr_matrix, eye, kron
from scipy.sparse.linalg import expm_multiply

from .moran_quartet_operator import (
    BACKGROUND_STATES, PAIR_BACKGROUND_STATES, PATTERNS, residual_vectors,
)
from .prescribed_frequency_quartet import prescribed_generator


@dataclass(frozen=True)
class RegimeBranch:
    duration: float = 0.0
    rate_low_to_high: float = 0.0
    rate_high_to_low: float = 0.0
    M: int = 10


@dataclass(frozen=True)
class TwoRegimeModel:
    p_low: float
    p_high: float
    rho: float

    def __post_init__(self):
        if not (0 <= self.p_low <= 1 and 0 <= self.p_high <= 1 and self.rho >= 0):
            raise ValueError("frequencies must be in [0,1] and rho nonnegative")


@dataclass(frozen=True)
class TwoRegimeSpeciesTreeResult:
    quartet_probabilities: np.ndarray
    early_topology_contributions: np.ndarray
    root_topology_contributions: np.ndarray
    root_entry: np.ndarray
    joint: np.ndarray | None
    pattern_probabilities: np.ndarray | None
    conditional_quartet_probabilities: np.ndarray | None
    diagnostics: dict[str, float]


def identity_inheritance() -> np.ndarray:
    return np.eye(2, dtype=float)


def _check_kernel(H: np.ndarray | None) -> np.ndarray:
    out = identity_inheritance() if H is None else np.asarray(H, dtype=float)
    if out.shape != (2, 2) or np.any(out < -1e-14) or not np.allclose(out.sum(axis=1), 1.0):
        raise ValueError("inheritance kernels must be 2x2 row-stochastic matrices")
    return out


def regime_generator(branch: RegimeBranch) -> np.ndarray:
    if branch.duration < 0 or branch.M <= 0 or branch.rate_low_to_high < 0 or branch.rate_high_to_low < 0:
        raise ValueError("invalid regime branch")
    return np.array([[-branch.rate_low_to_high, branch.rate_low_to_high],
                     [branch.rate_high_to_low, -branch.rate_high_to_low]], dtype=float)


def regime_operator(model: TwoRegimeModel, branch: RegimeBranch, n: int) -> csr_matrix:
    """Return R kron I + blockdiag(G_n(p_L)^T,G_n(p_H)^T)."""
    if n not in (1, 2, 3, 4):
        raise ValueError("n must be 1, 2, 3, or 4")
    d = {1: 2, 2: 6, 3: 11, 4: 19}[n]
    R = csr_matrix(regime_generator(branch))
    G = [csr_matrix(prescribed_generator(n, p, branch.M, model.rho).T)
         for p in (model.p_low, model.p_high)]
    return (kron(R, eye(d, format="csr"), format="csr") +
            block_diag(G, format="csr"))


def _lift(model: TwoRegimeModel, branch: RegimeBranch, n: int, value: np.ndarray) -> np.ndarray:
    value = np.asarray(value, dtype=float)
    d = {1: 2, 2: 6, 3: 11, 4: 19}[n]
    if value.shape != (2, d):
        raise ValueError(f"value must have shape (2,{d})")
    if branch.duration == 0:
        return value.copy()
    return np.asarray(expm_multiply(regime_operator(model, branch, n) * branch.duration, value.reshape(-1))).reshape(2, d)


def _terminal_value(model: TwoRegimeModel, bit: int | None = None) -> np.ndarray:
    out = np.zeros((2, 2), dtype=float)
    for e, p in enumerate((model.p_low, model.p_high)):
        if bit is None:
            out[e] = (1.0 - p, p)
        else:
            out[e, bit] = p if bit else 1.0 - p
    return out


def _terminal_at_parent(model, terminal_branch, H, bit=None):
    return _check_kernel(H) @ _lift(model, terminal_branch, 1, _terminal_value(model, bit))


def _pair_at_parent(model, terminal_branches, pair_branch, H1, H2, bits=None):
    u1 = _terminal_at_parent(model, terminal_branches[0], H1, None if bits is None else bits[0])
    u2 = _terminal_at_parent(model, terminal_branches[1], H2, None if bits is None else bits[1])
    young = np.zeros((2, 6), dtype=float)
    for index, (b1, b2) in enumerate(PAIR_BACKGROUND_STATES):
        young[:, index] = u1[:, b1] * u2[:, b2]
    return _lift(model, pair_branch, 2, young)


def _root_result(model, root_branch, F, initial_regime):
    pi = np.asarray(initial_regime, dtype=float)
    if pi.shape != (2,) or np.any(pi < 0) or not np.isclose(pi.sum(), 1):
        raise ValueError("initial_regime must be a length-2 probability vector")
    propagated = _lift(model, root_branch, 4, F)
    initial = np.column_stack([np.kron(pi, h) for h in residual_vectors()])
    return np.asarray(initial.T @ propagated.reshape(-1))


def _young_regime_distribution(model, branch, pi):
    return np.asarray(expm(regime_generator(branch).T * branch.duration) @ np.asarray(pi, dtype=float))


def balanced_two_regime_quartet_probabilities(
    *, model: TwoRegimeModel, root_branch: RegimeBranch,
    left_pair_branch: RegimeBranch, right_pair_branch: RegimeBranch,
    terminal_branches: tuple[RegimeBranch, RegimeBranch, RegimeBranch, RegimeBranch],
    initial_regime: np.ndarray, H_root_left: np.ndarray | None = None,
    H_root_right: np.ndarray | None = None, H_left_1: np.ndarray | None = None,
    H_left_2: np.ndarray | None = None, H_right_3: np.ndarray | None = None,
    H_right_4: np.ndarray | None = None, with_joint: bool = False,
) -> TwoRegimeSpeciesTreeResult:
    """Exact unconditional balanced-tree result for the regime model."""
    Hrl, Hrr = _check_kernel(H_root_left), _check_kernel(H_root_right)
    Hl1, Hl2 = _check_kernel(H_left_1), _check_kernel(H_left_2)
    Hr3, Hr4 = _check_kernel(H_right_3), _check_kernel(H_right_4)
    left_old = _pair_at_parent(model, terminal_branches[:2], left_pair_branch, Hl1, Hl2)
    right_old = _pair_at_parent(model, terminal_branches[2:], right_pair_branch, Hr3, Hr4)
    left = Hrl @ left_old
    right = Hrr @ right_old
    F = np.zeros((2, 19), dtype=float)
    early_by_regime = np.zeros((2, 3), dtype=float)
    for e in range(2):
        sl, sr = left[e, :4].sum(), right[e, :4].sum()
        early_by_regime[e, 0] = 1.0 - sl * sr
        for il, zl in enumerate(PAIR_BACKGROUND_STATES[:4]):
            for ir, zr in enumerate(PAIR_BACKGROUND_STATES[:4]):
                F[e, BACKGROUND_STATES.index(zl + zr)] = left[e, il] * right[e, ir]
    root = _root_result(model, root_branch, F, initial_regime)
    pi_young = _young_regime_distribution(model, root_branch, initial_regime)
    early = pi_young @ early_by_regime
    q = early + root
    joint = patterns = conditional = None
    if with_joint:
        joint = np.zeros((16, 3), dtype=float)
        for zi, pattern in enumerate(PATTERNS):
            lz = Hrl @ _pair_at_parent(model, terminal_branches[:2], left_pair_branch, Hl1, Hl2, tuple(map(int, pattern[:2])))
            rz = Hrr @ _pair_at_parent(model, terminal_branches[2:], right_pair_branch, Hr3, Hr4, tuple(map(int, pattern[2:])))
            Fz = np.zeros((2, 19), dtype=float); ez = np.zeros((2, 3), dtype=float)
            for e in range(2):
                ez[e, 0] = lz[e].sum() * rz[e].sum() - lz[e, :4].sum() * rz[e, :4].sum()
                for il, zl in enumerate(PAIR_BACKGROUND_STATES[:4]):
                    for ir, zr in enumerate(PAIR_BACKGROUND_STATES[:4]):
                        Fz[e, BACKGROUND_STATES.index(zl + zr)] = lz[e, il] * rz[e, ir]
            joint[zi] = pi_young @ ez + _root_result(model, root_branch, Fz, initial_regime)
        patterns = joint.sum(axis=1)
        conditional = np.divide(joint, patterns[:, None], out=np.full_like(joint, np.nan), where=patterns[:, None] > 0)
    return TwoRegimeSpeciesTreeResult(q, early, root, F, joint, patterns, conditional,
                                      {"probability_sum_error": float(q.sum() - 1), "root_entry_mass": float(F.sum())})


def unbalanced_two_regime_quartet_probabilities(
    *, model: TwoRegimeModel, root_branch: RegimeBranch,
    intermediate_branch: RegimeBranch, sister_branch: RegimeBranch,
    terminal_branches: tuple[RegimeBranch, RegimeBranch, RegimeBranch, RegimeBranch],
    initial_regime: np.ndarray, H_intermediate_pair: np.ndarray | None = None,
    H_intermediate_3: np.ndarray | None = None, H_root_intermediate: np.ndarray | None = None,
    H_root_4: np.ndarray | None = None, H_pair_1: np.ndarray | None = None,
    H_pair_2: np.ndarray | None = None, with_joint: bool = False,
) -> TwoRegimeSpeciesTreeResult:
    """Exact unconditional pectinate-tree result for the regime model."""
    Hip, Hi3 = _check_kernel(H_intermediate_pair), _check_kernel(H_intermediate_3)
    Hri, Hr4 = _check_kernel(H_root_intermediate), _check_kernel(H_root_4)
    Hp1, Hp2 = _check_kernel(H_pair_1), _check_kernel(H_pair_2)

    def evaluate(bits=None):
        pair = _pair_at_parent(model, terminal_branches[:2], sister_branch, Hp1, Hp2, None if bits is None else bits[:2])
        u3 = _terminal_at_parent(model, terminal_branches[2], Hi3, None if bits is None else bits[2])
        pair = Hip @ pair
        young = np.zeros((2, 11), dtype=float)
        # A sister-pair topology may be resolved before taxon 3 enters.  For
        # joint terminal-pattern probabilities, retain taxon-3 sampling mass
        # even though it is absent from the unresolved genealogy state.
        young[:, 8] = (pair[:, 4] + pair[:, 5]) * u3.sum(axis=1)
        for index, (b1, b2) in enumerate(PAIR_BACKGROUND_STATES[:4]):
            for b3 in (0, 1):
                young[:, index * 2 + b3] = pair[:, index] * u3[:, b3]
        triple_old = _lift(model, intermediate_branch, 3, young)
        triple = Hri @ triple_old
        u4 = _terminal_at_parent(model, terminal_branches[3], Hr4, None if bits is None else bits[3])
        F = np.zeros((2, 19), dtype=float)
        for e in range(2):
            for state in product((0, 1), repeat=3):
                index = sum(bit << (2 - j) for j, bit in enumerate(state))
                for b4 in (0, 1):
                    F[e, BACKGROUND_STATES.index(state + (b4,))] = triple[e, index] * u4[e, b4]
        # If a terminal pattern is requested, its taxon-4 sampling mass still
        # belongs to the joint event even when the topology was resolved
        # before lineage 4 entered.  The unconditional u4 sum is one, so this
        # factor is invisible in the marginal recursion.
        early = triple[:, 8:].copy() * u4.sum(axis=1)[:, None]
        return early, F

    early_by_regime, F = evaluate()
    root = _root_result(model, root_branch, F, initial_regime)
    pi_young = _young_regime_distribution(model, root_branch, initial_regime)
    early = pi_young @ early_by_regime
    q = early + root
    joint = patterns = conditional = None
    if with_joint:
        joint = np.zeros((16, 3), dtype=float)
        for zi, pattern in enumerate(PATTERNS):
            ez, Fz = evaluate(tuple(map(int, pattern)))
            joint[zi] = pi_young @ ez + _root_result(model, root_branch, Fz, initial_regime)
        patterns = joint.sum(axis=1)
        conditional = np.divide(joint, patterns[:, None], out=np.full_like(joint, np.nan), where=patterns[:, None] > 0)
    return TwoRegimeSpeciesTreeResult(q, early, root, F, joint, patterns, conditional,
                                      {"probability_sum_error": float(q.sum() - 1), "root_entry_mass": float(F.sum())})


__all__ = [
    "RegimeBranch", "TwoRegimeModel", "TwoRegimeSpeciesTreeResult",
    "identity_inheritance", "regime_generator", "regime_operator",
    "balanced_two_regime_quartet_probabilities",
    "unbalanced_two_regime_quartet_probabilities",
]
