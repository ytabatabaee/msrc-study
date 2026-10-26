"""Exact quartet kernel for a separate finite-regime frequency model.

This is not a reduction of the Moran model.  It defines a different model in
which the frequency regime is a two-state CTMC with fixed frequencies p_L and
p_H.  The resulting shared-population operator has 2*19 states.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import block_diag, csr_matrix, eye, kron
from scipy.sparse.linalg import expm_multiply

from .moran_quartet_operator import BACKGROUND_STATES, residual_vectors
from .prescribed_frequency_quartet import prescribed_generator


@dataclass(frozen=True)
class RegimeSwitchingQuartetResult:
    joint: np.ndarray
    quartet_probabilities: np.ndarray
    pattern_probabilities: np.ndarray
    conditional_quartet_probabilities: np.ndarray


def regime_switching_shared_population_quartet_probabilities(
    *, p_low: float, p_high: float, rate_low_to_high: float,
    rate_high_to_low: float, M: int, rho: float, T: float,
    initial_regime: np.ndarray, initial_count_distribution: np.ndarray | None = None,
) -> RegimeSwitchingQuartetResult:
    """Evaluate the exact two-regime shared-population quartet law.

    The regime CTMC is independent of the genealogy conditional on its state.
    The model is initialized at an older boundary; terminal sampling uses the
    frequency of the regime at the younger boundary through the joint
    operator.  ``initial_count_distribution`` is intentionally unsupported in
    this alternative model and is retained only as a guard against confusing
    it with the Moran solver.
    """
    if initial_count_distribution is not None:
        raise ValueError("the regime model has no Moran count coordinate")
    if not (0 <= p_low <= 1 and 0 <= p_high <= 1 and M > 0 and rho >= 0 and T >= 0):
        raise ValueError("invalid regime-model parameters")
    pi = np.asarray(initial_regime, dtype=float)
    if pi.shape != (2,) or np.any(pi < 0) or not np.isclose(pi.sum(), 1):
        raise ValueError("initial_regime must be a probability vector of length 2")
    if rate_low_to_high < 0 or rate_high_to_low < 0:
        raise ValueError("regime rates must be nonnegative")
    R = np.array([[-rate_low_to_high, rate_low_to_high],
                  [rate_high_to_low, -rate_high_to_low]], dtype=float)
    G = [prescribed_generator(4, p, M, rho) for p in (p_low, p_high)]
    # Match the existing Moran dual convention: the forward frequency-process
    # generator is used as R, while the backward genealogy blocks are G^T.
    A = kron(csr_matrix(R), eye(19, format="csr"), format="csr")
    A += block_diag([csr_matrix(g.T) for g in G], format="csr")
    terminal = np.zeros((38, 16), dtype=float)
    for regime, p in enumerate((p_low, p_high)):
        for state, z in enumerate(BACKGROUND_STATES):
            terminal[regime * 19 + state, state] = p ** sum(z) * (1.0 - p) ** (4 - sum(z))
    initial = np.column_stack([np.kron(pi, h) for h in residual_vectors()])
    propagated = expm_multiply(A * T, terminal) if T else terminal
    joint = np.asarray((initial.T @ propagated).T)
    patterns = joint.sum(axis=1)
    q = joint.sum(axis=0)
    conditional = np.divide(joint, patterns[:, None], out=np.full_like(joint, np.nan), where=patterns[:, None] > 0)
    return RegimeSwitchingQuartetResult(joint, q, patterns, conditional)


__all__ = ["RegimeSwitchingQuartetResult", "regime_switching_shared_population_quartet_probabilities"]
