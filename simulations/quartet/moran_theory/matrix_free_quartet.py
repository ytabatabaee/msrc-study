"""Matrix-free exact propagation for finite-M Moran quartet kernels.

This module does not reduce the mathematical state space.  It applies the
same count-by-genealogy operator as ``moran_quartet_operator`` without
materializing its sparse block matrix.  It is useful when the exact
count-indexed state is too large for convenient sparse assembly.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse.linalg import LinearOperator, expm_multiply

try:
    from .moran_quartet_operator import (
        BACKGROUND_STATES,
        G4,
        genealogy_generator,
        moran_generator,
        residual_vectors,
        _validate_pi,
    )
except ImportError:  # pragma: no cover
    from moran_quartet_operator import (
        BACKGROUND_STATES,
        G4,
        genealogy_generator,
        moran_generator,
        residual_vectors,
        _validate_pi,
    )


@dataclass(frozen=True)
class MatrixFreeQuartetResult:
    """Result matching ``MoranQuartetResult`` for the shared-population case."""

    joint: np.ndarray
    quartet_probabilities: np.ndarray
    pattern_probabilities: np.ndarray
    conditional_quartet_probabilities: np.ndarray


def matrix_free_operator(
    M: int,
    rho: float,
    n: int,
    selection: float = 0.0,
) -> LinearOperator:
    """Return the exact finite-M operator as a ``LinearOperator``.

    The operator is exactly

    ``Q_M kron I_d + blockdiag(G_n(k).T)``

    but stores only the Moran generator and the ``M+1`` small genealogy
    blocks.  Its dimension remains ``(M+1)*d``.
    """
    if n not in {1, 2, 3, 4}:
        raise ValueError("n must be 1, 2, 3, or 4")
    if M <= 0 or rho < 0:
        raise ValueError("M must be positive and rho must be nonnegative")
    d = {1: 2, 2: 6, 3: 11, 4: 19}[n]
    q = moran_generator(M, selection)
    blocks = [genealogy_generator(n, k, M, rho).transpose().tocsr() for k in range(M + 1)]
    size = (M + 1) * d

    def apply(values: np.ndarray) -> np.ndarray:
        values = np.asarray(values, dtype=float)
        if values.ndim == 1:
            rows = values.reshape(M + 1, d)
            out = q @ rows
            for k, block in enumerate(blocks):
                out[k] += block @ rows[k]
            return np.asarray(out).reshape(-1)
        if values.ndim == 2:
            columns = values.shape[1]
            rows = values.reshape(M + 1, d, columns)
            out = np.empty_like(rows)
            for column in range(columns):
                out[:, :, column] = q @ rows[:, :, column]
            for k, block in enumerate(blocks):
                out[k] += block @ rows[k]
            return out.reshape(size, columns)
        raise ValueError("operator input must be a vector or a two-dimensional matrix")

    def apply_transpose(values: np.ndarray) -> np.ndarray:
        values = np.asarray(values, dtype=float)
        if values.ndim == 1:
            rows = values.reshape(M + 1, d)
            out = q.transpose() @ rows
            for k, block in enumerate(blocks):
                out[k] += block.transpose() @ rows[k]
            return np.asarray(out).reshape(-1)
        if values.ndim == 2:
            columns = values.shape[1]
            rows = values.reshape(M + 1, d, columns)
            out = np.empty_like(rows)
            for column in range(columns):
                out[:, :, column] = q.transpose() @ rows[:, :, column]
            for k, block in enumerate(blocks):
                out[k] += block.transpose() @ rows[k]
            return out.reshape(size, columns)
        raise ValueError("operator input must be a vector or a two-dimensional matrix")

    trace = float(q.diagonal().sum() * d)
    trace += float(sum(block.diagonal().sum() for block in blocks))

    operator = LinearOperator(
        shape=(size, size), dtype=float,
        matvec=apply, matmat=apply,
        rmatvec=apply_transpose, rmatmat=apply_transpose,
    )
    operator._trace = trace
    return operator


def matrix_free_dual_lift(
    M: int,
    rho: float,
    T: float,
    values: np.ndarray,
    n: int,
    *,
    selection: float = 0.0,
) -> np.ndarray:
    """Apply the exact dual exponential without assembling a block matrix."""
    d = {1: 2, 2: 6, 3: 11, 4: 19}.get(n)
    if d is None:
        raise ValueError("n must be 1, 2, 3, or 4")
    array = np.asarray(values, dtype=float)
    if array.shape[0] != (M + 1) * d:
        raise ValueError("values has the wrong joint-state dimension")
    operator = matrix_free_operator(M, rho, n, selection)
    # Supplying the trace avoids SciPy's stochastic trace estimate and makes
    # repeated exact evaluations deterministic.
    trace = float(operator._trace) * float(T)
    return np.asarray(expm_multiply(operator * float(T), array, traceA=trace))


def matrix_free_shared_population_quartet_probabilities(
    M: int,
    rho: float,
    T: float,
    pi_old: np.ndarray,
    *,
    selection: float = 0.0,
) -> MatrixFreeQuartetResult:
    """Compute the exact shared-population unconditional quartet law.

    This is the matrix-free counterpart of
    ``integrated_quartet_probabilities``.  It retains all ``M+1`` count states
    and is therefore an implementation optimization, not an M-independent
    state reduction.
    """
    if T < 0:
        raise ValueError("T must be nonnegative")
    pi = _validate_pi(pi_old, M)
    terminal = np.zeros(((M + 1) * 19, 16), dtype=float)
    for k in range(M + 1):
        p = k / M
        for index, z in enumerate(BACKGROUND_STATES):
            terminal[k * 19 + index, index] = p ** sum(z) * (1.0 - p) ** (4 - sum(z))
    propagated = matrix_free_dual_lift(M, rho, T, terminal, 4, selection=selection)
    initial = np.column_stack([np.kron(pi, h) for h in residual_vectors()])
    joint = (initial.T @ propagated).T
    pattern = joint.sum(axis=1)
    quartet = joint.sum(axis=0)
    conditional = np.divide(
        joint, pattern[:, None], out=np.full_like(joint, np.nan),
        where=pattern[:, None] > 0,
    )
    return MatrixFreeQuartetResult(joint, quartet, pattern, conditional)


__all__ = [
    "MatrixFreeQuartetResult",
    "matrix_free_operator",
    "matrix_free_dual_lift",
    "matrix_free_shared_population_quartet_probabilities",
]
