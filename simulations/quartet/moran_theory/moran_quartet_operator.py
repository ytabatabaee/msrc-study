"""Exact finite-M Moran-integrated shared-population quartet kernel.

The Moran count runs forward from the older boundary to the younger boundary.
The structured genealogy is conditional on that forward path and runs backward
from the younger boundary.  The transpose blocks in ``joint_operator`` are a
dual/Feynman--Kac construction for averaging those two processes; they are not
a literal forward biological genealogy CTMC.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
from scipy.sparse import block_diag, csr_matrix, eye, kron, lil_matrix
from scipy.sparse.linalg import expm_multiply

try:  # scipy is used for the sparse operator; this import is just for checks.
    from scipy.linalg import expm
except ImportError:  # pragma: no cover
    expm = None


TOPOLOGY_NAMES = ("12|34", "13|24", "14|23")
PATTERNS = tuple("".join(str(x) for x in z) for z in product((0, 1), repeat=4))
BACKGROUND_STATES = tuple(product((0, 1), repeat=4))
BACKGROUND_INDEX = {z: i for i, z in enumerate(BACKGROUND_STATES)}
TOPOLOGY_PAIRS = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))


@dataclass(frozen=True)
class MoranRates:
    q_plus: float
    q_minus: float


@dataclass(frozen=True)
class MoranQuartetResult:
    joint: np.ndarray
    quartet_probabilities: np.ndarray
    pattern_probabilities: np.ndarray
    conditional_quartet_probabilities: np.ndarray


def moran_rates(k: int, M: int, selection: float = 0.0) -> MoranRates:
    """Match ``msrcsim.moran.moran_rates`` on the generation-equivalent scale."""
    if M <= 0:
        raise ValueError("M must be positive")
    if not 0 <= k <= M:
        raise ValueError("k must be between 0 and M")
    if 1.0 + float(selection) <= 0.0:
        raise ValueError("Moran selection requires 1 + s_M > 0")
    if k in (0, M):
        return MoranRates(0.0, 0.0)
    p = k / M
    p_sel = p if selection == 0.0 else p * (1.0 + selection) / (1.0 + selection * p)
    scale = M / 2.0
    return MoranRates(scale * p_sel * (1.0 - p), scale * (1.0 - p_sel) * p)


def moran_generator(M: int, selection: float = 0.0) -> csr_matrix:
    """Return the row generator on count states ``K=0,...,M``."""
    Q = lil_matrix((M + 1, M + 1), dtype=float)
    for k in range(1, M):
        rates = moran_rates(k, M, selection)
        Q[k, k + 1] = rates.q_plus
        Q[k, k - 1] = rates.q_minus
        Q[k, k] = -rates.q_plus - rates.q_minus
    return Q.tocsr()


def G4(k: int, M: int, rho: float) -> csr_matrix:
    """Return the 19-state conditional genealogy generator for count ``k``."""
    if M <= 0 or not 0 <= k <= M or rho < 0:
        raise ValueError("require M>0, 0<=k<=M, and rho>=0")
    p = k / M
    c1 = 1.0 / max(k, 1)
    c0 = 1.0 / max(M - k, 1)
    S = lil_matrix((16, 16), dtype=float)
    R = np.zeros((16, 3), dtype=float)
    for row, z in enumerate(BACKGROUND_STATES):
        exit_rate = 0.0
        for i, bit in enumerate(z):
            rate = rho * p if bit == 0 else rho * (1.0 - p)
            if rate:
                zz = list(z)
                zz[i] = 1 - bit
                S[row, BACKGROUND_INDEX[tuple(zz)]] += rate
                exit_rate += rate
        for topology, pairs in enumerate(TOPOLOGY_PAIRS):
            absorption = 0.0
            for i, j in pairs:
                if z[i] == z[j]:
                    absorption += c1 if z[i] == 1 else c0
            R[row, topology] = absorption
            exit_rate += absorption
        S[row, row] = -exit_rate
    out = lil_matrix((19, 19), dtype=float)
    out[:16, :16] = S
    out[:16, 16:] = csr_matrix(R)
    return out.tocsr()


def joint_operator(
    M: int, rho: float, selection: float = 0.0, *, moran_enabled: bool = True
) -> csr_matrix:
    """Return ``Q_M kron I_19 + blockdiag(G4(k).T)``.

    The transpose blocks implement the dual representation needed because the
    Moran path is forward in biological time while the genealogy is backward.
    """
    Q = moran_generator(M, selection) if moran_enabled else csr_matrix((M + 1, M + 1))
    blocks = [G4(k, M, rho).transpose().tocsr() for k in range(M + 1)]
    return (kron(Q, eye(19, format="csr"), format="csr") + block_diag(blocks, format="csr")).tocsr()


def beta(k: int, M: int, z: tuple[int, int, int, int]) -> float:
    """Terminal independent sampling weight conditional on ``K_T=k``."""
    p = k / M
    n1 = sum(z)
    return p**n1 * (1.0 - p) ** (4 - n1)


def terminal_vectors(M: int) -> np.ndarray:
    """Return an ``(M+1)*19 by 16`` matrix whose columns are ``b_z``."""
    out = np.zeros(((M + 1) * 19, 16), dtype=float)
    for k in range(M + 1):
        for a, z in enumerate(BACKGROUND_STATES):
            out[k * 19 + a, a] = beta(k, M, z)
    return out


def residual_vectors() -> np.ndarray:
    """Return the three 19-vectors ``h_r`` used at the older boundary."""
    out = np.full((3, 19), 1.0 / 3.0, dtype=float)
    for r in range(3):
        out[r, 16:] = 0.0
        out[r, 16 + r] = 1.0
    return out


def _validate_pi(pi_old: np.ndarray, M: int) -> np.ndarray:
    pi = np.asarray(pi_old, dtype=float)
    if pi.shape != (M + 1,) or np.any(pi < -1e-14) or not np.isclose(pi.sum(), 1.0):
        raise ValueError("pi_old must be a nonnegative probability vector of length M+1")
    return pi


def integrated_quartet_probabilities(
    M: int,
    rho: float,
    T: float,
    pi_old: np.ndarray,
    *,
    selection: float = 0.0,
    moran_enabled: bool = True,
) -> MoranQuartetResult:
    """Evaluate ``J[z,r] = a_r exp(A4 T) b_z`` by sparse expm multiplication."""
    if T < 0:
        raise ValueError("T must be nonnegative")
    pi = _validate_pi(pi_old, M)
    A = joint_operator(M, rho, selection, moran_enabled=moran_enabled)
    initial = np.column_stack([np.kron(pi, h) for h in residual_vectors()])
    terminal = terminal_vectors(M)
    propagated = expm_multiply(A * float(T), terminal)
    joint = (initial.T @ propagated).T
    pattern = joint.sum(axis=1)
    q = joint.sum(axis=0)
    conditional = np.divide(joint, pattern[:, None], out=np.full_like(joint, np.nan), where=pattern[:, None] > 0)
    return MoranQuartetResult(joint, q, pattern, conditional)


def pattern_probabilities_from_moran(M: int, T: float, pi_old: np.ndarray) -> np.ndarray:
    """Direct terminal pattern law ``pi_old exp(Q_M T) beta_z``."""
    from scipy.sparse.linalg import expm_multiply as sparse_expm_multiply

    pi = _validate_pi(pi_old, M)
    pi_T = sparse_expm_multiply(moran_generator(M).transpose() * float(T), pi)
    return np.array([
        sum(pi_T[k] * beta(k, M, tuple(int(x) for x in z)) for k in range(M + 1))
        for z in PATTERNS
    ])


def path_joint_probabilities(
    ks: tuple[int, ...] | list[int],
    durations: tuple[float, ...] | list[float],
    M: int,
    rho: float,
) -> np.ndarray:
    """Evaluate a manually specified forward path using the transpose dual.

    For ``ks=(k1,k2,...)``, this applies
    ``exp(G4(k1).T*d1) exp(G4(k2).T*d2) ...`` to each older-boundary
    residual vector.  It is the transpose of the simulator convention
    ``U_back=exp(G4(kn)dn)...exp(G4(k1)d1)``.
    """
    if len(ks) != len(durations) or not ks:
        raise ValueError("ks and durations must have the same nonzero length")
    if any(d < 0 for d in durations):
        raise ValueError("durations must be nonnegative")
    W = residual_vectors().T
    # Matrix multiplication gives U_back.T = exp(G1.T) exp(G2.T) ...;
    # acting on a column therefore applies the factors in reverse order.
    for k, duration in reversed(tuple(zip(ks, durations))):
        W = expm_multiply(G4(int(k), M, rho).transpose() * float(duration), W)
    final_k = int(ks[-1])
    B = np.array([[beta(final_k, M, z) if a == BACKGROUND_INDEX[z] else 0.0
                   for z in BACKGROUND_STATES] for a in range(19)])
    return B.T @ W


def fixed_frequency_result(M: int, rho: float, T: float) -> MoranQuartetResult:
    """Convenience special case: ``Q_M=0`` and fixed ``K=M/2``."""
    if M % 2:
        raise ValueError("M must be even for fixed p=1/2")
    pi = np.zeros(M + 1)
    pi[M // 2] = 1.0
    return integrated_quartet_probabilities(M, rho, T, pi, moran_enabled=False)
