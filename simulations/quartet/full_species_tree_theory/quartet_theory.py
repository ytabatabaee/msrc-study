"""Symmetric two-background four-taxon MSRC species-tree theorem.

For S=((1,2),(3,4)), the calculation separates the two sister-population
branches from the final shared structured population.  The latter is exactly
the H_m kernel used by msrc-sim's analytic.py.
"""

from __future__ import annotations

from itertools import product

import numpy as np
from scipy.linalg import expm

TOPOLOGY_NAMES = ("12|34", "13|24", "14|23")
PATTERNS = tuple("".join(x) for x in product("01", repeat=4))
PAIR_STATES = ("00", "01", "10", "11")
PAIR_TO_TOPOLOGY = {
    frozenset((0, 1)): 0, frozenset((2, 3)): 0,
    frozenset((0, 2)): 1, frozenset((1, 3)): 1,
    frozenset((0, 3)): 2, frozenset((1, 2)): 2,
}


def U(t: float, m: float) -> np.ndarray:
    """One-lineage background transition matrix exp([[−m,m],[m,−m]]t)."""
    if t < 0 or m < 0:
        raise ValueError("t and m must be non-negative")
    e = np.exp(-2.0 * m * t)
    return 0.5 * np.array([[1.0 + e, 1.0 - e], [1.0 - e, 1.0 + e]])


def pair_generator(m: float, lam: float) -> np.ndarray:
    """The four-state transient sister-pair generator K."""
    if m < 0 or lam < 0:
        raise ValueError("m and lambda must be non-negative")
    return np.array([
        [-(lam + 2 * m), m, m, 0],
        [m, -2 * m, 0, m],
        [m, 0, -2 * m, m],
        [0, m, m, -(lam + 2 * m)],
    ], dtype=float)


def pair_kernel_closed_form(t: float, m: float, lam: float) -> np.ndarray:
    """Closed form for exp(Kt), using the symmetric/antisymmetric modes."""
    if t < 0:
        raise ValueError("t must be non-negative")
    delta = np.sqrt(lam * lam + 16.0 * m * m)
    if delta == 0.0:
        return np.eye(4)
    A = np.array([[-(lam + 2 * m), 2 * m], [2 * m, -2 * m]], dtype=float)
    mu_plus = (-(lam + 4 * m) + delta) / 2.0
    mu_minus = (-(lam + 4 * m) - delta) / 2.0
    identity = np.eye(2)
    # Spectral interpolation is algebraically the same as the cosh/sinh
    # expression but remains stable when m is large.
    E = (
        np.exp(mu_plus * t) * (A - mu_minus * identity)
        - np.exp(mu_minus * t) * (A - mu_plus * identity)
    ) / delta
    D = np.exp(-(lam + 2 * m) * t)
    mixed_difference = np.exp(-2.0 * m * t)
    # Row-vector ordering is 00, 01, 10, 11. Recover the four coordinates
    # from S=x00+x11, M=x01+x10, D=x00-x11, E=x01-x10.
    T = np.zeros((4, 4), dtype=float)
    for j in range(4):
        initial_S = 1.0 if j in (0, 3) else 0.0
        initial_M = 1.0 if j in (1, 2) else 0.0
        initial_D = 1.0 if j == 0 else (-1.0 if j == 3 else 0.0)
        initial_E = 1.0 if j == 1 else (-1.0 if j == 2 else 0.0)
        S, M = np.array([initial_S, initial_M]) @ E
        D_j = initial_D * D
        E_j = initial_E * mixed_difference
        T[j, :] = [(S + D_j) / 2, (M + E_j) / 2, (M - E_j) / 2, (S - D_j) / 2]
    return T


def pair_survivor_distribution(
    z_pair: tuple[int, int] | str,
    *,
    terminal_lengths: tuple[float, float] = (0.0, 0.0),
    shared_length: float,
    m: float,
    lam: float,
) -> np.ndarray:
    """Unnormalized four-state distribution after a sister branch."""
    bits = tuple(int(x) for x in z_pair)
    if len(bits) != 2 or set(bits) - {0, 1}:
        raise ValueError("z_pair must contain two binary states")
    if shared_length < 0 or min(terminal_lengths) < 0:
        raise ValueError("branch lengths must be non-negative")
    left = U(terminal_lengths[0], m)[bits[0]]
    right = U(terminal_lengths[1], m)[bits[1]]
    entering = np.kron(left, right)
    return entering @ pair_kernel_closed_form(shared_length, m, lam)


def _hm_generator(m: float, lam: float):
    Q = np.zeros((16, 16), dtype=float)
    R = np.zeros((16, 3), dtype=float)
    states = list(product((0, 1), repeat=4))
    index = {z: i for i, z in enumerate(states)}
    for row, z in enumerate(states):
        exit_rate = 0.0
        for i, bit in enumerate(z):
            zz = list(z)
            zz[i] = 1 - bit
            Q[row, index[tuple(zz)]] += m
            exit_rate += m
        for a in range(4):
            for b in range(a + 1, 4):
                if z[a] == z[b]:
                    R[row, PAIR_TO_TOPOLOGY[frozenset((a, b))]] += lam
                    exit_rate += lam
        Q[row, row] -= exit_rate
    return Q, R, states


def hm_kernel(shared_length: float, m: float, lam: float) -> np.ndarray:
    """Return H_m, including 1/3 resolution of survival mass after T."""
    if shared_length < 0:
        raise ValueError("shared_length must be non-negative")
    Q, R, _ = _hm_generator(m, lam)
    B = np.zeros((19, 19), dtype=float)
    B[:16, :16] = Q
    B[:16, 16:] = R
    E = expm(B * shared_length)
    return E[:16, 16:] + E[:16, :16] @ np.ones((16, 1)) @ np.full((1, 3), 1.0 / 3.0)


def conditional_quartet_probabilities(
    z: str,
    *,
    terminal_lengths: tuple[float, float, float, float] = (0.0, 0.0, 0.0, 0.0),
    sister_lengths: tuple[float, float] = (0.0, 0.0),
    shared_length: float = 0.0,
    m: float,
    lam: float,
) -> np.ndarray:
    """Return q_Z=(1-s_L s_R)e1 + (ell_Z kron r_Z)H_m."""
    if len(z) != 4 or set(z) - set("01"):
        raise ValueError("z must be a four-character 0/1 pattern")
    if len(terminal_lengths) != 4 or len(sister_lengths) != 2:
        raise ValueError("terminal_lengths must have four and sister_lengths two entries")
    ell = pair_survivor_distribution(
        z[:2], terminal_lengths=terminal_lengths[:2],
        shared_length=sister_lengths[0], m=m, lam=lam,
    )
    right = pair_survivor_distribution(
        z[2:], terminal_lengths=terminal_lengths[2:],
        shared_length=sister_lengths[1], m=m, lam=lam,
    )
    s_left, s_right = ell.sum(), right.sum()
    beta = np.kron(ell, right)
    q = (1.0 - s_left * s_right) * np.array([1.0, 0.0, 0.0])
    q += beta @ hm_kernel(shared_length, m, lam)
    return q


def all_pattern_probabilities(**kwargs) -> dict[str, np.ndarray]:
    return {z: conditional_quartet_probabilities(z, **kwargs) for z in PATTERNS}


def pattern_classes() -> dict[str, tuple[str, str]]:
    classes = {}
    for z in PATTERNS:
        ones = z.count("1")
        cls = "4:0" if ones in (0, 4) else ("3:1" if ones in (1, 3) else "2:2")
        split = "Q1-associated" if z in {"0000", "1111"} else ""
        if cls == "2:2":
            split = {"0011": "W1", "1100": "W1", "0101": "W2", "1010": "W2", "0110": "W3", "1001": "W3"}[z]
        classes[z] = (cls, split)
    return classes


def unconditional_quartet_probabilities(pattern_weights: dict[str, float], **kwargs) -> np.ndarray:
    missing = set(PATTERNS) - set(pattern_weights)
    if missing:
        raise ValueError(f"missing pattern weights: {sorted(missing)}")
    return sum(pattern_weights[z] * conditional_quartet_probabilities(z, **kwargs) for z in PATTERNS)


def msc_quartet_probabilities(tau: float) -> np.ndarray:
    e = np.exp(-tau)
    return np.array([1.0 - 2.0 * e / 3.0, e / 3.0, e / 3.0])
