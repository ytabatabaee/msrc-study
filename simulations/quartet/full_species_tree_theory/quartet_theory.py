"""Exact four-taxon MSRC quartet probabilities on ((1,2),(3,4)).

The implementation follows the event rules in msrc-sim's
``structured_coalescent.py``.  It conditions on the four sampled background
states and uses transient states with four surviving sampled lineages.  A
coalescence is absorbing because its first pair determines the quartet.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
from scipy.linalg import expm

TOPOLOGY_NAMES = ("12|34", "13|24", "14|23")
TAXA = ("1", "2", "3", "4")
PATTERNS = tuple("".join(bits) for bits in product("01", repeat=4))
PAIR_TO_TOPOLOGY = {
    frozenset((0, 1)): 0,
    frozenset((2, 3)): 0,
    frozenset((0, 2)): 1,
    frozenset((1, 3)): 1,
    frozenset((0, 3)): 2,
    frozenset((1, 2)): 2,
}


@dataclass(frozen=True)
class Epoch:
    """A finite epoch with fixed population assignments and rates."""

    duration: float
    populations: tuple[str, str, str, str]
    p: float
    rho: float
    ne: float


def _validate_p(p: float) -> None:
    if not 0.0 < p < 1.0:
        raise ValueError("p must lie strictly between 0 and 1")


def _state_index(bits: tuple[int, int, int, int]) -> int:
    return sum(int(bit) << (3 - i) for i, bit in enumerate(bits))


def _rates(p: float, rho: float, ne: float) -> tuple[float, float, float, float]:
    _validate_p(p)
    if ne <= 0 or rho < 0:
        raise ValueError("ne must be positive and rho must be non-negative")
    return rho * p, rho * (1.0 - p), 1.0 / (2.0 * ne * (1.0 - p)), 1.0 / (2.0 * ne * p)


def transient_generator(populations: tuple[str, str, str, str], p: float, rho: float, ne: float):
    """Return (Q, R) for the 16 background states in one epoch.

    Q is the transient generator and R has one column for each quartet
    topology.  Rows sum to zero after R is appended to Q.
    """
    m01, m10, lambda0, lambda1 = _rates(p, rho, ne)
    Q = np.zeros((16, 16), dtype=float)
    R = np.zeros((16, 3), dtype=float)
    for bits in product((0, 1), repeat=4):
        i = _state_index(bits)
        exit_rate = 0.0
        for k, bit in enumerate(bits):
            rate = m01 if bit == 0 else m10
            if rate:
                switched = list(bits)
                switched[k] = 1 - bit
                Q[i, _state_index(tuple(switched))] += rate
                exit_rate += rate
        for a in range(4):
            for b in range(a + 1, 4):
                if populations[a] == populations[b] and bits[a] == bits[b]:
                    coal_rate = lambda0 if bits[a] == 0 else lambda1
                    topology = PAIR_TO_TOPOLOGY[frozenset((a, b))]
                    R[i, topology] += coal_rate
                    exit_rate += coal_rate
        Q[i, i] -= exit_rate
    return Q, R


def _finite_propagation(alpha: np.ndarray, epoch: Epoch):
    Q, R = transient_generator(epoch.populations, epoch.p, epoch.rho, epoch.ne)
    augmented = np.zeros((19, 19), dtype=float)
    augmented[:16, :16] = Q
    augmented[:16, 16:] = R
    extended = np.zeros(19, dtype=float)
    extended[:16] = alpha
    result = extended @ expm(augmented * epoch.duration)
    return result[:16], result[16:]


def _root_absorption(alpha: np.ndarray, p: float, rho: float, ne: float) -> np.ndarray:
    Q, R = transient_generator(("root",) * 4, p, rho, ne)
    # Integral_0^infinity exp(Qt) R dt = (-Q)^(-1) R.
    return alpha @ np.linalg.solve(-Q, R)


def conditional_quartet_probabilities(
    z: str | tuple[int, int, int, int],
    *,
    sister_length: float,
    root_length_before_mrca: float,
    rho: float,
    p: float,
    ne: float,
    terminal_p: float | None = None,
    sister_p: float | None = None,
    root_p: float | None = None,
    terminal_ne: float | None = None,
    sister_ne: float | None = None,
    root_ne: float | None = None,
) -> np.ndarray:
    """Compute P(Q | Z=z,S) for the balanced four-taxon species tree.

    ``sister_length`` is the length of the two terminal branches before the
    sister populations form. ``root_length_before_mrca`` is the length of the
    two-population ancestral epoch before both sisters enter the root. Rates
    are piecewise constant within these three epochs.
    """
    if isinstance(z, str):
        if len(z) != 4 or set(z) - set("01"):
            raise ValueError("z must be a four-character 0/1 pattern")
        bits = tuple(int(x) for x in z)
    else:
        bits = tuple(int(x) for x in z)
        if len(bits) != 4 or set(bits) - {0, 1}:
            raise ValueError("z must contain four binary states")
    if sister_length < 0 or root_length_before_mrca < 0:
        raise ValueError("epoch lengths must be non-negative")
    terminal_p = p if terminal_p is None else terminal_p
    sister_p = p if sister_p is None else sister_p
    root_p = p if root_p is None else root_p
    terminal_ne = ne if terminal_ne is None else terminal_ne
    sister_ne = ne if sister_ne is None else sister_ne
    root_ne = ne if root_ne is None else root_ne
    alpha = np.zeros(16)
    alpha[_state_index(bits)] = 1.0
    q = np.zeros(3)
    alpha, absorbed = _finite_propagation(
        alpha, Epoch(sister_length, ("1", "2", "3", "4"), terminal_p, rho, terminal_ne)
    )
    q += absorbed
    alpha, absorbed = _finite_propagation(
        alpha, Epoch(root_length_before_mrca, ("12", "12", "34", "34"), sister_p, rho, sister_ne)
    )
    q += absorbed
    q += _root_absorption(alpha, root_p, rho, root_ne)
    return q / q.sum()


def all_pattern_probabilities(**kwargs) -> dict[str, np.ndarray]:
    """Return conditional quartet probabilities for all 16 terminal patterns."""
    return {z: conditional_quartet_probabilities(z, **kwargs) for z in PATTERNS}


def pattern_probability_mixture(pattern_probabilities: dict[str, float], **kwargs) -> np.ndarray:
    """Mix the 16 conditional distributions with P(Z=z)."""
    missing = set(PATTERNS) - set(pattern_probabilities)
    if missing:
        raise ValueError(f"missing pattern probabilities: {sorted(missing)}")
    q = sum(float(pattern_probabilities[z]) * conditional_quartet_probabilities(z, **kwargs) for z in PATTERNS)
    return np.asarray(q, dtype=float)


def msc_quartet_probabilities(internal_branch_coalescent_length: float) -> np.ndarray:
    """Ordinary four-taxon MSC probabilities for ((1,2),(3,4))."""
    if internal_branch_coalescent_length < 0:
        raise ValueError("internal branch length must be non-negative")
    no_coal = np.exp(-internal_branch_coalescent_length)
    return np.array([1.0 - (2.0 / 3.0) * no_coal, no_coal / 3.0, no_coal / 3.0])


def shared_population_kernel(
    z: str | tuple[int, int, int, int],
    *,
    duration: float,
    m01: float,
    m10: float,
    lambda0: float,
    lambda1: float,
) -> np.ndarray:
    """Finite shared-population kernel used by the earlier H_m calculation.

    If no coalescence occurs during the finite interval, the remaining
    panmictic genealogy resolves each quartet with probability 1/3.  This is
    exactly the convention in ``msrcsim.analytic.compute_Hm``.
    """
    if isinstance(z, str):
        bits = tuple(int(x) for x in z)
    else:
        bits = tuple(int(x) for x in z)
    if len(bits) != 4 or set(bits) - {0, 1}:
        raise ValueError("z must contain four binary states")
    if duration < 0 or min(m01, m10, lambda0, lambda1) < 0:
        raise ValueError("duration and rates must be non-negative")
    Q, R = transient_generator(("shared",) * 4, 0.5, 0.0, 1.0)
    # Replace the biological parameterized rates by the requested rates.
    Q.fill(0.0)
    R.fill(0.0)
    for state in product((0, 1), repeat=4):
        i = _state_index(state)
        exit_rate = 0.0
        for k, bit in enumerate(state):
            rate = m01 if bit == 0 else m10
            switched = list(state)
            switched[k] = 1 - bit
            if rate:
                Q[i, _state_index(tuple(switched))] += rate
                exit_rate += rate
        for a in range(4):
            for b in range(a + 1, 4):
                if state[a] == state[b]:
                    rate = lambda0 if state[a] == 0 else lambda1
                    R[i, PAIR_TO_TOPOLOGY[frozenset((a, b))]] += rate
                    exit_rate += rate
        Q[i, i] -= exit_rate
    augmented = np.zeros((19, 19), dtype=float)
    augmented[:16, :16] = Q
    augmented[:16, 16:] = R
    result = np.zeros(19)
    result[_state_index(bits)] = 1.0
    result = result @ expm(augmented * duration)
    return result[16:] + result[:16].sum() / 3.0
