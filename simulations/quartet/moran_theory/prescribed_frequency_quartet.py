"""Exact quartet probabilities conditional on prescribed frequency histories.

The frequency trajectories are supplied externally.  This module therefore
uses only the small background/genealogy generators and does not integrate a
Moran count process.  It is exact conditional on the supplied trajectories,
not an unconditional substitute for the finite-M Moran solver.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Iterable

import numpy as np
from scipy.linalg import expm

from .moran_quartet_operator import (
    BACKGROUND_STATES,
    PAIR_BACKGROUND_STATES,
    PATTERNS,
    residual_vectors,
)


@dataclass(frozen=True)
class PrescribedTrajectory:
    """Piecewise-constant frequency history, ordered from young to old.

    Each segment is ``(p, duration)``.  The first segment supplies the
    terminal sampling frequency.  Branch boundaries retain lineage background
    labels; the supplied trajectories are already the histories to condition
    on, so no daughter-count inheritance kernel is applied here.
    """

    segments: tuple[tuple[float, float], ...]

    def __post_init__(self) -> None:
        if not self.segments:
            raise ValueError("a prescribed trajectory needs at least one segment")
        for p, duration in self.segments:
            if not 0.0 <= float(p) <= 1.0:
                raise ValueError("frequency values must lie in [0, 1]")
            if float(duration) < 0.0:
                raise ValueError("segment durations must be nonnegative")

    @property
    def terminal_frequency(self) -> float:
        return float(self.segments[0][0])

    @classmethod
    def constant(cls, p: float, duration: float) -> "PrescribedTrajectory":
        return cls(((float(p), float(duration)),))


@dataclass(frozen=True)
class PrescribedQuartetResult:
    quartet_probabilities: np.ndarray
    early_topology_contributions: np.ndarray
    root_topology_contributions: np.ndarray
    root_entry: np.ndarray
    joint: np.ndarray | None
    pattern_probabilities: np.ndarray | None
    conditional_quartet_probabilities: np.ndarray | None
    diagnostics: dict[str, float]


def _rates(p: float, M: int) -> tuple[float, float]:
    if M <= 0:
        raise ValueError("M must be positive")
    c1 = 1.0 / (M * max(float(p), 1.0 / M))
    c0 = 1.0 / (M * max(1.0 - float(p), 1.0 / M))
    return c0, c1


def prescribed_generator(n: int, p: float, M: int, rho: float) -> np.ndarray:
    """Return the row generator G_n for a prescribed frequency p."""
    if n not in {1, 2, 3, 4}:
        raise ValueError("n must be 1, 2, 3, or 4")
    if not 0.0 <= float(p) <= 1.0 or rho < 0.0:
        raise ValueError("p must be in [0,1] and rho must be nonnegative")
    c0, c1 = _rates(float(p), M)
    states = {
        1: tuple(product((0, 1), repeat=1)),
        2: tuple(product((0, 1), repeat=2)),
        3: tuple(product((0, 1), repeat=3)),
        4: tuple(product((0, 1), repeat=4)),
    }[n]
    index = {state: i for i, state in enumerate(states)}
    transient = len(states)
    absorbing = 0 if n == 1 else (2 if n == 2 else (3 if n in {3, 4} else 0))
    out = np.zeros((transient + absorbing, transient + absorbing), dtype=float)

    for row, state in enumerate(states):
        exit_rate = 0.0
        for i, bit in enumerate(state):
            rate = rho * float(p) if bit == 0 else rho * (1.0 - float(p))
            if rate:
                switched = list(state)
                switched[i] = 1 - bit
                out[row, index[tuple(switched)]] += rate
                exit_rate += rate

        if n == 2:
            if state == (0, 0):
                out[row, transient] = c0
                exit_rate += c0
            elif state == (1, 1):
                out[row, transient + 1] = c1
                exit_rate += c1
        elif n in {3, 4}:
            pairs = (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))) if n == 4 else (((0, 1),), ((0, 2),), ((1, 2),))
            for topology, pair_set in enumerate(pairs):
                for i, j in pair_set:
                    if state[i] == state[j]:
                        out[row, transient + topology] += c1 if state[i] else c0
                        exit_rate += c1 if state[i] else c0
        out[row, row] = -exit_rate

    if n == 2:
        out[transient, transient + 1] = rho * float(p)
        out[transient, transient] = -out[transient, transient + 1]
        out[transient + 1, transient] = rho * (1.0 - float(p))
        out[transient + 1, transient + 1] = -out[transient + 1, transient]
    return out


def _propagate_row(row: np.ndarray, trajectory: PrescribedTrajectory, n: int, M: int, rho: float) -> np.ndarray:
    value = np.asarray(row, dtype=float)
    for p, duration in trajectory.segments:
        if duration:
            value = value @ expm(prescribed_generator(n, p, M, rho) * float(duration))
    return value


def _terminal_row(trajectory: PrescribedTrajectory, bit: int | None) -> np.ndarray:
    p = trajectory.terminal_frequency
    if bit is None:
        return np.array([1.0 - p, p], dtype=float)
    if bit not in (0, 1):
        raise ValueError("terminal bits must be 0 or 1")
    if (bit == 0 and p == 1.0) or (bit == 1 and p == 0.0):
        return np.zeros(2, dtype=float)
    out = np.zeros(2, dtype=float)
    out[bit] = 1.0
    return out


def _lineage_at_parent(trajectory: PrescribedTrajectory, M: int, rho: float, bit: int | None) -> np.ndarray:
    return _propagate_row(_terminal_row(trajectory, bit), trajectory, 1, M, rho)


def _pair_at_parent(
    terminal_trajectories: tuple[PrescribedTrajectory, PrescribedTrajectory],
    sister_trajectory: PrescribedTrajectory,
    M_terminal: tuple[int, int], M_sister: int, rho: float,
    bits: tuple[int, int] | None,
) -> np.ndarray:
    u1 = _lineage_at_parent(terminal_trajectories[0], M_terminal[0], rho, None if bits is None else bits[0])
    u2 = _lineage_at_parent(terminal_trajectories[1], M_terminal[1], rho, None if bits is None else bits[1])
    young = np.zeros(6, dtype=float)
    for b1, b2 in PAIR_BACKGROUND_STATES:
        young[2 * b1 + b2] = u1[b1] * u2[b2]
    return _propagate_row(young, sister_trajectory, 2, M_sister, rho)


def _root_four_lineage(
    F: np.ndarray, root_trajectory: PrescribedTrajectory, M_root: int, rho: float,
) -> np.ndarray:
    return _propagate_row(F, root_trajectory, 4, M_root, rho)


def _pattern_mass(terminal_trajectories: tuple[PrescribedTrajectory, ...], pattern: str) -> float:
    mass = 1.0
    for trajectory, bit in zip(terminal_trajectories, pattern):
        p = trajectory.terminal_frequency
        mass *= p if bit == "1" else 1.0 - p
    return float(mass)


def balanced_prescribed_quartet_probabilities(
    *, M_root: int, M_left: int, M_right: int,
    M_terminal: tuple[int, int, int, int], rho: float,
    terminal_trajectories: tuple[PrescribedTrajectory, PrescribedTrajectory, PrescribedTrajectory, PrescribedTrajectory],
    sister_trajectories: tuple[PrescribedTrajectory, PrescribedTrajectory],
    root_trajectory: PrescribedTrajectory,
    with_joint: bool = False,
) -> PrescribedQuartetResult:
    """Evaluate the exact prescribed-history balanced-tree recursion."""
    left_by = {
        bits: _pair_at_parent(terminal_trajectories[:2], sister_trajectories[0], M_terminal[:2], M_left, rho, bits)
        for bits in PAIR_BACKGROUND_STATES
    }
    right_by = {
        bits: _pair_at_parent(terminal_trajectories[2:], sister_trajectories[1], M_terminal[2:], M_right, rho, bits)
        for bits in PAIR_BACKGROUND_STATES
    }
    left = _pair_at_parent(
        terminal_trajectories[:2], sister_trajectories[0], M_terminal[:2], M_left, rho, None,
    )
    right = _pair_at_parent(
        terminal_trajectories[2:], sister_trajectories[1], M_terminal[2:], M_right, rho, None,
    )
    s_left, s_right = left[:4].sum(), right[:4].sum()
    early = np.array([1.0 - s_left * s_right, 0.0, 0.0])
    F = np.zeros(19, dtype=float)
    for il, zl in enumerate(PAIR_BACKGROUND_STATES):
        for ir, zr in enumerate(PAIR_BACKGROUND_STATES):
            state = BACKGROUND_STATES.index(zl + zr)
            F[state] = left[il] * right[ir]
    root_state = _root_four_lineage(F, root_trajectory, M_root, rho)
    root = np.array([root_state @ h for h in residual_vectors()])
    q = early + root
    joint = patterns = conditional = None
    if with_joint:
        joint = np.zeros((len(PATTERNS), 3), dtype=float)
        for zi, pattern in enumerate(PATTERNS):
            mass = _pattern_mass(terminal_trajectories, pattern)
            if mass == 0.0:
                continue
            left_z = left_by[tuple(int(x) for x in pattern[:2])]
            right_z = right_by[tuple(int(x) for x in pattern[2:])]
            Fz = np.zeros(19, dtype=float)
            for il, zl in enumerate(PAIR_BACKGROUND_STATES):
                for ir, zr in enumerate(PAIR_BACKGROUND_STATES):
                    Fz[BACKGROUND_STATES.index(zl + zr)] = left_z[il] * right_z[ir]
            early_z = np.array([1.0 - left_z[:4].sum() * right_z[:4].sum(), 0.0, 0.0])
            root_z = _root_four_lineage(Fz, root_trajectory, M_root, rho)
            joint[zi] = mass * (
                early_z + np.array([root_z @ h for h in residual_vectors()])
            )
        patterns = joint.sum(axis=1)
        conditional = np.divide(joint, patterns[:, None], out=np.full_like(joint, np.nan), where=patterns[:, None] > 0)
    diagnostics = {
        "probability_sum_error": float(q.sum() - 1.0),
        "early_root_sum_error": float(early.sum() + root.sum() - q.sum()),
        "root_entry_mass": float(F.sum()),
    }
    return PrescribedQuartetResult(q, early, root, F, joint, patterns, conditional, diagnostics)


def unbalanced_prescribed_quartet_probabilities(
    *, M_root: int, M_intermediate: int, M_sister: int,
    M_terminal: tuple[int, int, int, int], rho: float,
    terminal_trajectories: tuple[PrescribedTrajectory, PrescribedTrajectory, PrescribedTrajectory, PrescribedTrajectory],
    sister_trajectory: PrescribedTrajectory,
    intermediate_trajectory: PrescribedTrajectory,
    root_trajectory: PrescribedTrajectory,
    with_joint: bool = False,
) -> PrescribedQuartetResult:
    """Evaluate the exact prescribed-history pectinate recursion."""
    def evaluate(pattern: str | None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        pair = _pair_at_parent(terminal_trajectories[:2], sister_trajectory, M_terminal[:2], M_sister, rho, None if pattern is None else (int(pattern[0]), int(pattern[1])))
        u3 = _lineage_at_parent(terminal_trajectories[2], M_terminal[2], rho, None if pattern is None else int(pattern[2]))
        young = np.zeros(11, dtype=float)
        young[8] = pair[4] + pair[5]
        for index, state in enumerate(product((0, 1), repeat=3)):
            pair_index = PAIR_BACKGROUND_STATES.index(state[:2])
            young[index] = pair[pair_index] * u3[state[2]]
        triple = _propagate_row(young, intermediate_trajectory, 3, M_intermediate, rho)
        u4 = _lineage_at_parent(terminal_trajectories[3], M_terminal[3], rho, None if pattern is None else int(pattern[3]))
        F = np.zeros(19, dtype=float)
        for state, index in zip(product((0, 1), repeat=3), range(8)):
            for b4 in (0, 1):
                F[BACKGROUND_STATES.index(state + (b4,))] = triple[index] * u4[b4]
        return triple[8:].copy(), F, triple

    early, F, _ = evaluate(None)
    root_state = _root_four_lineage(F, root_trajectory, M_root, rho)
    root = np.array([root_state @ h for h in residual_vectors()])
    q = early + root
    joint = patterns = conditional = None
    if with_joint:
        joint = np.zeros((len(PATTERNS), 3), dtype=float)
        for zi, pattern in enumerate(PATTERNS):
            mass = _pattern_mass(terminal_trajectories, pattern)
            if mass == 0.0:
                continue
            early_z, Fz, _ = evaluate(pattern)
            root_z = _root_four_lineage(Fz, root_trajectory, M_root, rho)
            joint[zi] = mass * (
                early_z + np.array([root_z @ h for h in residual_vectors()])
            )
        patterns = joint.sum(axis=1)
        conditional = np.divide(joint, patterns[:, None], out=np.full_like(joint, np.nan), where=patterns[:, None] > 0)
    diagnostics = {
        "probability_sum_error": float(q.sum() - 1.0),
        "early_root_sum_error": float(early.sum() + root.sum() - q.sum()),
        "root_entry_mass": float(F.sum()),
    }
    return PrescribedQuartetResult(q, early, root, F, joint, patterns, conditional, diagnostics)


__all__ = [
    "PrescribedTrajectory", "PrescribedQuartetResult", "prescribed_generator",
    "balanced_prescribed_quartet_probabilities",
    "unbalanced_prescribed_quartet_probabilities",
]
