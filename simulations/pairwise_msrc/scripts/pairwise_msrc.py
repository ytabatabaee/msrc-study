"""Exact two-lineage MSRC expectations for piecewise-constant histories.

The backward event generator mirrors msrcsim.structured_coalescent. Time is
in the same generation-equivalent units as the supplied Ne and rho values.
Background A is simulator state 0; background B is state 1 and has frequency p.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.linalg import expm

STATE_INDEX = {"AA": 0, "AB": 1, "BB": 2}


@dataclass(frozen=True)
class Segment:
    duration: float
    frequency_b: float
    effective_population_size: float
    rho: float


def transient_generator(p: float, ne: float, rho: float,
                        finite_copy_floor: bool = True) -> np.ndarray:
    """Row generator on transient states (AA, AB, BB).

    The 1/(2Ne) frequency floor is the legacy simulator convention. In
    particular, at p=0 an explicitly conditioned BB pair has coalescence rate
    1 per generation-equivalent time; this is a regularization, not the literal
    biological limit for an absent background.
    """
    p, ne, rho = float(p), float(ne), float(rho)
    if not (0.0 <= p <= 1.0):
        raise ValueError("frequency_b must lie in [0, 1]")
    if ne <= 0 or rho < 0:
        raise ValueError("Ne must be positive and rho nonnegative")
    floor = 1.0 / (2.0 * ne) if finite_copy_floor else 0.0
    c_a = 1.0 / (2.0 * ne * max(1.0 - p, floor))
    c_b = 1.0 / (2.0 * ne * max(p, floor))
    q_a_to_b = rho * p
    q_b_to_a = rho * (1.0 - p)
    return np.asarray([
        [-(2.0 * q_a_to_b + c_a), 2.0 * q_a_to_b, 0.0],
        [q_b_to_a, -(q_a_to_b + q_b_to_a), q_a_to_b],
        [0.0, 2.0 * q_b_to_a, -(2.0 * q_b_to_a + c_b)],
    ], dtype=float)


def _ordinary_wait(segments: Sequence[tuple[float, float | None]],
                   ancestral_ne: float) -> float:
    """Expected wait from origin through ordinary MSC segments and tail.

    A finite segment is (duration, Ne). Ne=None denotes a no-coalescence
    interval, such as two lineages in distinct populations. The final ancestor
    has constant Ne and extends indefinitely.
    """
    if ancestral_ne <= 0:
        raise ValueError("ancestral_ne must be positive")
    survival = 1.0
    wait = 0.0
    for duration, ne in segments:
        duration = float(duration)
        if duration < 0:
            raise ValueError("ancestral segment durations must be nonnegative")
        if duration == 0:
            continue
        rate = 0.0 if ne is None else 1.0 / (2.0 * float(ne))
        if rate == 0.0:
            wait += survival * duration
        else:
            wait += survival * (-np.expm1(-rate * duration)) / rate
            survival *= np.exp(-rate * duration)
    return float(wait + survival * 2.0 * float(ancestral_ne))


def expected_pair_times(
    segments: Sequence[Segment],
    post_origin_segments: Sequence[tuple[float, float | None]] = (),
    ancestral_ne: float = 10_000.0,
    finite_copy_floor: bool = True,
) -> dict[str, float]:
    """Return expected T_AA, T_AB, T_BB from the present to coalescence.

    The structured phase is propagated exactly per constant segment with an
    augmented 4x4 matrix exponential. The fourth coordinate accumulates
    transient survival time. At the finite structured-origin boundary, labeled
    states collapse into the ordinary MSC described by post_origin_segments
    and an infinite constant-Ne ancestral tail.
    """
    probability = np.eye(3)  # columns are initial states
    survival_integral = np.zeros(3)
    structured_age = 0.0
    ones = np.ones(3)
    for seg in segments:
        duration = float(seg.duration)
        if duration < 0:
            raise ValueError("structured segment durations must be nonnegative")
        if duration == 0:
            continue
        q = transient_generator(seg.frequency_b, seg.effective_population_size,
                                 seg.rho, finite_copy_floor)
        augmented = np.zeros((4, 4), dtype=float)
        augmented[:3, :3] = q.T
        augmented[3, :3] = ones
        transition_integral = expm(augmented * duration)
        survival_integral += transition_integral[3, :3] @ probability
        probability = transition_integral[:3, :3] @ probability
        structured_age += duration
    survival_at_origin = ones @ probability
    ordinary_residual = _ordinary_wait(post_origin_segments, ancestral_ne)
    # The survival integral already includes the elapsed structured-phase
    # contribution for lineages that survive to the boundary; add only their
    # expected *residual* ordinary-MSC waiting time, not the boundary age again.
    values = survival_integral + survival_at_origin * ordinary_residual
    return {state: float(values[i]) for state, i in STATE_INDEX.items()}


def ordinary_pair_mean_time(
    pre_origin_segments: Sequence[Segment],
    post_origin_segments: Sequence[tuple[float, float | None]] = (),
    ancestral_ne: float = 10_000.0,
) -> float:
    """Matched unstructured MSC mean for the same demography and durations."""
    parts = [(s.duration, s.effective_population_size) for s in pre_origin_segments]
    parts.extend(post_origin_segments)
    return _ordinary_wait(parts, ancestral_ne)


def predicted_statistics(times: dict[str, float], mutation_rate: float) -> dict[str, float]:
    """Infinite-sites moments and ratio of expected statistics.

    ``R_ratio_of_expectations`` is not E[2 dXY/(pi_A+pi_B)]; it is the
    large-window target formed from expected diversity and divergence.
    """
    mu = float(mutation_rate)
    if mu < 0:
        raise ValueError("mutation_rate must be nonnegative")
    pi_a = 2.0 * mu * times["AA"]
    pi_b = 2.0 * mu * times["BB"]
    d_ab = 2.0 * mu * times["AB"]
    denominator = pi_a + pi_b
    delta = d_ab - denominator / 2.0
    ratio = 2.0 * d_ab / denominator if denominator > 0 else float("nan")
    return {
        "E_pi_A": pi_a,
        "E_pi_B": pi_b,
        "E_d_AB": d_ab,
        "E_Delta_AB": delta,
        "R_ratio_of_expectations": ratio,
        "FST_ratio_of_expectations": delta / d_ab if d_ab > 0 else float("nan"),
    }


def history_segments(history, branch_id: str, ne: float, rho: float,
                     duration: float, max_segments: int = 1_000_000) -> list[Segment]:
    """Convert an msrc-sim piecewise-constant frequency path to segments."""
    t = 0.0
    out: list[Segment] = []
    while t < duration - 1e-12:
        p = float(history.frequency_at(branch_id, t))
        boundary = float(history.next_frequency_boundary(branch_id, t))
        end = min(float(duration), boundary)
        if not np.isfinite(end) or end <= t + 1e-12:
            end = float(duration)
        out.append(Segment(end - t, p, float(ne), float(rho)))
        t = end
        if len(out) > max_segments:
            raise RuntimeError("frequency history generated too many segments")
    return out
