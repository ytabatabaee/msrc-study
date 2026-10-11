#!/usr/bin/env python3
"""Explicit SimPhy-to-msrc-sim population-size conversion.

SimPhy D1 logs label -sp as a haploid effective population size N_H.  The
ordinary diploid coalescent in msrc-sim uses a pair rate 1/(2 N_e), so the
continuous mapping is N_e = N_H / 2.  Finite Wright-Fisher histories require
an integer diploid N_e and therefore use a separately reported approximation.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PopulationMapping:
    simphy_haploid_ne: int
    msrc_diploid_ne_exact: float
    msrc_diploid_ne_integer: int
    simphy_pair_coalescence_rate: float
    msrc_pair_coalescence_rate: float
    relative_rate_error: float


def map_simphy_haploid_to_msrc(simphy_haploid_ne: int) -> PopulationMapping:
    """Map one SimPhy D1 population size and quantify integer approximation."""
    n_h = int(simphy_haploid_ne)
    if n_h <= 0 or n_h != simphy_haploid_ne:
        raise ValueError("SimPhy haploid effective population size must be a positive integer")
    exact = n_h / 2.0
    candidates = {max(1, int(exact // 1)), max(1, int(exact // 1) + 1)}
    source_rate = 1.0 / n_h
    integer = min(candidates, key=lambda n: abs((1.0 / (2.0 * n)) - source_rate))
    msrc_rate = 1.0 / (2.0 * integer)
    return PopulationMapping(
        simphy_haploid_ne=n_h,
        msrc_diploid_ne_exact=exact,
        msrc_diploid_ne_integer=integer,
        simphy_pair_coalescence_rate=source_rate,
        msrc_pair_coalescence_rate=msrc_rate,
        relative_rate_error=(msrc_rate - source_rate) / source_rate,
    )


def generations_to_coalescent_units(generations: float, mapping: PopulationMapping, *, integer: bool = False) -> float:
    """Convert a branch measured in generations using the chosen convention."""
    ne = mapping.msrc_diploid_ne_integer if integer else mapping.msrc_diploid_ne_exact
    return float(generations) / (2.0 * ne)


def coalescent_units_to_generations(cu: float, mapping: PopulationMapping, *, integer: bool = False) -> float:
    ne = mapping.msrc_diploid_ne_integer if integer else mapping.msrc_diploid_ne_exact
    return float(cu) * (2.0 * ne)
