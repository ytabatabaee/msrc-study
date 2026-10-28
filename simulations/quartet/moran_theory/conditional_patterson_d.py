"""Branch-length Patterson's D summaries for four-taxon MSRC genealogies.

The statistic is computed from the explicit first two coalescences in a
complete rooted genealogy.  For taxa (P1, P2, P3, O), ABBA is supported by
the P2|P3 internal branch (topology P1O|P2P3), while BABA is supported by the
P1|P3 internal branch (topology P2O|P1P3).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np


TAXA = ("P1", "P2", "P3", "O")
Q1 = frozenset(("P1", "P2"))
Q2 = frozenset(("P1", "P3"))
Q3 = frozenset(("P1", "O"))


@dataclass(frozen=True)
class GenealogyBranchSummary:
    topology: int
    internal_branch_length: float


@dataclass(frozen=True)
class PattersonSummary:
    n_genealogies: int
    topology_counts: np.ndarray
    topology_frequencies: np.ndarray
    topology_branch_sums: np.ndarray
    mean_internal_branch_by_topology: np.ndarray
    abba_sum: float
    baba_sum: float
    d_patterson: float
    d_topology: float
    d_branch: float
    informative_genealogies: int


def _descendant_set(text: str) -> frozenset[str]:
    return frozenset(part for part in text.split(",") if part)


def classify_genealogy(result, taxa: tuple[str, str, str, str] = TAXA) -> GenealogyBranchSummary:
    """Classify a complete four-taxon genealogy and extract its internal edge.

    The first coalescence identifies one side of the unrooted quartet.  The
    interval between the first and second coalescences is the branch on which
    an infinite-sites mutation produces the corresponding two-two pattern.
    """
    if len(result.coalescences) < 2 or len(result.coalescence_times) < 2:
        raise ValueError("a four-taxon genealogy must contain two coalescences")
    first = result.coalescences[0]
    pair = _descendant_set(first.descendant_set_1) | _descendant_set(first.descendant_set_2)
    if len(pair) != 2:
        raise ValueError(f"first coalescence is not a two-taxon pair: {sorted(pair)}")
    complement = frozenset(taxa) - pair
    split = pair if len(pair) == 2 else complement
    if split == Q1 or complement == Q1:
        topology = 0
    elif split == Q2 or complement == Q2:
        topology = 1
    elif split == Q3 or complement == Q3:
        topology = 2
    else:
        raise ValueError(f"unrecognized quartet split: {sorted(pair)}")
    t1, t2 = map(float, result.coalescence_times[:2])
    if t2 < t1:
        raise ValueError("coalescence times are not nondecreasing")
    return GenealogyBranchSummary(topology, t2 - t1)


def summarize_genealogies(results: Iterable, taxa: tuple[str, str, str, str] = TAXA) -> PattersonSummary:
    counts = np.zeros(3, dtype=int)
    branch_sums = np.zeros(3, dtype=float)
    informative = 0
    for result in results:
        summary = classify_genealogy(result, taxa)
        counts[summary.topology] += 1
        branch_sums[summary.topology] += summary.internal_branch_length
        if summary.topology in (1, 2) and summary.internal_branch_length > 0:
            informative += 1
    n = int(counts.sum())
    freqs = counts / n if n else np.full(3, np.nan)
    means = np.divide(branch_sums, counts, out=np.full(3, np.nan), where=counts > 0)
    # ABBA = P2,P3 branch (Q3); BABA = P1,P3 branch (Q2).
    abba = float(branch_sums[2])
    baba = float(branch_sums[1])
    denom = abba + baba
    d = (abba - baba) / denom if denom > 0 else np.nan
    discordant = counts[1] + counts[2]
    d_top = (counts[1] - counts[2]) / discordant if discordant else np.nan
    branch_denom = branch_sums[1] + branch_sums[2]
    d_branch = (branch_sums[1] - branch_sums[2]) / branch_denom if branch_denom > 0 else np.nan
    return PattersonSummary(n, counts, freqs, branch_sums, means, abba, baba, d,
                            d_top, d_branch, informative)


def simulate_infinite_sites_counts(results: Iterable, sequence_replicates: int, mutation_rate: float,
                                   rng: np.random.Generator,
                                   taxa: tuple[str, str, str, str] = TAXA) -> tuple[int, int]:
    """Simulate only ABBA/BABA mutations under infinite sites.

    Mutations on terminal or other non-two-two branches cannot contribute to
    ABBA/BABA, so the two relevant internal edges are sufficient.
    """
    abba = baba = 0
    for result in results:
        summary = classify_genealogy(result, taxa)
        n = int(sequence_replicates)
        mutations = int(rng.poisson(max(0.0, mutation_rate * summary.internal_branch_length) * n))
        if summary.topology == 2:
            abba += mutations
        elif summary.topology == 1:
            baba += mutations
    return abba, baba


def ratio_d(abba: float, baba: float) -> float:
    denominator = float(abba) + float(baba)
    return (float(abba) - float(baba)) / denominator if denominator > 0 else float("nan")


def bootstrap_summary(results: list, n_bootstrap: int, rng: np.random.Generator,
                      taxa: tuple[str, str, str, str] = TAXA) -> dict[str, float]:
    """Genealogy-resampling uncertainty; independent history uncertainty is separate."""
    if not results or n_bootstrap <= 0:
        return {"D_patterson_se": float("nan"), "D_topology_se": float("nan"),
                "D_branch_se": float("nan")}
    values = np.empty((n_bootstrap, 3), dtype=float)
    for i in range(n_bootstrap):
        sample = [results[int(j)] for j in rng.integers(0, len(results), len(results))]
        summary = summarize_genealogies(sample, taxa)
        values[i] = (summary.d_patterson, summary.d_topology, summary.d_branch)
    return {
        "D_patterson_se": float(np.nanstd(values[:, 0], ddof=1)),
        "D_topology_se": float(np.nanstd(values[:, 1], ddof=1)),
        "D_branch_se": float(np.nanstd(values[:, 2], ddof=1)),
    }


__all__ = ["GenealogyBranchSummary", "PattersonSummary", "classify_genealogy",
           "summarize_genealogies", "simulate_infinite_sites_counts", "ratio_d",
           "bootstrap_summary"]
