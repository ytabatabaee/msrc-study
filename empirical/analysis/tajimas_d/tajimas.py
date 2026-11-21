"""Classical Tajima's D and neutral-site summaries for fixed haplotype samples.

Input is one allele-count vector per completely called position. Variant-only
input is allowed only when it is a complete record of segregating sites for the
sample panel. Callable-base denominators are an independent input and are never
inferred from variant records.
"""
from __future__ import annotations

import math
from typing import Iterable, Optional, Sequence


def allele_counts(genotypes: Sequence[Optional[Sequence[int]]], ploidy: int):
    """Return allele counts for complete genotypes; reject partial/mixed ploidy.

    Each genotype is a sequence of allele indices or None. Missing genotypes
    yield None so callers can apply a fixed complete-case sample rule.
    """
    if ploidy not in (1, 2):
        raise ValueError("only haploid or diploid genotypes are supported")
    counts = {}
    for gt in genotypes:
        if gt is None or len(gt) == 0 or any(a is None for a in gt):
            return None
        if len(gt) != ploidy:
            raise ValueError("genotype ploidy differs from declared ploidy")
        for a in gt:
            if not isinstance(a, int) or a < 0:
                raise ValueError("allele indices must be nonnegative integers")
            counts[a] = counts.get(a, 0) + 1
    return [counts[k] for k in sorted(counts)]


def constants(n: int):
    """Return a1, a2, e1, e2 for the classical constant-n variance."""
    if n < 4:
        raise ValueError("classical Tajima variance requires at least four haplotypes")
    a1 = sum(1.0 / i for i in range(1, n))
    a2 = sum(1.0 / (i * i) for i in range(1, n))
    b1 = (n + 1.0) / (3.0 * (n - 1.0))
    b2 = 2.0 * (n * n + n + 3.0) / (9.0 * n * (n - 1.0))
    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2.0) / (a1 * n) + a2 / (a1 * a1)
    e1 = c1 / a1
    e2 = c2 / (a1 * a1 + a2)
    return a1, a2, e1, e2


def tajima_from_counts(site_counts: Iterable[Sequence[int]], n: int,
                       callable_bases: Optional[int] = None):
    """Compute S, pi_total, theta_W, and classical Tajima's D.

    `site_counts` contains fully called positions; monomorphic positions may be
    present but are ignored in S. `pi_total` is the sum over sites of unbiased
    average pairwise differences, in differences per sequence (not per base).
    `callable_bases` is used only for separately labeled per-base summaries.
    Undefined D is returned as None with an explicit reason.
    """
    if n < 1:
        raise ValueError("n must be positive")
    counts_rows = [tuple(int(x) for x in row) for row in site_counts]
    if any(any(x < 0 for x in row) for row in counts_rows):
        raise ValueError("allele counts cannot be negative")
    if any(sum(row) != n for row in counts_rows):
        raise ValueError("every retained site must have the same fixed haplotype count n")
    if callable_bases is not None and callable_bases <= 0:
        raise ValueError("callable_bases must be positive")
    S = 0
    pi_total = 0.0
    for row in counts_rows:
        nonzero = [x for x in row if x]
        if len(nonzero) > 1:
            S += 1
        pi_total += (n * n - sum(x * x for x in row)) / (n * (n - 1.0)) if n > 1 else 0.0
    theta_w = None
    var = None
    D = None
    reason = ""
    if n <= 3:
        reason = "n_haplotypes<=3"
    else:
        a1, a2, e1, e2 = constants(n)
        theta_w = S / a1
        var = e1 * S + e2 * S * (S - 1)
        if S == 0:
            reason = "no_segregating_sites"
        elif var <= 0 or not math.isfinite(var):
            reason = "nonpositive_or_nonfinite_classical_variance"
        else:
            D = (pi_total - theta_w) / math.sqrt(var)
    return {
        "n_haplotypes": n,
        "n_sites_examined": len(counts_rows),
        "S": S,
        "pi_total": pi_total,
        "theta_W": theta_w,
        "tajimas_D": D,
        "classical_variance": var,
        "callable_bases": callable_bases,
        "pi_per_callable_base": pi_total / callable_bases if callable_bases and n > 1 else None,
        "thetaW_per_callable_base": theta_w / callable_bases if callable_bases and theta_w is not None else None,
        "D_undefined_reason": reason,
    }


def fixed_complete_columns(site_genotypes: Sequence[Sequence[Optional[Sequence[int]]]],
                           sample_indices: Sequence[int], ploidy: int):
    """Find fixed sample subset called at every site; return (indices, counts).

    Genotypes are site-major. Samples absent/missing/heterozygous-invalid under
    a caller's ploidy policy are excluded from the fixed subset for this window.
    """
    complete = []
    for idx in sample_indices:
        ok = True
        for row in site_genotypes:
            if idx >= len(row):
                ok = False; break
            gt = row[idx]
            if gt is None or len(gt) != ploidy or any(a is None or a < 0 for a in gt):
                ok = False; break
        if ok:
            complete.append(idx)
    counts = []
    for row in site_genotypes:
        c = allele_counts([row[i] for i in complete], ploidy)
        # If a fixed subset could not be found, preserve an empty result.
        if c is None:
            return complete, []
        counts.append(c)
    n = len(complete) * ploidy
    if counts and any(sum(c) != n for c in counts):
        raise RuntimeError("fixed complete subset did not yield constant n")
    return complete, counts


def weighted_delta(D_pooled, D_X, D_Y, n_X, n_Y):
    """Return descriptive unweighted and sample-size-weighted D contrasts."""
    if any(x is None for x in (D_pooled, D_X, D_Y)):
        return None, None
    unweighted = D_pooled - (D_X + D_Y) / 2.0
    weighted = D_pooled - (n_X * D_X + n_Y * D_Y) / (n_X + n_Y)
    return unweighted, weighted
