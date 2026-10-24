"""Pure numerical Patterson's D helpers for Anopheles 2La."""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

@dataclass(frozen=True)
class DResult:
    abba: float
    baba: float
    denominator: float
    d: float | None
    n_examined: int
    n_pass: int
    n_informative: int
    exclusions: dict[str, int]
    mean_callable_chromosomes: dict[str, float]
    outgroup_fixed_ancestral: int

def _counts(gt: np.ndarray, indices: list[int]) -> tuple[np.ndarray, np.ndarray]:
    sub = gt[:, indices, :]
    called = (sub >= 0).all(axis=2)
    chromosomes = np.where(called, sub, 0).sum(axis=(1, 2)).astype(float)
    n_called = called.sum(axis=1).astype(float) * 2
    return chromosomes, n_called

def calculate_window(
    gt: np.ndarray,
    groups: dict[str, list[int]],
    outgroup: list[int],
    *,
    max_site_missing: float = 0.25,
    min_callable_chromosomes: int = 2,
    min_mac: int = 0,
) -> DResult:
    """Calculate ratio-of-sums D from variants x samples x ploidy calls.

    Allele codes 0 and 1 are REF/ALT. The outgroup must be callable and fixed
    for one allele at a site; B is the allele different from that state.
    Multiallelic, ambiguous, missing, and low-MAC sites are excluded before
    calculating either ABBA or BABA.
    """
    gt = np.asarray(gt)
    if gt.ndim != 3 or gt.shape[2] != 2:
        raise ValueError("expected variants x samples x diploid-ploidy calls")
    n_sites = gt.shape[0]
    exclusions = {k: 0 for k in ("multiallelic", "site_missing", "insufficient_group_calls", "outgroup_missing", "outgroup_polymorphic", "outgroup_nonbinary", "low_mac")}
    abba = baba = 0.0
    informative = passed = fixed = 0
    call_totals = {k: 0.0 for k in ("P1", "P2", "P3", "O")}
    for site in range(n_sites):
        row = gt[site]
        called = row >= 0
        if np.any(row[called] > 1):
            exclusions["multiallelic"] += 1; continue
        if called.mean() == 0 or 1.0 - float(called.all(axis=1).mean()) > max_site_missing:
            exclusions["site_missing"] += 1; continue
        group_counts = {}
        group_called = {}
        for name, idx in (("P1", groups["P1"]), ("P2", groups["P2"]), ("P3", groups["P3"]), ("O", outgroup)):
            sub = row[idx]
            ok = (sub >= 0).all(axis=1)
            group_called[name] = int(ok.sum() * 2)
            group_counts[name] = float(np.where(ok, sub, 0).sum())
            call_totals[name] += group_called[name]
        if any(group_called[k] < min_callable_chromosomes for k in group_called):
            exclusions["insufficient_group_calls"] += 1; continue
        out_alleles = row[outgroup][(row[outgroup] >= 0)].astype(int)
        if out_alleles.size == 0:
            exclusions["outgroup_missing"] += 1; continue
        if np.any(out_alleles > 1):
            exclusions["outgroup_nonbinary"] += 1; continue
        if np.any(out_alleles != out_alleles[0]):
            exclusions["outgroup_polymorphic"] += 1; continue
        ancestral = int(out_alleles[0]); fixed += 1
        total_alt = sum(group_counts.values())
        total_called = sum(group_called.values())
        mac = min(total_alt, total_called - total_alt)
        if mac < min_mac:
            exclusions["low_mac"] += 1; continue
        passed += 1
        p = {name: (group_counts[name] / group_called[name]) for name in ("P1", "P2", "P3")}
        if ancestral == 1:
            p = {name: 1.0 - value for name, value in p.items()}
        a = (1.0 - p["P1"]) * p["P2"] * p["P3"]
        b = p["P1"] * (1.0 - p["P2"]) * p["P3"]
        abba += a; baba += b
        if a + b > 0: informative += 1
    denom = abba + baba
    d = (abba - baba) / denom if denom > 0 else None
    means = {name: call_totals[name] / n_sites if n_sites else 0.0 for name in call_totals}
    return DResult(abba, baba, denom, d, n_sites, passed, informative, exclusions, means, fixed)

def jackknife(block_values: list[tuple[float, float, int]], alpha: float = 0.05) -> dict[str, float | int | None]:
    """Delete-one-block jackknife for block ABBA/BABA sums."""
    usable = [(a, b) for a, b, n in block_values if a + b > 0 and n > 0]
    total_a, total_b = sum(a for a, _, _ in block_values), sum(b for _, b, _ in block_values)
    total_d = (total_a - total_b) / (total_a + total_b) if total_a + total_b else None
    n_blocks = len(usable)
    if n_blocks < 2 or total_d is None:
        return {"pooled_d": total_d, "se": None, "z": None, "ci_low": None, "ci_high": None, "n_blocks": n_blocks}
    estimates = []
    for a, b, _ in block_values:
        aa, bb = total_a - a, total_b - b
        estimates.append((aa - bb) / (aa + bb) if aa + bb else np.nan)
    estimates = np.asarray(estimates, dtype=float); estimates = estimates[np.isfinite(estimates)]
    if len(estimates) < 2:
        return {"pooled_d": total_d, "se": None, "z": None, "ci_low": None, "ci_high": None, "n_blocks": n_blocks}
    mean = float(estimates.mean())
    se = float(math.sqrt((len(estimates) - 1) / len(estimates) * np.sum((estimates - mean) ** 2)))
    return {"pooled_d": total_d, "se": se, "z": total_d / se if se else None, "ci_low": total_d - 1.96 * se if se else None, "ci_high": total_d + 1.96 * se if se else None, "n_blocks": n_blocks}

def block_id(position: int, size: int) -> str:
    return str(((int(position) - 1) // size) * size + 1)
