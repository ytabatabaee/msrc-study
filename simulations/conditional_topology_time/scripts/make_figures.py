#!/usr/bin/env python3
"""Make figures from the conditional topology/time result tables."""
from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "simulations/conditional_topology_time"
OUT = BASE / "figures"


def read(path: Path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def save(fig, name: str):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=220, bbox_inches="tight")


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    result = BASE / "results/validation"
    marginal = read(result / "marginal_summary.tsv")
    pattern_probs = read(result / "pattern_probabilities.tsv")
    times = read(result / "conditional_times.tsv")
    quartets = read(result / "conditional_quartets.tsv")
    decomp = read(result / "marginal_decomposition.tsv")
    colors = {"all": "#555555", "concordant_2_2": "#3b82f6", "discordant_2_2": "#ef4444", "3_1": "#10b981", "uniform": "#a855f7"}

    # Figure A: category-specific deviations from the marginal MSRC mean,
    # with the matched MSC mean shown separately as an unconditional comparator.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, M, rho in zip(axes, (20, 40), (0.0, 0.5)):
        cell = next(r for r in marginal if r["M"] == str(M) and r["rho"] == str(rho))
        cell_times = [r for r in times if r["M"] == str(M) and r["rho"] == str(rho) and r["statistic"] == "mean_pairwise_time"]
        cell_probs = [r for r in pattern_probs if r["M"] == str(M) and r["rho"] == str(rho)]
        overall = float(cell["msrc_mean_pairwise_time"])
        labels = ["all"] + ["concordant_2_2", "discordant_2_2", "3_1", "uniform"]
        vals = [0.0]
        lo = [0.0]
        hi = [0.0]
        for label in labels[1:]:
            sub = [r for r in cell_times if r["pattern_category"] == label]
            # Weight patterns by their observed natural frequencies within the
            # category, so this is a descriptive category conditional mean.
            prob_by_pattern = {r["pattern"]: float(r["empirical_probability"]) for r in cell_probs}
            weights = np.asarray([prob_by_pattern.get(r["pattern"], 0.0) for r in sub], dtype=float)
            if sub and weights.sum() > 0:
                means = np.asarray([float(r["mean"]) for r in sub])
                lows = np.asarray([float(r["bootstrap_lo"]) for r in sub])
                highs = np.asarray([float(r["bootstrap_hi"]) for r in sub])
                vals.append(float(np.average(means, weights=weights)) - overall)
                lo.append(float(np.average(lows, weights=weights)) - overall)
                hi.append(float(np.average(highs, weights=weights)) - overall)
            else:
                vals.append(math.nan); lo.append(math.nan); hi.append(math.nan)
        x = np.arange(len(labels))
        ax.bar(x, vals, color=[colors.get(x, "#555") for x in labels], alpha=.85)
        ax.vlines(x[1:], lo[1:], hi[1:], color="black", lw=1)
        ax.axhline(float(cell["msc_mean_pairwise_time"]) - overall, color="black", ls="--", lw=1, label="MSC minus MSRC")
        ax.axhline(0, color="#999", lw=.7)
        ax.set_xticks(x, ["all", "2:2 C", "2:2 D", "3:1", "uniform"], rotation=35, ha="right")
        ax.set_title(f"M={M}, rho={rho}"); ax.set_ylabel("mean-time difference from marginal MSRC")
    axes[0].legend(fontsize=8); fig.tight_layout(); save(fig, "conditional_vs_marginal_effects"); plt.close(fig)

    # Figure B: theory versus simulation for informative patterns.
    patterns = ["0011", "1010", "0101", "0110"]
    fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharey=True)
    for ax, pattern in zip(axes.flat, patterns):
        rows = [r for r in quartets if r["pattern"] == pattern and r["M"] == "40" and r["rho"] == "0.5"]
        x = np.arange(3); width = .36
        sim = [float(r["empirical_probability"]) for r in rows]
        th = [float(r["theory_probability"]) for r in rows]
        se = [float(r["empirical_mc_se"]) for r in rows]
        ax.bar(x - width / 2, sim, width, yerr=se, capsize=2, label="simulation")
        ax.bar(x + width / 2, th, width, label="theory")
        ax.set_title(pattern); ax.set_xticks(x, ["Q1", "Q2", "Q3"]); ax.set_ylim(0, .7)
    axes[0, 0].set_ylabel("P(Q | Z)"); axes[1, 0].set_ylabel("P(Q | Z)"); axes[0, 0].legend(fontsize=8)
    fig.suptitle("Conditional quartet probabilities: M=40, rho=0.5")
    fig.tight_layout(); save(fig, "conditional_quartet_probabilities"); plt.close(fig)

    # Figure C: conditional quantile bands for natural pattern strata.
    fig, ax = plt.subplots(figsize=(9, 5))
    for pattern in ["0011", "1010", "0101", "0110"]:
        rows = [r for r in times if r["pattern"] == pattern and r["M"] == "40" and r["rho"] == "0.5" and r["statistic"] == "mean_pairwise_time"]
        if not rows: continue
        r = rows[0]; x = float(r["mean"]); yerr = [[x - float(r["q05"])], [float(r["q95"]) - x]]
        ax.errorbar([pattern], [x], yerr=yerr, fmt="o", capsize=4, label=pattern)
    ax.set_ylabel("mean pairwise coalescence time"); ax.set_xlabel("terminal arrangement pattern"); ax.set_title("Conditional genealogy-time distributions: M=40, rho=0.5")
    fig.tight_layout(); save(fig, "conditional_coalescence_time_distributions"); plt.close(fig)

    # Figure D: pattern contributions to the marginal Q1 and time mean.
    rows = [r for r in decomp if r["M"] == "40" and r["rho"] == "0.5"]
    patterns = [f"{i:04b}" for i in range(16)]
    q1 = [sum(float(r["marginal_contribution"]) for r in rows if r["pattern"] == p and r["quantity"] == "quartet_probability" and r["topology"] == "12|34") for p in patterns]
    tm = [sum(float(r["marginal_contribution"]) for r in rows if r["pattern"] == p and r["quantity"] == "mean_pairwise_time") for p in patterns]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].bar(patterns, q1, color=[colors.get(next(r["pattern_category"] for r in rows if r["pattern"] == p), "#555") for p in patterns]); axes[0].set_ylabel("contribution to P(Q1)"); axes[0].tick_params(axis="x", rotation=70)
    axes[1].bar(patterns, tm, color=[colors.get(next(r["pattern_category"] for r in rows if r["pattern"] == p), "#555") for p in patterns]); axes[1].set_ylabel("contribution to E[mean pairwise time]"); axes[1].tick_params(axis="x", rotation=70)
    fig.suptitle("Natural pattern-weighted marginalization: M=40, rho=0.5"); fig.tight_layout(); save(fig, "marginalization_decomposition"); plt.close(fig)


if __name__ == "__main__":
    main()
