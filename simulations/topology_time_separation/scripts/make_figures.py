#!/usr/bin/env python3
"""Create the four topology/time-separation figures from saved TSV outputs."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "simulations/topology_time_separation"


def read(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def save(fig, name: str) -> None:
    out = BASE / "figures"
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(out / f"{name}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def figure_a() -> None:
    rows = read(BASE / "results/discovery/summary_metrics.tsv")
    fig, ax = plt.subplots(figsize=(7.2, 5.4), constrained_layout=True)
    rhos = sorted({f(r, "rho") for r in rows})
    markers = {0.0: "o", 0.1: "s", 0.5: "^"}
    for rho in rhos:
        z = [r for r in rows if f(r, "rho") == rho]
        ax.scatter([f(r, "delta_Q") for r in z], [f(r, "delta_T_mean_pairwise") for r in z], s=28, marker=markers.get(rho, "o"), label=f"rho={rho:g}", alpha=.72)
    ax.axvline(0, color="black", lw=.6); ax.axhline(0, color="black", lw=.6)
    ax.set_xlabel("topology separation, ΔQ (total variation)")
    ax.set_ylabel("standardized mean-pairwise-time separation, ΔT")
    ax.set_title("Discovery grid: topology versus genealogical-time separation")
    ax.legend(frameon=False, title="effective switching intensity")
    save(fig, "topology_vs_time_separation")


def figure_b() -> None:
    rows = read(BASE / "results/validation_long/time_metrics.tsv")
    stats = ["mean_pairwise_time", "internal_branch_length", "R_T"]
    labels = ["mean pairwise time", "internal quartet branch", "R_T"]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), constrained_layout=True)
    for ax, stat, label in zip(axes, stats, labels):
        z = [r for r in rows if r["statistic"] == stat]
        positions = []
        for i, r in enumerate(z):
            x = i + (0 if r["model"] == "msc" else .34)
            ax.errorbar(x, f(r, "mean"), yerr=[[f(r, "mean") - f(r, "q05")], [f(r, "q95") - f(r, "mean")]], fmt="o", ms=4, capsize=2, color="#4c78a8" if r["model"] == "msc" else "#e45756")
            positions.append(x)
        ax.set_title(label)
        ax.set_xticks([])
        ax.grid(axis="y", alpha=.2)
    axes[0].set_ylabel("value; whiskers are 5th–95th percentiles")
    axes[0].plot([], [], "o", color="#4c78a8", label="MSC")
    axes[0].plot([], [], "o", color="#e45756", label="MSRC")
    axes[0].legend(frameon=False)
    fig.suptitle("Validation scenarios: true genealogical-time summaries")
    save(fig, "genealogical_time_distributions")


def figure_c() -> None:
    rows = read(BASE / "results/validation_strong/theory_comparison.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), sharex=True, constrained_layout=True)
    for i, topo in enumerate(("q1", "q2", "q3")):
        x = np.arange(len(rows))
        axes[i].errorbar(x, [f(r, f"analytic_{topo}") for r in rows], fmt="o", label="theory", color="#4c78a8")
        axes[i].scatter(x, [f(r, f"empirical_{topo}") for r in rows], label="simulation", color="#e45756", s=18)
        axes[i].set_title(topo)
        axes[i].set_ylim(0, .5)
        axes[i].grid(axis="y", alpha=.2)
    axes[0].set_ylabel("quartet probability")
    axes[1].set_xlabel("validation cell index")
    axes[0].legend(frameon=False)
    fig.suptitle("Finite-M theory versus direct Moran MSRC simulation")
    save(fig, "quartet_theory_comparison")


def figure_d() -> None:
    rows = read(BASE / "results/validation_strong/rate_robustness.tsv")
    controls = [r["control"] for r in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)
    axes[0].errorbar(x, [f(r, "mean") for r in rows], yerr=[f(r, "sd") for r in rows], fmt="o", color="#4c78a8")
    axes[1].errorbar(x, [f(r, "mean_R_T") for r in rows], yerr=[f(r, "sd_R_T") for r in rows], fmt="o", color="#e45756")
    for ax, ylabel in zip(axes, ("mean scaled pair distance / true time", "relative contrast R_T")):
        ax.set_xticks(x, controls, rotation=35, ha="right")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", alpha=.2)
    fig.suptitle("Mutation-rate controls on true genealogical summaries")
    save(fig, "mutation_rate_sensitivity")


if __name__ == "__main__":
    figure_a(); figure_b(); figure_c(); figure_d()
