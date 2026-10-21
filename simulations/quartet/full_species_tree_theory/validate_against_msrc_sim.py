"""Direct msrc-sim validation for the symmetric species-tree theorem.

The simulator is run with p=1/2, so rho=2m and lambda=1/N. A long root
extension is used for direct comparisons so the simulator's terminal fallback
has negligible probability; finite-root-duration curves are evaluated exactly
by quartet_theory.py.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from quartet_theory import (
    PATTERNS,
    TOPOLOGY_NAMES,
    conditional_quartet_probabilities,
    conditional_scalar,
    topology_reversal_condition,
)


class FixedHistory:
    def frequency_at(self, branch_id: str, age: float) -> float:
        return 0.5

    def next_frequency_boundary(self, branch_id: str, age: float) -> float:
        return float("inf")

    def terminal_frequency(self, taxon: str) -> float:
        return 0.5


def simulate_one_pattern(z: str, *, replicates: int, m: float, lam: float, terminal_length: float, sister_length: float, root_extension: float, seed: int):
    from msrcsim.species_tree import SpeciesTree
    from msrcsim.structured_coalescent import simulate_genealogy

    ne = 1.0 / lam
    tree = SpeciesTree(
        f"((1:{terminal_length},2:{terminal_length}):{sister_length},"
        f"(3:{terminal_length},4:{terminal_length}):{sister_length});",
        int(round(ne)), root_extension,
    )
    rng = np.random.default_rng(seed)
    history = FixedHistory()
    sampled = dict(zip(tree.taxa, map(int, z)))
    counts = np.zeros(3, dtype=int)
    for locus in range(replicates):
        result = simulate_genealogy(locus, tree, history, sampled, 2.0 * m, 1.0, rng, record_events=False)
        counts[result.topology_index] += 1
    return counts / counts.sum(), counts


def exact(z: str, *, m: float, lam: float, terminal_length: float, sister_length: float, shared_length: float):
    return conditional_quartet_probabilities(
        z, terminal_lengths=(terminal_length,) * 4,
        sister_lengths=(sister_length,) * 2,
        shared_length=shared_length, m=m, lam=lam,
    )


def run(args):
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for m in args.m:
        for sister_length in args.sister_length:
            for z_index, z in enumerate(PATTERNS):
                theory = exact(z, m=m, lam=args.lam, terminal_length=args.terminal_length,
                               sister_length=sister_length, shared_length=args.root_duration)
                empirical, counts = simulate_one_pattern(
                    z, replicates=args.replicates, m=m, lam=args.lam,
                    terminal_length=args.terminal_length, sister_length=sister_length,
                    root_extension=args.root_extension,
                    seed=args.seed + z_index + int(1000 * m) + int(10000 * sister_length),
                )
                for j, name in enumerate(TOPOLOGY_NAMES):
                    rows.append({
                        "m": m, "lambda": args.lam, "sister_length": sister_length,
                        "root_duration": args.root_duration, "pattern": z,
                        "topology": name, "theory": theory[j], "empirical": empirical[j],
                        "count": counts[j], "replicates": args.replicates,
                    })
    with (out / "theory_vs_msrc_sim.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    representative = {"4:0": "0000", "3:1": "0001", "2:2 W2": "0101"}
    m_grid = np.logspace(-4, 3, 100)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True)
    for ax, (label, z) in zip(axes, representative.items()):
        curves = np.array([exact(z, m=m, lam=args.lam, terminal_length=args.terminal_length,
                                sister_length=args.sister_length[0], shared_length=args.root_duration) for m in m_grid])
        for j, name in enumerate(TOPOLOGY_NAMES):
            ax.semilogx(m_grid, curves[:, j], label=name)
        ax.set_title(f"{label}: Z={z}")
        ax.set_xlabel("m")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel("quartet probability")
    axes[-1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "quartet_probabilities_vs_m.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.5, 4))
    for sister_length in args.sister_length:
        q = [exact("0101", m=m, lam=args.lam, terminal_length=args.terminal_length,
                   sister_length=sister_length, shared_length=args.root_duration)[0] for m in m_grid]
        ax.semilogx(m_grid, q, label=f"T={sister_length:g}")
    tau = args.lam * sum(args.sister_length) / 2.0
    ax.axhline(1 - 2 * np.exp(-tau) / 3, color="black", linestyle="--", label="MSC limit")
    ax.set(xlabel="m", ylabel="P(12|34 | Z=0101, S)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "msc_limit_vs_m.png", dpi=180)
    plt.close(fig)

    errors = []
    for m in m_grid:
        errors.append(max(
            np.max(np.abs(exact(z, m=m, lam=args.lam, terminal_length=args.terminal_length,
                                sister_length=args.sister_length[0], shared_length=args.root_duration)
                          - conditional_scalar(z, args.sister_length[0], m, args.lam, args.root_duration)))
            for z in PATTERNS
        ))
    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.loglog(m_grid, np.maximum(errors, 1e-18))
    ax.set(xlabel="m", ylabel="max |matrix - scalar|", title="Scalar/matrix agreement")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "scalar_vs_matrix_validation.png", dpi=180)
    plt.close(fig)

    m_phase = np.logspace(-4, 2, 160)
    T_phase = np.linspace(0.0, 3.0, 160)
    reversal = np.array([
        [topology_reversal_condition(T, m, args.lam, args.root_duration) for m in m_phase]
        for T in T_phase
    ])
    fig, ax = plt.subplots(figsize=(5.5, 4))
    image = ax.pcolormesh(m_phase, T_phase, reversal, shading="auto", cmap="RdYlBu_r", vmin=0, vmax=1)
    ax.set_xscale("log")
    ax.set(xlabel="m", ylabel="sister branch length T", title="Conditional topology-reversal region")
    fig.colorbar(image, ax=ax, label="q2 > q1")
    fig.tight_layout()
    fig.savefig(out / "conditional_topology_reversal_region.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--m", type=float, nargs="+", default=[0.0, 0.01, 0.1, 1.0, 10.0, 100.0])
    parser.add_argument("--sister-length", type=float, nargs="+", default=[0.1, 0.3, 0.6])
    parser.add_argument("--terminal-length", type=float, default=0.0)
    parser.add_argument("--root-duration", type=float, default=20.0)
    parser.add_argument("--root-extension", type=float, default=20.0)
    parser.add_argument("--lam", type=float, default=0.1)
    parser.add_argument("--seed", type=int, default=20261019)
    parser.add_argument("--output", type=Path, default=Path("validation_output"))
    run(parser.parse_args())
