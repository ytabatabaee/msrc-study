"""Compare the exact evaluator with direct msrc-sim genealogies.

Run with PYTHONPATH pointing at the authoritative msrc-sim checkout, for
example:

  PYTHONPATH=/Users/ytabatabaee/Desktop/msrc-sim \
    .venv/bin/python validate_against_msrc_sim.py --replicates 5000
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from quartet_theory import PATTERNS, TOPOLOGY_NAMES, conditional_quartet_probabilities


class FixedHistory:
    def __init__(self, p: float):
        self.p = p

    def frequency_at(self, branch_id: str, age: float) -> float:
        return self.p

    def next_frequency_boundary(self, branch_id: str, age: float) -> float:
        return float("inf")

    def terminal_frequency(self, taxon: str) -> float:
        return self.p


def simulate_one_pattern(z: str, *, replicates: int, rho: float, p: float, ne: int, root_length: float, seed: int):
    from msrcsim.species_tree import SpeciesTree
    from msrcsim.structured_coalescent import simulate_genealogy

    tree = SpeciesTree(f"((1:0.2,2:0.2):{root_length},(3:0.2,4:0.2):{root_length});", ne, 20.0)
    history = FixedHistory(p)
    rng = np.random.default_rng(seed)
    counts = np.zeros(3, dtype=int)
    sampled = dict(zip(tree.taxa, map(int, z)))
    for locus in range(replicates):
        result = simulate_genealogy(locus, tree, history, sampled, rho, 1.0, rng, record_events=False)
        counts[result.topology_index] += 1
    return counts / counts.sum(), counts


def run(args):
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for rho in args.rho:
        for p in args.p:
            for root_length in args.root_length:
                for z_index, z in enumerate(PATTERNS):
                    theory = conditional_quartet_probabilities(
                        z, sister_length=0.2, root_length_before_mrca=root_length,
                        rho=rho, p=p, ne=args.ne,
                    )
                    empirical, counts = simulate_one_pattern(
                        z, replicates=args.replicates, rho=rho, p=p,
                        ne=args.ne, root_length=root_length,
                        seed=args.seed + z_index + int(1000 * rho) + int(10000 * p) + int(100000 * root_length),
                    )
                    for j, name in enumerate(TOPOLOGY_NAMES):
                        rows.append({
                            "rho": rho, "p": p, "root_length": root_length, "pattern": z, "topology": name,
                            "theory": theory[j], "empirical": empirical[j],
                            "count": counts[j], "replicates": args.replicates,
                        })
    table = out / "theory_vs_msrc_sim.csv"
    with table.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for rho in args.rho:
        sub = [r for r in rows if r["rho"] == rho and r["p"] == args.p[0] and r["root_length"] == args.root_length[0] and r["topology"] == "13|24"]
        ax.scatter([r["pattern"] for r in sub], [r["empirical"] for r in sub], s=12, label=f"simulation rho={rho:g}")
        ax.plot([r["pattern"] for r in sub], [r["theory"] for r in sub], linewidth=1.2, label=f"theory rho={rho:g}")
    ax.set(xlabel="sampled pattern Z", ylabel="P(13|24 | Z, S)")
    ax.tick_params(axis="x", rotation=90)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(out / "theory_vs_msrc_sim.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.5, 4))
    rho_grid = np.logspace(-2, 5, 60)
    msc = np.array([1 - (2 / 3) * np.exp(-args.root_length[0] / args.ne) for _ in rho_grid])
    q1 = [conditional_quartet_probabilities("0101", sister_length=0.2, root_length_before_mrca=args.root_length[0], rho=r, p=0.5, ne=args.ne)[0] for r in rho_grid]
    ax.semilogx(rho_grid, q1, label="exact MSRC")
    ax.axhline(msc[0], color="black", linestyle="--", label="MSC limit")
    ax.set(xlabel=r"switching rate $\rho$", ylabel="P(12|34 | Z=0101, S)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "convergence_to_msc.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=5000)
    parser.add_argument("--rho", type=float, nargs="+", default=[0.1, 1.0, 10.0, 1000.0])
    parser.add_argument("--p", type=float, nargs="+", default=[0.25, 0.5, 0.7])
    parser.add_argument("--root-length", type=float, nargs="+", default=[0.1, 0.3, 0.6])
    parser.add_argument("--ne", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20261019)
    parser.add_argument("--output", type=Path, default=Path("validation_output"))
    run(parser.parse_args())
