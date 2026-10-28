"""Small conditional MSRC Patterson's D pilot using the existing msrc-sim.

This fixes a predeclared terminal arrangement pattern by rejection sampling
history plus terminal sampling, then holds each accepted history fixed while
simulating independent genealogies. It does not select histories by D.
"""

from __future__ import annotations

import argparse
import csv
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

from .conditional_patterson_d import (
    bootstrap_summary,
    ratio_d,
    simulate_infinite_sites_counts,
    summarize_genealogies,
)


TREE_NEWICK = "(((P1:0,P2:0)N12:1,P3:1)N123:1,O:2)ROOT;"
TARGET_PATTERN = "0101"


def _rss_mib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / (1024.0 * 1024.0) if sys.platform == "darwin" else value / 1024.0


def _load_msrc(msrc_sim: Path):
    sys.path.insert(0, str(msrc_sim))
    from msrcsim.moran import simulate_frequency_history
    from msrcsim.rearrangement import Rearrangement
    from msrcsim.species_tree import SpeciesTree
    from msrcsim.structured_coalescent import simulate_genealogy, simulate_msc_genealogy
    return SpeciesTree, Rearrangement, simulate_frequency_history, simulate_genealogy, simulate_msc_genealogy


def _accept_history(tree, history, rng, target: str):
    sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
    pattern = "".join(str(sampled[taxon]) for taxon in ("P1", "P2", "P3", "O"))
    return sampled if pattern == target else None


def run_pilot(msrc_sim: Path, output: Path, *, histories: int = 3, genealogies: int = 200,
              seed: int = 20261028, M: int = 20, rho: float = 0.2):
    SpeciesTree, Rearrangement, simulate_history, simulate_genealogy, simulate_msc = _load_msrc(msrc_sim)
    output.mkdir(parents=True, exist_ok=True)
    tree = SpeciesTree(TREE_NEWICK, M // 2, 2.0)
    rearrangement = Rearrangement("quartet", "inversion", "ROOT", 0, M // 2, 0.0)
    rng = np.random.default_rng(seed)
    history_rows = []
    gene_rows = []
    accepted = 0
    attempts = 0
    t0 = time.perf_counter()
    while accepted < histories:
        attempts += 1
        history = simulate_history(tree, rearrangement, rng)
        sampled = _accept_history(tree, history, rng, TARGET_PATTERN)
        if sampled is None:
            continue
        accepted += 1
        attempts_for_history = attempts
        results = [simulate_genealogy(accepted, tree, history, sampled, rho, 1.0, rng, record_events=True)
                   for _ in range(genealogies)]
        summary = summarize_genealogies(results)
        abba, baba = simulate_infinite_sites_counts(results, 100, 1.0e-3, rng)
        uncertainty = bootstrap_summary(results, 200, rng)
        history_rows.append({
            "condition": "conditional_msrc", "history": accepted, "target_pattern": TARGET_PATTERN,
            "history_attempts": attempts_for_history, "history_acceptance_fraction": 1.0 / attempts_for_history,
            "M": M, "rho": rho, "genealogies": genealogies,
            "q1": summary.topology_frequencies[0], "q2": summary.topology_frequencies[1],
            "q3": summary.topology_frequencies[2], "mean_b_q1": summary.mean_internal_branch_by_topology[0],
            "mean_b_q2": summary.mean_internal_branch_by_topology[1], "mean_b_q3": summary.mean_internal_branch_by_topology[2],
            "D_patterson_branch": summary.d_patterson, "D_topology_q2_minus_q3": summary.d_topology,
            "D_branch_q2_minus_q3": summary.d_branch, "D_sites": ratio_d(abba, baba),
            **uncertainty,
            "ABBA_branch_sum": summary.abba_sum, "BABA_branch_sum": summary.baba_sum,
            "ABBA_sites": abba, "BABA_sites": baba, "P_Z_positive": True,
        })
        attempts = 0

    # Independent standard MSC baseline on the same rooted tree.
    msc_results = [simulate_msc(0, tree, rng, record_events=True) for _ in range(genealogies * histories)]
    msc_summary = summarize_genealogies(msc_results)
    msc_abba, msc_baba = simulate_infinite_sites_counts(msc_results, 100, 1.0e-3, rng)
    msc_uncertainty = bootstrap_summary(msc_results, 200, rng)
    history_rows.append({
        "condition": "standard_msc", "history": 0, "target_pattern": "",
        "history_attempts": 0, "history_acceptance_fraction": np.nan, "M": M, "rho": 0.0,
        "genealogies": len(msc_results), "q1": msc_summary.topology_frequencies[0],
        "q2": msc_summary.topology_frequencies[1], "q3": msc_summary.topology_frequencies[2],
        "mean_b_q1": msc_summary.mean_internal_branch_by_topology[0],
        "mean_b_q2": msc_summary.mean_internal_branch_by_topology[1],
        "mean_b_q3": msc_summary.mean_internal_branch_by_topology[2],
        "D_patterson_branch": msc_summary.d_patterson,
        "D_topology_q2_minus_q3": msc_summary.d_topology,
        "D_branch_q2_minus_q3": msc_summary.d_branch,
        "D_sites": ratio_d(msc_abba, msc_baba), "ABBA_branch_sum": msc_summary.abba_sum,
        "BABA_branch_sum": msc_summary.baba_sum, "ABBA_sites": msc_abba,
        "BABA_sites": msc_baba, "P_Z_positive": True,
        **msc_uncertainty,
    })
    elapsed = time.perf_counter() - t0
    fields = list(history_rows[0])
    with (output / "pilot_results.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(history_rows)
    (output / "pilot_manifest.json").write_text(json.dumps({
        "tree_newick": TREE_NEWICK, "target_pattern": TARGET_PATTERN, "M": M,
        "rho": rho, "histories": histories, "genealogies_per_history": genealogies,
        "seed": seed, "elapsed_seconds": elapsed, "peak_rss_mib": _rss_mib(),
        "msrc_sim": str(msrc_sim), "model": "conditional_history_rejection_then_fixed_history_genealogies",
    }, indent=2) + "\n")
    return history_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--msrc-sim", type=Path, default=Path("/Users/ytabatabaee/Desktop/msrc-sim"))
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results" / "conditional_patterson_d")
    parser.add_argument("--histories", type=int, default=3)
    parser.add_argument("--genealogies", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20261028)
    parser.add_argument("--M", type=int, default=20)
    parser.add_argument("--rho", type=float, default=0.2)
    args = parser.parse_args()
    print(run_pilot(args.msrc_sim, args.output, histories=args.histories,
                    genealogies=args.genealogies, seed=args.seed, M=args.M, rho=args.rho))


if __name__ == "__main__":  # pragma: no cover
    main()
