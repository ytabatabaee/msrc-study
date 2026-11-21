#!/usr/bin/env python3
"""Validate exact pair expectations against the authoritative msrc-sim API."""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pairwise_msrc import expected_pair_times, history_segments  # noqa: E402


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def bootstrap_ci(x: np.ndarray, rng: np.random.Generator, reps: int = 4000) -> tuple[float, float]:
    if len(x) < 2:
        return float(x.mean()), float(x.mean())
    draws = rng.integers(0, len(x), size=(reps, len(x)))
    means = x[draws].mean(axis=1)
    return tuple(np.quantile(means, [0.025, 0.975]).astype(float))


def run(msrc_path: Path, output: Path, histories_per_cell: int, seed: int,
        M: int, tau_coal: float, ancestral_ne: int) -> None:
    if str(msrc_path) not in sys.path:
        sys.path.insert(0, str(msrc_path))
    from msrcsim.moran import simulate_frequency_history as simulate_moran
    from msrcsim.wright_fisher import simulate_frequency_history as simulate_wf
    from msrcsim.rearrangement import Rearrangement
    from msrcsim.species_tree import SpeciesTree
    from msrcsim.structured_coalescent import simulate_genealogy

    ne = M // 2
    tau = tau_coal * M
    base_seed = np.random.SeedSequence(seed)
    children = iter(base_seed.spawn(100_000))
    rng_boot = np.random.default_rng(seed + 11)
    rows: list[dict] = []
    summaries: list[dict] = []
    pair_states = {"AA": (0, 0), "AB": (0, 1), "BB": (1, 1)}
    cells = [(model, p0, selection) for model in ("moran", "wright_fisher")
             for p0 in (0.25, 0.5, 0.75) for selection in (0.0, 0.1)]

    for cell_i, (model, p0, selection) in enumerate(cells):
        tree = SpeciesTree(
            f"((1:0,2:0)NODE:{tau})ROOT;", ne, 1.0,
            {"NODE": {"effective_population_size": ne, "selection_coefficient": selection},
             "ROOT": {"effective_population_size": ancestral_ne}},
        )
        initial_count = max(1, min(M - 1, int(round(M * p0))))
        rearrangement = Rearrangement(
            f"pair-{cell_i}", "inversion", "NODE", 0, initial_count, selection
        )
        history_generator = simulate_moran if model == "moran" else simulate_wf
        cell_rows = []
        # rho is in the legacy simulator's backward switching-rate convention.
        # gamma=2Ne*rho is dimensionless; one value tests histories rather than
        # selecting a biological p from the empirical arrangement sample counts.
        gamma = 1.0
        rho = gamma / (2.0 * ne)
        for history_id in range(histories_per_cell):
            hseed = next(children).generate_state(1, dtype=np.uint32)[0]
            hrng = np.random.default_rng(int(hseed))
            history = history_generator(tree, rearrangement, hrng)
            segments = history_segments(history, "NODE", ne, rho, tau)
            exact = expected_pair_times(segments, ancestral_ne=ancestral_ne)
            terminal_p = float(history.frequency_at("NODE", 0.0))
            end_status = next((r.status for r in history.by_branch["NODE"] if r.absolute_age <= 1e-12), "unknown")
            row = {
                "model": model, "initial_frequency_B": p0, "selection_coefficient": selection,
                "M_chromosome_copies": M, "Ne": ne, "tau_generations": tau,
                "tau_over_2Ne": tau_coal, "gamma_2Ne_rho": gamma, "rho": rho,
                "history_id": history_id, "terminal_frequency_B": terminal_p,
                "terminal_status": end_status,
            }
            grng = np.random.default_rng(int(next(children).generate_state(1, dtype=np.uint32)[0]))
            for name, states in pair_states.items():
                row[f"analytic_T_{name}"] = exact[name]
                sample = dict(zip(tree.taxa, states))
                genealogy = simulate_genealogy(
                    history_id, tree, history, sample, rho, 1.0, grng, record_events=False
                )
                row[f"simulated_T_{name}"] = float(genealogy.coalescence_times[0])
            rows.append(row)
            cell_rows.append(row)

        arr = cell_rows
        for name in pair_states:
            analytic = np.asarray([r[f"analytic_T_{name}"] for r in arr], float)
            simulated = np.asarray([r[f"simulated_T_{name}"] for r in arr], float)
            lo, hi = bootstrap_ci(analytic, rng_boot)
            diff = simulated - analytic
            diff_lo, diff_hi = bootstrap_ci(diff, rng_boot)
            summaries.append({
                "model": model, "initial_frequency_B": p0, "selection_coefficient": selection,
                "M_chromosome_copies": M, "Ne": ne, "tau_over_2Ne": tau_coal,
                "gamma_2Ne_rho": gamma, "pair_state": name,
                "n_independent_frequency_histories": len(arr),
                "mean_analytic_conditional_T": analytic.mean(),
                "history_bootstrap_95_low": lo, "history_bootstrap_95_high": hi,
                "mean_simulated_genealogy_T": simulated.mean(),
                "simulated_genealogy_cluster_se": simulated.std(ddof=1) / np.sqrt(len(simulated)),
                "mean_sim_minus_analytic": diff.mean(),
                "paired_history_bootstrap_diff_95_low": diff_lo,
                "paired_history_bootstrap_diff_95_high": diff_hi,
                "fixed_or_lost_frequency_histories": sum(r["terminal_status"] in {"fixed", "lost"} for r in arr),
                "uncertainty_unit": "independent forward frequency history; one genealogy per state/history",
            })

    write_tsv(output / "frequency_history_validation.tsv", rows)
    write_tsv(output / "frequency_history_validation_summary.tsv", summaries)
    manifest = {
        "purpose": "history-conditional pair coalescence expectations vs direct legacy genealogy events",
        "msrc_sim_path": str(msrc_path),
        "msrc_sim_commit": "1658db065a3b42a2fc1d0bebda874eecc4da6d3b",
        "models": ["continuous-time Moran", "discrete-generation Wright-Fisher"],
        "histories_per_cell": histories_per_cell,
        "independent_uncertainty_unit": "forward frequency history, not genealogy draws sharing a history",
        "genealogies_per_history": 3,
        "cells": len(cells), "seed": seed,
        "initial_frequency_B": [0.25, 0.5, 0.75],
        "selection_coefficients": [0.0, 0.1],
        "rho_mapping": "gamma=2Ne*rho=1; this is a simulation grid value, not inferred from cM/Mb",
        "structured_origin_age": tau,
        "ordinary_ancestor_Ne": ancestral_ne,
        "fixed_loss_conditioning": "AA/AB/BB genealogies are explicitly conditioned even if that terminal state is absent; simulator finite-copy floor retained",
    }
    (output / "frequency_history_validation_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"validated {len(rows)} history/state rows across {len(cells)} frequency-history cells")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--msrc-sim", type=Path, default=ROOT.parent / "msrc-sim")
    p.add_argument("--output", type=Path, default=ROOT / "simulations/pairwise_msrc/results")
    p.add_argument("--histories", type=int, default=200)
    p.add_argument("--seed", type=int, default=20261121)
    p.add_argument("--M", type=int, default=20)
    p.add_argument("--tau-coal", type=float, default=1.0)
    p.add_argument("--ancestral-Ne", type=int, default=10)
    a = p.parse_args()
    run(a.msrc_sim, a.output, a.histories, a.seed, a.M, a.tau_coal, a.ancestral_Ne)


if __name__ == "__main__":
    main()
