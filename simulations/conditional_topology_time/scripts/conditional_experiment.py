#!/usr/bin/env python3
"""Conditional topology/time summaries using the authoritative MSRC simulator."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MSRC = ROOT.parent / "msrc-sim"
SIGNAL_SCRIPT = ROOT / "simulations/topology_time_separation/scripts"
if str(SIGNAL_SCRIPT) not in sys.path:
    sys.path.insert(0, str(SIGNAL_SCRIPT))
from signal_separation import PAIRS, TOPOLOGY_NAMES, _msc_commit, parse_newick, topology_from_split, tree_metrics  # noqa: E402

PATTERNS = tuple(f"{i:04b}" for i in range(16))
CONCORDANT_22 = {"0011", "1100"}
DISCORDANT_22 = {"0101", "1010", "0110", "1001"}
UNIFORM = {"0000", "1111"}
THREE_ONE = {p for p in PATTERNS if p.count("1") in {1, 3}}
TIME_STATISTICS = (
    "T12", "T13", "T14", "T23", "T24", "T34", "mean_pairwise_time",
    "internal_branch_length", "total_tree_length", "R_T", "root_coalescence_time",
)


def category(pattern: str) -> str:
    if pattern in CONCORDANT_22:
        return "concordant_2_2"
    if pattern in DISCORDANT_22:
        return "discordant_2_2"
    if pattern in UNIFORM:
        return "uniform"
    if pattern in THREE_ONE:
        return "3_1"
    raise ValueError(f"unknown pattern {pattern}")


def bootstrap_mean(values: np.ndarray, rng: np.random.Generator, n_boot: int = 400) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return math.nan, math.nan
    if len(values) == 1:
        return float(values[0]), float(values[0])
    draws = rng.integers(0, len(values), size=(n_boot, len(values)))
    means = values[draws].mean(axis=1)
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0])
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def imports(msrc_path: Path):
    if str(msrc_path) not in sys.path:
        sys.path.insert(0, str(msrc_path))
    from msrcsim.moran import simulate_frequency_history
    from msrcsim.rearrangement import Rearrangement
    from msrcsim.species_tree import SpeciesTree
    from msrcsim.structured_coalescent import simulate_genealogy, simulate_msc_genealogy
    return simulate_frequency_history, Rearrangement, SpeciesTree, simulate_genealogy, simulate_msc_genealogy


def simulate_cell(msrc_path: Path, cell: dict[str, object], histories: int, seed: int) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    simulate_frequency_history, Rearrangement, SpeciesTree, simulate_genealogy, simulate_msc_genealogy = imports(msrc_path)
    M = int(cell["M"])
    tree = SpeciesTree(
        f"((1:0,2:0)AB:{cell['sister_length']},(3:0,4:0)CD:{cell['sister_length']})ROOT;",
        M // 2,
        float(cell["root_length"]),
    )
    rearrangement = Rearrangement("conditional", "inversion", "ROOT", 0, int(round(M * float(cell["initial_fraction"]))), 0.0)
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    history_rows: list[dict[str, object]] = []
    for history_id in range(histories):
        history = simulate_frequency_history(tree, rearrangement, rng)
        sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
        pattern = "".join(str(sampled[taxon]) for taxon in tree.taxa)
        history_rows.append({"history_id": history_id, "pattern": pattern, "terminal_frequency": history.terminal_frequency("1")})
        result = simulate_genealogy(history_id, tree, history, sampled, float(cell["rho"]), 1.0, rng, record_events=False)
        metrics = tree_metrics(result.newick)
        row = {k: v for k, v in metrics.items() if k not in {"pair_distances", "newick"}}
        row.update({
            "model": "msrc", "history_id": history_id, "genealogy_id": 0,
            "pattern": pattern, "pattern_category": category(pattern),
            "topology_index": int(result.topology_index), "topology": topology_from_split(result.topology_index),
            "root_coalescence_time": max(result.coalescence_times),
        })
        rows.append({**cell, **row})
    # MSC is an unconditional comparator only. It has no intrinsic Z variable.
    for history_id in range(histories):
        result = simulate_msc_genealogy(history_id, tree, rng, record_events=False)
        metrics = tree_metrics(result.newick)
        row = {k: v for k, v in metrics.items() if k not in {"pair_distances", "newick"}}
        row.update({
            "model": "msc", "history_id": history_id, "genealogy_id": 0,
            "pattern": "", "pattern_category": "undefined_msc_Z",
            "topology_index": int(result.topology_index), "topology": topology_from_split(result.topology_index),
            "root_coalescence_time": max(result.coalescence_times),
        })
        rows.append({**cell, **row})
    return rows, history_rows


def theory_cell(msrc_path: Path, cell: dict[str, object]) -> dict[str, object]:
    # Import lazily so tests of classification and aggregation do not require
    # SciPy; actual theory comparisons use the configured analysis environment.
    theory_path = ROOT / "simulations/quartet/moran_theory"
    if str(theory_path) not in sys.path:
        sys.path.insert(0, str(theory_path))
    from full_species_tree import full_species_tree_quartet_probabilities
    M = int(cell["M"])
    pi = np.zeros(M + 1)
    pi[int(round(M * float(cell["initial_fraction"])))] = 1.0
    result = full_species_tree_quartet_probabilities(
        M_root=M, M_left=M, M_right=M, M_terminal=(M, M, M, M),
        rho=float(cell["rho"]), root_length=float(cell["root_length"]),
        sister_lengths=(float(cell["sister_length"]), float(cell["sister_length"])),
        terminal_lengths=(0.0, 0.0, 0.0, 0.0), pi_old=pi,
    )
    return {"joint": np.asarray(result.joint), "pattern": np.asarray(result.pattern_probabilities), "conditional": np.asarray(result.conditional_quartet_probabilities), "q": np.asarray(result.quartet_probabilities)}


def summarize_patterns(rows: list[dict[str, object]], cell: dict[str, object], theory: dict[str, object], seed: int, min_report_n: int = 20) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    msrc = [r for r in rows if r["model"] == "msrc"]
    rng = np.random.default_rng(seed)
    by_pattern = {p: [r for r in msrc if r["pattern"] == p] for p in PATTERNS}
    conditional_quartets = []
    conditional_times = []
    topology_times = []
    pattern_probs = []
    for zi, pattern in enumerate(PATTERNS):
        sub = by_pattern[pattern]
        n = len(sub)
        empirical_p = n / len(msrc)
        theoretical_p = float(theory["pattern"][zi])
        pattern_probs.append({**cell, "pattern": pattern, "pattern_category": category(pattern), "n_histories": n, "empirical_probability": empirical_p, "empirical_mc_se": math.sqrt(empirical_p * (1 - empirical_p) / len(msrc)), "theory_probability": theoretical_p, "absolute_error": abs(empirical_p - theoretical_p), "adequate_for_conditional_summary": str(n >= min_report_n).lower()})
        counts = np.bincount([int(r["topology_index"]) for r in sub], minlength=3) if sub else np.zeros(3, dtype=int)
        q = counts / n if n else np.full(3, math.nan)
        theory_q = theory["conditional"][zi]
        for qi, topology in enumerate(TOPOLOGY_NAMES):
            conditional_quartets.append({**cell, "pattern": pattern, "pattern_category": category(pattern), "topology": topology, "empirical_probability": float(q[qi]), "empirical_mc_se": math.sqrt(float(q[qi]) * (1 - float(q[qi])) / n) if n else math.nan, "theory_probability": float(theory_q[qi]), "absolute_error": abs(float(q[qi]) - float(theory_q[qi])) if n else math.nan, "n_histories": n, "adequate_for_conditional_summary": str(n >= min_report_n).lower()})
        for stat in TIME_STATISTICS:
            vals = np.asarray([float(r[stat if stat != "internal_branch_length" else f"internal_branch_length_q{int(r['topology_index'])}"]) for r in sub], dtype=float)
            lo, hi = bootstrap_mean(vals, rng)
            conditional_times.append({**cell, "pattern": pattern, "pattern_category": category(pattern), "statistic": stat, "mean": float(np.nanmean(vals)) if len(vals) else math.nan, "sd": float(np.nanstd(vals, ddof=1)) if len(vals) > 1 else math.nan, "q05": float(np.nanquantile(vals, .05)) if len(vals) else math.nan, "q50": float(np.nanquantile(vals, .50)) if len(vals) else math.nan, "q95": float(np.nanquantile(vals, .95)) if len(vals) else math.nan, "bootstrap_lo": lo, "bootstrap_hi": hi, "n_histories": n, "adequate_for_conditional_summary": str(n >= min_report_n).lower()})
        for qi, topology in enumerate(TOPOLOGY_NAMES):
            subq = [r for r in sub if int(r["topology_index"]) == qi]
            for stat in ("mean_pairwise_time", "internal_branch_length", "total_tree_length", "R_T", "root_coalescence_time"):
                key = stat if stat != "internal_branch_length" else f"internal_branch_length_q{qi}"
                vals = np.asarray([float(r[key]) for r in subq], dtype=float)
                lo, hi = bootstrap_mean(vals, rng)
                topology_times.append({**cell, "pattern": pattern, "pattern_category": category(pattern), "topology": topology, "statistic": stat, "mean": float(np.nanmean(vals)) if len(vals) else math.nan, "sd": float(np.nanstd(vals, ddof=1)) if len(vals) > 1 else math.nan, "q05": float(np.nanquantile(vals, .05)) if len(vals) else math.nan, "q50": float(np.nanquantile(vals, .50)) if len(vals) else math.nan, "q95": float(np.nanquantile(vals, .95)) if len(vals) else math.nan, "bootstrap_lo": lo, "bootstrap_hi": hi, "n_histories": len(subq), "adequate_for_conditional_summary": str(len(subq) >= min_report_n).lower()})
    return pattern_probs, conditional_quartets, conditional_times, topology_times


def marginal_decomposition(rows: list[dict[str, object]], cell: dict[str, object], theory: dict[str, object]) -> list[dict[str, object]]:
    msrc = [r for r in rows if r["model"] == "msrc"]
    out = []
    for zi, pattern in enumerate(PATTERNS):
        sub = [r for r in msrc if r["pattern"] == pattern]
        p = len(sub) / len(msrc)
        for qi, topology in enumerate(TOPOLOGY_NAMES):
            qz = sum(int(r["topology_index"]) == qi for r in sub) / len(sub) if sub else math.nan
            out.append({**cell, "pattern": pattern, "pattern_category": category(pattern), "quantity": "quartet_probability", "topology": topology, "pattern_probability": p, "conditional_value": qz, "marginal_contribution": p * qz if math.isfinite(qz) else math.nan, "theory_pattern_probability": float(theory["pattern"][zi]), "theory_conditional_value": float(theory["conditional"][zi, qi]), "theory_marginal_contribution": float(theory["joint"][zi, qi])})
        vals = np.asarray([float(r["mean_pairwise_time"]) for r in sub], dtype=float)
        mean_value = float(np.mean(vals)) if len(vals) else math.nan
        out.append({**cell, "pattern": pattern, "pattern_category": category(pattern), "quantity": "mean_pairwise_time", "topology": "all", "pattern_probability": p, "conditional_value": mean_value, "marginal_contribution": p * mean_value if math.isfinite(mean_value) else math.nan, "theory_pattern_probability": float(theory["pattern"][zi]), "theory_conditional_value": math.nan, "theory_marginal_contribution": math.nan})
    return out


def marginal_summary(rows: list[dict[str, object]], cell: dict[str, object]) -> dict[str, object]:
    out = dict(cell)
    for model in ("msrc", "msc"):
        sub = [r for r in rows if r["model"] == model]
        counts = np.bincount([int(r["topology_index"]) for r in sub], minlength=3) / len(sub)
        out.update({f"{model}_q{qi + 1}": float(counts[qi]) for qi in range(3)})
        values = np.asarray([float(r["mean_pairwise_time"]) for r in sub])
        out[f"{model}_mean_pairwise_time"] = float(np.mean(values))
        out[f"{model}_mean_pairwise_time_sd"] = float(np.std(values, ddof=1))
        out[f"{model}_mean_pairwise_time_mc_se"] = float(np.std(values, ddof=1) / math.sqrt(len(values)))
        out[f"{model}_n"] = len(sub)
    out["delta_Q_msrc_minus_msc"] = float(sum(abs(out[f"msrc_q{i}"] - out[f"msc_q{i}"]) for i in (1, 2, 3)) / 2)
    out["delta_mean_pairwise_time_msrc_minus_msc"] = out["msrc_mean_pairwise_time"] - out["msc_mean_pairwise_time"]
    out["delta_mean_pairwise_time_standardized"] = abs(out["delta_mean_pairwise_time_msrc_minus_msc"]) / out["msc_mean_pairwise_time_sd"] if out["msc_mean_pairwise_time_sd"] else math.nan
    # The following reconstructions use natural MSRC pattern frequencies and
    # are an internal identity check, not a second estimator.
    msrc = [r for r in rows if r["model"] == "msrc"]
    q_recon = [0.0, 0.0, 0.0]
    t_recon = 0.0
    for pattern in PATTERNS:
        sub = [r for r in msrc if r["pattern"] == pattern]
        p = len(sub) / len(msrc)
        for qi in range(3): q_recon[qi] += p * sum(int(r["topology_index"]) == qi for r in sub) / len(sub) if sub else 0.0
        t_recon += p * np.mean([float(r["mean_pairwise_time"]) for r in sub]) if sub else 0.0
    out.update({f"reconstructed_msrc_q{i + 1}": q_recon[i] for i in range(3)})
    out["reconstructed_msrc_mean_pairwise_time"] = float(t_recon)
    out["reconstruction_max_abs_error"] = max([abs(out[f"reconstructed_msrc_q{i}"] - out[f"msrc_q{i}"]) for i in (1, 2, 3)] + [abs(out["reconstructed_msrc_mean_pairwise_time"] - out["msrc_mean_pairwise_time"])])
    return out


def run_stage(output: Path, msrc_path: Path, grid: list[dict[str, object]], histories: int, seed: int, stage: str) -> None:
    output.mkdir(parents=True, exist_ok=True)
    all_pattern, all_cond_q, all_cond_t, all_top_t, all_decomp, all_theory, all_marginal = [], [], [], [], [], [], []
    for index, cell in enumerate(grid):
        rows, history_rows = simulate_cell(msrc_path, cell, histories, seed + index * 1009)
        theory = theory_cell(msrc_path, cell)
        p, cq, ct, tt = summarize_patterns(rows, cell, theory, seed + index * 17)
        all_pattern.extend(p); all_cond_q.extend(cq); all_cond_t.extend(ct); all_top_t.extend(tt); all_decomp.extend(marginal_decomposition(rows, cell, theory))
        all_marginal.append(marginal_summary(rows, cell))
        for zi, pattern in enumerate(PATTERNS):
            sub = [r for r in rows if r["model"] == "msrc" and r["pattern"] == pattern]
            for qi, topology in enumerate(TOPOLOGY_NAMES):
                empirical = sum(int(r["topology_index"]) == qi for r in sub) / len(sub) if sub else math.nan
                joint = sum(int(r["topology_index"]) == qi for r in sub) / histories
                all_theory.append({**cell, "pattern": pattern, "topology": topology, "n_histories": len(sub), "empirical_joint": joint, "empirical_joint_mc_se": math.sqrt(joint * (1 - joint) / histories), "empirical_pattern": len(sub) / histories, "empirical_conditional": empirical, "empirical_conditional_mc_se": math.sqrt(empirical * (1 - empirical) / len(sub)) if sub else math.nan, "theory_joint": float(theory["joint"][zi, qi]), "theory_pattern": float(theory["pattern"][zi]), "theory_conditional": float(theory["conditional"][zi, qi]), "joint_abs_error": abs(joint - float(theory["joint"][zi, qi])), "conditional_abs_error": abs(empirical - float(theory["conditional"][zi, qi])) if sub else math.nan})
    write_tsv(output / "pattern_probabilities.tsv", all_pattern)
    write_tsv(output / "conditional_quartets.tsv", all_cond_q)
    write_tsv(output / "conditional_times.tsv", all_cond_t)
    write_tsv(output / "topology_conditioned_times.tsv", all_top_t)
    write_tsv(output / "marginal_decomposition.tsv", all_decomp)
    write_tsv(output / "theory_comparison.tsv", all_theory)
    write_tsv(output / "marginal_summary.tsv", all_marginal)
    (output / "parameter_grid.tsv").write_text("M\tinitial_fraction\trho\tsister_length\troot_length\n" + "\n".join("\t".join(str(cell[k]) for k in ("M", "initial_fraction", "rho", "sister_length", "root_length")) for cell in grid) + "\n")
    manifest = {"stage": stage, "experiment": "conditional topology and coalescence-time signals under MSRC", "msrc_sim_path": str(msrc_path), "msrc_sim_commit": _msc_commit(msrc_path), "histories_per_cell": histories, "genes_per_history": 1, "seed": seed, "independent_unit": "forward Moran frequency history", "grid_cells": len(grid), "theory": "simulations/quartet/moran_theory/full_species_tree.py; P(Q,Z), P(Z), P(Q|Z)", "msc_conditional_Z": "undefined; ordinary MSC is used only as an unconditional comparator", "pattern_order": list(PATTERNS), "conditional_min_report_n": 20}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--msrc-sim", type=Path, default=DEFAULT_MSRC)
    p.add_argument("--output", type=Path, default=ROOT / "simulations/conditional_topology_time/results/discovery")
    p.add_argument("--histories", type=int, default=500)
    p.add_argument("--seed", type=int, default=20261121)
    p.add_argument("--M", nargs="+", type=int, default=[20, 40])
    p.add_argument("--initial-fraction", nargs="+", type=float, default=[0.25, 0.5, 0.75])
    p.add_argument("--rho", nargs="+", type=float, default=[0.0, 0.5])
    p.add_argument("--sister-length", nargs="+", type=float, default=[0.25])
    p.add_argument("--root-length", nargs="+", type=float, default=[0.5])
    p.add_argument("--stage", default="discovery")
    args = p.parse_args()
    grid = [{"M": M, "initial_fraction": f, "rho": rho, "sister_length": sister, "root_length": root} for M in args.M for f in args.initial_fraction for rho in args.rho for sister in args.sister_length for root in args.root_length]
    run_stage(args.output, args.msrc_sim, grid, args.histories, args.seed, args.stage)


if __name__ == "__main__":
    main()
