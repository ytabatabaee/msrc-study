#!/usr/bin/env python3
"""Topology/time separation experiment using the authoritative msrc-sim API.

The script records true genealogical quantities from simulator Newick output;
it does not infer times from sequence branch lengths and does not implement a
second genealogy simulator.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_MSRC = ROOT.parent / "msrc-sim"
TREE_NEWICK = "((1:0,2:0)AB:{sister},(3:0,4:0)CD:{sister})ROOT;"
TOPOLOGY_NAMES = ("12|34", "13|24", "14|23")
PAIRS = (("1", "2"), ("1", "3"), ("1", "4"), ("2", "3"), ("2", "4"), ("3", "4"))


@dataclass
class Node:
    children: list["Node"]
    name: str = ""
    length: float = 0.0
    parent: "Node | None" = None


def parse_newick(text: str) -> Node:
    text = text.strip().rstrip(";")
    pos = 0

    def node() -> Node:
        nonlocal pos
        if text[pos] == "(":
            pos += 1
            children = [node()]
            while text[pos] == ",":
                pos += 1
                children.append(node())
            if text[pos] != ")":
                raise ValueError("malformed Newick")
            pos += 1
            name = label()
            length = branch_length()
            out = Node(children, name, length)
            for child in children:
                child.parent = out
            return out
        name = label()
        return Node([], name, branch_length())

    def label() -> str:
        nonlocal pos
        start = pos
        while pos < len(text) and text[pos] not in ":,()":
            pos += 1
        return text[start:pos].strip()

    def branch_length() -> float:
        nonlocal pos
        if pos >= len(text) or text[pos] != ":":
            return 0.0
        pos += 1
        start = pos
        while pos < len(text) and text[pos] not in ",()":
            pos += 1
        return float(text[start:pos])

    root = node()
    if pos != len(text):
        raise ValueError(f"unexpected Newick suffix: {text[pos:]}")
    return root


def tree_metrics(newick: str, edge_multipliers: dict[frozenset[str], float] | None = None) -> dict[str, object]:
    root = parse_newick(newick)
    depths: dict[int, float] = {}
    leaves: dict[int, frozenset[str]] = {}
    edges: list[tuple[int, float, frozenset[str]]] = []

    def visit(n: Node, depth: float) -> frozenset[str]:
        depths[id(n)] = depth
        if not n.children:
            leaves[id(n)] = frozenset((n.name,))
            return leaves[id(n)]
        desc = frozenset().union(*(visit(c, depth + c.length) for c in n.children))
        leaves[id(n)] = desc
        for c in n.children:
            edges.append((id(c), c.length, leaves[id(c)]))
        return desc

    visit(root, 0.0)
    root_height = max(depths.values())
    nodes = []
    def collect(n: Node) -> None:
        nodes.append(n)
        for c in n.children:
            collect(c)
    collect(root)

    def mrca(a: str, b: str) -> Node:
        candidates = [n for n in nodes if a in leaves[id(n)] and b in leaves[id(n)]]
        return max(candidates, key=lambda n: depths[id(n)])

    def pair_distance(a: str, b: str) -> float:
        na, nb = next(n for n in nodes if n.name == a), next(n for n in nodes if n.name == b)
        ancestors_a: dict[int, float] = {id(na): 0.0}
        cur = na
        while cur.parent is not None:
            ancestors_a[id(cur.parent)] = ancestors_a[id(cur)] + cur.length * (edge_multipliers or {}).get(leaves[id(cur)], 1.0)
            cur = cur.parent
        cur = nb
        d_b = 0.0
        while id(cur) not in ancestors_a:
            if cur.parent is None:
                raise ValueError("disconnected Newick")
            d_b += cur.length * (edge_multipliers or {}).get(leaves[id(cur)], 1.0)
            cur = cur.parent
        return ancestors_a[id(cur)] + d_b

    pair_times = {f"T{a}{b}": root_height - depths[id(mrca(a, b))] for a, b in PAIRS}
    distances = {(a, b): pair_distance(a, b) for a, b in PAIRS}
    d = distances
    if all(x in d for x in PAIRS):
        q = 0
        # The simulator already supplies topology_index; this is only for
        # the internal-branch metric and uses the standard four-point formula.
        internal = (
            (d[("1", "3")] + d[("2", "4")] - d[("1", "2")] - d[("3", "4")]) / 2.0
        )
        # Keep all three split-specific lengths for validation/debugging.
        internal_by_topology = (
            internal,
            (d[("1", "2")] + d[("3", "4")] - d[("1", "3")] - d[("2", "4")]) / 2.0,
            (d[("1", "2")] + d[("3", "4")] - d[("1", "4")] - d[("2", "3")]) / 2.0,
        )
    else:
        internal_by_topology = (math.nan, math.nan, math.nan)
    mean_pair = float(np.mean(list(pair_times.values())))
    denominator = pair_times["T13"] + pair_times["T24"]
    ratio = (denominator / (pair_times["T12"] + pair_times["T34"])) if pair_times["T12"] + pair_times["T34"] > 0 else math.nan
    return {
        **pair_times,
        "mean_pairwise_time": mean_pair,
        "R_T": ratio,
        "total_tree_length": float(sum(n.length for n in nodes if n.parent is not None)),
        "internal_branch_length_q0": internal_by_topology[0],
        "internal_branch_length_q1": internal_by_topology[1],
        "internal_branch_length_q2": internal_by_topology[2],
        "pair_distances": distances,
        "newick": newick,
    }


def topology_from_split(result_topology: int) -> str:
    if int(result_topology) not in range(3):
        raise ValueError(f"invalid topology index {result_topology}")
    return TOPOLOGY_NAMES[int(result_topology)]


def _msc_commit(msrc_path: Path) -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=msrc_path, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "unknown"


def _imports(msrc_path: Path):
    if str(msrc_path) not in sys.path:
        sys.path.insert(0, str(msrc_path))
    from msrcsim.moran import simulate_frequency_history
    from msrcsim.rearrangement import Rearrangement
    from msrcsim.species_tree import SpeciesTree
    from msrcsim.structured_coalescent import simulate_genealogy, simulate_msc_genealogy
    return simulate_frequency_history, Rearrangement, SpeciesTree, simulate_genealogy, simulate_msc_genealogy


def _one_genealogy(result, model: str, pattern: str | None, history_id: int, gene_id: int) -> dict[str, object]:
    metrics = tree_metrics(result.newick)
    row = {k: v for k, v in metrics.items() if k not in {"pair_distances", "newick"}}
    row.update({"model": model, "history_id": history_id, "genealogy_id": gene_id, "topology_index": int(result.topology_index), "topology": topology_from_split(result.topology_index), "pattern": pattern or "", "n_coalescences": len(result.coalescences), "last_coalescence_time": max(result.coalescence_times) if result.coalescence_times else math.nan})
    return row


def _simulate_cell(msrc_path: Path, M: int, initial_fraction: float, rho: float, sister_length: float, root_length: float, histories: int, genes_per_history: int, seed: int, record_trees: bool = False):
    simulate_frequency_history, Rearrangement, SpeciesTree, simulate_genealogy, simulate_msc_genealogy = _imports(msrc_path)
    tree = SpeciesTree(TREE_NEWICK.format(sister=sister_length), M // 2, root_length)
    rr = Rearrangement("topology_time", "inversion", "ROOT", 0, int(round(M * initial_fraction)), 0.0)
    rng = np.random.default_rng(seed)
    rows: list[dict[str, object]] = []
    history_pattern_rows: list[dict[str, object]] = []
    for history_id in range(histories):
        history = simulate_frequency_history(tree, rr, rng)
        sampled = {taxon: int(rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
        pattern = "".join(str(sampled[taxon]) for taxon in tree.taxa)
        history_pattern_rows.append({"history_id": history_id, "pattern": pattern, "terminal_frequency": history.terminal_frequency("1")})
        for gene_id in range(genes_per_history):
            result = simulate_genealogy(gene_id, tree, history, sampled, rho, 1.0, rng, record_events=False)
            row = _one_genealogy(result, "msrc", pattern, history_id, gene_id)
            if record_trees:
                row["newick"] = result.newick
            rows.append(row)
    # MSC has one genealogy per independent history; its ``history_id`` is the
    # independent sampling unit for the same summary tables.
    for history_id in range(histories * genes_per_history):
        result = simulate_msc_genealogy(history_id, tree, rng, record_events=False)
        rows.append(_one_genealogy(result, "msc", "", history_id, 0))
    return rows, history_pattern_rows


def _theory(msrc_path: Path, M: int, initial_fraction: float, rho: float, sister_length: float, root_length: float) -> dict[str, object]:
    sys.path.insert(0, str(ROOT / "simulations" / "quartet" / "moran_theory"))
    from full_species_tree import full_species_tree_quartet_probabilities
    pi = np.zeros(M + 1)
    pi[int(round(M * initial_fraction))] = 1.0
    result = full_species_tree_quartet_probabilities(
        M_root=M, M_left=M, M_right=M, M_terminal=(M, M, M, M), rho=rho,
        root_length=root_length, sister_lengths=(sister_length, sister_length),
        terminal_lengths=(0.0,) * 4, pi_old=pi,
    )
    return {"q1": float(result.quartet_probabilities[0]), "q2": float(result.quartet_probabilities[1]), "q3": float(result.quartet_probabilities[2]), "joint": result.joint, "pattern": result.pattern_probabilities}


def _distribution_distance(a: np.ndarray, b: np.ndarray) -> float:
    a = np.sort(np.asarray(a, dtype=float)); b = np.sort(np.asarray(b, dtype=float))
    return float(np.mean(np.abs(a - b)))


def summarize_cell(rows: list[dict[str, object]], M: int, initial_fraction: float, rho: float, sister_length: float, root_length: float, theory: dict[str, object] | None) -> tuple[dict, list[dict], list[dict]]:
    out = {"M": M, "initial_fraction": initial_fraction, "rho": rho, "sister_length": sister_length, "root_length": root_length, "n_rows": len(rows), "n_msrc_histories": len({r["history_id"] for r in rows if r["model"] == "msrc"}), "n_msrc_genealogies": sum(r["model"] == "msrc" for r in rows), "n_msc_genealogies": sum(r["model"] == "msc" for r in rows)}
    topo_rows = []
    time_rows = []
    for model in ("msc", "msrc"):
        sub = [r for r in rows if r["model"] == model]
        counts = np.bincount([int(r["topology_index"]) for r in sub], minlength=3) / len(sub)
        mean_time = np.asarray([float(r["mean_pairwise_time"]) for r in sub])
        internal = np.asarray([float(r[f"internal_branch_length_q{int(r['topology_index'])}"]) for r in sub])
        total = np.asarray([float(r["total_tree_length"]) for r in sub])
        ratio = np.asarray([float(r["R_T"]) for r in sub])
        for k, name in enumerate(TOPOLOGY_NAMES):
            topo_rows.append({"M": M, "initial_fraction": initial_fraction, "rho": rho, "sister_length": sister_length, "root_length": root_length, "model": model, "topology": name, "probability": counts[k], "n": len(sub), "mc_se": math.sqrt(counts[k] * (1 - counts[k]) / len(sub))})
        for name, values in (("mean_pairwise_time", mean_time), ("internal_branch_length", internal), ("total_tree_length", total), ("R_T", ratio)):
            time_rows.append({"M": M, "initial_fraction": initial_fraction, "rho": rho, "sister_length": sister_length, "root_length": root_length, "model": model, "statistic": name, "mean": float(np.nanmean(values)), "sd": float(np.nanstd(values, ddof=1)), "q05": float(np.nanquantile(values, .05)), "q50": float(np.nanquantile(values, .5)), "q95": float(np.nanquantile(values, .95)), "n": int(np.isfinite(values).sum())})
        out[f"{model}_q1"], out[f"{model}_q2"], out[f"{model}_q3"] = map(float, counts)
        out[f"{model}_mean_pairwise_time"] = float(np.mean(mean_time))
        out[f"{model}_mean_pairwise_time_sd"] = float(np.std(mean_time, ddof=1))
        out[f"{model}_internal_branch_mean"] = float(np.mean(internal))
        out[f"{model}_total_tree_length_mean"] = float(np.mean(total))
    out["delta_Q"] = float(.5 * np.sum(np.abs(np.array([out["msrc_q1"], out["msrc_q2"], out["msrc_q3"]]) - np.array([out["msc_q1"], out["msc_q2"], out["msc_q3"]]))))
    out["delta_T_mean_pairwise"] = float(abs(out["msrc_mean_pairwise_time"] - out["msc_mean_pairwise_time"]) / out["msc_mean_pairwise_time_sd"]) if out["msc_mean_pairwise_time_sd"] > 0 else math.nan
    n_q = min(out["n_msrc_genealogies"], out["n_msc_genealogies"])
    q_se = .5 * sum(math.sqrt((out[f"msrc_q{k}"] * (1 - out[f"msrc_q{k}"]) + out[f"msc_q{k}"] * (1 - out[f"msc_q{k}"])) / n_q) for k in (1, 2, 3))
    out["delta_Q_mc_se_approx"] = float(q_se)
    out["delta_T_mc_se_approx"] = float(math.sqrt(out["msrc_mean_pairwise_time_sd"] ** 2 / n_q + out["msc_mean_pairwise_time_sd"] ** 2 / n_q) / out["msc_mean_pairwise_time_sd"]) if out["msc_mean_pairwise_time_sd"] > 0 else math.nan
    if theory is not None:
        out.update({"theory_q1": theory["q1"], "theory_q2": theory["q2"], "theory_q3": theory["q3"], "theory_abs_error_q1": abs(out["msrc_q1"] - theory["q1"]), "theory_abs_error_q2": abs(out["msrc_q2"] - theory["q2"]), "theory_abs_error_q3": abs(out["msrc_q3"] - theory["q3"])})
    return out, topo_rows, time_rows


def _rate_rows(rows: list[dict[str, object]], seed: int) -> list[dict[str, object]]:
    """Apply synthetic rates to stored pairwise distances for one validation cell."""
    rng = np.random.default_rng(seed)
    # The table compares true time with summaries under common, window-wide,
    # and branch-specific multipliers. Pairwise ratio behavior is exact under
    # the common/window-wide controls and generally changes under branch rates.
    records = []
    selected = [r for r in rows if r["model"] == "msrc"]
    for control in ("true_time", "constant_rate", "window_rate", "branch_rate"):
        values = []
        ratios = []
        for r in selected:
            if control == "true_time":
                vals = [float(r["T12"]), float(r["T13"]), float(r["T14"]), float(r["T23"]), float(r["T24"]), float(r["T34"])]
                ratio = float(r["R_T"])
            else:
                if "newick" not in r:
                    raise ValueError("rate-robustness rows require stored true genealogies")
                root = parse_newick(str(r["newick"]))
                edges: list[tuple[Node, frozenset[str]]] = []
                def collect_edges(n: Node) -> frozenset[str]:
                    if not n.children:
                        return frozenset((n.name,))
                    desc = frozenset().union(*(collect_edges(c) for c in n.children))
                    for child in n.children:
                        child_desc = collect_desc(child)
                        edges.append((child, child_desc))
                    return desc
                def collect_desc(n: Node) -> frozenset[str]:
                    return frozenset((n.name,)) if not n.children else frozenset().union(*(collect_desc(c) for c in n.children))
                collect_edges(root)
                factor = 1.0 if control == "constant_rate" else float(rng.lognormal(0.0, .5)) if control == "window_rate" else None
                multipliers = {desc: (factor if factor is not None else float(rng.lognormal(0.0, .35))) for _edge, desc in edges}
                distances = tree_metrics(str(r["newick"]), multipliers)["pair_distances"]
                vals = [float(distances[p]) for p in PAIRS]
                ratio = (vals[1] + vals[4]) / (vals[0] + vals[5]) if vals[0] + vals[5] > 0 else math.nan
            values.append(float(np.mean(vals))); ratios.append(ratio)
        records.append({"control": control, "statistic": "mean_pairwise_time_or_scaled_distance", "mean": float(np.mean(values)), "sd": float(np.std(values, ddof=1)), "mean_R_T": float(np.nanmean(ratios)), "sd_R_T": float(np.nanstd(ratios, ddof=1)), "n": len(values)})
    return records


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--msrc-sim", type=Path, default=DEFAULT_MSRC)
    parser.add_argument("--output", type=Path, default=ROOT / "simulations/topology_time_separation/results/discovery")
    parser.add_argument("--histories", type=int, default=100)
    parser.add_argument("--genes-per-history", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20261021)
    parser.add_argument("--M", type=int, nargs="+", default=[20, 40])
    parser.add_argument("--initial-fraction", type=float, nargs="+", default=[0.25, 0.5, 0.75])
    parser.add_argument("--rho", type=float, nargs="+", default=[0.0, 0.1, 0.5])
    parser.add_argument("--sister-length", type=float, nargs="+", default=[0.0, 0.25])
    parser.add_argument("--root-length", type=float, nargs="+", default=[0.25, 0.5])
    parser.add_argument("--validation", action="store_true")
    args = parser.parse_args()
    grid = [{"M": M, "initial_fraction": f, "rho": rho, "sister_length": sister, "root_length": root} for M in args.M for f in args.initial_fraction for rho in args.rho for sister in args.sister_length for root in args.root_length]
    args.output.mkdir(parents=True, exist_ok=True)
    topo_all: list[dict] = []; time_all: list[dict] = []; metrics: list[dict] = []; theory_rows: list[dict] = []; rate_rows: list[dict] = []
    for index, cell in enumerate(grid):
        rows, _patterns = _simulate_cell(args.msrc_sim, **cell, histories=args.histories, genes_per_history=args.genes_per_history, seed=args.seed + index * 1009, record_trees=(index == 0))
        theory = _theory(args.msrc_sim, **cell)
        summary, topo, times = summarize_cell(rows, **cell, theory=theory)
        metrics.append(summary); topo_all.extend(topo); time_all.extend(times)
        theory_rows.append({"M": cell["M"], "initial_fraction": cell["initial_fraction"], "rho": cell["rho"], "sister_length": cell["sister_length"], "root_length": cell["root_length"], "analytic_q1": theory["q1"], "analytic_q2": theory["q2"], "analytic_q3": theory["q3"], "empirical_q1": summary["msrc_q1"], "empirical_q2": summary["msrc_q2"], "empirical_q3": summary["msrc_q3"], "absolute_error_q1": summary["theory_abs_error_q1"], "absolute_error_q2": summary["theory_abs_error_q2"], "absolute_error_q3": summary["theory_abs_error_q3"], "scope": "theory_applies"})
        if index == 0:
            rate_rows.extend(_rate_rows(rows, args.seed + 991))
    write_tsv(args.output / "topology_metrics.tsv", topo_all); write_tsv(args.output / "time_metrics.tsv", time_all); write_tsv(args.output / "rate_robustness.tsv", rate_rows); write_tsv(args.output / "theory_comparison.tsv", theory_rows); write_tsv(args.output / "summary_metrics.tsv", metrics)
    (args.output / "parameter_grid.tsv").write_text("M\tinitial_fraction\trho\tsister_length\troot_length\n" + "\n".join("\t".join(str(cell[k]) for k in ("M", "initial_fraction", "rho", "sister_length", "root_length")) for cell in grid) + "\n")
    manifest = {"experiment": "MSRC topology-coalescence-time signal separation", "msrc_sim_path": str(args.msrc_sim), "msrc_sim_commit": _msc_commit(args.msrc_sim), "population_process": "authoritative Moran process", "time_observables": "true genealogy Newick/coalescence output", "histories_per_cell": args.histories, "genes_per_history": args.genes_per_history, "seed": args.seed, "validation_mode": args.validation, "grid_cells": len(grid), "theory_scope": "full finite-M quartet topology/joint law; no analytic time distribution", "command": " ".join(sys.argv)}
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"cells": len(grid), "histories_per_cell": args.histories, "genes_per_history": args.genes_per_history, "output": str(args.output)}))


if __name__ == "__main__":
    main()
