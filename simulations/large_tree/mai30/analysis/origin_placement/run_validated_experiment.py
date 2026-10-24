#!/usr/bin/env python3
"""Run resumable A/B Mai30 quartet validation with the authoritative msrc-sim.

Each compressed history file is an atomic checkpoint.  The forward history is
simulated once and reused across all rho values in that file.  Only compact
genealogy/topology diagnostics are retained; this script does not create a
full 30-taxon gene-tree collection.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import yaml

MSRC_ROOT = Path("/Users/ytabatabaee/Desktop/msrc-sim")
ROOT = Path(__file__).resolve().parents[5]
MAI = ROOT / "simulations/large_tree/mai30"
OUT = MAI / "analysis/origin_placement"
sys.path.insert(0, str(MSRC_ROOT))
sys.path.insert(0, str(MAI / "analysis/scripts"))

from compute_original_baseline import parse_newick, quartet_key  # noqa: E402
from msrcsim import SpeciesTree  # noqa: E402
from msrcsim.population_process import simulate_population_history  # noqa: E402
from msrcsim.rearrangement import Rearrangement  # noqa: E402
from msrcsim.structured_coalescent import simulate_genealogy  # noqa: E402
from msrcsim.wright_fisher import FrequencyHistory, FrequencyRecord  # noqa: E402

TREES = [
    ("outgroup.0.species.1.5.genes.1.5", "007"),
    ("outgroup.0.species.1.5.genes.1.5", "042"),
    ("outgroup.0.species.1.5.genes.1.5", "074"),
    ("outgroup.0.species.5.genes.5", "005"),
]
RHOS = [0.0, 0.0001, 0.001, 0.002, 0.01, 0.02]
SCENARIOS = ("A_existing_internal", "B_deep_internal")
SEED = 20261022
SIMULATOR_COMMIT = subprocess.check_output(["git", "-C", str(MSRC_ROOT), "rev-parse", "HEAD"], text=True).strip()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def config_for(condition: str, replicate: str) -> dict:
    path = MAI / "datasets/configs/mai30_msrc_pilot" / f"{condition}__replicate_{replicate}__msrc_no_suppression__freq_0.5.yaml"
    return yaml.safe_load(path.read_text())


def atomic_gzip_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent))
    os.close(fd)
    try:
        with gzip.open(name, "wt", encoding="utf-8") as handle:
            json.dump(payload, handle, separators=(",", ":"))
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def load_history(path: Path) -> dict:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def history_from_payload(payload: dict) -> FrequencyHistory:
    records = [FrequencyRecord(**row) for row in payload["records"]]
    by_branch: dict[str, list[FrequencyRecord]] = {}
    for record in records:
        by_branch.setdefault(record.branch_id, []).append(record)
    return FrequencyHistory(records=records, by_branch=by_branch,
                            structured_origin_age=payload["structured_origin_age"],
                            structured_origin_branch=payload["structured_origin_branch"])


def topology_labels(gene_tree, quartet: tuple[str, ...], species_key) -> tuple[str, str]:
    actual = quartet_key(gene_tree, quartet)
    labels = {}
    for i, pairs in enumerate(((0, 1), (0, 2), (0, 3)), 1):
        left = frozenset((quartet[pairs[0]], quartet[pairs[1]]))
        right = frozenset(quartet) - left
        labels[left if tuple(sorted(left)) <= tuple(sorted(right)) else right] = f"Q{i}"
    return labels.get(actual, "unresolved"), labels.get(species_key, "unresolved")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--histories", type=int, default=100)
    parser.add_argument("--rhos", nargs="+", type=float, default=RHOS)
    parser.add_argument("--tree-indices", nargs="+", type=int, default=[0, 1, 2, 3])
    parser.add_argument("--scenarios", nargs="+", choices=SCENARIOS, default=list(SCENARIOS))
    parser.add_argument("--output-dir", type=Path, default=OUT / "validated_100")
    parser.add_argument("--panel", type=Path, default=OUT / "selected_quartets.tsv")
    args = parser.parse_args()
    if any(rho not in RHOS for rho in args.rhos):
        raise SystemExit(f"rho values must be from {RHOS}")
    panel_rows = read_tsv(args.panel)
    scenario_rows = read_tsv(OUT / "origin_scenarios.tsv")
    panel = {}
    for row in panel_rows:
        if row["scenario"] in args.scenarios:
            key = (row["condition"], row["replicate"], row["scenario"])
            panel.setdefault(key, []).append(tuple(row["quartet"].split(",")))
    run_rows = []
    for tree_index, (condition, replicate) in enumerate(TREES):
        if tree_index not in args.tree_indices:
            continue
        config = config_for(condition, replicate)
        for scenario_index, scenario_name in enumerate(args.scenarios):
            scenario = next(r for r in scenario_rows if r["condition"] == condition and r["replicate"] == replicate and r["scenario"] == scenario_name)
            key = (condition, replicate, scenario_name)
            quartets = panel[key]
            tree = SpeciesTree(config["species_tree"]["newick"], config["species_tree"]["default_effective_population_size"], float(scenario["root_extension"]), ultrametric_tolerance=1e-6)
            species_tree = parse_newick(config["species_tree"]["newick"])
            species_keys = {",".join(q): quartet_key(species_tree, q) for q in quartets}
            rr = config["rearrangement"]
            rearrangement = Rearrangement("validated_origin", rr["type"], scenario["origin_branch"], int(scenario["origin_time_from_branch_start"]), int(rr["initial_copy_count"]), 0.0)
            config_dir = args.output_dir / f"{condition}__replicate_{replicate}__{scenario_name}"
            config_dir.mkdir(parents=True, exist_ok=True)
            start = time.perf_counter()
            completed = 0
            for history_index in range(args.histories):
                path = config_dir / f"history_{history_index:04d}.json.gz"
                expected = {"condition": condition, "replicate": replicate, "scenario": scenario_name,
                            "history": history_index, "rho_grid": args.rhos, "simulator_commit": SIMULATOR_COMMIT}
                if path.exists():
                    saved = load_history(path)
                    if all(saved.get(k) == v for k, v in expected.items()):
                        completed += 1
                        continue
                    raise RuntimeError(f"checkpoint metadata mismatch: {path}")
                seed = SEED + tree_index * 10_000_000 + scenario_index * 100_000 + history_index
                history = simulate_population_history(tree, rearrangement, np.random.default_rng(seed), config)
                sampled_rng = np.random.default_rng(seed + 777)
                sampled = {taxon: int(sampled_rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
                origin_records = history.by_branch[scenario["origin_branch"]]
                genealogy_rows = []
                for rho_index, rho in enumerate(args.rhos):
                    genealogy_seed = SEED + 1_000_000_000 + tree_index * 10_000_000 + scenario_index * 100_000 + history_index * 100 + rho_index
                    result = simulate_genealogy(history_index, tree, history, sampled, rho, 1.0, np.random.default_rng(genealogy_seed), True)
                    gene_tree = parse_newick(result.newick)
                    switches = sum(event.event_type == "switch" for event in result.events)
                    coals = sum(event.event_type == "coalescence" for event in result.events)
                    for quartet in quartets:
                        qkey = ",".join(quartet)
                        observed_topology, species_topology = topology_labels(gene_tree, quartet, species_keys[qkey])
                        genealogy_rows.append({"rho": rho, "quartet": qkey, "pattern": "".join(str(sampled[t]) for t in quartet),
                                               "topology": observed_topology, "species_topology": species_topology,
                                               "n_switches": switches, "n_coalescences": coals})
                payload = {**expected, "seed": seed, "simulator_root": str(MSRC_ROOT),
                           "origin_branch": scenario["origin_branch"], "origin_time_from_branch_start": int(scenario["origin_time_from_branch_start"]),
                           "root_extension": float(scenario["root_extension"]), "sampled_arrangements": sampled,
                           "origin_polymorphic_fraction": float(sum(0 < r.copy_count_A1 < r.population_chromosomes for r in origin_records) / len(origin_records)),
                           "origin_end_polymorphic": int(0 < origin_records[-1].copy_count_A1 < origin_records[-1].population_chromosomes),
                           "ever_fixed": int(any(r.status == "fixed" for r in origin_records)),
                           "ever_lost": int(any(r.status == "lost" for r in origin_records)),
                           "terminal_polymorphic": int(any(0 < history.terminal_frequency(t) < 1 for t in tree.taxa)),
                           "structured_origin_age": history.structured_origin_age, "structured_origin_branch": history.structured_origin_branch,
                           "records": [asdict(r) for r in history.records], "genealogies": genealogy_rows}
                atomic_gzip_json(path, payload)
                completed += 1
            run_rows.append({"condition": condition, "replicate": replicate, "scenario": scenario_name, "histories_requested": args.histories,
                             "histories_completed": completed, "rho_count": len(args.rhos), "panel_quartets": len(quartets),
                             "seconds": time.perf_counter() - start, "status": "complete" if completed == args.histories else "partial",
                             "simulator_commit": SIMULATOR_COMMIT, "output_dir": str(config_dir)})
    manifest = {"status": "complete" if all(r["status"] == "complete" for r in run_rows) else "partial",
                "histories_requested": args.histories, "rhos": args.rhos, "scenarios": args.scenarios,
                "trees": [TREES[i] for i in args.tree_indices], "forward_histories_reused_across_rho": True,
                "genealogies_per_history_per_rho": 1, "simulator_commit": SIMULATOR_COMMIT,
                "seed": SEED, "generated_utc": time.time(), "raw_checkpoints": str(args.output_dir)}
    (args.output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with (args.output_dir / "forward_histories_manifest.tsv").open("w", newline="") as handle:
        fields = list(run_rows[0])
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(run_rows)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
