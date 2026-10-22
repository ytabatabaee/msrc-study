#!/usr/bin/env python3
"""Repository-maintained Mai30 stage adapters.

All simulator imports are inside ``forward``/``backward``/``panel``.  This
module is therefore suitable for static inspection and schema tests without a
local simulator checkout.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from statistics import NormalDist
from typing import Any

from mai30_workflow import (
    BASE_SEED,
    EXPECTED_SIMULATOR_COMMIT,
    HERE,
    MAI,
    MANIFEST,
    RHOS,
    TREE_SPECS,
    atomic_json,
    atomic_write,
    backward_dir,
    config_metadata,
    config_path,
    config_key,
    forward_dir,
    genealogy_seed,
    history_path,
    history_payload,
    kind_for,
    leaf_names,
    load_simulator,
    load_yaml,
    output_manifest,
    output_root,
    parse_splits,
    quartet_topology,
    read_gzip_json,
    rebuild_history,
    scenario_row,
    seed_for,
    selected_row,
    sha256_file,
    sha256_text,
    task_dir,
    task_manifest,
    validate_history_payload,
    validate_newick,
    write_gzip_json,
    require_live_compute_allocation,
)


def checkpoint_rows(row: dict[str, str], root: Path) -> list[tuple[int, Path, dict[str, Any]]]:
    result = []
    for history_id in range(int(row["history_count"])):
        path = history_path(row, history_id, root)
        if not path.exists():
            continue
        payload = read_gzip_json(path)
        validate_history_payload(payload, row, history_id)
        result.append((history_id, path, payload))
    return result


def make_tree(row: dict[str, str], api: dict[str, Any]) -> tuple[Any, dict[str, Any], dict[str, str]]:
    config = load_yaml(config_path(row["condition"], row["replicate"]))
    scenario = scenario_row(row["condition"], row["replicate"], row["scenario"])
    tree = api["SpeciesTree"](
        config["species_tree"]["newick"],
        config["species_tree"]["default_effective_population_size"],
        float(scenario["root_extension"]),
        ultrametric_tolerance=1e-6,
    )
    return tree, config, scenario


def forward(row: dict[str, str], task_id: int, root: Path) -> None:
    api = load_simulator()
    tree, config, scenario = make_tree(row, api)
    metadata = config_metadata(row)
    destination = forward_dir(row, root)
    expected = int(row["history_count"])
    complete = checkpoint_rows(row, root)
    if len(complete) == expected:
        output_manifest(task_manifest("forward" if row["stage"] == "forward" else "benchmark_forward", task_id, root), row["stage"], task_id, row, [x[1] for x in complete], expected, len(complete), [str(x[0]) for x in complete], extra={"checkpoint_dir": str(destination), **metadata})
        return
    if complete and len(complete) < expected:
        # Existing checkpoints are valid and can be resumed.  Invalid files
        # are rejected by checkpoint_rows rather than silently replaced.
        pass
    destination.mkdir(parents=True, exist_ok=True)
    tree_index, scenario_index = (TREE_SPECS.index((row["condition"], row["replicate"])), 0 if row["scenario"].startswith("A") else 1)
    rr = config["rearrangement"]
    rearrangement = api["Rearrangement"](
        "validated_origin", rr["type"], scenario["origin_branch"],
        int(scenario["origin_time_from_branch_start"]), int(rr["initial_copy_count"]), 0.0,
    )
    for history_id in range(expected):
        path = history_path(row, history_id, root)
        if path.exists():
            payload = read_gzip_json(path)
            validate_history_payload(payload, row, history_id)
            continue
        seed = seed_for(tree_index, scenario_index, history_id)
        import numpy as np  # imported only inside a compute-stage adapter
        history = api["simulate_population_history"](tree, rearrangement, np.random.default_rng(seed), config)
        sampled_rng = np.random.default_rng(seed + 5_000_000)
        sampled = {taxon: int(sampled_rng.random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
        payload = history_payload(row, history_id, tree, history, sampled, seed, metadata)
        write_gzip_json(path, payload)
        # Read back the compressed checkpoint before it can count as complete.
        validate_history_payload(read_gzip_json(path), row, history_id)
    complete = checkpoint_rows(row, root)
    if len(complete) != expected:
        raise RuntimeError(f"forward checkpoint count is {len(complete)}, expected {expected}")
    output_manifest(task_manifest(row["stage"], task_id, root), row["stage"], task_id, row, [x[1] for x in complete], expected, len(complete), [str(x[0]) for x in complete], extra={"checkpoint_dir": str(destination), **metadata})


def forward_persistence(row: dict[str, str], task_id: int, root: Path) -> None:
    rows_for_config = [selected_row("forward", i) for i in range(8)]
    output = task_dir("persistence", task_id, root) / "forward_persistence.tsv"
    lines = ["condition\treplicate\tscenario\thistory_id\tterminal_persistence\tjoint_persistence\n"]
    for forward_row in rows_for_config:
        for history_id, _, payload in checkpoint_rows(forward_row, root):
            frequencies = payload.get("terminal_frequencies", {})
            persistence = sum(float(v) for v in frequencies.values()) / max(1, len(frequencies))
            joint = math.prod(float(v) for v in frequencies.values()) if frequencies else 0.0
            lines.append(f"{forward_row['condition']}\t{forward_row['replicate']}\t{forward_row['scenario']}\t{history_id}\t{persistence:.12g}\t{joint:.12g}\n")
    if len(lines) != 1 + 8 * int(row["history_count"]):
        raise RuntimeError("forward persistence did not find all validated origin scenarios")
    atomic_write(output, "".join(lines).encode())
    output_manifest(task_manifest("persistence", task_id, root), row["stage"], task_id, row, [output], len(lines) - 1, len(lines) - 1, columns=lines[0].rstrip().split("\t"))


def panel(row: dict[str, str], task_id: int, root: Path) -> None:
    api = load_simulator()
    output = task_dir("panel", task_id, root) / "selected_quartets.tsv"
    rows_out: list[str] = ["condition\treplicate\tscenario\tquartet_id\ttaxon_a\ttaxon_b\ttaxon_c\ttaxon_d\texposure_count\tselection_rule\tspecies_topology\n"]

    def ancestors(node: Any) -> list[Any]:
        result = []
        while node is not None:
            result.append(node)
            node = node.parent
        return result

    def mrca(tree: Any, names: set[str]) -> Any:
        nodes = {node.name: node for node in tree.nodes()}
        candidates = ancestors(nodes[next(iter(names))])
        for candidate in candidates:
            if all(candidate in ancestors(nodes[name]) for name in names):
                return candidate
        raise RuntimeError("species-tree MRCA not found")

    def split(tree: Any, quartet: tuple[str, ...]) -> tuple[frozenset[str], frozenset[str]]:
        for left in itertools.combinations(quartet, 2):
            right = tuple(name for name in quartet if name not in left)
            if mrca(tree, set(left)) is not mrca(tree, set(right)):
                a, b = frozenset(left), frozenset(right)
                return (a, b) if tuple(sorted(a)) <= tuple(sorted(b)) else (b, a)
        raise RuntimeError(f"unresolved species quartet: {quartet}")

    for condition, replicate in TREE_SPECS:
        panel_row = {"condition": condition, "replicate": replicate, "scenario": "A_existing_internal", "rho": "", "history_count": "0", "stage": "panel"}
        tree, _, _ = make_tree(panel_row, api)
        for scenario_name in ("A_existing_internal", "B_deep_internal"):
            scenario = scenario_row(condition, replicate, scenario_name)
            descendants = set(tree.descendants(scenario["origin_branch"])) if scenario["origin_branch"] != tree.root.name else set(tree.taxa)
            bins: dict[int, list[dict[str, Any]]] = {2: [], 3: [], 4: []}
            for quartet in itertools.combinations(sorted(tree.taxa), 4):
                exposed = len(descendants & set(quartet))
                if exposed not in bins:
                    continue
                left, right = split(tree, quartet)
                all_mrca = mrca(tree, set(quartet))
                pair_mrca = max(mrca(tree, set(left)).age, mrca(tree, set(right)).age)
                bins[exposed].append({"quartet": quartet, "exposure_count": exposed, "internal_length": float(all_mrca.age - pair_mrca)})
            for exposure in (2, 3, 4):
                candidates = sorted(bins[exposure], key=lambda item: (item["internal_length"], item["quartet"]))
                chosen = candidates[:2] + candidates[-2:]
                seen: set[tuple[str, ...]] = set()
                for rank, item in enumerate(chosen, 1):
                    quartet = item["quartet"]
                    if quartet in seen:
                        continue
                    seen.add(quartet)
                    qid = f"{condition}_{replicate}_{scenario_name}_{exposure}_{rank}"
                    topology = quartet_topology(load_yaml(config_path(condition, replicate))["species_tree"]["newick"], quartet)
                    rule = f"geometry-only exposure={exposure}; {'short' if rank <= 2 else 'long'} induced internal branch"
                    rows_out.append("\t".join([condition, replicate, scenario_name, qid, *quartet, str(exposure), rule, topology]) + "\n")
    atomic_write(output, "".join(rows_out).encode())
    output_manifest(task_manifest("panel", task_id, root), row["stage"], task_id, row, [output], len(rows_out) - 1, len(rows_out) - 1, columns=rows_out[0].rstrip().split("\t"))


def backward(row: dict[str, str], task_id: int, root: Path) -> None:
    api = load_simulator()
    tree, config, _ = make_tree(row, api)
    checkpoints = checkpoint_rows(row, root)
    expected = int(row["history_count"])
    if len(checkpoints) != expected:
        raise RuntimeError(f"backward stage requires {expected} validated forward histories, found {len(checkpoints)}")
    destination = backward_dir(row, root)
    output = destination / "genealogies.tsv"
    if output.exists():
        validate_genealogy_file(output, tree.taxa, expected)
        output_manifest(task_manifest(row["stage"], task_id, root), row["stage"], task_id, row, [output], expected, expected, [str(i) for i in range(expected)], columns=GENEALOGY_COLUMNS)
        return
    destination.mkdir(parents=True, exist_ok=True)
    tree_index, scenario_index = (TREE_SPECS.index((row["condition"], row["replicate"])), 0 if row["scenario"].startswith("A") else 1)
    rho_index = RHOS.index(row["rho"])
    rr = config["rearrangement"]
    lines = ["\t".join(GENEALOGY_COLUMNS) + "\n"]
    for history_id, _, payload in checkpoints:
        history = rebuild_history(payload, row, tree)
        result = api["simulate_genealogy"](
            history_id, tree, history, payload["sampled"], float(row["rho"]), 1.0,
            __import__("numpy").random.default_rng(genealogy_seed(tree_index, scenario_index, rho_index, history_id)), True,
        )
        validate_newick(result.newick, tree.taxa)
        events = getattr(result, "events", [])
        switches = sum(1 for event in events if getattr(event, "event_type", getattr(event, "type", "")) == "switch")
        coalescences = sum(1 for event in events if getattr(event, "event_type", getattr(event, "type", "")) == "coalescence")
        fields = [str(history_id), str(history_id), row["rho"], row["condition"], row["replicate"], row["scenario"], str(genealogy_seed(tree_index, scenario_index, rho_index, history_id)), str(switches), str(coalescences), result.newick.rstrip("\n")]
        lines.append("\t".join(fields) + "\n")
    atomic_write(output, "".join(lines).encode())
    validate_genealogy_file(output, tree.taxa, expected)
    output_manifest(task_manifest(row["stage"], task_id, root), row["stage"], task_id, row, [output], expected, expected, [str(i) for i in range(expected)], columns=GENEALOGY_COLUMNS, extra={"taxa_sha256": sha256_text("\n".join(sorted(tree.taxa))), **config_metadata(row)})


GENEALOGY_COLUMNS = ["genealogy_id", "history_id", "rho", "condition", "replicate", "scenario", "seed", "switches", "coalescences", "newick"]


def validate_genealogy_file(path: Path, taxa: Any, expected: int) -> None:
    with path.open() as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if len(records) != expected or [r["genealogy_id"] for r in records] != [str(i) for i in range(expected)]:
        raise RuntimeError(f"genealogy count/IDs invalid in {path}")
    for record in records:
        validate_newick(record["newick"], taxa)


def panel_records(root: Path) -> list[dict[str, str]]:
    path = task_dir("panel", 0, root) / "selected_quartets.tsv"
    if not path.exists():
        raise RuntimeError("validated quartet panel is missing")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def analysis(row: dict[str, str], task_id: int, root: Path) -> None:
    backward_rows = [r for r in [selected_row("backward", i) for i in range(48)]]
    panel_data = panel_records(root)
    output = task_dir("analysis", task_id, root)
    probability = output / "quartet_probabilities.tsv"
    margins = output / "quartet_margins.tsv"
    rows_probability: list[list[str]] = []
    rows_margin: list[list[str]] = []
    for backward_row in backward_rows:
        path = backward_dir(backward_row, root) / "genealogies.tsv"
        if not path.exists():
            raise RuntimeError(f"missing validated genealogy collection: {path}")
        with path.open(newline="") as handle:
            genes = list(csv.DictReader(handle, delimiter="\t"))
        if len(genes) != int(backward_row["history_count"]):
            raise RuntimeError(f"genealogy count mismatch for {config_key(backward_row)}")
        relevant = [q for q in panel_data if q["condition"] == backward_row["condition"] and q["replicate"] == backward_row["replicate"]]
        for quartet in relevant:
            labels = tuple(quartet[x] for x in ("taxon_a", "taxon_b", "taxon_c", "taxon_d"))
            counts: dict[str, int] = {}
            for gene in genes:
                topology = quartet_topology(gene["newick"], labels)
                counts[topology] = counts.get(topology, 0) + 1
            total = len(genes)
            observed = quartet["species_topology"]
            p = counts.get(observed, 0) / total
            alternatives = sorted((v / total for k, v in counts.items() if k not in {observed, "unresolved", "missing"}), reverse=True)
            competitor = alternatives[0] if alternatives else 0.0
            margin = p - competitor
            rows_probability.append([backward_row["condition"], backward_row["replicate"], backward_row["scenario"], backward_row["rho"], quartet["quartet_id"], str(total), f"{p:.12g}", f"{counts.get('unresolved', 0) / total:.12g}"])
            rows_margin.append([backward_row["condition"], backward_row["replicate"], backward_row["scenario"], backward_row["rho"], quartet["quartet_id"], str(total), f"{margin:.12g}"])
    probability_header = ["condition", "replicate", "scenario", "rho", "quartet_id", "history_count", "p_species", "p_unresolved"]
    margin_header = ["condition", "replicate", "scenario", "rho", "quartet_id", "history_count", "margin"]
    output.mkdir(parents=True, exist_ok=True)
    atomic_write(probability, ("\t".join(probability_header) + "\n" + "\n".join("\t".join(row) for row in rows_probability) + "\n").encode())
    atomic_write(margins, ("\t".join(margin_header) + "\n" + "\n".join("\t".join(row) for row in rows_margin) + "\n").encode())
    output_manifest(task_manifest("analysis", task_id, root), row["stage"], task_id, row, [probability, margins], len(rows_probability), len(rows_probability), columns=probability_header, extra={"margin_rows": len(rows_margin)})


def validation(row: dict[str, str], task_id: int, root: Path) -> None:
    source = task_dir("analysis", 0, root) / "quartet_margins.tsv"
    if not source.exists():
        raise RuntimeError("quartet analysis output is missing")
    with source.open(newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if not records:
        raise RuntimeError("quartet analysis output is empty")
    output = task_dir("validation", task_id, root) / "quartet_reversal_validation.tsv"
    header = ["condition", "replicate", "scenario", "rho", "quartet_id", "margin", "supported_reversal"]
    rows_out = [header]
    for record in records:
        margin = float(record["margin"])
        rows_out.append([record["condition"], record["replicate"], record["scenario"], record["rho"], record["quartet_id"], record["margin"], str(margin < 0).lower()])
    atomic_write(output, ("\t".join(header) + "\n" + "\n".join("\t".join(x) for x in rows_out[1:]) + "\n").encode())
    output_manifest(task_manifest("validation", task_id, root), row["stage"], task_id, row, [output], len(records), len(records), columns=header, extra={"gate": "descriptive_only; independent validation required"})


def astral(row: dict[str, str], task_id: int, root: Path) -> None:
    astral_bin = os.environ.get("ASTRAL")
    if not astral_bin:
        raise RuntimeError("ASTRAL must be supplied in the Bridges2 environment for the optional stage")
    destination = task_dir("astral", task_id, root)
    input_path = destination / "all_gene_trees.tre"
    trees: list[str] = []
    for i in range(48):
        backward_row = selected_row("backward", i)
        path = backward_dir(backward_row, root) / "genealogies.tsv"
        if not path.exists():
            raise RuntimeError(f"missing validated ASTRAL input: {path}")
        with path.open(newline="") as handle:
            for gene in csv.DictReader(handle, delimiter="\t"):
                trees.append(gene["newick"])
    if len(trees) != 48 * 100:
        raise RuntimeError(f"ASTRAL input has {len(trees)} trees; expected 4800")
    atomic_write(input_path, ("\n".join(trees) + "\n").encode())
    output = destination / "species_tree.tre"
    completed = subprocess.run([astral_bin, "-i", str(input_path), "-o", str(output)], cwd=str(destination))
    if completed.returncode != 0 or not output.exists() or not output.read_text().strip():
        raise RuntimeError("ASTRAL failed or did not create a species tree")
    validate_newick(output.read_text().strip(), leaf_names(trees[0]))
    output_manifest(task_manifest("astral", task_id, root), row["stage"], task_id, row, [input_path, output], len(trees), len(trees), columns=["newick"])


DISPATCH = {
    "forward": forward,
    "benchmark_forward": forward,
    "persistence": forward_persistence,
    "panel": panel,
    "backward": backward,
    "benchmark_backward": backward,
    "analysis": analysis,
    "validation": validation,
    "astral": astral,
}


def main(stage: str | None = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", default=stage, choices=tuple(DISPATCH))
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--output-dir", type=Path, default=output_root())
    args = parser.parse_args()
    require_live_compute_allocation()
    row = selected_row(args.stage, args.task_id)
    DISPATCH[args.stage](row, args.task_id, args.output_dir.resolve())


if __name__ == "__main__":
    main()
