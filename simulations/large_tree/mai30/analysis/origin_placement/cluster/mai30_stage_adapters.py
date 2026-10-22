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
    induced_quartet_split,
    kind_for,
    leaf_names,
    load_simulator,
    load_yaml,
    output_manifest,
    output_root,
    parse_splits,
    canonical_quartet_topologies,
    quartet_topology,
    read_gzip_json,
    rebuild_history,
    scenario_row,
    boundary_observations,
    sampled_pattern,
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
        payload.update({
            "origin_branch": scenario["origin_branch"],
            "origin_time_from_branch_start": int(scenario["origin_time_from_branch_start"]),
            "root_extension": float(scenario["root_extension"]),
            "rearrangement_type": rr["type"],
            "initial_copy_count": int(rr["initial_copy_count"]),
            "species_tree_newick": config["species_tree"]["newick"],
        })
        write_gzip_json(path, payload)
        # Read back the compressed checkpoint before it can count as complete.
        validate_history_payload(read_gzip_json(path), row, history_id)
    complete = checkpoint_rows(row, root)
    if len(complete) != expected:
        raise RuntimeError(f"forward checkpoint count is {len(complete)}, expected {expected}")
    output_manifest(task_manifest(row["stage"], task_id, root), row["stage"], task_id, row, [x[1] for x in complete], expected, len(complete), [str(x[0]) for x in complete], extra={"checkpoint_dir": str(destination), **metadata})


def forward_persistence(row: dict[str, str], task_id: int, root: Path) -> None:
    rows_for_config = [selected_row("forward", i) for i in range(8)]
    output_dir = task_dir("persistence", task_id, root)
    output = output_dir / "forward_persistence.tsv"
    joint_output = output_dir / "joint_ancestral_persistence.tsv"
    pattern_output = output_dir / "terminal_arrangement_patterns.tsv"
    lines = ["condition\treplicate\tscenario\tboundary_id\tbranch_id\tparent_branch_id\tboundary_role\tabsolute_age\tpopulation_chromosomes\tmean_arrangement_frequency\tprobability_polymorphic\tprobability_fixed\tprobability_lost\thistory_count\n"]
    joint_lines = ["condition\treplicate\tscenario\tboundary_set\tprobability_joint_persistence\thistory_count\n"]
    pattern_lines = ["condition\treplicate\tscenario\tpattern\tcount\tprobability\thistory_count\ttaxa_order\n"]
    for forward_row in rows_for_config:
        histories = checkpoint_rows(forward_row, root)
        expected = int(forward_row["history_count"])
        if len(histories) != expected:
            raise RuntimeError(f"persistence requires {expected} validated histories for {config_key(forward_row)}")
        by_boundary: dict[str, list[dict[str, Any]]] = {}
        patterns: dict[str, int] = {}
        taxa_order = None
        origin_branch = None
        for history_id, _, payload in histories:
            observations = boundary_observations(payload)
            if not observations:
                raise RuntimeError(f"checkpoint {history_id} has no recorded population boundaries")
            for observation in observations:
                by_boundary.setdefault(observation["boundary_id"], []).append(observation)
            patterns[sampled_pattern(payload)] = patterns.get(sampled_pattern(payload), 0) + 1
            taxa_order = ",".join(sorted(payload["taxa"]))
            origin_branch = payload.get("structured_origin_branch")
        for boundary_id, observations in sorted(by_boundary.items()):
            if len(observations) != expected:
                raise RuntimeError(f"boundary {boundary_id} is not present in every history")
            first = observations[0]
            if any((x["population_chromosomes"], x["branch_id"], x["parent_branch_id"], x["role"]) != (first["population_chromosomes"], first["branch_id"], first["parent_branch_id"], first["role"]) for x in observations):
                raise RuntimeError(f"boundary metadata differs across histories: {boundary_id}")
            mean_frequency = sum(x["arrangement_frequency"] for x in observations) / expected
            probability_polymorphic = sum(x["polymorphic"] for x in observations) / expected
            probability_fixed = sum(x["fixed"] for x in observations) / expected
            probability_lost = sum(x["lost"] for x in observations) / expected
            lines.append("\t".join([
                forward_row["condition"], forward_row["replicate"], forward_row["scenario"], boundary_id,
                first["branch_id"], str(first["parent_branch_id"] or ""), first["role"], f"{first['absolute_age']:.12g}",
                str(first["population_chromosomes"]), f"{mean_frequency:.12g}", f"{probability_polymorphic:.12g}",
                f"{probability_fixed:.12g}", f"{probability_lost:.12g}", str(expected),
            ]) + "\n")
        if origin_branch is None:
            raise RuntimeError("checkpoint lacks structured origin branch")
        joint_set = (f"{origin_branch}:start", f"{origin_branch}:end")
        joint_values = []
        for history_id, _, payload in histories:
            states = {x["boundary_id"]: x for x in boundary_observations(payload)}
            if any(boundary not in states for boundary in joint_set):
                raise RuntimeError(f"origin boundary set is incomplete in history {history_id}")
            joint_values.append(int(all(states[boundary]["polymorphic"] for boundary in joint_set)))
        joint_lines.append("\t".join([forward_row["condition"], forward_row["replicate"], forward_row["scenario"], ",".join(joint_set), f"{sum(joint_values) / expected:.12g}", str(expected)]) + "\n")
        for pattern, count in sorted(patterns.items()):
            pattern_lines.append("\t".join([forward_row["condition"], forward_row["replicate"], forward_row["scenario"], pattern, str(count), f"{count / expected:.12g}", str(expected), taxa_order or ""]) + "\n")
    output_dir.mkdir(parents=True, exist_ok=True)
    atomic_write(output, "".join(lines).encode())
    atomic_write(joint_output, "".join(joint_lines).encode())
    atomic_write(pattern_output, "".join(pattern_lines).encode())
    files = [output, joint_output, pattern_output]
    total_rows = len(lines) + len(joint_lines) + len(pattern_lines) - 3
    output_manifest(task_manifest("persistence", task_id, root), row["stage"], task_id, row, files, total_rows, total_rows, columns=lines[0].rstrip().split("\t"), extra={"persistence_definition": "I_v=1{0<K_v<M_v} at recorded branch start/end boundaries; joint set is origin branch start and end", "terminal_patterns_are_sampled_histories": True})


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

    for condition, replicate in TREE_SPECS:
        panel_row = {"condition": condition, "replicate": replicate, "scenario": "A_existing_internal", "rho": "", "history_count": "0", "stage": "panel"}
        tree, config, _ = make_tree(panel_row, api)
        species_newick = config["species_tree"]["newick"]
        for scenario_name in ("A_existing_internal", "B_deep_internal"):
            scenario = scenario_row(condition, replicate, scenario_name)
            descendants = set(tree.descendants(scenario["origin_branch"])) if scenario["origin_branch"] != tree.root.name else set(tree.taxa)
            bins: dict[int, list[dict[str, Any]]] = {2: [], 3: [], 4: []}
            for quartet in itertools.combinations(sorted(tree.taxa), 4):
                exposed = len(descendants & set(quartet))
                if exposed not in bins:
                    continue
                topology_split = induced_quartet_split(species_newick, quartet)
                if topology_split is None:
                    continue
                left, right = (set(topology_split[0]), set(topology_split[1]))
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
                    topology = quartet_topology(species_newick, quartet)
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
        validate_genealogy_file(output, tree.taxa, expected, row)
        output_manifest(task_manifest(row["stage"], task_id, root), row["stage"], task_id, row, [output], expected, expected, [str(i) for i in range(expected)], columns=GENEALOGY_COLUMNS, extra=config_metadata(row))
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
    validate_genealogy_file(output, tree.taxa, expected, row)
    output_manifest(task_manifest(row["stage"], task_id, root), row["stage"], task_id, row, [output], expected, expected, [str(i) for i in range(expected)], columns=GENEALOGY_COLUMNS, extra={"taxa_sha256": sha256_text("\n".join(sorted(tree.taxa))), **config_metadata(row)})


GENEALOGY_COLUMNS = ["genealogy_id", "history_id", "rho", "condition", "replicate", "scenario", "seed", "switches", "coalescences", "newick"]


def validate_genealogy_file(path: Path, taxa: Any, expected: int, row: dict[str, str] | None = None) -> None:
    with path.open() as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if len(records) != expected or [r["genealogy_id"] for r in records] != [str(i) for i in range(expected)]:
        raise RuntimeError(f"genealogy count/IDs invalid in {path}")
    for record in records:
        if record["history_id"] != record["genealogy_id"]:
            raise RuntimeError(f"genealogy/history ID mismatch in {path}")
        if row is not None and (record["rho"] != row["rho"] or record["condition"] != row["condition"] or record["replicate"] != row["replicate"] or record["scenario"] != row["scenario"]):
            raise RuntimeError(f"genealogy metadata mismatch in {path}")
        validate_newick(record["newick"], taxa)


def panel_records(root: Path) -> list[dict[str, str]]:
    path = task_dir("panel", 0, root) / "selected_quartets.tsv"
    if not path.exists():
        raise RuntimeError("validated quartet panel is missing")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def analysis(row: dict[str, str], task_id: int, root: Path) -> None:
    backward_rows = [selected_row("backward", i) for i in range(48)]
    panel_data = panel_records(root)
    output = task_dir("analysis", task_id, root)
    probability = output / "quartet_probabilities.tsv"
    margins = output / "quartet_margins.tsv"
    observations = output / "quartet_history_observations.tsv"
    probability_header = ["condition", "replicate", "scenario", "rho", "quartet_id", "history_count", "p_species", "p_discordant_1", "p_discordant_2", "p_unresolved", "probability_sum"]
    margin_header = ["condition", "replicate", "scenario", "rho", "quartet_id", "history_count", "margin_species_d1", "margin_species_d2", "margin_d1_d2", "margin_d2_d1"]
    observation_header = ["condition", "replicate", "scenario", "rho", "quartet_id", "history_id", "topology", "species_topology"]
    rows_probability: list[list[str]] = []
    rows_margin: list[list[str]] = []
    rows_observations: list[list[str]] = []
    for backward_row in backward_rows:
        path = backward_dir(backward_row, root) / "genealogies.tsv"
        if not path.exists():
            raise RuntimeError(f"missing validated genealogy collection: {path}")
        with path.open(newline="") as handle:
            genes = list(csv.DictReader(handle, delimiter="\t"))
        if len(genes) != int(backward_row["history_count"]):
            raise RuntimeError(f"genealogy count mismatch for {config_key(backward_row)}")
        relevant = [q for q in panel_data if q["condition"] == backward_row["condition"] and q["replicate"] == backward_row["replicate"] and q["scenario"] == backward_row["scenario"]]
        for quartet in relevant:
            labels = tuple(quartet[x] for x in ("taxon_a", "taxon_b", "taxon_c", "taxon_d"))
            all_topologies = canonical_quartet_topologies(labels)
            species_topology = quartet["species_topology"]
            if species_topology not in all_topologies:
                raise RuntimeError(f"panel has invalid species topology for {quartet['quartet_id']}: {species_topology}")
            discordant = tuple(topology for topology in all_topologies if topology != species_topology)
            counts = {topology: 0 for topology in all_topologies}
            unresolved = 0
            for gene in genes:
                topology = quartet_topology(gene["newick"], labels)
                if topology in counts:
                    counts[topology] += 1
                else:
                    unresolved += 1
                    topology = "unresolved"
                rows_observations.append([backward_row["condition"], backward_row["replicate"], backward_row["scenario"], backward_row["rho"], quartet["quartet_id"], gene["history_id"], topology, species_topology])
            total = len(genes)
            probabilities = {key: value / total for key, value in counts.items()}
            p_species = probabilities[species_topology]
            p_d1, p_d2 = (probabilities[discordant[0]], probabilities[discordant[1]])
            probability_sum = p_species + p_d1 + p_d2 + unresolved / total
            if abs(probability_sum - 1.0) > 1e-10:
                raise RuntimeError(f"quartet probability normalization failed for {quartet['quartet_id']}")
            rows_probability.append([backward_row["condition"], backward_row["replicate"], backward_row["scenario"], backward_row["rho"], quartet["quartet_id"], str(total), f"{p_species:.12g}", f"{p_d1:.12g}", f"{p_d2:.12g}", f"{unresolved / total:.12g}", f"{probability_sum:.12g}"])
            rows_margin.append([backward_row["condition"], backward_row["replicate"], backward_row["scenario"], backward_row["rho"], quartet["quartet_id"], str(total), f"{p_species - p_d1:.12g}", f"{p_species - p_d2:.12g}", f"{p_d1 - p_d2:.12g}", f"{p_d2 - p_d1:.12g}"])
    output.mkdir(parents=True, exist_ok=True)
    atomic_write(probability, ("\t".join(probability_header) + "\n" + "\n".join("\t".join(item) for item in rows_probability) + "\n").encode())
    atomic_write(margins, ("\t".join(margin_header) + "\n" + "\n".join("\t".join(item) for item in rows_margin) + "\n").encode())
    atomic_write(observations, ("\t".join(observation_header) + "\n" + "\n".join("\t".join(item) for item in rows_observations) + "\n").encode())
    output_manifest(task_manifest("analysis", task_id, root), row["stage"], task_id, row, [probability, margins, observations], len(rows_probability) + len(rows_observations), len(rows_probability) + len(rows_observations), columns=probability_header, extra={"margin_rows": len(rows_margin), "history_observation_rows": len(rows_observations), "independent_replication_unit": "forward_history"})


def validation(row: dict[str, str], task_id: int, root: Path) -> None:
    analysis_dir = task_dir("analysis", 0, root)
    source = analysis_dir / "quartet_margins.tsv"
    observation_file = analysis_dir / "quartet_history_observations.tsv"
    if not source.exists() or not observation_file.exists():
        raise RuntimeError("quartet analysis outputs are incomplete")
    with source.open(newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    with observation_file.open(newline="") as handle:
        observations = list(csv.DictReader(handle, delimiter="\t"))
    if not records or not observations:
        raise RuntimeError("quartet analysis output is empty")
    by_key: dict[tuple[str, str, str, str, str], list[dict[str, str]]] = {}
    for observation in observations:
        key = tuple(observation[name] for name in ("condition", "replicate", "scenario", "rho", "quartet_id"))
        by_key.setdefault(key, []).append(observation)
    comparisons = max(1, len(records) * 4)
    critical_z = NormalDist().inv_cdf(1.0 - 0.05 / (2.0 * comparisons))
    output = task_dir("validation", task_id, root) / "quartet_reversal_validation.tsv"
    header = ["condition", "replicate", "scenario", "rho", "quartet_id", "margin_species_d1", "margin_species_d2", "margin_d1_d2", "margin_d2_d1", "se_species_d1", "se_species_d2", "se_d1_d2", "se_d2_d1", "critical_z", "supported_d1", "supported_d2", "supported_reversal", "exploratory_negative_species_margin", "validation_scope"]
    rows_out = [header]
    for record in records:
        key = tuple(record[name] for name in ("condition", "replicate", "scenario", "rho", "quartet_id"))
        group = by_key.get(key, [])
        if len(group) != int(record["history_count"]):
            raise RuntimeError(f"history observation count mismatch for {key}")
        species = group[0]["species_topology"]
        topologies = canonical_quartet_topologies((species.split("::")[0].split("|")[0], species.split("::")[0].split("|")[1], species.split("::")[1].split("|")[0], species.split("::")[1].split("|")[1]))
        discordant = tuple(topology for topology in topologies if topology != species)
        if len(discordant) != 2:
            raise RuntimeError(f"invalid species topology in observations for {key}")

        def margin_and_se(left: str, right: str) -> tuple[float, float]:
            values = [int(item["topology"] == left) - int(item["topology"] == right) for item in group]
            mean = sum(values) / len(values)
            variance = sum((value - mean) ** 2 for value in values) / max(1, len(values) - 1)
            return mean, math.sqrt(variance / len(values))

        ms1, se1 = margin_and_se(species, discordant[0])
        ms2, se2 = margin_and_se(species, discordant[1])
        md1, sed1 = margin_and_se(discordant[0], discordant[1])
        md2, sed2 = margin_and_se(discordant[1], discordant[0])
        supported_d1 = ms1 < 0 and md1 > 0 and ms1 + critical_z * se1 < 0 and md1 - critical_z * sed1 > 0
        supported_d2 = ms2 < 0 and md2 > 0 and ms2 + critical_z * se2 < 0 and md2 - critical_z * sed2 > 0
        rows_out.append([record["condition"], record["replicate"], record["scenario"], record["rho"], record["quartet_id"], f"{ms1:.12g}", f"{ms2:.12g}", f"{md1:.12g}", f"{md2:.12g}", f"{se1:.12g}", f"{se2:.12g}", f"{sed1:.12g}", f"{sed2:.12g}", f"{critical_z:.12g}", str(supported_d1).lower(), str(supported_d2).lower(), str((supported_d1 or supported_d2)).lower(), str(ms1 < 0 or ms2 < 0).lower(), "independent_history_multiple_testing_gate"])
    atomic_write(output, ("\t".join(header) + "\n" + "\n".join("\t".join(x) for x in rows_out[1:]) + "\n").encode())
    output_manifest(task_manifest("validation", task_id, root), row["stage"], task_id, row, [output], len(records), len(records), columns=header, extra={"gate": "two pairwise dominant-discordant comparisons with Bonferroni-adjusted normal intervals", "exploratory_and_independent_outputs_separate": True})


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
