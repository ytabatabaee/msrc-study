#!/usr/bin/env python3
"""Run a compact origin-placement discovery experiment with msrc-sim.

The script uses the authoritative simulator API.  One forward Wright--Fisher
history is generated per tree/scenario/history and reused across all rho
values, while genealogy RNG streams remain treatment-specific.  It writes
aggregated diagnostics only; no large genealogy collections are retained.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
import random
import sys
from collections import Counter, defaultdict
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


TREES = [
    ("outgroup.0.species.1.5.genes.1.5", "007"),
    ("outgroup.0.species.1.5.genes.1.5", "042"),
    ("outgroup.0.species.1.5.genes.1.5", "074"),
    ("outgroup.0.species.5.genes.5", "005"),
]
RHO_GRID = [0.0, 0.0001, 0.001, 0.002, 0.01, 0.02, 0.1, 1.0]
N_HISTORIES = 1
N_LOCI = 1
PANEL_SIZE = 120
SEED = 20261022
GENEALOGY_SCENARIOS = {"A_existing_internal", "B_deep_internal", "C_root_adjacent"}


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        if not rows:
            handle.write("")
            return
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def load_base(condition: str, replicate: str) -> dict:
    path = MAI / "datasets/configs/mai30_msrc_pilot" / f"{condition}__replicate_{replicate}__msrc_no_suppression__freq_0.5.yaml"
    return yaml.safe_load(path.read_text())


def scenario_rows(condition: str, replicate: str, config: dict, tree: SpeciesTree) -> list[dict]:
    rr = config["rearrangement"]
    original_branch = rr["origin_branch"]
    original_time = int(rr["origin_time_from_branch_start"])
    candidates = []
    for branch in tree.branches.values():
        if branch.branch_id == "ROOT":
            continue
        descendants = tree.descendants(branch.child_node)
        duration = branch.older_age - branch.younger_age
        if 4 <= len(descendants) <= len(tree.taxa) - 2 and duration >= 2:
            candidates.append((branch.older_age, len(descendants), branch.branch_id, duration))
    if not candidates:
        raise ValueError(f"no deep internal branch for {condition}/{replicate}")
    _, _, deep_branch, deep_duration = max(candidates)
    root_duration = tree.root_extension
    root_branch = tree.root.name
    rows = [
        {"scenario": "A_existing_internal", "origin_branch": original_branch, "origin_time_from_branch_start": original_time, "root_extension": tree.root_extension, "rule": "original pilot configuration"},
        {"scenario": "B_deep_internal", "origin_branch": deep_branch, "origin_time_from_branch_start": max(1, int(deep_duration - min(50, max(1, deep_duration - 1)))), "root_extension": tree.root_extension, "rule": "oldest non-root internal branch, with approximately 50 generations of structured exposure when possible"},
        {"scenario": "C_root_adjacent", "origin_branch": root_branch, "origin_time_from_branch_start": max(1, int(root_duration - 50)), "root_extension": tree.root_extension, "rule": "root-extension origin 50 generations before root speciation"},
        {"scenario": "D_older_root_extended", "origin_branch": root_branch, "origin_time_from_branch_start": 2, "root_extension": 2 * tree.root_extension, "rule": "root-extension origin near older boundary with doubled root extension"},
    ]
    for row in rows:
        branch = tree.branches[row["origin_branch"]]
        origin_age = branch.older_age - row["origin_time_from_branch_start"]
        is_root = row["origin_branch"] == root_branch
        row.update({
            "condition": condition,
            "replicate": replicate,
            "n_taxa": len(tree.taxa),
            "descendant_taxa": len(tree.descendants(branch.child_node)) if not is_root else len(tree.taxa),
            "origin_younger_age": branch.younger_age,
            "origin_older_age": branch.older_age if not is_root else tree.root.age + row["root_extension"],
            "origin_age": origin_age if not is_root else tree.root.age + row["root_extension"] - row["origin_time_from_branch_start"],
            "structured_duration_to_root_or_branch_end": row["origin_time_from_branch_start"],
        })
    return rows


def choose_panel(taxa: tuple[str, ...]) -> list[tuple[str, str, str, str]]:
    combos = list(itertools.combinations(taxa, 4))
    rng = random.Random(SEED)
    rng.shuffle(combos)
    return sorted(combos[:PANEL_SIZE])


def executable_rhos(scenario: str) -> list[float]:
    # rho=1 is retained in parameter_grid.tsv as a requested feasibility
    # point.  The doubled-root scenario is intentionally capped at rho=.02;
    # its long structured interval makes rho=.1 and 1.0 an event-count stress
    # test rather than a useful discovery setting.
    return [rho for rho in RHO_GRID if rho <= 0.02]


def classify_pattern(pattern: str, quartet: tuple[str, ...], species_side: frozenset[str]) -> str:
    n = pattern.count("1")
    if n in (0, 4):
        return "uniform_4_0"
    if n in (1, 3):
        return "3_1"
    ones = {taxon for taxon, bit in zip(quartet, pattern) if bit == "1"}
    qset = set(quartet)
    side = set(species_side) & qset
    return "concordant_2_2" if ones in (side, qset - side) else "discordant_2_2"


def exposure(tree: SpeciesTree, scenario: dict) -> dict:
    descendants = set(tree.descendants(scenario["origin_branch"]))
    n = len(tree.taxa)
    total = math.comb(n, 4)
    exposed = 0
    for q in itertools.combinations(tree.taxa, 4):
        if len(descendants & set(q)) >= 2:
            exposed += 1
    return {"condition": scenario["condition"], "replicate": scenario["replicate"], "scenario": scenario["scenario"], "n_taxa": n, "descendant_taxa": len(descendants), "total_quartets": total, "quartets_with_at_least_two_descendants": exposed, "fraction_structurally_exposed": exposed / total}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    scenario_rows_all = []
    exposure_rows = []
    for condition, replicate in TREES:
        config = load_base(condition, replicate)
        tree = SpeciesTree(config["species_tree"]["newick"], config["species_tree"]["default_effective_population_size"], config["species_tree"]["root_extension"], ultrametric_tolerance=1e-6)
        rows = scenario_rows(condition, replicate, config, tree)
        for row in rows:
            scenario_rows_all.append(row)
            exposure_rows.append(exposure(tree, row))
    write_tsv(OUT / "origin_scenarios.tsv", scenario_rows_all)
    write_tsv(OUT / "quartet_exposure.tsv", exposure_rows)

    parameter_rows = []
    for scenario in scenario_rows_all:
        for rho in RHO_GRID:
            status = "executed" if scenario["scenario"] in GENEALOGY_SCENARIOS and rho in executable_rhos(scenario["scenario"]) else "planned_not_executed"
            parameter_rows.append({**{k: scenario[k] for k in ("condition", "replicate", "scenario", "origin_branch", "origin_age", "descendant_taxa")}, "rho": rho, "n_histories": N_HISTORIES, "loci_per_history": N_LOCI, "status": status})
    write_tsv(OUT / "parameter_grid.tsv", parameter_rows)

    persistence = []
    pattern_counts = Counter()
    quartet_counts = Counter()
    history_quartet_counts = Counter()
    switching = Counter()
    scenario_by_key = {(r["condition"], r["replicate"], r["scenario"]): r for r in scenario_rows_all}

    for tree_index, (condition, replicate) in enumerate(TREES):
        base = load_base(condition, replicate)
        base_n = int(base["species_tree"]["default_effective_population_size"])
        original_tree = SpeciesTree(base["species_tree"]["newick"], base_n, base["species_tree"]["root_extension"], ultrametric_tolerance=1e-6)
        panel = choose_panel(original_tree.taxa)
        reference = {q: quartet_key(parse_newick(base["species_tree"]["newick"]), q) for q in panel}
        for scenario_index, scenario in enumerate([r for r in scenario_rows_all if r["condition"] == condition and r["replicate"] == replicate]):
            if scenario["scenario"] not in GENEALOGY_SCENARIOS:
                persistence.append({"condition": condition, "replicate": replicate, "scenario": scenario["scenario"], "history": "", "origin_branch": scenario["origin_branch"], "origin_age": scenario["origin_age"], "status": "not_executed_long_root_stress_case"})
                continue
            tree = SpeciesTree(base["species_tree"]["newick"], base_n, scenario["root_extension"], ultrametric_tolerance=1e-6)
            rr = Rearrangement("origin_discovery", "inversion", scenario["origin_branch"], int(scenario["origin_time_from_branch_start"]), int(base["rearrangement"]["initial_copy_count"]), 0.0)
            for history_index in range(N_HISTORIES):
                seed = SEED + tree_index * 1_000_000 + scenario_index * 10_000 + history_index
                history = simulate_population_history(tree, rr, np.random.default_rng(seed), base)
                origin_records = history.by_branch[scenario["origin_branch"]]
                persistence.append({
                    "condition": condition, "replicate": replicate, "scenario": scenario["scenario"], "history": history_index,
                    "origin_branch": scenario["origin_branch"], "origin_age": scenario["origin_age"],
                    "origin_branch_polymorphic_fraction": sum(0 < r.copy_count_A1 < r.population_chromosomes for r in origin_records) / len(origin_records),
                    "origin_branch_end_polymorphic": int(0 < origin_records[-1].copy_count_A1 < origin_records[-1].population_chromosomes),
                    "origin_branch_ever_fixed": int(any(r.status == "fixed" for r in origin_records)),
                    "origin_branch_ever_lost": int(any(r.status == "lost" for r in origin_records)),
                    "root_terminal_polymorphic": int(0 < history.frequency_at(tree.root.name, tree.root.age + tree.root_extension) < 1),
                    "status": "executed",
                })
                sampled = {taxon: int(np.random.default_rng(seed + 77).random() < history.terminal_frequency(taxon)) for taxon in tree.taxa}
                for rho_index, rho in enumerate(executable_rhos(scenario["scenario"])):
                    for locus in range(N_LOCI):
                        rng = np.random.default_rng(seed + 1000 + rho_index * 100 + locus)
                        result = simulate_genealogy(locus, tree, history, sampled, rho, 1.0, rng, record_events=True)
                        gene = parse_newick(result.newick)
                        switching[(condition, replicate, scenario["scenario"], rho, history_index)] += sum(e.event_type == "switch" for e in result.events)
                        for quartet in panel:
                            pattern = "".join(str(sampled[t]) for t in quartet)
                            category = classify_pattern(pattern, quartet, reference[quartet])
                            actual = quartet_key(gene, quartet)
                            truth = reference[quartet]
                            qlabel = "species_tree" if actual == truth else "discordant"
                            pattern_counts[(condition, replicate, scenario["scenario"], rho, pattern, category)] += 1
                            quartet_counts[(condition, replicate, scenario["scenario"], rho, ",".join(quartet), qlabel)] += 1
                            history_quartet_counts[(condition, replicate, scenario["scenario"], rho, history_index, ",".join(quartet), qlabel)] += 1

    write_tsv(OUT / "persistence_results.tsv", persistence)

    pattern_rows = []
    for key, count in sorted(pattern_counts.items()):
        condition, replicate, scenario, rho, pattern, category = key
        denom = sum(v for (c, r, s, rr, _p, _cat), v in pattern_counts.items() if (c, r, s, rr) == (condition, replicate, scenario, rho))
        pattern_rows.append({"condition": condition, "replicate": replicate, "scenario": scenario, "rho": rho, "pattern": pattern, "category": category, "count": count, "probability": count / denom})
    write_tsv(OUT / "arrangement_patterns.tsv", pattern_rows)

    margin_rows = []
    reversal_rows = []
    for key in sorted({k[:5] for k in quartet_counts}):
        condition, replicate, scenario, rho, quartet = key
        counts = {label: value for (c, r, s, rr, q, label), value in quartet_counts.items() if (c, r, s, rr, q) == key}
        total = sum(counts.values())
        qs = {label: counts.get(label, 0) / total for label in ("species_tree", "discordant")}
        margin = qs["species_tree"] - qs["discordant"]
        history_margins = []
        for h in range(N_HISTORIES):
            hc = {label: value for (c, r, s, rr, hh, q, label), value in history_quartet_counts.items() if (c, r, s, rr, hh, q) == (*key, h)}
            ht = sum(hc.values())
            if ht:
                history_margins.append(hc.get("species_tree", 0) / ht - hc.get("discordant", 0) / ht)
        se = float(np.std(history_margins, ddof=1) / math.sqrt(len(history_margins))) if len(history_margins) > 1 else float("nan")
        margin_rows.append({"condition": condition, "replicate": replicate, "scenario": scenario, "rho": rho, "quartet": quartet, "p_species_tree": qs["species_tree"], "p_discordant": qs["discordant"], "quartet_margin": margin, "history_margin_se": se, "n_genealogies": total})
    write_tsv(OUT / "quartet_margins.tsv", margin_rows)

    for key in sorted({(r["condition"], r["replicate"], r["scenario"], r["rho"]) for r in margin_rows}):
        rows = [r for r in margin_rows if (r["condition"], r["replicate"], r["scenario"], r["rho"]) == key]
        negative = [r for r in rows if float(r["quartet_margin"]) < 0]
        supported = [r for r in rows if math.isfinite(float(r["history_margin_se"])) and float(r["quartet_margin"]) + 1.96 * float(r["history_margin_se"]) < 0]
        reversal_rows.append({"condition": key[0], "replicate": key[1], "scenario": key[2], "rho": key[3], "n_panel_quartets": len(rows), "fraction_negative_margin": len(negative) / len(rows), "fraction_approx_supported_reversal": len(supported) / len(rows), "mean_margin": statistics_mean(float(r["quartet_margin"]) for r in rows), "status": "discovery_only"})
    write_tsv(OUT / "quartet_reversal_results.tsv", reversal_rows)

    switch_rows = []
    for key, value in sorted(switching.items()):
        condition, replicate, scenario, rho, history = key
        switch_rows.append({"condition": condition, "replicate": replicate, "scenario": scenario, "rho": rho, "history": history, "n_genealogies": N_LOCI, "mean_switches": value / N_LOCI})
    write_tsv(OUT / "switching_results.tsv", switch_rows)
    write_tsv(OUT / "astral_results.tsv", [{"status": "not_executed", "reason": "discovery run did not meet a prespecified supported-reversal or large-margin-shift criterion; no full gene-tree collections generated"}])
    (OUT / "manifest.json").write_text(json.dumps({"status": "discovery_executed", "simulator": str(MSRC_ROOT), "trees": TREES, "rho_grid": RHO_GRID, "executed_rho_rule": "rho<=0.02; rho=.1 and 1.0 retained as planned stress points", "histories_per_configuration": N_HISTORIES, "loci_per_history": N_LOCI, "panel_size": PANEL_SIZE, "all_27405_quartets_structural_exposure": True, "all_27405_quartets_genealogy_analysis": False, "forward_histories_reused_across_rho": True, "validation_stage": "not_executed", "rng_seed": SEED}, indent=2) + "\n")


def statistics_mean(values):
    values = list(values)
    return sum(values) / len(values) if values else float("nan")


if __name__ == "__main__":
    main()
