#!/usr/bin/env python3
"""Diagnose weak Mai30 MSRC effects from retained pilot outputs.

This is a retrospective audit.  It does not simulate new genealogies.  The
pilot retains per-history frequency records, terminal states, and per-locus
switch counts, which are sufficient for the primary switching, persistence,
sorting, and conditional-quartet summaries.  Lineage-level event ordering and
branch-specific backward switch counts were not written by the pilot and are
reported as unavailable rather than reconstructed.
"""

from __future__ import annotations

import csv
import itertools
import json
import random
import statistics
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
MAI = ROOT / "simulations/large_tree/mai30"
PILOT = MAI / "datasets/processed/mai30_scaled_wf_runs"
OUT = MAI / "analysis/mechanism_diagnostics"
sys_path = MAI / "analysis/scripts"
import sys
sys.path.insert(0, str(sys_path))
from compute_original_baseline import parse_newick, quartet_key  # noqa: E402


GROUPS = {
    "msrc_no_suppression": 1.0,
    "msrc_moderate_suppression": 0.5,
    "msrc_strong_suppression": 0.1,
}
BASE_RATE = 0.02
PANEL_SEED = 20261022


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def group_parts(group: str) -> tuple[str, str, str]:
    parts = group.split("__")
    condition = parts[0]
    replicate = next(x.removeprefix("replicate_") for x in parts if x.startswith("replicate_"))
    treatment = next(x for x in parts if x.startswith("msrc_"))
    return condition, replicate, treatment


def make_panel(taxa: list[str], n: int = 12) -> list[tuple[str, str, str, str]]:
    combos = list(itertools.combinations(taxa, 4))
    rng = random.Random(PANEL_SEED)
    rng.shuffle(combos)
    return sorted(combos[:n])


def state_for(tree_line: str, quartet: tuple[str, str, str, str]) -> str:
    side = quartet_key(parse_newick(tree_line), quartet)
    a = "".join(sorted(side))
    b = "".join(sorted(set(quartet) - set(side)))
    return "|".join(sorted((a, b)))


def class_for(pattern: str, quartet: tuple[str, str, str, str], reference_set: frozenset[str]) -> str:
    n = pattern.count("1")
    if n in (0, 4):
        return "uniform_4_0"
    if n in (1, 3):
        return "3_1"
    ones = {taxon for taxon, bit in zip(quartet, pattern) if bit == "1"}
    quartet_set = set(quartet)
    reference_side = set(reference_set) & quartet_set
    # A 2:2 pattern is concordant when its two arrangement-1 taxa occupy one
    # side of the quartet's actual species-tree split.  This keeps the
    # classification correct for arbitrary taxon panels rather than assuming
    # the tuple order is always (1,2,3,4).
    if ones == reference_side or ones == quartet_set - reference_side:
        return "concordant_2_2"
    return "discordant_2_2"


def terminal_pattern(states: dict[str, int], quartet: tuple[str, str, str, str]) -> str:
    return "".join(str(states[taxon]) for taxon in quartet)


def forward_integrals(records: list[dict[str, str]], rho: float) -> tuple[float, float, float]:
    """Return branch duration, rho*integral(p), rho*integral(1-p)."""
    ordered = sorted(records, key=lambda r: float(r["absolute_age"]), reverse=True)
    duration = 0.0
    h01 = 0.0
    h10 = 0.0
    for current, nxt in zip(ordered, ordered[1:]):
        dt = max(0.0, float(current["absolute_age"]) - float(nxt["absolute_age"]))
        p = float(current["frequency_A1"])
        duration += dt
        h01 += rho * p * dt
        h10 += rho * (1.0 - p) * dt
    return duration, h01, h10


def main() -> None:
    manifest = read_tsv(MAI / "datasets/processed/mai30_scaled_wf_pilot/pilot_manifest_freq_0.5.tsv")
    switching_rows = []
    persistence_rows = []
    pattern_counts: Counter = Counter()
    conditional_counts: Counter = Counter()
    decomposition_counts: Counter = Counter()
    panel_by_taxa: dict[tuple[str, ...], list[tuple[str, str, str, str]]] = {}
    group_history_count = Counter()

    for summary in manifest:
        group = summary["group"]
        condition, replicate, treatment = group_parts(group)
        fraction = GROUPS[treatment]
        rho = BASE_RATE * fraction
        group_root = PILOT / group
        histories = sorted((group_root / "history_runs").glob("history_*"))
        if not histories:
            histories = sorted((PILOT / group).glob("history_*"))
        # The pilot's combined directory contains no history subdirectories;
        # use the authoritative run root for per-history diagnostics.
        run_group = PILOT.parent / "mai30_scaled_wf_runs" / group
        histories = sorted(run_group.glob("history_*"))
        if len(histories) != 25:
            raise ValueError(f"{group}: expected 25 retained histories, found {len(histories)}")
        species_tree_text = (histories[0] / "species_tree.nwk").read_text()
        species_tree = parse_newick(species_tree_text)
        for history in histories:
            group_history_count[group] += 1
            states = {r["taxon"]: int(r["sampled_state"]) for r in read_tsv(history / "terminal_states.tsv")}
            blocks = read_tsv(history / "genealogy_blocks.tsv")
            taxa = tuple(sorted(states))
            panel = panel_by_taxa.setdefault(taxa, make_panel(list(taxa)))
            switch_counts = [int(r["num_switch_events"]) for r in blocks if r["local_model"] == "msrc"]
            switching_rows.append({
                "group": group, "condition": condition, "replicate": replicate,
                "treatment": treatment, "rho": rho, "history": history.name,
                "n_genealogies": len(switch_counts),
                "mean_switches_per_genealogy": statistics.mean(switch_counts),
                "zero_switch_fraction": sum(x == 0 for x in switch_counts) / len(switch_counts),
                "max_switches": max(switch_counts),
                "switch_count_source": "genealogy_blocks.tsv; per-locus total only",
                "branch_specific_switches": "unavailable in retained pilot",
                "coalescence_order_relative_to_switch": "unavailable in retained pilot",
            })

            records = read_tsv(history / "arrangement_history.tsv")
            by_branch = defaultdict(list)
            for record in records:
                by_branch[record["branch_id"]].append(record)
            for branch, branch_records in sorted(by_branch.items()):
                duration, h01, h10 = forward_integrals(branch_records, rho)
                starts = [r for r in branch_records if r["is_branch_start"] == "True"]
                ends = [r for r in branch_records if r["is_branch_end"] == "True"]
                start_poly = sum(0 < int(r["copy_count_A1"]) < int(r["population_chromosomes"]) for r in starts)
                end_poly = sum(0 < int(r["copy_count_A1"]) < int(r["population_chromosomes"]) for r in ends)
                poly_records = sum(0 < int(r["copy_count_A1"]) < int(r["population_chromosomes"]) for r in branch_records)
                persistence_rows.append({
                    "group": group, "condition": condition, "replicate": replicate,
                    "treatment": treatment, "rho": rho, "history": history.name,
                    "branch_id": branch, "branch_duration": duration,
                    "polymorphic_at_branch_start": start_poly,
                    "polymorphic_at_branch_end": end_poly,
                    "fraction_recorded_generations_polymorphic": poly_records / len(branch_records),
                    "integrated_hazard_0_to_1": h01,
                    "integrated_hazard_1_to_0": h10,
                })

            reference_states = {quartet: quartet_key(species_tree, quartet) for quartet in panel}
            for block in blocks:
                parsed_gene = parse_newick(block["newick"])
                for quartet in panel:
                    actual_set = quartet_key(parsed_gene, quartet)
                    reference_set = reference_states[quartet]
                    actual = "|".join(sorted(("".join(sorted(actual_set)), "".join(sorted(set(quartet) - set(actual_set))))))
                    reference = "|".join(sorted(("".join(sorted(reference_set)), "".join(sorted(set(quartet) - set(reference_set))))))
                    pattern = terminal_pattern(states, quartet)
                    category = class_for(pattern, quartet, reference_set)
                    key = (group, ",".join(quartet), pattern, category)
                    pattern_counts[key] += 1
                    qlabel = "species_tree" if actual == reference else actual
                    conditional_counts[(group, ",".join(quartet), pattern, category, qlabel)] += 1
                    decomposition_counts[(group, ",".join(quartet), pattern, category, qlabel)] += 1

    write_tsv(OUT / "switching_diagnostics.tsv", switching_rows)
    write_tsv(OUT / "persistence_probabilities.tsv", persistence_rows)

    pattern_rows = []
    for (group, quartet, pattern, category), count in sorted(pattern_counts.items()):
        total = sum(v for (g, q, p, c), v in pattern_counts.items() if g == group and q == quartet)
        pattern_rows.append({"group": group, "quartet": quartet, "pattern": pattern, "category": category, "count": count, "probability": count / total, "unit": "genealogy; terminal state is history-shared"})
    write_tsv(OUT / "arrangement_pattern_probabilities.tsv", pattern_rows)

    cond_rows = []
    for key, count in sorted(conditional_counts.items()):
        group, quartet, pattern, category, qlabel = key
        denom = sum(v for (g, q, p, c, _), v in conditional_counts.items() if (g, q, p, c) == (group, quartet, pattern, category))
        cond_rows.append({"group": group, "quartet": quartet, "pattern": pattern, "category": category, "topology": qlabel, "count": count, "conditional_probability": count / denom, "n_pattern": denom})
    write_tsv(OUT / "conditional_quartet_probabilities.tsv", cond_rows)

    decomp_rows = []
    for (group, quartet, pattern, category, qlabel), count in sorted(decomposition_counts.items()):
        total = sum(v for (g, q, _p, _c, _label), v in decomposition_counts.items() if g == group and q == quartet)
        pz = sum(v for (g, q, p, c, _label), v in decomposition_counts.items() if (g, q, p, c) == (group, quartet, pattern, category)) / total
        decomp_rows.append({"group": group, "quartet": quartet, "pattern": pattern, "category": category, "topology": qlabel, "joint_probability": count / total, "pattern_probability": pz, "conditional_probability": (count / total) / pz if pz else ""})
    write_tsv(OUT / "marginal_quartet_decomposition.tsv", decomp_rows)

    # The expanded rate study is a predeclared design artifact only.  It is
    # deliberately not labeled as executed data.
    grid = []
    for condition, replicate in [("outgroup.0.species.1.5.genes.1.5", "007"), ("outgroup.0.species.1.5.genes.1.5", "042"), ("outgroup.0.species.5.genes.5", "005"), ("outgroup.0.species.1.5.genes.1.5", "074")]:
        for rho_value in [0.0, 0.0002, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1]:
            grid.append({"condition": condition, "replicate": replicate, "rho": rho_value, "structured_duration_source": "existing tree/root extension", "status": "planned_not_executed", "rationale": "log grid includes current rates 0.002, 0.01, 0.02 and higher-hazard regimes"})
    write_tsv(OUT / "switching_rate_grid.tsv", grid)
    write_tsv(OUT / "parameter_sensitivity.tsv", [{"parameter": x, "status": "planned_not_executed", "reason": "requires new authoritative simulations"} for x in ["initial_frequency", "structured_duration", "population_size", "origin_placement", "internal_branch_length"]])

    # Existing ASTRAL pilot results are copied into a clearly labeled audit
    # table, not treated as results from the planned rate grid.
    astral = read_tsv(MAI / "analysis/results/mai30_scaled_wf_astral.tsv")
    for row in astral:
        row["analysis_status"] = "completed_pilot_only"
        row["rate_grid_status"] = "not_available"
    write_tsv(OUT / "astral_results.tsv", astral)
    (OUT / "manifest.json").write_text(json.dumps({"status": "retrospective_pilot_diagnostics_plus_planned_grid", "pilot_groups": len(manifest), "pilot_histories_per_group": 25, "pilot_genealogies_per_group": 500, "base_rate": BASE_RATE, "pattern_classification": "2:2 patterns classified against each quartet's induced species-tree split", "simulator_outputs_used": ["arrangement_history.tsv", "terminal_states.tsv", "genealogy_blocks.tsv"], "unavailable_from_retained_outputs": ["lineage-level switch time and branch", "coalescence-before/after-switch ordering", "backward time spent in each background"], "expanded_rate_grid_executed": False}, indent=2) + "\n")
    print(f"wrote diagnostics under {OUT}")


if __name__ == "__main__":
    main()
