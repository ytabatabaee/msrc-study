from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
import run_pilot_diagnostics as diag  # noqa: E402


def read(name: str) -> list[dict[str, str]]:
    with (HERE / name).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_effective_switching_rates_are_the_declared_mapping():
    assert diag.BASE_RATE * diag.GROUPS["msrc_no_suppression"] == 0.02
    assert diag.BASE_RATE * diag.GROUPS["msrc_moderate_suppression"] == 0.01
    assert diag.BASE_RATE * diag.GROUPS["msrc_strong_suppression"] == 0.002


def test_actual_switching_is_ordered_across_pilot_treatments():
    rows = read("switching_diagnostics.tsv")
    by_treatment = defaultdict(list)
    for row in rows:
        by_treatment[row["treatment"]].append(float(row["mean_switches_per_genealogy"]))
    means = {key: sum(values) / len(values) for key, values in by_treatment.items()}
    assert means["msrc_no_suppression"] > means["msrc_moderate_suppression"] > means["msrc_strong_suppression"]
    assert len(rows) == 300


def test_pattern_probabilities_sum_to_one_per_group_and_quartet():
    rows = read("arrangement_pattern_probabilities.tsv")
    totals = defaultdict(float)
    for row in rows:
        totals[(row["group"], row["quartet"])] += float(row["probability"])
    assert totals
    assert all(abs(value - 1.0) < 1e-12 for value in totals.values())


def test_conditional_topology_probabilities_sum_to_one_when_pattern_observed():
    rows = read("conditional_quartet_probabilities.tsv")
    totals = defaultdict(float)
    for row in rows:
        totals[(row["group"], row["quartet"], row["pattern"])] += float(row["conditional_probability"])
    assert totals
    assert all(abs(value - 1.0) < 1e-12 for value in totals.values())


def test_marginal_decomposition_is_normalized_and_reconstructs_pattern_weights():
    rows = read("marginal_quartet_decomposition.tsv")
    joint_totals = defaultdict(float)
    pattern_weights = {}
    for row in rows:
        key = (row["group"], row["quartet"])
        joint_totals[key] += float(row["joint_probability"])
        pattern_weights.setdefault((row["group"], row["quartet"], row["pattern"]), float(row["pattern_probability"]))
    assert all(abs(value - 1.0) < 1e-12 for value in joint_totals.values())
    by_pattern = defaultdict(float)
    for (group, quartet, pattern), value in pattern_weights.items():
        by_pattern[(group, quartet)] += value
    assert all(abs(value - 1.0) < 1e-12 for value in by_pattern.values())


def test_persistence_probabilities_are_valid():
    rows = read("persistence_probabilities.tsv")
    assert len(rows) == 17700
    for row in rows:
        for field in ("polymorphic_at_branch_start", "polymorphic_at_branch_end", "fraction_recorded_generations_polymorphic"):
            value = float(row[field])
            assert 0.0 <= value <= 1.0


def test_expanded_rate_grid_is_design_only_and_staged_manifest_is_unchanged():
    grid = read("switching_rate_grid.tsv")
    assert grid and {row["status"] for row in grid} == {"planned_not_executed"}
    staged = Path(__file__).resolve().parents[2] / "staged_production/staged_manifest.tsv"
    assert staged.exists()
    assert sum(1 for _ in staged.open()) == 91


def test_species_split_classification_does_not_depend_on_tuple_order():
    quartet = ("a", "b", "c", "d")
    split = frozenset({"a", "b"})
    assert diag.class_for("1100", quartet, split) == "concordant_2_2"
    assert diag.class_for("0011", quartet, split) == "concordant_2_2"
    assert diag.class_for("1010", quartet, split) == "discordant_2_2"
    permuted = ("a", "c", "b", "d")
    assert diag.class_for("1010", permuted, split) == "concordant_2_2"
