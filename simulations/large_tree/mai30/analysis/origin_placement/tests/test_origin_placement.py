from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path


HERE = Path(__file__).resolve().parents[1]


def read(name):
    with (HERE / name).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_four_scenarios_and_all_quartet_exposure_rows():
    scenarios = read("origin_scenarios.tsv")
    exposure = read("quartet_exposure.tsv")
    assert len(scenarios) == 16
    assert {row["scenario"] for row in scenarios} == {"A_existing_internal", "B_deep_internal", "C_root_adjacent", "D_older_root_extended"}
    assert len(exposure) == 16
    assert {int(row["total_quartets"]) for row in exposure} == {27405}
    assert all(0 <= float(row["fraction_structurally_exposed"]) <= 1 for row in exposure)


def test_root_scenarios_expose_all_quartets():
    rows = [row for row in read("quartet_exposure.tsv") if row["scenario"].startswith(("C_", "D_"))]
    assert rows and all(int(row["quartets_with_at_least_two_descendants"]) == 27405 for row in rows)


def test_forward_persistence_rows_are_explicitly_unconditional():
    rows = read("persistence_results.tsv")
    assert len(rows) == 16
    executed = [row for row in rows if row["status"] == "executed_forward_history"]
    assert len(executed) == 16
    for row in executed:
        assert 0 <= float(row["origin_branch_polymorphic_fraction"]) <= 1
        assert row["history"] == "0"


def test_genealogy_and_astral_outputs_do_not_claim_execution():
    for name in ("quartet_margins.tsv", "quartet_reversal_results.tsv", "astral_results.tsv"):
        rows = read(name)
        assert rows and all(row["status"] == "not_executed" for row in rows)


def test_rho_stress_points_are_preserved_as_planned():
    rows = read("parameter_grid.tsv")
    assert {float(row["rho"]) for row in rows} == {0.0, 0.0001, 0.001, 0.002, 0.01, 0.02, 0.1, 1.0}
    assert any(float(row["rho"]) == 1.0 and row["status"] == "planned_not_executed" for row in rows)
