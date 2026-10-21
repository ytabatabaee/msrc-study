from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "datasets/processed/mai30_production_manifest.tsv"


def test_production_manifest_has_three_conditions_and_900_tasks():
    with MANIFEST.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 900
    assert len({row["task_id"] for row in rows}) == 900
    assert Counter(row["condition"] for row in rows) == {
        "outgroup.0.species.0.15.genes.0.15": 300,
        "outgroup.0.species.1.5.genes.1.5": 300,
        "outgroup.0.species.5.genes.5": 300,
    }
    assert Counter(row["treatment"] for row in rows) == {
        "msrc_no_suppression": 300,
        "msrc_moderate_suppression": 300,
        "msrc_strong_suppression": 300,
    }


def test_manifest_paths_and_parameters_are_reproducible():
    with MANIFEST.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert all(not row["species_tree_path"].startswith("/") for row in rows)
    assert all(row["initial_frequency"] == "0.5" for row in rows)
    assert all(row["scale_factor"] == "0.001" for row in rows)
    assert all(row["n_histories"] == "25" and row["loci_per_history"] == "20" for row in rows)
    assert all(row["status"] == "planned" for row in rows)
