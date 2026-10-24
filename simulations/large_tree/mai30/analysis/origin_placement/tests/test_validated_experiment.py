from __future__ import annotations

import csv
import gzip
import json
from pathlib import Path


HERE = Path(__file__).resolve().parents[1]
RUN = HERE / "validated_100"


def rows(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_validated_panel_is_geometry_only_and_stratified():
    panel = rows(HERE / "selected_quartets.tsv")
    assert len(panel) == 96
    assert {row["scenario"] for row in panel} == {"A_existing_internal", "B_deep_internal"}
    counts = {}
    for row in panel:
        key = (row["condition"], row["replicate"], row["scenario"])
        counts[key] = counts.get(key, 0) + 1
        assert len(row["quartet"].split(",")) == 4
        assert len(set(row["quartet"].split(","))) == 4
    assert set(counts.values()) == {12}


def test_checkpoint_provenance_and_atomic_payload_contract():
    files = sorted(RUN.glob("*/history_*.json.gz"))
    assert files
    with gzip.open(files[0], "rt", encoding="utf-8") as handle:
        payload = json.load(handle)
    assert payload["simulator_commit"] == "bfba58d627f6f467af649cdb0b37f969349e0aab"
    assert payload["rho_grid"] == [0.0, 0.0001, 0.001, 0.002, 0.01, 0.02]
    assert len(payload["genealogies"]) == 12 * 6
    assert payload["records"]
    assert payload["sampled_arrangements"]


def test_incomplete_configurations_are_not_promoted_to_complete():
    manifest = rows(RUN / "forward_histories_manifest.tsv")
    assert manifest
    for row in manifest:
        completed = int(row["histories_completed"])
        requested = int(row["histories_requested"])
        assert row["status"] == ("complete" if completed == requested else "partial")


def test_reversal_gate_is_stricter_than_nominal_gate():
    reversal = rows(RUN / "quartet_reversal_validation.tsv")
    assert reversal
    assert all(row["multiple_comparison_supported"] in {"0", "1"} for row in reversal)
    # A nominal row cannot be silently treated as a validated reversal.
    assert not any(row["multiple_comparison_supported"] == "1" for row in reversal)
