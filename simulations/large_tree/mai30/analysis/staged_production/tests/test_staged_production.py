from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[3]
STAGED = ROOT / "analysis/staged_production"


def rows(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_staged_selection_has_ten_replicates_per_condition():
    selected = rows(STAGED / "selected_replicates.tsv")
    assert len(selected) == 30
    counts = {}
    for row in selected:
        counts[row["condition"]] = counts.get(row["condition"], 0) + 1
    assert set(counts.values()) == {10}
    assert len({(r["condition"], r["replicate"]) for r in selected}) == 30


def test_staged_manifest_has_three_treatments_per_selected_replicate():
    selected = {(r["condition"], r["replicate"]) for r in rows(STAGED / "selected_replicates.tsv")}
    manifest = rows(STAGED / "staged_manifest.tsv")
    assert len(manifest) == 90
    assert {(r["condition"], r["replicate"]) for r in manifest} == selected
    for key in selected:
        treatments = [r["treatment"] for r in manifest if (r["condition"], r["replicate"]) == key]
        assert sorted(treatments) == ["msrc_moderate_suppression", "msrc_no_suppression", "msrc_strong_suppression"]
    assert {r["output_root"] for r in manifest} == {"datasets/processed/mai30_staged_runs"}


def test_scaled_tree_is_not_raw_generation_tree():
    script_dir = ROOT / "analysis/scripts"
    import sys
    sys.path.insert(0, str(script_dir))
    from run_astral_scaled_pilot import parse_lengths_correct

    condition = "outgroup.0.species.1.5.genes.1.5"
    rep = "042"
    raw = ROOT / "datasets/raw" / condition / rep / "s_tree.trees"
    scaled_run = next((ROOT / "datasets/processed/mai30_scaled_wf_runs").glob(
        f"{condition}__replicate_{rep}__*__freq_0.5/history_000"
    ))
    raw_lengths = parse_lengths_correct(raw.read_text())[2]
    scaled_lengths = parse_lengths_correct((scaled_run / "species_tree.nwk").read_text())[2]
    common = set(raw_lengths) & set(scaled_lengths)
    assert common
    ratios = [raw_lengths[s] / scaled_lengths[s] for s in common if scaled_lengths[s] > 0]
    assert ratios
    assert all(ratio == pytest.approx(1000.0, rel=1e-9) for ratio in ratios)

    config = json.loads((scaled_run / "simulation_config.json").read_text())
    assert config["species_tree"]["default_effective_population_size"] < 100000


def test_corrected_pilot_branch_table_has_scaled_lengths():
    branch_rows = rows(STAGED / "results/pilot_branch_lengths.tsv")
    rep042 = [float(r["true_culength"]) for r in branch_rows if r["replicate"] == "042"]
    assert rep042
    assert max(rep042) < 200.0
    assert max(rep042) > 100.0
