from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "empirical/house_mouse_t_complex/results/stage6_distance_discontinuity"
DIST = ROOT / "empirical/house_mouse_t_complex/results/stage4_distance/stage4_window_distance.tsv"


def read(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_recomputed_distance_and_relative_track_match_frozen_values():
    diagnostics = {row["locus_id"]: row for row in read(OUT / "window_tree_diagnostics.tsv")}
    source = {row["locus_id"]: row for row in read(DIST)}
    assert len(diagnostics) == len(source) == 4046
    assert max(abs(float(row["B_difference"])) for row in diagnostics.values()) == 0.0
    for locus, row in diagnostics.items():
        assert float(row["A_relative"]) == float(source[locus]["A_relative"])
        assert abs(float(row["A_distance"]) - float(source[locus]["A_distance"])) < 1e-15


def test_coordinate_order_and_change_point_are_data_driven():
    diagnostics = read(OUT / "window_tree_diagnostics.tsv")
    starts = [int(row["start_bp"]) for row in diagnostics]
    assert starts == sorted(starts)
    top = read(OUT / "change_point_diagnostics.tsv")[0]
    assert abs(float(top["coordinate_bp"]) - 14737499.5) < 1
    assert float(top["fold_right_over_left"]) > 30


def test_branch_validity_and_outlier_sensitivity_are_recorded():
    diagnostics = read(OUT / "window_tree_diagnostics.tsv")
    assert all(int(row["negative_length_branches"]) == 0 for row in diagnostics)
    assert all(int(row["nonfinite_branches"]) == 0 for row in diagnostics)
    downstream = [row for row in diagnostics if float(row["midpoint_bp"]) >= 14737499.5]
    assert all(row["max_terminal_tip"] == "Mmd_FRA4_18B.variant" for row in downstream)
    comparison = {row["region"]: row for row in read(OUT / "upstream_downstream_comparison.tsv")}
    assert float(comparison["downstream"]["B_excluding_Mmd_FRA4_18B_mean"]) < 0.01

