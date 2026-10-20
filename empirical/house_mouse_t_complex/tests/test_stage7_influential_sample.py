import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "empirical/house_mouse_t_complex/results/stage7_influential_sample"


def read(name):
    with (OUT / name).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_stage4_reproduction_and_consistent_loo_tracks():
    rows = read("window_original_vs_leave_one_out.tsv")
    assert len(rows) == 4046
    assert max(abs(float(r["B_all"]) - float(r["B_source_stage4"])) for r in rows) == 0
    assert max(abs(float(r["A_distance_original"]) - float(r["A_source_stage4"])) for r in rows) < 1e-12
    assert max(abs(float(r["A_relative_original"]) - float(r["A_relative_source_stage4"])) for r in rows) < 1e-12
    for row in rows[::101]:
        assert abs(float(row["A_relative_leave_one_out"]) - float(row["A_distance_leave_one_out"]) / float(row["B_leave_one_out"])) < 1e-12


def test_sample_influence_and_branch_validity_are_recorded():
    audit = read("sample_branch_length_audit.tsv")
    assert len(audit) == 4046
    assert all(row["max_terminal_tip"] == "Mmd_FRA4_18B.variant" for row in audit if float(row["midpoint_bp"]) >= 14737499.5)
    assert all(int(row["n_negative_branches"]) == 0 for row in audit)
    assert all(int(row["n_nonfinite_branches"]) == 0 for row in audit)
    summary = {(r["region"], r["track"]): r for r in read("influence_summary.tsv")}
    assert float(summary[("downstream", "B_all")]["mean"]) > 0.3
    assert float(summary[("downstream", "B_leave_one_out")]["mean"]) < 0.01


def test_leave_one_out_quartets_are_normalized():
    rows = read("quartet_leave_one_out_validation.tsv")
    assert len(rows) == 4046
    for row in rows[::97]:
        assert abs(float(row["sum_original"]) - 1) < 1e-12
        assert abs(float(row["sum_leave_one_out"]) - 1) < 1e-12
        assert float(row["max_original_source_abs_diff"]) < 1e-9
        assert math.isfinite(float(row["D_leave_one_out"]))


def test_boundary_sensitivity_keeps_single_and_multi_models_separate():
    rows = read("boundary_detection_sensitivity.tsv")
    assert {row["model"] for row in rows} == {"single_interval", "multi_interval_0_to_4"}
    assert any(row["method"] == "original_B" and float(row["iou_union"]) > 0.3 for row in rows)
    assert any(row["method"] == "leave_one_out_B" for row in rows)
