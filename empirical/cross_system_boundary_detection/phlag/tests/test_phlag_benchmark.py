import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts/14_phlag_cross_system_benchmark.py"
spec = importlib.util.spec_from_file_location("phlag_benchmark", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def test_parse_finite_and_invalid_outputs(tmp_path):
    p = tmp_path / "out.txt"
    p.write_text("# header\n0,1,0\n0.2,0.7,0.1\n")
    state, posterior, status = mod.parse_phlag_output(p)
    assert status == "valid"
    assert state == [0, 1, 0]
    assert posterior == [0.2, 0.7, 0.1]
    p.write_text("0,0\nnan,nan\n")
    _, _, status = mod.parse_phlag_output(p)
    assert status == "nonfinite_posterior"


def test_taxon_mapping_detects_mismatch():
    assert mod.read_tree_tips("(a:1,b:1);") == {"a", "b"}
    assert mod.read_tree_tips("(a:1,c:1);") != {"a", "b"}


def test_union_branch_aggregation_uses_max_posterior():
    rows = [
        {"window_id": "w1", "posterior_anomaly": "0.2", "predicted_state": "0"},
        {"window_id": "w1", "posterior_anomaly": "0.8", "predicted_state": "1"},
        {"window_id": "w2", "posterior_anomaly": "nan", "predicted_state": "0"},
    ]
    result = mod.aggregate_branch_predictions(rows)
    assert result["w1"]["posterior_anomaly_max"] == 0.8
    assert result["w1"]["predicted_state_union"] == 1
    assert result["w2"]["predicted_state_union"] == 0


def test_interval_iou_and_numeric_classification():
    assert mod.interval_iou((0, 10), (5, 15)) == 5 / 15
    assert mod.interval_iou((0, 10), (20, 30)) == 0
    assert mod.classify_run([0, 0], [0.2, 0.3], "valid") == "valid_no_anomaly"
    assert mod.classify_run([0, 1], [float("nan"), 0.3], "nonfinite_posterior") == "degenerate_or_invalid"
