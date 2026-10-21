import importlib.util
from pathlib import Path

BASE = Path(__file__).parents[1]
SPEC = importlib.util.spec_from_file_location("stage13d", BASE / "scripts/complete_phlag_benchmark.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_five_reference_intervals_and_right_censoring():
    refs = MODULE.inversion_intervals()
    assert [r["name"] for r in refs] == ["Inv1", "Inv2", "Inv3", "Inv4", "Inv5"]
    assert MODULE.observable_reference(refs[-1]) == (19_401_000, MODULE.OBS_END)


def test_disconnected_union_iou_is_not_hull_iou():
    metrics = MODULE.score([(0, 9), (100, 109)], [(0, 9), (100, 109)])
    assert metrics["iou"] == 1.0
    broad = MODULE.score([(0, 109)], [(0, 9), (100, 109)])
    assert broad["iou"] < 1.0


def test_adjacent_and_overlapping_intervals_merge():
    assert MODULE.merge_intervals([(0, 9), (10, 19), (19, 30)], max_gap=0) == [(0, 30)]


def test_phlag_parser_rejects_length_mismatch(tmp_path):
    output = tmp_path / "phlag.txt"
    output.write_text("0,1\n0.1,0.9\n")
    try:
        MODULE.parse_phlag_output(output, 3)
    except ValueError as exc:
        assert "length mismatch" in str(exc)
    else:
        raise AssertionError("length mismatch must fail loudly")


def test_phlag_parser_accepts_finite_sequences(tmp_path):
    output = tmp_path / "phlag.txt"
    output.write_text("# header\n0,1,0\n0.1,0.9,0.2\n")
    states, probs, status = MODULE.parse_phlag_output(output, 3)
    assert states == ["0", "1", "0"]
    assert probs == [0.1, 0.9, 0.2]
    assert status == "valid_informative"


def test_manifest_paths_are_relative_after_score_only():
    manifest = BASE / "run_manifest.tsv"
    if not manifest.exists():
        return
    import csv
    with manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    assert all(not Path(row["output"]).is_absolute() for row in rows)
