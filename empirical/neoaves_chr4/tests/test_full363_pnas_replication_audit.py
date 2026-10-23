from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "empirical/neoaves_chr4/scripts/full363_pnas_replication_audit.py"
SPEC = importlib.util.spec_from_file_location("full363_pnas_replication_audit", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_root_independent_rf_ignores_child_order():
    left = MODULE.load_tree(b"((A,B),(C,D));")
    right = MODULE.load_tree(b"((D,C),(B,A));")
    assert MODULE.rf(left, right) == (0, 0, 0)


def test_canonical_split_detects_the_alternative():
    baseline = MODULE.load_tree(b"((A,B),(C,D));")
    alternative = MODULE.load_tree(b"((A,C),(B,D));")
    result = MODULE.rf(baseline, alternative)
    assert result == (2, 1, 1)
    assert MODULE.canonical_group_split({"A", "B"}, baseline["taxa"]) in baseline["splits"]
    assert MODULE.canonical_group_split({"A", "B"}, alternative["taxa"]) not in alternative["splits"]


def test_focal_state_reports_q1_q2_and_absent():
    roles = ({"A"}, {"B"}, {"C"})
    assert MODULE.focal_state(MODULE.load_tree(b"((A,B),(C,D));"), roles) == "q1"
    assert MODULE.focal_state(MODULE.load_tree(b"((A,C),(B,D));"), roles) == "q2"
    assert MODULE.focal_state(MODULE.load_tree(b"((A,D),(B,C));"), roles) == "q3"
    assert MODULE.focal_state(MODULE.load_tree(b"(A,B,C,D);"), roles) == "absent_or_unresolved"


def test_cached_audit_outputs_record_expected_reference_comparisons():
    path = ROOT / "empirical/neoaves_chr4/results/stage4e/full363_pnas_tree_comparison.tsv"
    rows = { (row["left"], row["right"]): row for row in MODULE.read_tsv(path) }
    assert rows[("T0", "PNAS_FULL")]["rf"] == "0"
    assert rows[("T_PNAS", "PNAS_OUTLIER_REMOVED")]["rf"] == "2"
    assert rows[("PNAS_FULL", "PNAS_OUTLIER_REMOVED")]["rf"] == "0"


def test_focal_branch_table_preserves_absent_splits():
    path = ROOT / "empirical/neoaves_chr4/results/stage4e/full363_focal_branch_lengths.tsv"
    rows = MODULE.read_tsv(path)
    by_key = {(r["treatment"], r["clade"], r["state"]): r for r in rows}
    assert by_key[("T0", "N62", "q1")]["present"] == "True"
    assert by_key[("T0", "N62", "q2")]["present"] == "False"
    assert by_key[("T_PNAS", "N62", "q2")]["present"] == "True"
    assert by_key[("T_PNAS", "N62", "q1")]["CU_length"] == "NA"
