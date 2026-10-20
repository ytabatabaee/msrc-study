"""Deterministic tests for the Atlantic-cod Stage-12b extension."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "empirical/atlantic_cod/scripts/12b_robust_boundary_validation.py"
spec = importlib.util.spec_from_file_location("stage12b", SCRIPT)
stage12b = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(stage12b)


def test_stage12_reproduction_is_exact() -> None:
    mod, rows = stage12b.load_stage12()
    result = stage12b.baseline_reproduction(mod, rows, stage12b.load_regions())
    assert len(result) == 20
    assert all(r["exact_within_1e-9"] for r in result)


def test_fast_max_matches_candidate_scan() -> None:
    x = np.arange(30, dtype=float).reshape(10, 3)
    candidates = stage12b.candidate_profile(x)
    expected = max(0.0, max(float(r["score"]) for r in candidates))
    assert abs(stage12b.max_score_fast(x) - expected) < 1e-10


def test_direction_constraint_rejects_inconsistent_shift() -> None:
    rows = [{"D": float(x), "A": float(-x), "topo_a_third": 0.0, "topo_b_third": 0.0} for x in range(10)]
    fit, flags, _ = stage12b.fit_method(rows, "C_direction")
    assert fit["status"] == "no_detection"
    assert sum(flags) == 0


def test_stage12b_has_separate_output_namespace() -> None:
    text = SCRIPT.read_text()
    assert "boundary_detection_stage12b" in text
    assert 'RESULTS = EMP / "results" / "boundary_detection_stage12b"' in text
