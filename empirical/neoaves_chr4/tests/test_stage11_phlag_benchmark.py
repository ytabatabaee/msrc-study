#!/usr/bin/env python3
import importlib.util
import sys
from pathlib import Path

SCRIPT = Path(__file__).parents[1] / "scripts/11_phlag_benchmark.py"
spec = importlib.util.spec_from_file_location("stage11", SCRIPT)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


def test_interval_set_metrics_merge_overlaps_and_preserve_gaps():
    m = mod.interval_metrics([(0, 10), (5, 15), (30, 40)], [(0, 15), (30, 35)])
    assert m["predicted_bp"] == 25
    assert m["reference_bp"] == 20
    assert m["intersection_bp"] == 20
    assert m["iou"] == 20 / 25


def test_robust_scale_is_deterministic():
    import numpy as np
    x = np.array([[1.0, 10.0], [2.0, 11.0], [3.0, 12.0]])
    a = mod.robust_scale(x)[0]
    b = mod.robust_scale(x)[0]
    assert np.allclose(a, b)


def test_stage6_columbea_rows_are_ordered_and_normalized():
    rows = mod.load_locus_rows()
    assert rows
    mids = [float(r["midpoint"]) for r in rows]
    assert mids == sorted(mids)
    assert max(abs(sum(float(r[k]) for k in ("q1_species_history", "q2_alternative_1", "q3_alternative_2")) - 1) for r in rows) < 1e-8


def test_phlag_and_detector_are_separate_methods():
    assert "Phlag" not in "stage6 QQS summaries"
    assert mod.PHLAG.name == "phlag"
