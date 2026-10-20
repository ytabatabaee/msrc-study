#!/usr/bin/env python3
"""Unit tests for the gap-aware Stage 12 interval detector."""
from pathlib import Path
import importlib.util
import unittest
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "empirical" / "atlantic_cod" / "scripts" / "12_boundary_detection.py"
spec = importlib.util.spec_from_file_location("stage12_boundary", SCRIPT)
stage12 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stage12)


class Stage12BoundaryTests(unittest.TestCase):
    def test_no_change_option_wins_for_constant_features(self):
        fit = stage12.interval_fit(np.zeros((20, 2)), min_segment=3, penalty_factor=1.0)
        self.assertEqual(fit["status"], "no_detection")
        self.assertEqual(fit["score"], 0.0)

    def test_internal_change_point(self):
        x = np.vstack([np.zeros((8, 1)), np.ones((8, 1)), np.zeros((8, 1))])
        fit = stage12.interval_fit(x, min_segment=3, penalty_factor=0.1)
        self.assertEqual(fit["status"], "internal")
        self.assertEqual((fit["start_index"], fit["end_index"]), (8, 15))

    def test_right_and_left_censoring(self):
        right = np.vstack([np.zeros((12, 1)), np.ones((8, 1))])
        left = np.vstack([np.ones((8, 1)), np.zeros((8, 1))])
        self.assertEqual(stage12.interval_fit(right, 3, 0.1)["status"], "right_censored")
        self.assertEqual(stage12.interval_fit(left, 3, 0.1)["status"], "left_censored")

    def test_genomic_gap_splits_window_runs(self):
        rows = [{"start": "1", "end": "100"}, {"start": "101", "end": "200"}, {"start": "1000", "end": "1100"}]
        self.assertEqual(stage12.runs(rows, [1, 1, 1]), [(1, 200, 0, 1), (1000, 1100, 2, 2)])

    def test_metrics_iou_f1_and_partial_window(self):
        rows = [{"start": "1", "end": "100"}, {"start": "101", "end": "200"}, {"start": "201", "end": "300"}]
        m = stage12.metrics(rows, [0, 1, 1], (150, 250), (1, 300))
        self.assertAlmostEqual(m["precision"], 1.0)
        self.assertAlmostEqual(m["recall"], 1.0)
        self.assertAlmostEqual(m["f1"], 1.0)
        self.assertGreater(m["interval_iou_full"], 0.0)


if __name__ == "__main__":
    unittest.main()
