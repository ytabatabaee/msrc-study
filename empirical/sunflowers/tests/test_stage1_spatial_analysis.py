"""Focused tests for fixed quartet scoring and spatial interval summaries."""
from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/stage1_spatial_analysis.py"
SPEC = importlib.util.spec_from_file_location("sunflower_stage1", SCRIPT)
stage1 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage1)


class QuartetTests(unittest.TestCase):
    def test_extracts_all_fixed_individual_quartets(self):
        tips = ["ANN1", "ANN2", "ARG1", "ARG2", "PET1", "PET2", "DEB1", "DEB2", "OTHER"]
        combos = stage1.quartet_combinations(tips)
        self.assertEqual(len(combos), 16)
        self.assertTrue(all(len(set(combo)) == 4 for combo in combos))
        self.assertTrue(all("OTHER" not in combo for combo in combos))
        par_combos=stage1.quartet_combinations(["ANN1","ARG1","PET1","PAR1","PAR2"],"par")
        self.assertEqual(len(par_combos),2)

    def test_extracts_each_induced_topology(self):
        # Terminal lengths 1 and internal length 2 yield a 2-unit quartet edge.
        splits = [
            (2, 4, 4, 4, 4, 2),  # AB | PD (q1)
            (4, 2, 4, 4, 2, 4),  # AP | BD (q2)
            (4, 4, 2, 2, 4, 4),  # AD | BP (q3)
        ]
        for expected, distances in enumerate(splits):
            split, internal, _ = stage1.score_quartet_distances(*distances)
            self.assertEqual(split, expected)
            self.assertAlmostEqual(internal, 2.0)

    def test_polytomy_is_unresolved(self):
        split, internal, _ = stage1.score_quartet_distances(2, 2, 2, 2, 2, 2)
        self.assertIsNone(split)
        self.assertEqual(internal, 0.0)

    def test_mean_pairwise_distance(self):
        self.assertAlmostEqual(stage1.mean_pairwise_distance([0.1, 0.3, 0.5]), 0.3)


class SpatialTests(unittest.TestCase):
    def test_genomic_overlap_unions_overlapping_parts(self):
        self.assertEqual(stage1.overlap_bp(5, 25, [(0, 10), (8, 20)]), 15)
        self.assertEqual(stage1.overlap_bp(20, 25, [(0, 10), (8, 20)]), 0)

    def test_adjacent_requires_nonnegative_gap_at_most_threshold(self):
        a = {"start": 0, "end": 10_000}
        self.assertTrue(stage1.adjacent(a, {"start": 25_000, "end": 35_000}))
        self.assertFalse(stage1.adjacent(a, {"start": 30_001, "end": 40_001}))
        self.assertFalse(stage1.adjacent(a, {"start": 9_000, "end": 20_000}))

    def test_run_detection_does_not_bridge_missing_intervals(self):
        rows = []
        for start, q2 in ((0, .8), (10_000, .8), (100_000, .8), (110_000, .7)):
            rows.append({"start": start, "end": start + 8_000,
                         "dominant_topology": "q2_ann_pet", "q1_ann_arg": .1,
                         "q2_ann_pet": q2, "q3_arg_pet": .1})
        runs = stage1.detect_runs(rows, "ann05.01", "interior")
        lengths = [r["n_adjacent_windows"] for r in runs if r["topology"] == "q2_ann_pet"]
        self.assertEqual(lengths, [2, 2])


if __name__ == "__main__":
    unittest.main()
