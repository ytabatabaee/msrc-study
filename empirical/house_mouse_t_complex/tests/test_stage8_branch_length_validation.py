#!/usr/bin/env python3
"""Regression checks for house-mouse Stage 8 branch validation."""
from __future__ import annotations

import csv
import math
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "empirical/house_mouse_t_complex/results/stage8_branch_length_validation"


def rows(path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


class Stage8BranchTests(unittest.TestCase):
    def test_all_source_windows_and_thresholds(self):
        z = rows(OUT / "branch_length_distribution.tsv")
        self.assertEqual(len(z), 4046)
        focal = [float(r["focal_terminal_branch"]) for r in z]
        self.assertGreater(max(focal), 50)
        self.assertEqual(sum(x >= 1 for x in focal), 2748)
        self.assertEqual(sum(x >= 5 for x in focal), 1678)
        self.assertEqual(sum(x >= 10 for x in focal), 309)
        self.assertEqual(sum(x >= 20 for x in focal), 212)
        self.assertEqual(sum(x >= 30 for x in focal), 83)

    def test_stage4_distance_reproduction(self):
        z = rows(OUT / "robust_distance_comparison.tsv")
        self.assertEqual(len(z), 4046)
        for r in z[::101]:
            self.assertAlmostEqual(float(r["B_all"]), float(r["stage4_mean_pair_distance_no_spretus"]), places=12)
            self.assertAlmostEqual(float(r["A_distance_original"]), float(r["stage4_A_distance"]), places=12)
            self.assertAlmostEqual(float(r["A_relative_original"]), float(r["stage4_A_relative"]), places=12)

    def test_source_filter_comparison_identifies_filter_specific_pattern(self):
        z = rows(OUT / "source_filter_comparison.tsv")
        data2 = [r for r in z if r["filtering_procedure"] == "Coverage_AlleleRatio_Filtered"]
        self.assertTrue(any(float(r["focal_terminal_branch"]) > 40 for r in data2))
        for label in ("PASS_SNPs", "Coverage_Filtered"):
            for r in z:
                if r["filtering_procedure"] == label and r["focal_terminal_branch"]:
                    self.assertLess(float(r["focal_terminal_branch"]), 1)

    def test_quartet_support_is_normalized(self):
        z = rows(OUT / "sample_branch_length_audit.tsv")
        for r in z[::101]:
            self.assertAlmostEqual(sum(float(r[f"q_{k}_original"]) for k in ("species", "t_alt", "other", "unresolved")), 1.0, places=10)
            self.assertAlmostEqual(sum(float(r[f"q_{k}_leave_one_out"]) for k in ("species", "t_alt", "other", "unresolved")), 1.0, places=10)

    def test_boundary_results_are_fixed_configuration(self):
        z = rows(OUT / "boundary_sensitivity.tsv")
        self.assertEqual(len(z), 8)
        self.assertTrue(all(r["model"] == "single_interval" for r in z))
        self.assertTrue(all("references post-fit" in r["configuration"] for r in z))


if __name__ == "__main__":
    unittest.main()
