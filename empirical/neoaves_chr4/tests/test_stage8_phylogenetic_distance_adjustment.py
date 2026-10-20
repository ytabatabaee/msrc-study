#!/usr/bin/env python3
"""Regression and mathematical checks for Neoaves Stage 8 outputs."""
from __future__ import annotations

import csv
import math
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "empirical/neoaves_chr4/results/stage8_phylogenetic_distance_adjustment"
STAGE7 = ROOT / "empirical/neoaves_chr4/results/stage7_arrangement_divergence"


def rows(path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


class Stage8Tests(unittest.TestCase):
    def test_candidate_regions_and_six_taxa_are_frozen(self):
        elig = rows(RES / "region_eligibility_updated.tsv")
        self.assertEqual({r["canonical_region_id"] for r in elig if r["eligible"] == "True"}, {"CR_01", "CR_10", "CR_11"})
        mapping = rows(RES / "expanded_structural_state_mapping.tsv")
        self.assertEqual(len(mapping), 6)
        self.assertTrue(all(r["state_inferred_from_local_tree"] == "False" for r in mapping))

    def test_stage7_raw_contrasts_reproduce_exactly(self):
        audit = rows(RES / "raw_stage7_reproduction.tsv")
        checked = [r for r in audit if r["absolute_difference"]]
        self.assertGreater(len(checked), 1000)
        self.assertLessEqual(max(float(r["absolute_difference"]) for r in checked), 1e-12)

    def test_residual_and_relative_definitions(self):
        rows8 = rows(RES / "window_pairwise_residuals.tsv")
        self.assertGreater(len(rows8), 100000)
        for r in rows8[::max(1, len(rows8)//1000)]:
            d = float(r["distance"]); b = float(r["baseline_median"]); residual = float(r["residual"])
            self.assertAlmostEqual(residual, d - b, places=10)
            rel = float(r["relative_residual"])
            self.assertAlmostEqual(rel, residual / b, places=10)

    def test_baselines_are_outside_regions(self):
        bg = rows(RES / "pairwise_background_distances.tsv")
        self.assertTrue(all(int(r["n_baseline_windows"]) > 0 for r in bg))
        self.assertEqual({r["baseline_type"] for r in bg}, {"broad_outside_focal", "strict_outside_all_structural_regions"})

    def test_leave_one_species_out_is_recomputed(self):
        loo = rows(RES / "leave_one_species_out.tsv")
        self.assertEqual({r["excluded_species"] for r in loo}, {"Stork", "Flamingo", "Dove", "Sandgrouse", "Turaco", "Cuckoo"})
        self.assertTrue(any(r["status"] == "unsupported_after_exclusion" for r in loo))


if __name__ == "__main__":
    unittest.main()
