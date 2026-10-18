from __future__ import annotations

import csv
import json
import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "empirical/atlantic_cod/scripts/10b_fit_conditional_msrc_quartets.py"
RESULTS = ROOT / "empirical/atlantic_cod/results"


def load_module():
    spec = importlib.util.spec_from_file_location("stage10b_test_module", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def read_tsv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class Stage10bTests(unittest.TestCase):
    def test_original_stage8b_scaffold_and_no_mrca_likelihood(self):
        module = load_module()
        manifest = json.loads((RESULTS / "stage10b_manifest.json").read_text())
        tree_key = "empirical/atlantic_cod/results/stage8b_baseline_original_times.nwk"
        self.assertEqual(manifest["inputs"][tree_key], module.sha256(module.TREE))
        self.assertTrue(manifest["scaffold"]["rooted"])
        self.assertTrue(manifest["no_stage7_mrca_in_likelihood"])
        self.assertNotIn("choose_root(", SCRIPT.read_text())


    def test_exact_kernel_and_limits(self):
        module = load_module()
        self.assertLess(module.validate_frozen_kernel(), 1e-12)
        self.assertAlmostEqual(sum(module.conditional_probabilities(0.2, 1.0, 0.5)), 1.0)
        self.assertAlmostEqual(module.A_m(0.5, 1.0, 0.0), 1.0 - __import__("math").exp(-1.0))
        self.assertLess(module.A_m(0.5, 1.0, 1e9), 1e-8)


    def test_eligibility_is_restricted_to_supported_configuration(self):
        rows = read_tsv(RESULTS / "stage10b_quartet_eligibility.tsv")
        for row in rows:
            if row["eligible"] == "true":
                self.assertEqual(row["exclusion_reason"], "eligible")
                self.assertIn(row["local_class"], {"arrangement", "baseline", "third"})
                self.assertEqual(len(row["state_pattern"]), 4)
                self.assertEqual(sorted(row["state_pattern"].count(x) for x in "01"), [2, 2])
                self.assertGreater(float(row["interval_t"]), 0)
        self.assertGreater(sum(row["eligible"] == "true" for row in rows), 0)


    def test_contiguous_cv_has_no_train_test_overlap_and_profiles_reoptimize(self):
        assignments = read_tsv(RESULTS / "stage10b_cv_assignment.tsv")
        self.assertEqual(len({(r["lg"], r["window_id"]) for r in assignments}), len(assignments))
        by_lg = {}
        for row in assignments:
            by_lg.setdefault(row["lg"], []).append(row)
        for rows in by_lg.values():
            ordered = sorted(rows, key=lambda r: int(r["start"]))
            folds = [int(r["fold"]) for r in ordered]
            self.assertEqual(folds, sorted(folds))
            self.assertTrue(set(folds) <= set(range(5)))
        predictions = read_tsv(RESULTS / "stage10b_cv_predictions.tsv")
        assigned = {(r["lg"], r["window_id"]): int(r["fold"]) for r in assignments}
        self.assertTrue(all(assigned[(r["lg"], r["window_id"])] == int(r["fold"]) for r in predictions))
        profiles = read_tsv(RESULTS / "stage10b_parameter_profiles.tsv")
        self.assertTrue(profiles and all(r["nuisance_reoptimized"] == "true" for r in profiles))


    def test_legacy_stage10_hashes_unchanged(self):
        manifest = json.loads((RESULTS / "stage10b_manifest.json").read_text())
        module = load_module()
        for relative, expected in manifest["legacy_stage10_hashes"].items():
            self.assertEqual(module.sha256(ROOT / relative), expected)
