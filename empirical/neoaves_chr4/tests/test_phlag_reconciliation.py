import csv
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical/neoaves_chr4"
RES = EMP / "results/stage11_phlag_benchmark"
MULTI = EMP / "results/stage11_original_multitaxon_detector"


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_reconcile_module():
    path = EMP / "scripts/reconcile_phlag_benchmark.py"
    spec = importlib.util.spec_from_file_location("reconcile_phlag_benchmark", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class PhlagReconciliationTests(unittest.TestCase):
    def test_validated_phlag_counts_and_finite_posteriors(self):
        manifest = read_tsv(RES / "matched_locus_manifest.tsv")
        trees = sum(1 for _ in (RES / "matched_inputs/Columbea.gtrees").open())
        predictions = read_tsv(RES / "window_level_predictions.tsv")
        validation = read_tsv(RES / "phlag_output_validation.tsv")[0]
        self.assertEqual(len(manifest), 5619)
        self.assertEqual(len(manifest), trees)
        self.assertEqual(len(predictions), len(manifest))
        self.assertEqual(int(validation["n_states"]), len(manifest))
        self.assertEqual(int(validation["n_posteriors"]), len(manifest))
        self.assertTrue(all(math.isfinite(float(row["posterior_anomaly"])) for row in predictions))
        self.assertTrue(all(0 <= float(row["posterior_anomaly"]) <= 1 for row in predictions))

    def test_seven_intervals_are_generated_from_flagged_manifest_rows(self):
        module = load_reconcile_module()
        rows = read_tsv(RES / "window_level_predictions.tsv")
        flagged_file = read_tsv(RES / "phlag_anomalous_windows.tsv")
        self.assertEqual(
            [(row["input_index"], row["start"], row["end"]) for row in flagged_file],
            [(row["input_index"], row["start"], row["end"]) for row in rows if row["state"] == "1"],
        )
        manifest = read_tsv(RES / "matched_locus_manifest.tsv")
        self.assertEqual(
            [(row["input_index"], row["start"], row["end"]) for row in rows],
            [(row["input_index"], row["start"], row["end"]) for row in manifest],
        )
        expected = []
        for row in rows:
            if row["state"] != "1":
                continue
            start, end = int(row["start"]), int(row["end"])
            if expected and start - expected[-1][1] <= 250000:
                expected[-1] = (expected[-1][0], max(expected[-1][1], end))
            else:
                expected.append((start, end))
        observed = module.parse_intervals(RES / "phlag_predicted_intervals.tsv")
        self.assertEqual(expected, observed)
        self.assertEqual(len(observed), 7)

    def test_stage11b_metadata_was_not_replaced_by_multitaxon_metadata(self):
        manifest = json.loads((RES / "manifest.json").read_text())
        self.assertEqual(manifest["stage"], "11B")
        self.assertEqual(manifest["n_input_trees"], 5619)
        self.assertEqual(manifest["n_manifest"], 5619)
        self.assertEqual(manifest["n_excluded"], 477)
        multi = read_tsv(MULTI / "our_detector_predictions.tsv")
        self.assertEqual({row["n_loci"] for row in multi}, {"6096"})
        self.assertEqual({row["method"] for row in multi}, {"topology", "topology_plus_B"})

    def test_method_domains_and_features_are_explicit(self):
        rows = read_tsv(RES / "interval_comparison.tsv")
        self.assertEqual(len(rows), 48)
        domains = {row["method"]: row["evaluation"] for row in rows if row["target"] == "all_canonical_regions_union"}
        self.assertEqual(domains["official_phlag_matched_four_taxon"], "matched_input_post_fit")
        self.assertEqual(domains["our_matched_four_taxon_topology"], "matched_input_post_fit")
        self.assertEqual(domains["our_original_multitaxon_topology"], "unmatched_input_post_fit")
        self.assertEqual(domains["our_original_multitaxon_topology_plus_B"], "unmatched_input_post_fit")
        self.assertIn("D2", (MULTI / "methods_text.md").read_text())
        self.assertIn("B_observed", (MULTI / "methods_text.md").read_text())
        quartet_header = (RES / "matched_four_taxon_quartets.tsv").read_text().splitlines()[0]
        self.assertNotIn("canonical_region", quartet_header)

    def test_regenerated_tables_match_source_predictions(self):
        module = load_reconcile_module()
        intervals = module.detector_intervals()
        refs = read_tsv(EMP / "results/canonical_structural_regions.tsv")
        observed = read_tsv(RES / "interval_comparison.tsv")
        by_key = {(row["method"], row["target"]): row for row in observed}
        reference_union = [(int(row["canonical_start"]), int(row["canonical_end"])) for row in refs]
        for method, predicted in intervals.items():
            expected = module.metrics(predicted, reference_union)
            row = by_key[(method, "all_canonical_regions_union")]
            self.assertEqual(int(row["n_predicted_intervals"]), expected["n_predicted_intervals"])
            self.assertEqual(int(row["intersection_bp"]), expected["intersection_bp"])
            self.assertEqual(float(row["iou"]), expected["iou"])


if __name__ == "__main__":
    unittest.main()
