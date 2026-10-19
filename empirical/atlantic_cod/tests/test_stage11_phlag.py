#!/usr/bin/env python3
"""Regression checks for the Atlantic-cod Phlag benchmark outputs."""
from pathlib import Path
import csv, hashlib, re, unittest

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical" / "atlantic_cod"
OUT = EMP / "results" / "phlag_benchmark"
SCRIPT = EMP / "scripts" / "11_phlag_benchmark.py"


class Stage11PhlagTests(unittest.TestCase):
    def test_rooted_species_copy_and_no_root_search(self):
        text = (OUT / "inputs" / "baseline_population_tree_12taxa_labeled.nwk").read_text()
        self.assertTrue(text.startswith("[&R]"))
        source = SCRIPT.read_text()
        self.assertNotIn("choose_root", source)
        self.assertNotIn("root-edge", source)

    def test_focal_catalog_and_branch_count(self):
        with (OUT / "inputs" / "focal_branch_catalog.tsv").open() as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        self.assertEqual(len(rows), 10)
        self.assertEqual(sum(r["phlag_eligible"] == "true" for r in rows), 9)
        self.assertEqual([r["label"] for r in rows if r["phlag_eligible"] == "false"], ["N009"])

    def test_phlag_inputs_withhold_truth_coordinates(self):
        for p in (OUT / "inputs").glob("*.gene_trees"):
            self.assertNotIn("inside_inversion", p.read_text())
            self.assertNotIn("inversion_id", p.read_text())
        for p in (OUT / "inputs").glob("*.window_mapping.tsv"):
            self.assertNotIn("inside_inversion", p.read_text())
            self.assertNotIn("overlaps_inversion", p.read_text())

    def test_independent_lg_runs_and_raw_outputs(self):
        commands = (OUT / "phlag_commands.txt").read_text().splitlines()
        self.assertEqual(len(commands), 36)
        self.assertEqual({re.search(r"-g .*?/inputs/(LG\d+).gene_trees", x).group(1) for x in commands}, {"LG01", "LG02", "LG07", "LG12"})
        self.assertEqual(len(list((OUT / "raw").glob("LG*__N*.txt"))), 36)

    def test_metrics_are_present_and_phlag_default_is_explicitly_conservative(self):
        with (OUT / "detector_metrics.tsv").open() as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        phlag = [r for r in rows if r["method"] == "Phlag union" and r["lg"] != "MEAN"]
        self.assertEqual(len(phlag), 4)
        self.assertTrue(all(float(r["f1"]) == 0.0 and float(r["interval_iou"]) == 0.0 for r in phlag))
        self.assertTrue((OUT / "report.md").exists())
        self.assertTrue((EMP / "figures" / "atlantic_cod_phlag_comparison.pdf").exists())

    def test_stage10_files_are_not_touched_by_stage11_script(self):
        source = SCRIPT.read_text()
        self.assertNotIn("write_text", source.split("STAGE10", 1)[0] if "STAGE10" in source else "")
        self.assertNotIn("stage10_", source)


if __name__ == "__main__":
    unittest.main()
