import csv
import importlib.util
import sys
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "11b_repair_phlag_benchmark.py"
spec = importlib.util.spec_from_file_location("stage11b_repair", SCRIPT)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)


class Stage11BTests(unittest.TestCase):
    def test_manifest_gene_tree_and_phlag_lengths(self):
        manifest = mod.read_tsv(mod.RES / "matched_locus_manifest.tsv")
        n_trees = sum(1 for _ in (mod.INP / "Columbea.gtrees").open())
        valid = mod.read_tsv(mod.RES / "phlag_output_validation.tsv")[0]
        self.assertEqual(len(manifest), n_trees)
        self.assertGreater(n_trees, 0)
        self.assertEqual(int(valid["n_states"]), n_trees)
        self.assertEqual(int(valid["n_posteriors"]), n_trees)

    def test_strict_output_parser_rejects_length_mismatch(self):
        p = Path("/tmp/stage11b_bad_phlag.txt")
        p.write_text("0,1\n0.1,0.2\n")
        with self.assertRaises(ValueError):
            mod.parse_phlag_output(p, 3)

    def test_physical_gap_merging(self):
        rows = [{"state": 1, "start": 0, "end": 10, "input_index": 0},
                {"state": 1, "start": 20, "end": 30, "input_index": 1},
                {"state": 1, "start": 100, "end": 110, "input_index": 2}]
        self.assertEqual(len(mod.merge_flagged(rows, 10)), 2)
        self.assertEqual(len(mod.merge_flagged(rows, 9)), 3)

    def test_manifest_coordinates_are_ordered_and_unique(self):
        rows = mod.read_tsv(mod.RES / "matched_locus_manifest.tsv")
        self.assertEqual(len({r["Gene"] for r in rows}), len(rows))
        starts = [int(r["start"]) for r in rows]
        self.assertEqual(starts, sorted(starts))
        self.assertTrue(all(r["all_selected_taxa_present"] == "true" for r in rows))


if __name__ == "__main__":
    unittest.main()
