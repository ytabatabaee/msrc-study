import importlib.util
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/multiresolution_tracks.py"
spec = importlib.util.spec_from_file_location("multires", SCRIPT)
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

class TestMultiresolutionTracks(unittest.TestCase):
    def test_boundary_assignment(self):
        self.assertEqual(mod.window_index(0, 100000), 0)
        self.assertEqual(mod.window_index(99999, 100000), 0)
        self.assertEqual(mod.window_index(100000, 100000), 100000)

    def test_pooled_arithmetic_and_zero(self):
        self.assertEqual(mod.d_from_counts(8, 2), 0.6)
        self.assertIsNone(mod.d_from_counts(0, 0))

    def test_empty_and_low_information_windows(self):
        rows=[{"design":"X","filter":"all","start":0,"end":10,"midpoint":5,"ABBA":0,"BABA":0,"informative":0}, {"design":"X","filter":"all","start":100000,"end":100010,"midpoint":100005,"ABBA":20,"BABA":10,"informative":30}]
        out=mod.aggregate(rows, 100000)
        zero=next(r for r in out if r["window_start"]==0)
        self.assertEqual(zero["Patterson_D"], "")
        self.assertEqual(zero["n_loci"], 1)
        low=next(r for r in out if r["window_start"]==100000)
        self.assertEqual(low["low_information"], 1)

    def test_region_overlap(self):
        self.assertEqual(mod.overlapping_regions(25030000, 25030001), "PNAS_1")
        self.assertEqual(mod.overlapping_regions(0, 100), "outside_all_PNAS")
        self.assertIn("PNAS_1", mod.overlapping_regions(24900000, 25100000))

    def test_reproducibility(self):
        rows=[{"design":"X","filter":"all","start":0,"end":10,"midpoint":5,"ABBA":3,"BABA":1,"informative":4}]
        self.assertEqual(mod.aggregate(rows,100000), mod.aggregate(rows,100000))

if __name__ == "__main__": unittest.main()
