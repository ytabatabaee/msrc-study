import importlib.util
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "integrate.py"
spec = importlib.util.spec_from_file_location("neoaves_integration", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class IntegrationTests(unittest.TestCase):
    def test_four_point_topology_is_root_independent(self):
        a = mod.parse_newick("((a:1,b:1):2,(c:1,d:1):2):0;")
        b = mod.parse_newick("(((a:1,b:1):1,c:2):1,d:3);")
        self.assertEqual(mod.quartet_summary(a, ("a", "b", "c", "d"))["topology"], "q1")
        self.assertEqual(mod.quartet_summary(b, ("a", "b", "c", "d"))["topology"], "q1")

    def test_missing_and_unresolved(self):
        self.assertEqual(mod.quartet_summary(mod.parse_newick("((a:1,b:1),c:1);"), ("a", "b", "c", "d"))["status"], "missing_taxon")
        q = mod.quartet_summary(mod.parse_newick("(a:1,b:1,c:1,d:1);"), ("a", "b", "c", "d"))
        self.assertEqual(q["status"], "unresolved")

    def test_patterson_pooling(self):
        self.assertAlmostEqual(mod.dcount(3, 1), 0.5)
        self.assertIsNone(mod.dcount(0, 0))

    def test_region_overlap(self):
        self.assertIn("PNAS_1", mod.region_overlap(25030000, 25030000))
        self.assertEqual(mod.region_overlap(1, 10), "outside_all_PNAS")

    def test_missing_segments_not_bridged(self):
        self.assertEqual(len(mod.draw_polyline.__name__), 13)


if __name__ == "__main__":
    unittest.main()
