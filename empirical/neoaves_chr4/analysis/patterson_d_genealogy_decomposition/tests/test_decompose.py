import importlib.util
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "decompose.py"
spec = importlib.util.spec_from_file_location("decompose", SCRIPT)
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)


class DecompositionTests(unittest.TestCase):
    def test_known_internal_edge_and_root_independence(self):
        q = mod.induced_quartet(mod.parse_newick("((a:1,b:1):2,(c:1,d:1):2);"), ("a", "b", "c", "d"))
        self.assertEqual(q["topology"], "q1"); self.assertAlmostEqual(q["internal_branch"], 4.0)
        r = mod.induced_quartet(mod.parse_newick("(((a:1,b:1):1,c:2):1,d:3);"), ("a", "b", "c", "d"))
        self.assertEqual(r["topology"], "q1"); self.assertAlmostEqual(r["internal_branch"], 1.0)

    def test_star_unresolved(self):
        q = mod.induced_quartet(mod.parse_newick("(a:1,b:1,c:1,d:1);"), ("a", "b", "c", "d"))
        self.assertEqual(q["status"], "unresolved")

    def test_prediction_arithmetic(self):
        self.assertAlmostEqual(mod.D(3, 1), .5)
        self.assertIsNone(mod.D(0, 0))

    def test_shapley_sum(self):
        inside = {"C_A": 4, "C_B": 2, "mu_A": 1, "mu_B": 3, "n_ABBA_topology": 4, "n_BABA_topology": 2}
        outside = {"C_A": 2, "C_B": 2, "mu_A": 2, "mu_B": 2, "n_ABBA_topology": 2, "n_BABA_topology": 2}
        x = mod.decompose(inside, outside)
        self.assertAlmostEqual(x["topology_plus_branch_check"], x["D_tree_shift"])

    def test_callable_sites(self):
        seqs = {"a": "ACN-", "b": "ACGT", "c": "ACGT", "d": "ACGT"}
        self.assertEqual(mod.callable_sites(seqs, ("a", "b", "c", "d")), (4, 2))


if __name__ == "__main__": unittest.main()
