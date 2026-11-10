import importlib.util
import math
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "finite_sites.py"
spec = importlib.util.spec_from_file_location("finite_sites", SCRIPT)
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)


class FiniteSitesTests(unittest.TestCase):
    def test_jc69_transition_matrix(self):
        p = mod.pattern_probability([0, 0, 0, 0, 0], "q1", "JC69", {"pi": [0.25]*4, "kappa": 2, "pairs": {}})
        self.assertAlmostEqual(sum(p), 0.0)

    def test_zero_branch_star_is_symmetric(self):
        x = mod.pattern_probability([0, 0, 0, 0, 0], "q1", "JC69", {"pi": [0.25]*4, "kappa": 2, "pairs": {}})
        self.assertAlmostEqual(x[0], x[1])

    def test_long_branches_remain_nonnegative(self):
        x = mod.pattern_probability([10, 10, 10, 10, 10], "q2", "JC69", {"pi": [0.25]*4, "kappa": 2, "pairs": {}})
        self.assertGreaterEqual(x[0], 0); self.assertGreaterEqual(x[1], 0)

    def test_swapping_q2_q3_changes_label_not_probability_validity(self):
        p = {"pi": [0.25]*4, "kappa": 2, "pairs": {}}
        q2 = mod.pattern_probability([.1,.2,.3,.4,.05], "q2", "JC69", p)
        q3 = mod.pattern_probability([.1,.2,.3,.4,.05], "q3", "JC69", p)
        self.assertTrue(all(v >= 0 for v in q2+q3))

    def test_q2_and_q3_have_expected_abba_baba_orientation(self):
        p = {"pi": [0.25]*4, "kappa": 2, "pairs": {}}
        # Short pendants and a resolved internal edge make the split-specific
        # site-pattern asymmetry unambiguous under JC69.
        q2 = mod.pattern_probability([.01,.01,.01,.01,.10], "q2", "JC69", p)
        q3 = mod.pattern_probability([.01,.01,.01,.01,.10], "q3", "JC69", p)
        self.assertGreater(q2[1], q2[0])  # BABA supports q2
        self.assertGreater(q3[0], q3[1])  # ABBA supports q3

    def test_gtr_and_gamma_are_finite(self):
        p = {"pi": [0.3, 0.2, 0.2, 0.3], "kappa": 3, "pairs": {("A", "C"): 2, ("A", "G"): 4, ("A", "T"): 1, ("C", "G"): 1, ("C", "T"): 3, ("G", "T"): 2}}
        x = mod.model_row([.1, .2, .3, .4, .05], "q1", "GTR_Gamma", p, 100)
        self.assertTrue(all(math.isfinite(v) and v >= 0 for v in x))

    def test_division_zero(self):
        self.assertIsNone(mod.div(0,0))


if __name__ == "__main__": unittest.main()
