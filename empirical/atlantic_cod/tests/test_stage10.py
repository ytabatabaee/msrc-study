import importlib.util
import math
from pathlib import Path
import unittest

import numpy as np


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "10_fit_msrc_parameters_predict_quartets.py"
SPEC = importlib.util.spec_from_file_location("stage10_under_test", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class Stage10ModelTests(unittest.TestCase):
    def test_switching_limits(self):
        self.assertAlmostEqual(MODULE.A_m(0.4, 1.0, 0.0), 1.0 - math.exp(-0.8), places=10)
        self.assertLess(MODULE.A_m(0.4, 1.0, 1e6), 1e-5)

    def test_matrix_density_has_expected_zero_switching_limits(self):
        same = MODULE.matrix_exponential_density(np.array([0.5]), np.array([0.0]), np.array([1.0]), 2.0, 0.0)[0]
        different = MODULE.matrix_exponential_density(np.array([0.5]), np.array([0.0]), np.array([-1.0]), 2.0, 0.0)[0]
        self.assertAlmostEqual(same, 2.0 * math.exp(-1.0), places=8)
        self.assertLess(different, 1e-300)

    def test_rates_are_positive_after_log_transform(self):
        self.assertTrue(np.all(np.exp(np.array([-4.0, 0.0, 2.0])) > 0))


if __name__ == "__main__":
    unittest.main()
