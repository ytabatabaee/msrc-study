#!/usr/bin/env python3
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from population_mapping import (  # noqa: E402
    coalescent_units_to_generations,
    generations_to_coalescent_units,
    map_simphy_haploid_to_msrc,
)


class PopulationMappingTests(unittest.TestCase):
    def test_even_haploid_size_is_exact(self):
        mapping = map_simphy_haploid_to_msrc(414468)
        self.assertEqual(mapping.msrc_diploid_ne_integer, 207234)
        self.assertEqual(mapping.relative_rate_error, 0.0)
        self.assertAlmostEqual(generations_to_coalescent_units(414468, mapping), 1.0)

    def test_odd_size_reports_finite_approximation(self):
        mapping = map_simphy_haploid_to_msrc(414467)
        self.assertEqual(mapping.msrc_diploid_ne_exact, 207233.5)
        self.assertIn(mapping.msrc_diploid_ne_integer, (207233, 207234))
        self.assertNotEqual(mapping.relative_rate_error, 0.0)
        self.assertAlmostEqual(
            mapping.simphy_pair_coalescence_rate,
            1.0 / mapping.simphy_haploid_ne,
        )

    def test_inverse_conversion(self):
        mapping = map_simphy_haploid_to_msrc(1000000)
        self.assertTrue(math.isclose(coalescent_units_to_generations(2.5, mapping), 2500000.0))


if __name__ == "__main__":
    unittest.main()
