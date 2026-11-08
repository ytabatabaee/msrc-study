import csv
import unittest
from pathlib import Path

from empirical.microbotryum.scripts.stage2_analysis import PAIR_NAMES, species_pair, cross_species_summary


class Stage2AnalysisTests(unittest.TestCase):
    def test_three_interpretable_species_pairs(self):
        self.assertEqual(set(PAIR_NAMES), {"cat|lag", "cat|sup", "lag|sup"})

    def test_species_pair_recognition_is_fixed(self):
        self.assertEqual(species_pair("Mv-cat-C212-A1|Mv-cat-C212-A2|Mv-sup-1065-A1|Mv-sup-1065-A2"), "cat|sup")
        self.assertEqual(species_pair("Mv-cat-C212-A1|Mv-lag-1253-A1|Mv-sup-1065-A1|Mv-sup-1065-A2"), None)

    def test_real_cross_pair_results_are_not_overstated(self):
        path = Path("empirical/microbotryum/results/all_quartet_assignments.tsv")
        if not path.exists():
            self.skipTest("Stage 1 assignments are absent")
        with path.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        summary = cross_species_summary(rows)
        black = {(row["stratum"], row["species_pair"]): row for row in summary}
        self.assertAlmostEqual(float(black[("black", "cat|sup")]["q2_mating_type"]), 117 / 134)
        self.assertAlmostEqual(float(black[("black", "cat|lag")]["q2_mating_type"]), 0.0)
        self.assertAlmostEqual(float(black[("black", "lag|sup")]["q2_mating_type"]), 0.0)


if __name__ == "__main__":
    unittest.main()
