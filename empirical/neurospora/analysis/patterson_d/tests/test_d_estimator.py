import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_patterson_d import site_contribution
from validate_outgroup_alignment import alignment_columns
from array import array


class PattersonEstimatorTests(unittest.TestCase):
    def test_pure_abba(self):
        r = site_contribution([ord("A")], [ord("G")], [ord("G")], ord("A"))
        self.assertEqual(r["status"], "polarizable")
        self.assertEqual((r["abba"], r["baba"]), (1.0, 0.0))

    def test_pure_baba(self):
        r = site_contribution([ord("G")], [ord("A")], [ord("G")], ord("A"))
        self.assertEqual((r["abba"], r["baba"]), (0.0, 1.0))

    def test_balanced_haploid_population_frequencies(self):
        r = site_contribution([ord("A"), ord("G")], [ord("A"), ord("G")], [ord("G")], ord("A"))
        self.assertAlmostEqual(r["abba"], 0.25)
        self.assertAlmostEqual(r["baba"], 0.25)

    def test_fixed_alt_outgroup_orients_reference_as_derived(self):
        # A is the reference allele, G the ALT allele, and the outgroup is fixed ALT.
        r = site_contribution([ord("A")], [ord("G")], [ord("G")], ord("G"))
        self.assertEqual(r["ancestral"], ord("G"))
        self.assertEqual(r["derived"], ord("A"))
        self.assertEqual(r["status"], "polarizable")
        self.assertEqual((r["abba"], r["baba"]), (0.0, 0.0))

    def test_missing_and_ambiguous_outgroup_are_excluded(self):
        self.assertEqual(site_contribution([65], [71], [71], ord("-"))["status"], "outgroup_gap")
        self.assertEqual(site_contribution([65], [71], [71], ord("N"))["status"], "outgroup_ambiguous")

    def test_missing_calls_and_multiallelic_sites_are_excluded(self):
        self.assertEqual(site_contribution([65, 78, 78], [71], [71], 65, 0.5)["status"], "insufficient_calls")
        self.assertEqual(site_contribution([65], [71], [67], 65)["status"], "multiallelic")
        self.assertEqual(site_contribution([65], [71], [71], 67)["status"], "outgroup_other")

    def test_complement_recoding_preserves_site_contributions(self):
        complement = {65: 84, 67: 71, 71: 67, 84: 65}
        original = site_contribution([65, 71], [71, 65], [71], 65)
        recoded = site_contribution([complement[x] for x in [65, 71]], [complement[x] for x in [71, 65]], [complement[71]], complement[65])
        self.assertAlmostEqual(original["abba"], recoded["abba"])
        self.assertAlmostEqual(original["baba"], recoded["baba"])

    def test_reference_position_maps_through_alignment_gaps(self):
        self.assertEqual(list(alignment_columns(b"A-CG", array("I", [1, 2, 3]))), [0, 2, 3])


if __name__ == "__main__":
    unittest.main()
