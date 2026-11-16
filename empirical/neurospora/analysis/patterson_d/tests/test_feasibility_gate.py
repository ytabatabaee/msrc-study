import csv
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / "empirical/neurospora/analysis/patterson_d"


def read_tsv(name):
    with (BASE / "results" / name).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class FeasibilityGateTests(unittest.TestCase):
    def test_archive_coordinate_audit_confirms_each_window_has_50_sites(self):
        audit = json.loads((BASE / "results/coordinate_audit.json").read_text())
        self.assertEqual(audit["windows"], 45303)
        self.assertEqual(audit["genotype_rows"], 2265369)
        self.assertEqual(audit["genotype_rows_per_window_histogram"], {"50": 45303})
        self.assertEqual(audit["rows_outside_window_intervals"], 219)

    def test_only_documented_sample_pairings_are_promoted(self):
        rows = {r["sample_id"]: r for r in read_tsv("cj_mating_type_resolution.tsv")}
        self.assertEqual(rows["CJ02"]["mating_type"], "a")
        self.assertEqual(rows["CJ03"]["mating_type"], "A")
        self.assertEqual(rows["CJ04"]["mating_type"], "a")
        self.assertEqual(rows["CJ15"]["mating_type"], "A")
        self.assertEqual(rows["CJ16"]["mating_type"], "a")
        unresolved = {s for s, r in rows.items() if r["mating_type"] == "unresolved"}
        self.assertEqual(unresolved, {"CJ07", "CJ11", "CJ12", "CJ18"})

    def test_outgroup_callability_is_measured_for_each_region(self):
        rows = read_tsv("outgroup_callability.tsv")
        self.assertEqual(len(rows), 3)
        self.assertEqual([int(r["N_discreta_callable_rows"]) for r in rows], [301014, 149665, 1474160])
        self.assertTrue(all("validated" in r["status"] for r in rows))

    def test_quartet_is_frozen_before_d_calculation(self):
        rows = read_tsv("../config/frozen_quartets.tsv")
        self.assertEqual(len(rows), 1)
        self.assertIn("frozen before D inspection", rows[0]["status"])


if __name__ == "__main__":
    unittest.main()
