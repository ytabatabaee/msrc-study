import csv
import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "empirical/neoaves_chr4/scripts/pnas48_unresolved317_audit.py"
spec = importlib.util.spec_from_file_location("pnas48_unresolved317_audit", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)
RESULTS = ROOT / "empirical/neoaves_chr4/results/stage4e"


class PNAS48Unresolved317Tests(unittest.TestCase):
    def test_leaf_and_polytomy_fixtures(self):
        tree = "(A:1,(B:1,C:1,D:1):1);"
        self.assertEqual(mod.taxa(tree), ("A", "B", "C", "D"))
        self.assertTrue(mod.has_polytomy(tree))
        self.assertEqual(mod.normalized("A-b"), "ab")
        self.assertEqual(mod.normalized_key(("A-b", "C")), ("ab", "c"))

    def test_global_matching_does_not_greedily_reuse_candidates(self):
        graph = {"row1": {"locus1"}, "row2": {"locus1", "locus2"}}
        matching = mod.hopcroft_karp(graph)
        self.assertEqual(len(matching), 2)
        self.assertEqual(set(matching.values()), {"locus1", "locus2"})

    def test_real_audit_has_exactly_317_rows_and_no_verified_mapping(self):
        csv.field_size_limit(20_000_000)
        with (RESULTS / "pnas48_unresolved317.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 317)
        self.assertTrue(all(row["classification"] != "verified mapping" for row in rows))
        self.assertEqual(sum(row["classification"] == "source discrepancy" for row in rows), 200)
        self.assertEqual(sum(row["classification"] == "treatment-equivalent mapping" for row in rows), 26)
        self.assertEqual(sum(row["classification"] == "unresolved" for row in rows), 91)

    def test_treatment_audit_preserves_partial_membership_status(self):
        with (RESULTS / "pnas48_treatment_membership_audit.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["treatment"] for row in rows}, set(mod.TREATMENTS))
        self.assertTrue(all(row["membership_verifiable"] == "Partial" for row in rows))
        self.assertTrue(all(int(row["potentially_affecting_membership_count"]) > 0 for row in rows))


if __name__ == "__main__":
    unittest.main()
