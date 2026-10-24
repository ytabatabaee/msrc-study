import importlib.util
from pathlib import Path
import unittest


SCRIPT = Path(__file__).parents[1] / "scripts" / "05_stage4e_analyze_j48_mixtures.py"
SPEC = importlib.util.spec_from_file_location("j48_mixture_analysis", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class J48MixtureAnalysisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.validation, cls.topology, cls.split_presence, cls.focal, cls.branches, cls.meta = MODULE.validate_and_analyze()

    def test_canonical_split_is_unrooted(self):
        taxa = {"a", "b", "c", "d"}
        self.assertEqual(MODULE.canonical_split({"a", "b"}, taxa), MODULE.canonical_split({"c", "d"}, taxa))

    def test_all_trees_have_48_taxa_and_cu_su_agree(self):
        rows = [row for row in self.validation if row["treatment"].startswith("J48_MIX")]
        self.assertEqual(len(rows), 13)
        self.assertTrue(all(row["n_taxa_cu"] == 48 and row["n_taxa_su"] == 48 for row in rows))
        self.assertTrue(all(row["cu_su_rf"] == 0 for row in rows))
        self.assertTrue(all(row["output_hashes_match"] == "true" for row in rows))

    def test_expected_rf_series(self):
        got = [(float(row["epsilon"]), int(row["rf_vs_NONCHR4"]), int(row["rf_vs_T0"])) for row in self.topology]
        expected = [(0.00, 0, 6), (0.02, 4, 6), (0.05, 2, 4), (0.10, 6, 0),
                    (0.15, 6, 0), (0.20, 6, 0), (0.25, 6, 0), (0.30, 6, 0),
                    (0.40, 6, 0), (0.50, 6, 0), (0.60, 6, 0), (0.75, 6, 0), (1.00, 6, 0)]
        self.assertEqual(got, expected)

    def test_novel_split_is_retained(self):
        by_epsilon = {float(row["epsilon"]): row for row in self.topology}
        self.assertEqual(int(by_epsilon[0.02]["n_novel_splits"]), 1)
        self.assertEqual(int(by_epsilon[0.05]["n_novel_splits"]), 0)

    def test_absent_branch_is_missing(self):
        row = next(row for row in self.split_presence if row["epsilon"] == 0.0 and row["relationship_label"] == "N62_q1")
        self.assertEqual(row["present_cu"], "false")
        self.assertEqual(row["CU_length"], "")
        self.assertEqual(row["localPP_cu"], "")

    def test_focal_states_identify_first_changes(self):
        states = {(float(row["epsilon"]), row["focal_relationship"]): row["observed_state_cu"] for row in self.focal}
        self.assertEqual(states[(0.00, "N62")], "q2")
        self.assertEqual(states[(0.02, "N62")], "q1")
        self.assertEqual(states[(0.05, "Columbea")], "q2")
        self.assertEqual(states[(0.10, "Columbea")], "q1")
        self.assertEqual(states[(0.10, "N61")], "q2")

    def test_natural_fraction_is_manifest_derived(self):
        self.assertEqual(self.meta["natural_count"], 6262)
        self.assertEqual(self.meta["total_count"], 63430)


if __name__ == "__main__":
    unittest.main()
