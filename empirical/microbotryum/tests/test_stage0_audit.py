import tempfile
import unittest
from pathlib import Path

from empirical.microbotryum.scripts.stage0_audit import classify_quartet, parse_newick, split_newicks


class Stage0AuditTests(unittest.TestCase):
    def test_newick_records_and_branch_lengths(self):
        text = "((S1_a1:0.1,S1_a2:0.2):0.3,(S2_a1:0.1,S2_a2:0.2):0.3);"
        self.assertEqual(len(split_newicks(text)), 1)
        taxa, weighted, resolved = parse_newick(text)
        self.assertEqual(taxa, {"S1_a1", "S1_a2", "S2_a1", "S2_a2"})
        self.assertTrue(weighted)
        self.assertTrue(resolved)

    def test_quartet_labels_are_fixed(self):
        q = ("S1_a1", "S1_a2", "S2_a1", "S2_a2")
        self.assertEqual(classify_quartet("S1_a1,S1_a2|S2_a1,S2_a2", q), "q1")
        self.assertEqual(classify_quartet("S1_a1,S2_a1|S1_a2,S2_a2", q), "q2")
        self.assertEqual(classify_quartet("S1_a1,S2_a2|S1_a2,S2_a1", q), "q3")

    def test_missing_or_malformed_tree_is_not_silently_valid(self):
        with self.assertRaises(ValueError):
            parse_newick("((a,b),(c,d)")


if __name__ == "__main__":
    unittest.main()
