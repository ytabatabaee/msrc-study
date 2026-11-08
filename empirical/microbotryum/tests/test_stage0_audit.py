import unittest
from pathlib import Path

from empirical.microbotryum.scripts.stage0_audit import FOCAL, NewickParser, parse_gene_rows, parse_newick, split_for_quartet, topology_label


class Stage0AuditTests(unittest.TestCase):
    def test_newick_records_and_branch_lengths(self):
        text = "((S1_a1:0.1,S1_a2:0.2):0.3,(S2_a1:0.1,S2_a2:0.2):0.3);"
        _, taxa, weighted, resolved = parse_newick(text)
        self.assertEqual(taxa, {"S1_a1", "S1_a2", "S2_a1", "S2_a2"})
        self.assertTrue(weighted)
        self.assertTrue(resolved)

    def test_quartet_labels_are_fixed(self):
        q = tuple(FOCAL.values())
        tree = NewickParser("((Mv-cat-C212-A1,Mv-cat-C212-A2),(Mv-sup-1065-A1,Mv-sup-1065-A2));").parse()
        self.assertEqual(topology_label(split_for_quartet(tree, q)), "q1")
        tree = NewickParser("((Mv-cat-C212-A1,Mv-sup-1065-A1),(Mv-cat-C212-A2,Mv-sup-1065-A2));").parse()
        self.assertEqual(topology_label(split_for_quartet(tree, q)), "q2")
        tree = NewickParser("((Mv-cat-C212-A1,Mv-sup-1065-A2),(Mv-cat-C212-A2,Mv-sup-1065-A1));").parse()
        self.assertEqual(topology_label(split_for_quartet(tree, q)), "q3")

    def test_missing_or_malformed_tree_is_not_silently_valid(self):
        with self.assertRaises(ValueError):
            parse_newick("((a,b),(c,d)")

    def test_downloaded_inventory_and_stratum_identifiers(self):
        path = Path("data/microbotryum/raw/gene_trees")
        if not path.exists():
            self.skipTest("local Dryad gene_trees file is absent")
        rows = parse_gene_rows(path)
        self.assertEqual(len(rows), 260)
        self.assertTrue(all(len(row["taxa"]) == 6 for row in rows))
        self.assertTrue(all(row["stratum"] and row["gene_number"] for row in rows))
        self.assertEqual(len({row["stratum"] for row in rows}), 11)


if __name__ == "__main__":
    unittest.main()
