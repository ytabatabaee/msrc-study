import sys
import unittest
from pathlib import Path
from io import StringIO

from Bio import Phylo

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "empirical/primula/scripts"))
import stage0_audit as audit


class Stage0Tests(unittest.TestCase):
    def test_split_newick_respects_quotes_and_comments(self):
        records = audit.split_newick_records("('A; taxon':1,B:2)[&x;y];(C,D);")
        self.assertEqual(len(records), 2)

    def test_tree_quality_and_branch_lengths(self):
        info = audit.inspect_newick("((A:1,B:1):2,(C:1,D:1):2);")
        self.assertEqual(info["n_taxa"], 4)
        self.assertTrue(info["branch_lengths_complete"])
        self.assertEqual(info["polytomies"], 0)

    def test_quartet_split_and_missing_taxa(self):
        tree = Phylo.read(StringIO("((A:1,B:1):1,(C:1,D:1):1);"), "newick")
        self.assertEqual(audit.quartets(tree, ("A", "B", "C", "D")), (("A", "B"), ("C", "D")))
        self.assertIsNone(audit.quartets(tree, ("A", "B", "C", "E")))

    def test_patristic_distance(self):
        tree = Phylo.read(StringIO("((A:1,B:2):3,C:4);"), "newick")
        self.assertEqual(audit.patristic(tree, "A", "B"), 3.0)

    def test_coordinate_mapping_rejects_order_only_or_missing_mapping(self):
        rows = [{"tree_index": "1", "chromosome": "chr1", "start_bp": "1", "end_bp": "50000", "mapping_status": "unverified"}]
        with self.assertRaises(ValueError):
            audit.validate_explicit_mapping(rows, 1)
        rows[0]["mapping_status"] = "verified"
        audit.validate_explicit_mapping(rows, 1)


if __name__ == "__main__":
    unittest.main()
