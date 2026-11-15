import io
import math
import unittest

from Bio import Phylo

from empirical.neurospora.scripts.run_stage1 import classify_split, normalize_weights, summarise


class Stage1ParserTests(unittest.TestCase):
    def test_twisst_subtree_counts_normalize_to_frequencies(self):
        self.assertEqual(normalize_weights([8.0, 112.0, 0.0]), [1/15, 14/15, 0.0])
        self.assertTrue(all(map(lambda x: x != x, normalize_weights([float("nan")] * 3))))

    def test_quartet_summary_reports_cu_only_for_symmetric_discordants(self):
        symmetric = [{"q1": .8, "q2": .1, "q3": .1, "start": 1, "end": 2}]
        self.assertAlmostEqual(summarise(symmetric, "x")["MSC_CU_if_symmetric"], -math.log(.3))
        asymmetric = [{"q1": .2, "q2": .1, "q3": .7, "start": 1, "end": 2}]
        self.assertTrue(summarise(asymmetric, "x")["MSC_CU_if_symmetric"] != summarise(asymmetric, "x")["MSC_CU_if_symmetric"])

    def test_quartet_split_labels_match_published_topologies(self):
        q1 = Phylo.read(io.StringIO("((CJ03:1,CJ04:1):2,(Ncra:1,8817:1):2);"), "newick")
        q2 = Phylo.read(io.StringIO("((CJ03:1,Ncra:1):2,(CJ04:1,8817:1):2);"), "newick")
        q3 = Phylo.read(io.StringIO("((CJ03:1,8817:1):2,(CJ04:1,Ncra:1):2);"), "newick")
        self.assertEqual(classify_split(q1, "CJ03", "CJ04", "Ncra", "8817")[0], 0)
        self.assertEqual(classify_split(q2, "CJ03", "CJ04", "Ncra", "8817")[0], 1)
        self.assertEqual(classify_split(q3, "CJ03", "CJ04", "Ncra", "8817")[0], 2)

    def test_incomplete_quartet_is_not_classified(self):
        tree = Phylo.read(io.StringIO("((CJ03,CJ04),8817);"), "newick")
        self.assertIsNone(classify_split(tree, "CJ03", "CJ04", "Ncra", "8817"))


if __name__ == "__main__":
    unittest.main()
