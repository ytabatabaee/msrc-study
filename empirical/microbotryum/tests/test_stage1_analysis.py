import itertools
import unittest

from empirical.microbotryum.scripts.stage0_audit import NewickParser, split_for_quartet
from empirical.microbotryum.scripts.stage1_analysis import (
    QUARTET,
    TAXA_ORDER,
    all_assignments,
    meaningful_labels,
    primary_metrics,
    primary_split_labels,
)


class Stage1AnalysisTests(unittest.TestCase):
    def test_all_fifteen_quartets_for_six_taxa(self):
        self.assertEqual(len(list(itertools.combinations(TAXA_ORDER, 4))), 15)

    def test_root_independent_primary_topology(self):
        first = NewickParser("((Mv-cat-C212-A1:1,Mv-cat-C212-A2:1):2,(Mv-sup-1065-A1:1,Mv-sup-1065-A2:1):2);").parse()
        second = NewickParser("(Mv-cat-C212-A1:1,(Mv-cat-C212-A2:1,(Mv-sup-1065-A1:1,Mv-sup-1065-A2:1):2):2);").parse()
        self.assertEqual(primary_split_labels(split_for_quartet(first, QUARTET)), "q1_species")
        self.assertEqual(primary_split_labels(split_for_quartet(second, QUARTET)), "q1_species")

    def test_fixed_primary_labels_and_meaningful_quartet(self):
        meaningful, labels = meaningful_labels(tuple(QUARTET))
        self.assertTrue(meaningful)
        self.assertCountEqual(labels, ["q1_species", "q2_mating_type", "q3_third"])

    def test_all_quartet_summary_normalizes_to_one(self):
        tree = NewickParser("(((Mv-cat-C212-A1:1,Mv-cat-C212-A2:1):1,(Mv-lag-1253-A1:1,Mv-lag-1253-A2:1):1):1,(Mv-sup-1065-A1:1,Mv-sup-1065-A2:1):1);").parse()
        row = {"root": tree, "taxa": set(TAXA_ORDER), "stratum": "synthetic", "gene_label": "synthetic-g1", "gene_number": 1}
        assignments, summaries = all_assignments([row])
        self.assertEqual(len(assignments), 15)
        self.assertEqual(len(summaries), 15)
        for summary in summaries:
            self.assertAlmostEqual(float(summary["q1_or_resolution1"]) + float(summary["q2_or_resolution2"]) + float(summary["q3_or_resolution3"]), 1.0)

    def test_patristic_and_internal_branch_stats(self):
        tree = NewickParser("((Mv-cat-C212-A1:1,Mv-cat-C212-A2:1):2,(Mv-sup-1065-A1:1,Mv-sup-1065-A2:1):2);").parse()
        metrics = primary_metrics(tree)
        self.assertAlmostEqual(metrics["B"], 14 / 3)
        self.assertAlmostEqual(metrics["internal_branch"], 4.0)
        self.assertAlmostEqual(metrics["A_opposite_minus_same"], 0.0)


if __name__ == "__main__":
    unittest.main()
