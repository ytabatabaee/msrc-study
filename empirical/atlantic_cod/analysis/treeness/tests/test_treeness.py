import importlib.util
import math
import unittest
from pathlib import Path

SCRIPT=Path(__file__).resolve().parents[1]/'scripts/run_treeness.py'
spec=importlib.util.spec_from_file_location('run_treeness',SCRIPT)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)

class TreenessTests(unittest.TestCase):
    def test_manual_tree_with_zero_length_and_root_stem(self):
        # Internal edge sum = 2; all non-root edges = 2 + 1 + 0 + 3 + 4 = 10.
        # Root stem (9) is excluded.
        tree='((A:1,B:0):2,(C:3,D:4):0):9;'
        got=mod.treeness(tree,set('ABCD'))
        self.assertAlmostEqual(got['internal_branch_sum'],2.0)
        self.assertAlmostEqual(got['total_branch_sum'],10.0)
        self.assertAlmostEqual(got['treeness'],.2)
        self.assertEqual(got['n_zero_length_edges'],2)
        self.assertEqual(got['root_branch_length'],9)

    def test_star_tree_treeness_zero(self):
        got=mod.treeness('(A:1,B:2,C:3);',set('ABC'))
        self.assertEqual(got['treeness'],0)

    def test_missing_branch_invalid(self):
        with self.assertRaisesRegex(ValueError,'missing_branch_lengths'):
            mod.treeness('(A,B:1);',set('AB'))

    def test_negative_branch_invalid(self):
        with self.assertRaisesRegex(ValueError,'negative_branch_length'):
            mod.treeness('(A:-1,B:1);',set('AB'))

    def test_taxon_set_mismatch(self):
        with self.assertRaisesRegex(ValueError,'taxon_set_mismatch'):
            mod.treeness('(A:1,C:1);',set('AB'))

if __name__=='__main__':unittest.main()
