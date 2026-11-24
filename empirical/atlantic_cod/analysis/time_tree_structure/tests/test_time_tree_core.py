import importlib.util
import math
from pathlib import Path
import unittest

CORE=Path(__file__).resolve().parents[1]/'scripts'/'time_tree_core.py'
spec=importlib.util.spec_from_file_location('time_tree_core',CORE);core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)

class TimeTreeCoreTests(unittest.TestCase):
    def test_identical_trees_have_identical_mrca_vectors(self):
        nwk='((A:1,B:1):2,C:3);'
        a=core.extract_mrca_ages(nwk);b=core.extract_mrca_ages(nwk)
        self.assertTrue(a['valid']);self.assertEqual(a['pair_ages'],b['pair_ages'])
        self.assertAlmostEqual(a['pair_ages'][('A','B')],1)

    def test_distinct_rooted_topologies_and_all_signed_qs(self):
        a=core.extract_mrca_ages('((A:1,B:1):2,C:3);')
        b=core.extract_mrca_ages('((A:1,C:1):2,B:3);')
        c=core.extract_mrca_ages('((B:1,C:1):2,A:3);')
        self.assertNotEqual(a['pair_ages'],b['pair_ages'])
        ab=core.triplet_features(a['taxa'],a['pair_ages'],a['root_age'])[0]
        ac=core.triplet_features(b['taxa'],b['pair_ages'],b['root_age'])[0]
        bc=core.triplet_features(c['taxa'],c['pair_ages'],c['root_age'])[0]
        self.assertEqual(ab['topology'],'A+B|C');self.assertEqual(ab['resolution_topologies']['r1'],'A+B|C');self.assertGreater(ab['q_contrasts']['r1'],0)
        self.assertLess(ab['q_contrasts']['r2'],0);self.assertLess(ab['q_contrasts']['r3'],0)
        self.assertEqual(ac['topology'],'A+C|B');self.assertEqual(ac['resolution_topologies']['r2'],'A+C|B');self.assertGreater(ac['q_contrasts']['r2'],0)
        self.assertEqual(bc['topology'],'B+C|A');self.assertEqual(bc['resolution_topologies']['r3'],'B+C|A');self.assertGreater(bc['q_contrasts']['r3'],0)
        self.assertGreater(ab['separation_by_resolution']['r1'],0)
        self.assertEqual(ab['separation_by_resolution']['r2'],0)

    def test_uniform_rescaling_changes_absolute_not_relative_vectors(self):
        a=core.extract_mrca_ages('((A:1,B:1):2,C:3);')
        b=core.extract_mrca_ages('((A:7,B:7):14,C:21);')
        self.assertTrue(a['valid'] and b['valid'])
        self.assertAlmostEqual(b['root_age'],7*a['root_age'])
        for pair in a['pair_ages']:
            self.assertAlmostEqual(a['pair_ages'][pair]/a['root_age'],b['pair_ages'][pair]/b['root_age'])
            self.assertAlmostEqual(b['pair_ages'][pair],7*a['pair_ages'][pair])

    def test_polytomy_and_zero_internal_branch_stay_unresolved(self):
        f=core.extract_mrca_ages('(A:3,B:3,C:3);')
        self.assertTrue(f['valid'])
        t=core.triplet_features(f['taxa'],f['pair_ages'],f['root_age'])[0]
        self.assertTrue(t['unresolved']);self.assertEqual(t['topology'],'unresolved')
        self.assertEqual(t['temporal_separation'],0)

    def test_non_ultrametric_tree_is_not_repaired(self):
        f=core.extract_mrca_ages('((A:1,B:1):2,C:2);')
        self.assertFalse(f['valid']);self.assertIn('non_ultrametric',f['issue'])

    def test_missing_and_inconsistent_taxa_are_invalid(self):
        missing=core.extract_mrca_ages('((A:1,B):2,C:3);')
        wrong=core.extract_mrca_ages('((A:1,B:1):2,C:3);',expected_taxa=('A','B','D'))
        self.assertFalse(missing['valid']);self.assertIn('missing_branch_length',missing['issue'])
        self.assertFalse(wrong['valid']);self.assertIn('inconsistent_taxon_set',wrong['issue'])

    def test_zero_root_age_and_negative_branch_are_invalid(self):
        zero=core.extract_mrca_ages('(A:0,B:0,C:0);')
        neg=core.extract_mrca_ages('((A:-1,B:1):2,C:3);')
        self.assertFalse(zero['valid']);self.assertIn('zero_or_invalid_root_age',zero['issue'])
        self.assertFalse(neg['valid']);self.assertIn('negative_branch_length',neg['issue'])

    def test_root_stem_does_not_change_mrca_ages(self):
        a=core.extract_mrca_ages('((A:1,B:1):2,C:3);')
        b=core.extract_mrca_ages('((A:1,B:1):2,C:3):100;')
        self.assertTrue(b['valid']);self.assertAlmostEqual(a['root_age'],b['root_age'])
        self.assertEqual(a['pair_ages'],b['pair_ages'])

if __name__=='__main__':unittest.main()
