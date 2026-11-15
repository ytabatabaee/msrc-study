"""Synthetic checks against the exact estimator helpers used in the scan."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_patterson_d import outgroup_fixed_ref, site_contrib

class SyntheticCounts(unittest.TestCase):
    def test_pure_abba_baba(self):
        self.assertEqual(site_contrib(0,1,1),(1,0))
        self.assertEqual(site_contrib(1,0,1),(0,1))
    def test_balanced_and_uninformative(self):
        self.assertEqual(site_contrib(.5,.5,1),(.25,.25))
        self.assertEqual(site_contrib(0,0,0),(0,0))
    def test_missing_outgroup_is_not_polarized(self):
        self.assertFalse(outgroup_fixed_ref({'o1':0.0,'o2':None},['o1','o2']))
        self.assertTrue(outgroup_fixed_ref({'o1':0.0,'o2':0.0},['o1','o2']))
    def test_derived_outgroup_is_not_ancestral(self):
        self.assertFalse(outgroup_fixed_ref({'o1':0.0,'o2':1.0},['o1','o2']))

if __name__=='__main__':unittest.main()
