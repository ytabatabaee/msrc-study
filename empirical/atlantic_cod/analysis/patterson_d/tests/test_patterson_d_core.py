import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from patterson_d_core import (parse_gt,alt_frequency,fixed_outgroup_state,site_contributions,d_from_sums)

class TestPattersonD(unittest.TestCase):
    def test_pure_abba(self):
        a,b=site_contributions(0,1,1,0);self.assertEqual((a,b),(1,0));self.assertEqual(d_from_sums(a,b),1)
    def test_pure_baba(self):
        a,b=site_contributions(1,0,1,0);self.assertEqual((a,b),(0,1));self.assertEqual(d_from_sums(a,b),-1)
    def test_balanced(self):
        a,b=site_contributions(.5,.5,1,0);self.assertEqual(a,b);self.assertEqual(d_from_sums(a,b),0)
    def test_heterozygote_frequency(self):
        p,n=alt_frequency(['0/1','1|1','./.'],2);self.assertEqual(n,2);self.assertEqual(p,.75)
    def test_missing_population(self):
        p,n=alt_frequency(['./.','./.'],1);self.assertIsNone(p);self.assertEqual(n,0)
    def test_fixed_alt_outgroup_repolarizes(self):
        self.assertEqual(fixed_outgroup_state('1/1'),1)
        self.assertEqual(site_contributions(1,0,0,1),site_contributions(0,1,1,0))
    def test_polymorphic_outgroup_rejected(self):self.assertIsNone(fixed_outgroup_state('0/1'))
    def test_ref_alt_reversal_invariance(self):
        a,b=site_contributions(.2,.7,.9,0)
        aa,bb=site_contributions(.8,.3,.1,1)
        self.assertAlmostEqual(a,aa);self.assertAlmostEqual(b,bb)
    def test_invalid_and_multiallelic_genotypes(self):
        self.assertIsNone(parse_gt('1/2'));self.assertIsNone(parse_gt('./.'))

if __name__=='__main__':unittest.main()
