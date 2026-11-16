import importlib.util
import pathlib
import unittest

SCRIPT = pathlib.Path(__file__).parents[1] / 'scripts' / 'run_validation.py'
spec = importlib.util.spec_from_file_location('run_validation', SCRIPT)
rv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rv)

class PattersonDValidationMath(unittest.TestCase):
    def test_pure_abba(self):
        a,b = rv.contributions(0,1,1,0)
        self.assertEqual((a,b),(1,0))
        self.assertEqual(rv.dval(a,b),1)
    def test_pure_baba(self):
        a,b = rv.contributions(1,0,1,0)
        self.assertEqual((a,b),(0,1))
        self.assertEqual(rv.dval(a,b),-1)
    def test_balanced(self):
        a,b = rv.contributions(.5,.5,1,0)
        self.assertAlmostEqual(a,b)
        self.assertAlmostEqual(rv.dval(a,b),0)
    def test_fixed_alt_outgroup_repolarizes(self):
        abba,baba = rv.contributions(0,1,1,0)
        abba2,baba2 = rv.contributions(1,0,0,1)
        self.assertEqual((abba,baba),(abba2,baba2))
    def test_heterozygotes_enter_frequency(self):
        p,n=rv.frequency(['0/1','1/1'])
        self.assertEqual((p,n),(.75,2))
    def test_missing_does_not_enter_frequency(self):
        p,n=rv.frequency(['0/1','./.'])
        self.assertEqual((p,n),(.5,1))
    def test_heterozygous_outgroup_rejected_by_parser_caller(self):
        self.assertEqual(rv.gt_parse('0/1'),(0,1))
        self.assertIsNone(rv.gt_parse('./.'))
    def test_ref_alt_relabeling_invariance(self):
        base=rv.contributions(.2,.8,.7,0)
        reversed_labels=rv.contributions(.8,.2,.3,1)
        self.assertAlmostEqual(base[0],reversed_labels[0])
        self.assertAlmostEqual(base[1],reversed_labels[1])

if __name__=='__main__': unittest.main()
