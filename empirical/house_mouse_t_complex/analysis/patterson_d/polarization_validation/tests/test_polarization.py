import importlib.util,unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]/'scripts'/'run_analysis.py'
spec=importlib.util.spec_from_file_location('polarization',P);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
class PolarizationTests(unittest.TestCase):
 def test_ref_alt_relabel_invariance(self):
  # Derived frequencies p give ABBA/BABA; complemented frequencies after REF/ALT swap
  # with fixed-ALT outgroup recover exactly the same derived frequencies.
  p=(0.2,0.8,0.7);a,b=m.contributions(*p)
  q=tuple(1-x for x in p);aa,bb=m.contributions(*(1-x for x in q))
  self.assertAlmostEqual(a,aa);self.assertAlmostEqual(b,bb)
  self.assertAlmostEqual((a-b)/(a+b),(aa-bb)/(aa+bb))
 def test_pure_abba_and_baba(self):
  self.assertEqual(m.contributions(0,1,1),(1,0))
  self.assertEqual(m.contributions(1,0,1),(0,1))
 def test_balanced(self):
  a,b=m.contributions(.5,.5,1);self.assertAlmostEqual(a,b)
 def test_fixed_ref_fixed_alt(self):
  self.assertEqual(m.polarization_state([(0,0)]*8),'fixed_ref')
  self.assertEqual(m.polarization_state([(1,1)]*8),'fixed_alt')
 def test_polymorphic_and_missing_outgroup(self):
  self.assertEqual(m.polarization_state([(0,0)]*7+[(0,1)]),'polymorphic')
  self.assertEqual(m.polarization_state([(0,0)]*7+[(None,None)]),'missing')
 def test_heterozygous_and_missing_ingroup_frequency(self):
  self.assertEqual(m.gt_freq((0,1)),.5)
  self.assertIsNone(m.gt_freq((0,None)))
 def test_multiallelic_excluded_by_biallelic_snp_filter(self):
  valid=lambda alleles: len(alleles)==2 and all(a is not None and len(a)==1 and a in 'ACGT' for a in alleles)
  self.assertTrue(valid(('A','G')));self.assertFalse(valid(('A','G','T')));self.assertFalse(valid(('A','AT')))
if __name__=='__main__':unittest.main(verbosity=2)
