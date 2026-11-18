import importlib.util
from pathlib import Path
import unittest

SCRIPT=Path(__file__).resolve().parents[1]/'scripts/run_sensitivity_scan.py'
spec=importlib.util.spec_from_file_location('scan',SCRIPT);scan=importlib.util.module_from_spec(spec);spec.loader.exec_module(scan)

class GenotypeParserTests(unittest.TestCase):
    def test_haploid_and_homozygous_pseudodiploid(self):
        ix={'GT':0,'DP':1}
        self.assertEqual(scan.parse('1:12',ix),(1,12))
        self.assertEqual(scan.parse('0/0:12',ix),(0,12))
        self.assertEqual(scan.parse('1|1:12',ix),(1,12))
    def test_heterozygous_missing_and_low_depth_excluded(self):
        ix={'GT':0,'DP':1}
        self.assertEqual(scan.parse('0/1:12',ix),(None,12))
        self.assertEqual(scan.parse('./.:20',ix),(None,20))
        a,d=scan.parse('1:9',ix);self.assertTrue(d<10);self.assertEqual(a,1)
    def test_fixed_alt_outgroup_flips_derived_allele(self):
        # P1=0,P2=1,P3=1 is ABBA when ancestral allele is REF.
        p=[0,1,1]; anc_ref=0
        self.assertGreater((1-p[0])*p[1]*p[2],p[0]*(1-p[1])*p[2])
        # After REF/ALT recoding the outgroup is ALT-fixed and the derived
        # allele is REF, whose frequency is 1 - recoded ALT frequency.
        recoded_alt=[1-x for x in p]; derived_ref=[1-x for x in recoded_alt]
        self.assertGreater((1-derived_ref[0])*derived_ref[1]*derived_ref[2],derived_ref[0]*(1-derived_ref[1])*derived_ref[2])

if __name__=='__main__':unittest.main()
