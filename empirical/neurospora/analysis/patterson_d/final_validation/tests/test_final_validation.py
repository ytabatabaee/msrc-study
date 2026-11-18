import csv
import math
import sys
import unittest
from pathlib import Path

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0,str(SCRIPTS))
from run_final_validation import scalar  # noqa: E402


class IndependentScalarCounterTests(unittest.TestCase):
    def test_abba_orientation(self):
        self.assertEqual(scalar([65],[71],[71],65),('polarizable',1.0,0.0))

    def test_baba_orientation(self):
        self.assertEqual(scalar([71],[65],[71],65),('polarizable',0.0,1.0))

    def test_unequal_population_frequencies(self):
        status,a,b=scalar([65,65,71],[65,71],[71],65)
        self.assertEqual(status,'polarizable')
        self.assertAlmostEqual(a,1.0/3.0)
        self.assertAlmostEqual(b,1.0/6.0)

    def test_partial_calls_at_half_group_threshold(self):
        status,a,b=scalar([65,78,71,78],[65,71],[71],65)
        self.assertEqual(status,'polarizable')
        self.assertAlmostEqual(a,.25)
        self.assertAlmostEqual(b,.25)

    def test_group_below_half_callability_is_excluded(self):
        self.assertEqual(scalar([65,78,78,78],[65,71],[71],65),('insufficient_calls',0.0,0.0))

    def test_fixed_alt_outgroup_reorients_reference_allele(self):
        self.assertEqual(scalar([71],[65],[65],71),('polarizable',1.0,0.0))

    def test_gap_and_ambiguity_are_excluded(self):
        self.assertEqual(scalar([65],[71],[71],45),('outgroup_gap',0.0,0.0))
        self.assertEqual(scalar([65],[71],[71],78),('outgroup_ambiguous',0.0,0.0))

    def test_outgroup_not_in_biallelic_ingroups(self):
        self.assertEqual(scalar([65],[71],[71],67),('outgroup_other',0.0,0.0))

    def test_multistate_site_is_excluded(self):
        self.assertEqual(scalar([65],[71],[67],65),('multiallelic',0.0,0.0))

    def test_zero_contribution_at_noninformative_pattern(self):
        status,a,b=scalar([65],[65],[71],65)
        self.assertEqual(status,'polarizable')
        self.assertEqual((a,b),(0.0,0.0))


class RawSiteAgreementTests(unittest.TestCase):
    def test_stratified_raw_site_validation_is_exact(self):
        p=Path(__file__).resolve().parents[1]/'results/estimator_validation.tsv'
        with p.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertGreaterEqual(len(rows),200)
        self.assertEqual({r['linkage_group'] for r in rows},{'LGI','LGII','LGIII','LGIV','LGV','LGVI','LGVII'})
        self.assertTrue(all(float(r['max_abs_difference'])==0 for r in rows))
        self.assertTrue(all(r['scalar_status']==r['vector_status'] for r in rows))


if __name__=='__main__':unittest.main()
