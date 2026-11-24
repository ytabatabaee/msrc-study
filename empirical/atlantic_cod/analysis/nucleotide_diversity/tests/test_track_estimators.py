import importlib.util
from pathlib import Path
import csv
import unittest
SCRIPT=Path(__file__).resolve().parents[1]/'scripts/estimate_tracks.py'
spec=importlib.util.spec_from_file_location('track_estimator',SCRIPT)
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
BASE=Path(__file__).resolve().parents[1]
class TrackTests(unittest.TestCase):
    def test_manual_finite_sample_and_dxy(self):
        px,py,d=mod.site_stats(1,4,4,4)
        self.assertAlmostEqual(px,.5);self.assertAlmostEqual(py,0);self.assertAlmostEqual(d,.75)
        self.assertAlmostEqual(mod.hudson_fst(px,py,d),2/3)
    def test_low_sample_site_is_ineligible(self):
        self.assertIsNone(mod.site_stats(1,2,1,4))
    def test_missing_genotype_is_not_called(self):
        self.assertIsNone(mod.call_ok('./.:.:.:.:.', ['GT','AD','DP','GQ','PL']))
    def test_fst_zero_denominator_is_undefined(self):
        self.assertEqual(mod.hudson_fst(.1,.1,0),'NA')
    def test_classification_does_not_invent_heterozygotes(self):
        p=BASE/'results/individual_arrangement_assignments.tsv'
        if not p.exists():self.skipTest('run audit_and_report.py first')
        with p.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertEqual(set(r['arrangement_genotype'] for r in rows),{'XX','YY'})
        self.assertEqual(len(rows),96)
    def test_outside_cohort_comparison_is_population_grouped(self):
        p=BASE/'results/inside_outside_cohort_group_windows.tsv'
        if not p.exists():self.skipTest('run audit_and_report.py first')
        with p.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        outside=[r for r in rows if r['comparison_scope']=='outside']
        self.assertTrue(outside)
        self.assertTrue(all('population-group comparison only' in r['interpretation'] for r in outside))
        self.assertTrue(all(r['pi_X_group_SNP']!='NA' for r in outside))
        boundary=[r for r in rows if r['comparison_scope']=='boundary']
        self.assertTrue(all(r['pi_X_group_SNP']=='NA' for r in boundary))

    def test_outside_windows_remain_unassigned(self):
        p=BASE/'results/atlantic_cod_pi_dxy_fst_windows.tsv'
        if not p.exists():self.skipTest('run audit_and_report.py first')
        with p.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        self.assertTrue(all(r['pi_X_SNP']=='NA' for r in rows if r['inside_suppressed_recombination']=='false'))
if __name__=='__main__':unittest.main()
