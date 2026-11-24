import csv
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
def rows(path):
    with open(path,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
class DatasetOutputTests(unittest.TestCase):
    def test_fire_ant_has_all_windows_and_matched_contrast(self):
        p=ROOT/'empirical/fire_ants_chr16/analysis/tajima_d/results/fire_ant_chr16_tajima_d_windows.tsv'
        r=rows(p); self.assertEqual(len(r),96)
        self.assertTrue(all(x['matched_Delta_D_unweighted']!='' for x in r))
        self.assertTrue(all((x['n_SB'],x['n_Sb'],x['pooled_n_haplotypes'])==('78','30','108') for x in r))
        self.assertTrue(all(x['n_sites_matched_complete_108']!='' for x in r))
    def test_cod_does_not_invent_outside_snp_records(self):
        p=ROOT/'empirical/atlantic_cod/analysis/tajima_d/results/atlantic_cod_tajima_d_source_audit.tsv'
        r=rows(p); outside=[x for x in r if x['subset']=='outside']
        self.assertEqual(len(outside),1); self.assertEqual(outside[0]['n_records_used'],'0')
    def test_anopheles_full_arm_tracks_and_callable_denominators(self):
        p=ROOT/'empirical/anopheles_2la/analysis/tajima_d/results/anopheles_tajima_d_50kb.tsv'
        r=rows(p); self.assertGreater(len(r),9000)
        self.assertTrue(all(x['callable_bases']!='' or x['D_undefined_reason']!='' for x in r))
        self.assertTrue(any(x['pi_per_callable_base'] not in ('','NA') for x in r))
        self.assertTrue(all(int(x['n_haplotypes'])>=4 and int(x['n_haplotypes'])%2==0 for x in r))
if __name__=='__main__':unittest.main()
