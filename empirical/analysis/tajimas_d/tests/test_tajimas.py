import math
import unittest
from tajimas import allele_counts, constants, fixed_complete_columns, tajima_from_counts, tajima_from_sufficient_stats, weighted_delta

class TajimaTests(unittest.TestCase):
    def test_hand_calculated_biallelic_site(self):
        x=tajima_from_counts([[1,3]],4)
        a1,_,e1,_=constants(4)
        self.assertEqual(x['S'],1)
        self.assertAlmostEqual(x['pi_total'],0.5)
        self.assertAlmostEqual(x['theta_W'],1/a1)
        self.assertAlmostEqual(x['tajimas_D'],(0.5-1/a1)/math.sqrt(e1))

    def test_sufficient_statistics_match_site_vectors(self):
        x=tajima_from_counts([[1,3],[2,2],[4,0]],4)
        y=tajima_from_sufficient_stats(4,x['S'],x['pi_total'],x['n_sites_examined'])
        for key in ('S','pi_total','theta_W','tajimas_D','classical_variance'):
            self.assertAlmostEqual(x[key],y[key])

    def test_unbiased_pi_multiallelic(self):
        x=tajima_from_counts([[1,1,2],[4]],4)
        self.assertEqual(x['S'],1)
        self.assertAlmostEqual(x['pi_total'],(16-(1+1+4))/12)
        self.assertEqual(x['n_sites_examined'],2)

    def test_no_invariants_influence_D_but_countable(self):
        a=tajima_from_counts([[1,3],[2,2]],4)
        b=tajima_from_counts([[0,4],[1,3],[2,2]],4)
        self.assertAlmostEqual(a['tajimas_D'],b['tajimas_D'])
        self.assertEqual(b['S'],2)

    def test_small_n_and_no_snp_are_explicitly_undefined(self):
        self.assertEqual(tajima_from_counts([[1,2]],3)['D_undefined_reason'],'n_haplotypes<=3')
        self.assertEqual(tajima_from_counts([[4]],4)['D_undefined_reason'],'no_segregating_sites')

    def test_fixed_complete_columns_keep_constant_sample_count(self):
        rows=[[(0,), (1,), (0,)],[(1,), (None,), (0,)]]
        ix,counts=fixed_complete_columns(rows,[0,1,2],1)
        self.assertEqual(ix,[0,2]);self.assertEqual(counts,[[2],[1,1]])

    def test_ploidy_and_sample_size(self):
        self.assertEqual(allele_counts([(0,1),(1,1)],2),[1,3])
        self.assertEqual(allele_counts([(0,),(1,)],1),[1,1])
        with self.assertRaises(ValueError): allele_counts([(0,1)],1)

    def test_weighted_descriptive_contrasts(self):
        u,w=weighted_delta(1.0,0.2,0.6,8,10)
        self.assertAlmostEqual(u,0.6)
        self.assertAlmostEqual(w,1-(8*.2+10*.6)/18)

if __name__=='__main__':unittest.main()
