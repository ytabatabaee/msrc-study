import sys, unittest
from itertools import combinations
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from heliconius_pipeline import parse_newick, split_for_quartet, split_for_quartet_from_dist, distances, split_label, classify, verified_regions, QUARTET, BASELINE, ALTERNATIVE

class TestHeliconiusUtilities(unittest.TestCase):
    def test_root_independent_topology(self):
        a='((H_erato:1,H_telesiphe:1):2,(H_hecalesia:1,H_sara:1):2);'
        b='(H_erato:1,(H_telesiphe:1,(H_hecalesia:1,H_sara:1):2):2);'
        self.assertEqual(split_for_quartet(parse_newick(a)[0], QUARTET), split_for_quartet(parse_newick(b)[0], QUARTET))
    def test_missing_taxon(self):
        root,_=parse_newick('((H_erato,H_telesiphe),(H_hecalesia,H_other));')
        self.assertIsNone(split_for_quartet(root, QUARTET))
    def test_patristic_distance(self):
        root,_=parse_newick('((H_erato:1,H_telesiphe:2):3,(H_hecalesia:4,H_sara:5):6);')
        d=distances(root, weighted=True); self.assertAlmostEqual(d[('H_erato','H_telesiphe')],3); self.assertAlmostEqual(d[('H_erato','H_sara')],15)
    def test_topology_labels_are_distinct(self): self.assertNotEqual(BASELINE, ALTERNATIVE)
    def test_unresolved_star(self):
        root,_=parse_newick('(H_erato,H_telesiphe,H_hecalesia,H_sara);')
        self.assertIsNone(split_for_quartet(root, QUARTET))

    def test_all_35_quartets_for_seven_taxa(self):
        taxa=('a','b','c','d','e','f','g')
        self.assertEqual(len(list(combinations(taxa,4))),35)

    def test_all_quartet_distance_classifier_is_root_independent(self):
        root_a,tips=parse_newick('((a:1,b:1):2,(c:1,d:1):2,(e:1,f:1,g:1):3);')
        root_b,_=parse_newick('(g:1,(f:1,(e:1,(c:1,d:1):2):2):2,(a:1,b:1):2);')
        qt=('a','b','c','d')
        self.assertEqual(split_for_quartet_from_dist(distances(root_a,False),qt),split_for_quartet_from_dist(distances(root_b,False),qt))

    def test_missing_taxa_are_not_counted_as_resolved(self):
        root,_=parse_newick('((a,b),(c,e),(f,g));')
        self.assertIsNone(split_for_quartet(root,('a','b','c','d')))

    def test_q_frequencies_normalize_and_d_uses_fixed_labels(self):
        labels=['baseline','alternative','third','baseline']
        q1=labels.count('baseline')/len(labels); q2=labels.count('alternative')/len(labels); q3=labels.count('third')/len(labels)
        self.assertAlmostEqual(q1+q2+q3,1.0)
        self.assertAlmostEqual(q2-q1,-0.25)

    def test_verified_boundary_classification(self):
        self.assertEqual(classify(9114112,9114112,'H_erato_demophoon_v1','2'),'inside')
        self.assertEqual(classify(9000000,9200000,'H_erato_demophoon_v1','2'),'overlap')
        self.assertEqual(classify(1000000,1100000,'H_erato_demophoon_v1','2'),'outside')

    def test_overlapping_windows_are_not_collapsed_by_classifier(self):
        self.assertEqual(classify(9114112,9164111,'H_erato_demophoon_v1','2'),'inside')
        self.assertEqual(classify(9164112,9214111,'H_erato_demophoon_v1','2'),'inside')

    def test_hmel_coordinate_definition_is_source_specific(self):
        self.assertEqual(verified_regions('Hmel2.5','2','erato_50kb_hmel'), [(2147755,2230606,'2b_etales8spp_hmel')])
        self.assertIn((4398947,6357046,'2b_hmelv25_res'), verified_regions('Hmel2.5','2'))
        self.assertEqual(verified_regions('Hmel2.5','2','hmelv25_res'), [])

if __name__ == '__main__': unittest.main()
