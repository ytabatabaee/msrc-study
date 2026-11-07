import sys, unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from heliconius_pipeline import parse_newick, split_for_quartet, distances, QUARTET, BASELINE, ALTERNATIVE

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

if __name__ == '__main__': unittest.main()
