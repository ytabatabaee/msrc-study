import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import analyze_topologies as a


class TopologyMappingTests(unittest.TestCase):
    def test_published_topology_mapping_preserves_other_classes(self):
        self.assertEqual(a.SPECIES, ('topo3',))
        self.assertEqual(set(a.RECENT), {'topo9', 'topo12', 'topo13'})
        self.assertEqual(a.ANCIENT, ('topo10',))
        self.assertEqual(a.topology_class('topo3'), 'species')
        self.assertEqual(a.topology_class('topo9'), 'recent')
        self.assertEqual(a.topology_class('topo10'), 'ancient')
        self.assertEqual(a.topology_class('topo15'), 'other')


class CoordinateTests(unittest.TestCase):
    def test_coordinate_alignment_uses_keys(self):
        x = [{'scaffold':'LG1','start':'1','end':'100000'},
             {'scaffold':'LG1','start':'100001','end':'200000'}]
        y = list(reversed(x))
        self.assertEqual(a.validate_coordinate_keys((x,y)), {('LG1',1,100000),('LG1',100001,200000)})

    def test_coordinate_mismatch_and_duplicates_fail(self):
        x = [{'scaffold':'LG1','start':'1','end':'100000'}]
        y = [{'scaffold':'LG1','start':'100001','end':'200000'}]
        with self.assertRaises(ValueError): a.validate_coordinate_keys((x,y))
        with self.assertRaises(ValueError): a.validate_coordinate_keys((x,x+x))


class RunTests(unittest.TestCase):
    @staticmethod
    def row(chrom, start, end, value, informative=1):
        return {'scaffold':chrom,'start':start,'end':end,'q':value,'is_informative':informative}

    def test_runs_split_on_chromosome_and_physical_gap(self):
        rows=[self.row('LG1',1,100000,.8),self.row('LG1',100001,200000,.9),
              self.row('LG1',300001,400000,.9),self.row('LG2',1,100000,.9)]
        runs=a.physical_runs(rows,'q',.5)
        self.assertEqual([(r['chromosome'],r['start'],r['end'],r['n_windows']) for r in runs],
                         [('LG1',1,200000,2),('LG1',300001,400000,1),('LG2',1,100000,1)])

    def test_threshold_and_missing_window_break_runs(self):
        rows=[self.row('LG1',1,100000,.6),self.row('LG1',100001,200000,.9,0),
              self.row('LG1',200001,300000,.8),self.row('LG1',300001,400000,.4)]
        runs=a.physical_runs(rows,'q',.5)
        self.assertEqual([r['n_windows'] for r in runs],[1,1])

    def test_autocorrelation_does_not_bridge_coordinate_gap(self):
        rows=[self.row('LG1',1,100000,.1),self.row('LG1',100001,200000,.2),
              self.row('LG1',300001,400000,.9)]
        ac=a.autocorrelation(rows,'q',max_dist=500000)
        self.assertEqual(sum(r['n_pairs'] for r in ac),1)


if __name__ == '__main__':
    unittest.main()
