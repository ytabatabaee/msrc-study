import importlib.util,unittest
from pathlib import Path
P=Path(__file__).resolve().parents[1]/'scripts'/'run_analysis.py'
spec=importlib.util.spec_from_file_location('quartets',P);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
g=[{'a'},{'b'},{'c'},{'d'}];mp={x:x for x in 'abcd'}
class QuartetScoringTests(unittest.TestCase):
 def test_q1_split_and_unique_count(self):
  x=m.tree_stats('((a:1,b:1):1,(c:1,d:1):1);',g,mp)
  self.assertEqual(x['possible'],1);self.assertEqual(x['resolved'],1);self.assertEqual(x['unresolved'],0)
  self.assertEqual((x['q1'],x['q2'],x['q3']),(1.,0.,0.))
 def test_q2_q3_splits(self):
  x=m.tree_stats('((a:1,c:1):1,(b:1,d:1):1);',g,mp);self.assertEqual((x['q1'],x['q2'],x['q3']),(0.,1.,0.))
  x=m.tree_stats('((b:1,c:1):1,(a:1,d:1):1);',g,mp);self.assertEqual((x['q1'],x['q2'],x['q3']),(0.,0.,1.))
 def test_star_is_unresolved(self):
  x=m.tree_stats('(a:1,b:1,c:1,d:1);',g,mp);self.assertEqual(x['resolved'],0);self.assertEqual(x['unresolved'],1)
 def test_extra_leaves_do_not_multiply_quartet(self):
  gg=[{'a','a2'},{'b','b2'},{'c','c2'},{'d','d2'}];mm={x:x for x in ['a','a2','b','b2','c','c2','d','d2']}
  x=m.tree_stats('(((a:1,a2:1):1,(b:1,b2:1):1):1,((c:1,c2:1):1,(d:1,d2:1):1):1);',gg,mm)
  self.assertEqual(x['resolved'],x['possible']);self.assertAlmostEqual(sum([x['q1'],x['q2'],x['q3']]),1.0)
if __name__=='__main__':unittest.main(verbosity=2)
