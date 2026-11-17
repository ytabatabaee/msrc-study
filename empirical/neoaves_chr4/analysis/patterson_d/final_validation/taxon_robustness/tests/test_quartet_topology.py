import io
import sys
import unittest
from pathlib import Path
from Bio import Phylo
HERE=Path(__file__).resolve()
sys.path.insert(0,str(HERE.parents[1]/'scripts'))
from compare_exact_gene_trees import quartet

class ExactQuartetTests(unittest.TestCase):
    q=('A','B','C','O')
    def test_all_three_unrooted_splits(self):
        cases=[('((A,B),(C,O));','q1'),('((A,C),(B,O));','q2'),('((B,C),(A,O));','q3')]
        for nwk,want in cases:
            self.assertEqual(quartet(Phylo.read(io.StringIO(nwk),'newick'),self.q),want)
    def test_star_is_unresolved(self):
        self.assertIsNone(quartet(Phylo.read(io.StringIO('(A,B,C,O);'),'newick'),self.q))

if __name__=='__main__':unittest.main()
