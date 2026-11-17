import io
import sys
import unittest
from pathlib import Path
HERE=Path(__file__).resolve()
sys.path.insert(0,str(HERE.parents[1]/'scripts'))
from compute_taxon_sensitivity import count_independent,parse_fasta

class RawAlignmentReferenceTests(unittest.TestCase):
    def test_fasta_taxa_case_normalization_and_duplicate_detection(self):
        data=b'>P1 sample\na\n>P2\nC\n>P1 duplicate\nG\n'
        seqs,dups=parse_fasta(io.BytesIO(data))
        self.assertEqual(seqs,{'P1':'G','P2':'C'})
        self.assertEqual(dups,['P1'])

    def test_exact_orientation_ambiguity_multistate_and_tv(self):
        seqs={'P1':'AGACCA-','P2':'GACACAA','P3':'GGCCGAA','O':'AAAAAAA'}
        all_sites=count_independent(seqs,False)
        tv=count_independent(seqs,True)
        self.assertEqual((all_sites['ABBA'],all_sites['BABA']),(2,2))
        self.assertEqual(all_sites['informative_sites'],4)
        self.assertEqual((tv['ABBA'],tv['BABA']),(1,1))
        self.assertEqual(tv['informative_sites'],2)
        self.assertEqual(all_sites['ambiguous_sites'],1)
        self.assertEqual(all_sites['multistate_nonpattern_sites'],1)
        self.assertEqual(all_sites['alignment_length'],7)

    def test_alignment_length_equality_is_explicitly_detectable(self):
        seqs={'P1':'AC','P2':'A','P3':'AC','O':'AC'}
        self.assertEqual(len({len(x) for x in seqs.values()}),2)

if __name__=='__main__':unittest.main()
