import importlib.util
from pathlib import Path
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/patterson_d.py"
spec = importlib.util.spec_from_file_location("patterson_d", SCRIPT)
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)

class TestPattersonD(unittest.TestCase):
    def test_abba_baba(self):
        seqs = {"P1":"A", "P2":"G", "P3":"G", "O":"A"}
        self.assertEqual(mod.count_sites(seqs), (1, 0, 0, 1))
        seqs = {"P1":"G", "P2":"A", "P3":"G", "O":"A"}
        self.assertEqual(mod.count_sites(seqs), (0, 1, 0, 1))
    def test_filters(self):
        seqs = {"P1":"N", "P2":"A", "P3":"A", "O":"A"}
        self.assertEqual(mod.count_sites(seqs)[2], 1)
    def test_pooled_ratio(self):
        self.assertAlmostEqual(mod.patterson_d(12, 4), 0.5)
        self.assertIsNone(mod.patterson_d(0, 0))
    def test_transition_filter(self):
        seqs = {"P1":"A", "P2":"G", "P3":"G", "O":"A"}
        self.assertEqual(mod.count_sites(seqs, transversions_only=True)[3], 0)
        seqs = {"P1":"A", "P2":"C", "P3":"C", "O":"A"}
        self.assertEqual(mod.count_sites(seqs, transversions_only=True)[3], 1)
    def test_independent_reference_arithmetic(self):
        # A deliberately independent one-column reference implementation.
        columns = [("A","G","G","A"),("G","A","G","A"),("A","C","C","A"),("A","A","C","A")]
        abba=sum(p1==o and p2==p3 and p2!=o for p1,p2,p3,o in columns)
        baba=sum(p2==o and p1==p3 and p1!=o for p1,p2,p3,o in columns)
        self.assertEqual((abba,baba),(2,1))
        self.assertAlmostEqual(mod.patterson_d(abba,baba),1/3)

if __name__ == "__main__": unittest.main()
