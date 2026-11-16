import unittest
from pathlib import Path
import sys

SCRIPT=Path(__file__).resolve().parents[1]/"scripts"/"run_stage2.py"
sys.path.insert(0,str(SCRIPT.parent))
import run_stage2 as d
import csv

class Stage2EstimatorTests(unittest.TestCase):
    def test_genotype_haploid_diploid_and_ambiguous_calls(self):
        self.assertEqual(d.vcf_gt("1:12","GT:DP"),(1,"called"))
        self.assertEqual(d.vcf_gt("0/0:12","GT:DP"),(0,"called"))
        self.assertEqual(d.vcf_gt("1|1:12","GT:DP"),(1,"called"))
        self.assertEqual(d.vcf_gt("0/1:12","GT:DP"),(None,"ambiguous"))
        self.assertEqual(d.vcf_gt("./.:12","GT:DP"),(None,"missing"))
        self.assertEqual(d.vcf_gt("1:9","GT:DP"),(None,"missing"))

    def test_pure_abba_baba_balanced_and_fixed_alt_outgroup(self):
        def estimate(raw, anc, outn=4):
            c=d.init_counter(); n=tuple(len(x) for x in raw)
            d.add_site(c,n,outn,raw,[anc]*outn)
            return c
        abba=estimate([[0],[1],[1]],0)
        self.assertEqual((abba["abba"],abba["baba"]),(1.0,0.0))
        baba=estimate([[1],[0],[1]],0)
        self.assertEqual((baba["abba"],baba["baba"]),(0.0,1.0))
        balanced=estimate([[0,1],[0,1],[1,1]],0)
        self.assertEqual((balanced["abba"],balanced["baba"]),(.25,.25))
        fixed_alt=estimate([[1],[0],[0]],1)
        self.assertEqual((fixed_alt["abba"],fixed_alt["baba"]),(1.0,0.0))

    def test_ref_alt_reversal_is_invariant(self):
        # Recode each allele and reverse the fixed outgroup; derived frequencies stay identical.
        a=d.init_counter();b=d.init_counter()
        d.add_site(a,(2,2,2),4,[[0,0],[1,1],[1,0]],[0,0,0,0])
        d.add_site(b,(2,2,2),4,[[1,1],[0,0],[0,1]],[1,1,1,1])
        self.assertAlmostEqual(a["abba"],b["abba"])
        self.assertAlmostEqual(a["baba"],b["baba"])

    def test_all_local_tree_quartet_frequencies_are_normalized(self):
        p=d.BASE/"results/induced_quartet_topologies.tsv"
        if not p.exists(): self.skipTest("Stage 2 topology output not generated")
        with p.open() as f:
            for r in csv.DictReader(f,delimiter="\t"):
                if r["D_topology"]=="NA":continue
                self.assertAlmostEqual(sum(float(r[k]) for k in ("q1","q2","q3")),1.0,places=10)
                self.assertEqual(int(r["n_resolved_sample_quartets"]),int(r["n_possible_sample_quartets"]))

    def test_transversion_sensitivity_preserves_supergene_direction(self):
        p=d.BASE/"results/transversions_inside_outside.tsv"
        if not p.exists(): self.skipTest("Transversion sensitivity not generated")
        with p.open() as f:
            rs={r["quartet_id"]:r for r in csv.DictReader(f,delimiter="\t")}
        for q in ("Q1_primary_Argentina","Q2_primary_Argentina"):
            self.assertGreater(float(rs[q]["D_supergene"]),0.7)
            self.assertGreater(float(rs[q]["Delta_D_supergene_minus_chr16outside"]),0.7)

    def test_strict_callability_is_marked_underpowered(self):
        p=d.BASE/"results/strict_inside_outside.tsv"
        if not p.exists(): self.skipTest("Strict sensitivity not generated")
        with p.open() as f:
            rs={r["quartet_id"]:r for r in csv.DictReader(f,delimiter="\t")}
        self.assertEqual(int(rs["Q1_primary_Argentina"]["informative_supergene"]),0)
        self.assertEqual(int(rs["Q2_primary_Argentina"]["informative_supergene"]),2)
        self.assertEqual(int(rs["Q2_primary_Argentina"]["informative_chr16A_plus_chr16B"]),0)

if __name__=="__main__":unittest.main()
