import csv
import importlib.util
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
SCRIPT=ROOT/"empirical/fire_ants_chr16/scripts/10_add_inversion_annotations.py"
spec=importlib.util.spec_from_file_location("inversion_annotations",SCRIPT)
mod=importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

class InversionAnnotationTests(unittest.TestCase):
    @staticmethod
    def read_rows(path):
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle,delimiter="\t"))

    def test_three_validated_target_axis_records(self):
        self.assertEqual(len(mod.INVERSIONS),3)
        self.assertEqual({r["inversion_id"] for r in mod.INVERSIONS},{"In(16)1","In(16)2","In(16)3"})
        for r in mod.INVERSIONS:
            self.assertEqual(r["chromosome"],"chr16")
            self.assertGreater(r["start_bp"],0)
            self.assertGreater(r["end_bp"],r["start_bp"])
            self.assertIn("Stolle Si_gnGA",r["target_assembly"])
            self.assertIn("outer coordinate envelope",r["mapping_status"])
            self.assertTrue(r["breakpoint_certainty"])

    def test_coordinate_mapping_and_reported_assembled_lengths_are_separate(self):
        lengths={r["inversion_id"]:r["assembled_bp_approx"] for r in mod.INVERSIONS}
        self.assertEqual(lengths,{"In(16)1":9_590_000,"In(16)2":840_000,"In(16)3":1_200_000})
        for r in mod.INVERSIONS:
            envelope=r["end_bp"]-r["start_bp"]+1
            self.assertGreater(envelope,r["assembled_bp_approx"])
        genomic_order=sorted(mod.INVERSIONS,key=lambda r:r["start_bp"])
        self.assertEqual([r["inversion_id"] for r in genomic_order],["In(16)3","In(16)2","In(16)1"])
        output=self.read_rows(ROOT/"empirical/fire_ants_chr16/results/chr16_inversion_intervals.tsv")
        for row in output:
            self.assertEqual(int(row["length_bp"]),int(row["end_bp"])-int(row["start_bp"])+1)

    def test_inclusive_overlap_and_partial_window_logic(self):
        self.assertEqual(mod.overlap_len(100,200,200,300),1)
        self.assertEqual(mod.overlap_len(100,199,200,300),0)
        self.assertEqual(mod.overlap_len(100,250,150,300),101)

    def test_window_coverage_counts_and_zero_inversion(self):
        p=ROOT/"empirical/fire_ants_chr16/results/chr16_inversion_coverage_summary.tsv"
        rows=self.read_rows(p)
        by={(r["inversion_id"],r["quartet_id"]):r for r in rows}
        for q in ("Q1_primary_Argentina","Q2_primary_Argentina"):
            self.assertEqual(by[("In(16)1",q)]["source_midpoint_assigned_tree_windows"],"47")
            self.assertEqual(by[("In(16)2",q)]["source_midpoint_assigned_tree_windows"],"4")
            self.assertEqual(by[("In(16)3",q)]["source_midpoint_assigned_tree_windows"],"0")
            self.assertEqual(by[("In(16)3",q)]["informative_snps_in_source_midpoint_assigned_windows"],"0")
            self.assertEqual(by[("In(16)3",q)]["fraction_inversion_length_covered_by_analyzed_gene_bases"],"NA")

    def test_descriptive_summary_marks_zero_coverage_na(self):
        p=ROOT/"empirical/fire_ants_chr16/results/chr16_inversion_specific_descriptive_summary.tsv"
        rows=self.read_rows(p)
        for r in rows:
            if r["inversion_id"]=="In(16)3":
                self.assertEqual(r["patterson_D_pooled"],"NA")
                self.assertEqual(r["mean_q_species"],"NA")
                self.assertEqual(r["status"],"no assigned windows; NA")

    def test_frozen_results_and_region_definitions_unchanged(self):
        mod.verify_frozen_inputs()
        rows=self.read_rows(ROOT/"data/fire_ants_chr16/metadata/region_manifest.tsv")
        r=next(x for x in rows if x.get("region_id")=="chr16_supergene")
        self.assertEqual((int(r["coordinate_start"]),int(r["coordinate_end"])),mod.SUPERGENE)

    def test_versioned_figures_exist_and_have_structural_labels(self):
        stems=[ROOT/"empirical/fire_ants_chr16/figures/fire_ants_main_with_inversions",
               ROOT/"empirical/fire_ants_chr16/figures/fire_ants_chr16_topology_distance_inversions",
               ROOT/"empirical/fire_ants_chr16/analysis/patterson_d/figures/figure_B_D_and_TWISST_with_inversions"]
        for stem in stems:
            for ext in ("png","pdf","svg"):
                self.assertTrue(stem.with_suffix("."+ext).is_file())
        svg=stems[1].with_suffix(".svg").read_text()
        for label in ("In(16)1","In(16)2","In(16)3","Supergene analysis span"):
            self.assertIn(label,svg)

if __name__=="__main__": unittest.main()
