import csv
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / "empirical/fire_ants_chr16/analysis/patterson_d"

def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))

class Stage1AuditTests(unittest.TestCase):
    def test_vcf_manifest_integrity_and_tip_coverage(self):
        m = json.loads((BASE / "results/stage1_manifest.json").read_text())
        self.assertTrue(m["vcf_sha256_match"])
        self.assertEqual(m["vcf_sample_count"], 368)
        self.assertEqual(m["tree_tips_exactly_in_vcf"], 267)
        self.assertEqual(m["window_genes_in_supplementary_data_2"], 852)
        self.assertEqual(m["window_gene_scaffolds_in_vcf_header"], 41)
        self.assertTrue(m["record_scan_complete"])

    def test_scanned_records_have_exact_window_gene_assignments(self):
        s = json.loads((BASE / "results/vcf_scan_summary.json").read_text())
        self.assertEqual(s["vcf_variant_records"], 2650251)
        self.assertEqual(s["variant_records_in_213_tree_window_busco_genes"], 311890)
        self.assertEqual(s["assignment_scope"], "exact gene-body intervals; no 1-kb flanks or unsupported scaffold-to-chromosome coordinate lift")

    def test_frozen_quartets_are_rooted_and_membership_counts_reconcile(self):
        designs = rows(BASE / "config/frozen_quartets.tsv")
        members = rows(BASE / "config/frozen_sample_membership.tsv")
        self.assertEqual(len(designs), 8)
        for d in designs:
            self.assertEqual(d["baseline_sister_check"], "PASS")
            self.assertEqual(d["outgroup_root_check"], "PASS")
            sub = [r for r in members if r["quartet_id"] == d["quartet_id"]]
            for role in ("P1", "P2", "P3", "O"):
                self.assertEqual(sum(r["role"] == role for r in sub), int(d[f"n_{role}"]))
            ids = [r["sample_id"] for r in sub]
            self.assertEqual(len(ids), len(set(ids)))

    def test_primary_sample_membership_matches_stage2_audit_filters(self):
        audit = {r["vcf_sample_id"]:r for r in rows(BASE / "results/sample_mapping_audit.tsv")}
        primary = [r for r in rows(BASE / "config/frozen_sample_membership.tsv") if r["quartet_id"] in {"Q1_primary_Argentina", "Q2_primary_Argentina"}]
        for item in primary:
            r = audit[item["sample_id"]]
            self.assertEqual(r["local_tree_tip"], "1")
            self.assertEqual(r["random_allele_treatment"], "0")
            self.assertEqual(r["pooled_worker_sample"], "0")
            self.assertEqual(r["caste"], "Male")

if __name__ == "__main__":
    unittest.main()
