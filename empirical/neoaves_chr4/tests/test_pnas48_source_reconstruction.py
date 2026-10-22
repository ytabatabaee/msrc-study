import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "empirical/neoaves_chr4/scripts/pnas48_source_reconstruction.py"
spec = importlib.util.spec_from_file_location("pnas48_source_reconstruction", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class PNAS48SourceReconstructionTests(unittest.TestCase):
    def test_header_locus_ids_are_unique_and_complete(self):
        rows = mod.header_rows()
        ids = [row["locus_id"] for row in rows]
        self.assertEqual(len(rows), 63430)
        self.assertEqual(len(ids), len(set(ids)))

    def test_source_archive_checksums(self):
        self.assertEqual(mod.sha256(mod.STILLER / "63430.named.gene.trees.gz"), mod.EXPECTED["named_sha256"])
        self.assertEqual(mod.sha256(mod.STILLER / "63430.aLRT-0.95-collapsed.gene.trees.gz"), mod.EXPECTED["collapsed_sha256"])

    def test_frozen_outlier_mask_count(self):
        rows = mod.header_rows()
        excluded = {
            row["locus_id"] for row in rows
            if row["chromosome"] == "chr4" and (
                25030000 <= int(row["start"]) <= 32670000
                or 33510000 <= int(row["start"]) <= 34470000
                or 44130000 <= int(row["start"]) <= 56810000
            )
        }
        self.assertEqual(len(excluded), 1431)

    def test_existing_j48_manifest_is_untouched_by_source_audit(self):
        rows = mod.existing_j48_rows()
        self.assertEqual(len(rows), 7)
        self.assertTrue(all(row["sha256"] for row in rows))


if __name__ == "__main__":
    unittest.main()
