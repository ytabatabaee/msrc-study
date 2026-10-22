import importlib.util
from pathlib import Path

BASE = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location("pnas48_data_audit", BASE / "scripts/pnas48_data_audit.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def test_frozen_pnas48_mapping_has_48_unique_stiller_taxa():
    rows = audit.read_tsv(BASE / "data/stage4d/jarvis48_taxa.tsv")
    assert len(rows) == 48
    assert len({r["source_2014_taxon"] for r in rows}) == 48
    assert len({r["stiller2024_taxon"] for r in rows}) == 48


def test_pnas_mask_matches_published_predicate_for_every_manifest_locus():
    rows = audit.read_tsv(BASE / "results/stage4d_locus_manifest.tsv")
    assert len(rows) == 63430
    for row in rows:
        expected = audit.published_outlier(row["chromosome"], int(float(row["start"])))
        assert (row["published_outlier_region"] == "True") == expected


def test_existing_j48_inputs_are_preserved_and_manifest_hashes_match():
    for row in audit.read_tsv(BASE / "results/stage4e_input_manifest.tsv"):
        if row["scope"] != "jarvis48":
            continue
        path = audit.ROOT / row["relative_path"]
        assert path.exists()
        assert audit.sha256(path) == row["sha256"]


def test_audit_classifies_exact_published_input_as_unverified():
    report = (BASE / "results/stage4e/pnas48_data_audit.md").read_text()
    assert "reconstructed published-design reproduction" in report
    assert "No raw `.tre` was overwritten" in report
