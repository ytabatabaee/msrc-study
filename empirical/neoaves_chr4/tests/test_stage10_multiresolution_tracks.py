from __future__ import annotations

import csv
import importlib.util
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "empirical/neoaves_chr4/scripts/10_multiresolution_tracks.py"
OUT = ROOT / "empirical/neoaves_chr4/results/stage10_multiresolution_tracks"


def module():
    spec = importlib.util.spec_from_file_location("stage10", SCRIPT)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def read(path):
    with path.open() as h:
        return list(csv.DictReader(h, delimiter="\t"))


def test_all_requested_resolution_outputs_exist():
    for resolution in (100, 250, 500):
        assert (OUT / f"binned_quartet_support_{resolution}kb.tsv").exists()
        assert (OUT / f"binned_background_distance_{resolution}kb.tsv").exists()


def test_binned_support_sums_to_one():
    for resolution in (100, 250, 500):
        values = read(OUT / f"binned_quartet_support_{resolution}kb.tsv")
        assert max(abs(float(r["support_sum"]) - 1) for r in values) < 1e-9


def test_stage9_background_rows_are_preserved():
    raw = read(ROOT / "empirical/neoaves_chr4/results/stage9_full_taxon_background/background_distance_estimators.tsv")
    assert len(raw) == 6140
    assert all("B_observed" in r for r in raw)
    assert len(read(OUT / "binned_background_distance_250kb.tsv")) == 361


def test_bin_assignment_is_deterministic_and_row_order_invariant():
    m = module()
    support = read(ROOT / "empirical/neoaves_chr4/results/stage6_genomic_tracks/window_quartet_support.tsv")
    a = m.aggregate_support(support, 250_000)
    shuffled = support[:]
    random.Random(19).shuffle(shuffled)
    b = m.aggregate_support(shuffled, 250_000)
    key = lambda r: (r["bin_start"], r["clade"])
    da = {key(r): (r["n_loci"], r["q1"], r["q2"], r["q3"]) for r in a}
    db = {key(r): (r["n_loci"], r["q1"], r["q2"], r["q3"]) for r in b}
    assert da == db


def test_contrasts_are_recomputed_from_binned_support():
    for resolution in (100, 250, 500):
        for r in read(OUT / f"binned_quartet_support_{resolution}kb.tsv"):
            assert abs(float(r["D2"]) - (float(r["q2"]) - float(r["q1"]))) < 1e-12
            assert abs(float(r["D3"]) - (float(r["q3"]) - float(r["q1"]))) < 1e-12


def test_structural_region_coordinates_are_unchanged():
    source = read(ROOT / "data/neoaves_chr4/processed/canonical_structural_regions.tsv")
    comparison = read(OUT / "structural_region_comparison.tsv")
    assert {r["region"] for r in comparison} == {r["canonical_region_id"] for r in source}
