from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "empirical/house_mouse_t_complex/results/stage5_unified_tracks"
Q = ROOT / "empirical/house_mouse_t_complex/results/stage2_fixed_quartet_scan.tsv"
DIST = ROOT / "empirical/house_mouse_t_complex/results/stage4_distance/stage4_window_distance.tsv"
META = ROOT / "data/house_mouse_t_complex/processed/house_mouse_t_complex_ml_5kb_metadata.tsv"


def rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_all_tips_primary_track_and_normalization():
    out = rows(OUT / "window_unified_tracks.tsv")
    assert len(out) == 4046
    for row in out:
        q = [float(row[k]) for k in ("q_species", "q_t_alt", "q_other")]
        assert all(-1e-12 <= value <= 1 + 1e-12 for value in q)
        assert abs(float(row["D"]) - (float(row["q_t_alt"]) - float(row["q_species"]))) < 1e-9
        assert abs(float(row["D"]) + float(row["D_old_species_minus_alt"])) < 1e-9
        assert abs(sum(q) + float(row["q_unresolved"]) - 1.0) < 2e-8


def test_distance_values_and_coordinates_are_frozen():
    out = {r["locus_id"]: r for r in rows(OUT / "window_unified_tracks.tsv")}
    source = {r["locus_id"]: r for r in rows(DIST)}
    metadata = {r["locus_id"]: r for r in rows(META)}
    assert set(out) == set(source) == set(metadata)
    for locus, row in out.items():
        assert row["start_bp"] == source[locus]["start_bp"] == metadata[locus]["start_bp"]
        assert row["end_bp"] == source[locus]["end_bp"] == metadata[locus]["end_bp"]
        assert float(row["A_distance"]) == float(source[locus]["A_distance"])
        assert float(row["A_relative"]) == float(source[locus]["A_relative"])


def test_primary_source_is_all_tips_and_historical_sign_is_preserved():
    primary = [r for r in rows(Q) if r["treatment"] == "ALL_TIPS"]
    out = rows(OUT / "window_unified_tracks.tsv")
    assert len(primary) == len(out) == 4046
    assert {r["locus_id"] for r in primary} == {r["locus_id"] for r in out}
    assert all(float(r["D"]) == float(r["q_t_alt"]) - float(r["q_species"]) for r in out)


def test_250kb_output_has_observed_counts_and_no_interpolation():
    binned = rows(OUT / "binned_250kb_tracks.tsv")
    assert binned
    assert all(int(row["n_windows"]) > 0 for row in binned)
    assert sum(int(row["n_windows"]) for row in binned) == 4046
