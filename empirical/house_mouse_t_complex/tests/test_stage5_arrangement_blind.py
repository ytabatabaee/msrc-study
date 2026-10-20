from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "empirical/house_mouse_t_complex/results/stage5_unified_tracks"
DIST = ROOT / "empirical/house_mouse_t_complex/results/stage4_distance/stage4_window_distance.tsv"


def read(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_arrangement_blind_track_matches_frozen_stage4():
    unified = {row["locus_id"]: row for row in read(OUT / "window_unified_tracks.tsv")}
    source = {row["locus_id"]: row for row in read(DIST)}
    assert set(unified) == set(source)
    for locus, row in unified.items():
        assert float(row["mean_pair_distance_no_spretus"]) == float(source[locus]["mean_pair_distance_no_spretus"])
        assert float(row["A_distance"]) == float(source[locus]["A_distance"])
        assert float(row["A_relative"]) == float(source[locus]["A_relative"])


def test_blind_track_order_and_250kb_coverage():
    unified = read(OUT / "window_unified_tracks.tsv")
    starts = [int(row["start_bp"]) for row in unified]
    assert starts == sorted(starts)
    assert len(unified) == 4046
    binned = read(OUT / "binned_250kb_tracks_arrangement_blind.tsv")
    assert sum(int(row["n_observed_windows"]) for row in binned) == len(unified)
    assert all(0 < float(row["coverage_fraction"]) <= 1 for row in binned)
    binned_100 = read(OUT / "binned_100kb_tracks_arrangement_blind.tsv")
    assert sum(int(row["n_observed_windows"]) for row in binned_100) == len(unified)
    assert all(0 < float(row["coverage_fraction"]) <= 1 for row in binned_100)


def test_blind_outputs_exist_and_use_frozen_units():
    for suffix in ("5kb", "100kb", "250kb"):
        assert (ROOT / f"empirical/house_mouse_t_complex/figures/house_mouse_t_complex_topology_distance_blind_{suffix}.pdf").exists()
        assert (ROOT / f"empirical/house_mouse_t_complex/figures/house_mouse_t_complex_topology_distance_blind_{suffix}.png").exists()
    methods = (OUT / "arrangement_blind_methods.md").read_text()
    assert "substitution units" in methods
    assert "coalescence times" not in methods
