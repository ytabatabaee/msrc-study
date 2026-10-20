from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "empirical/house_mouse_t_complex/results/stage5_unified_tracks"


def rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_all_resolutions_cover_the_same_observed_windows():
    for size in (25, 50, 100, 250):
        path = OUT / (f"binned_{size}kb_tracks.tsv" if size != 250 else "multiresolution_binned_250kb_tracks.tsv")
        records = rows(path)
        assert sum(int(row["n_observed_windows"]) for row in records) == 4046
        assert all(0 < float(row["coverage_fraction"]) <= 1 for row in records)
        assert all(0 < float(row["observed_window_fraction"]) <= 1 for row in records)


def test_common_bin_boundaries_and_pooled_support_normalization():
    starts = {}
    for size in (25, 50, 100, 250):
        path = OUT / (f"binned_{size}kb_tracks.tsv" if size != 250 else "multiresolution_binned_250kb_tracks.tsv")
        records = rows(path)
        starts[size] = {(row["bin_start_bp"], row["bin_end_bp"]) for row in records}
        for row in records:
            support = sum(float(row[key]) for key in ("q_species", "q_t_alt", "q_other", "q_unresolved"))
            assert abs(support - 1.0) < 2e-8
            assert abs(float(row["D"]) - (float(row["q_t_alt"]) - float(row["q_species"]))) < 1e-10
    assert all(starts[size] for size in starts)
    assert starts[25] != starts[50]


def test_multiresolution_outputs_and_no_boundary_detection_side_effect():
    summary = rows(OUT / "multiresolution_summary.tsv")
    assert {int(row["resolution_bp"]) for row in summary} == {25_000, 50_000, 100_000, 250_000}
    assert all(int(row["n_incomplete_bins"]) >= 0 for row in summary)
    manifest = (OUT / "multiresolution_manifest.json").read_text()
    assert '"boundary_detection_rerun": false' in manifest
