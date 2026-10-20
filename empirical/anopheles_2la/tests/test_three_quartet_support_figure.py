from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "empirical/anopheles_2la/results/three_quartet_support"
TOPO = ROOT / "empirical/anopheles_2la/results/stage2r_quartet_topology_signal.tsv"
DIST = ROOT / "empirical/anopheles_2la/results/stage12_distance/stage12_window_distance.tsv"


def read(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_three_supports_and_contrast_are_consistent():
    rows = read(OUT / "window_three_quartet_support.tsv")
    validation = {row["window_id"]: row for row in read(OUT / "D_consistency_validation.tsv")}
    assert len(rows) == 988
    valid = [row for row in rows if row["n_informative_quartets"] != "0"]
    assert len(valid) == 978
    for row in valid:
        q = [float(row[key]) for key in ("q_species", "q_arrangement", "q_third")]
        assert all(0 <= x <= 1 for x in q)
        assert abs(sum(q) - 1) < 2e-8
        assert abs(float(row["D_reconstructed"]) - (float(row["q_arrangement"]) - float(row["q_species"]))) < 1e-12
        assert abs(float(validation[row["window_id"]]["D_difference"])) < 2e-9


def test_distance_tracks_and_coordinates_match_frozen_stage12():
    output = {row["window_id"]: row for row in read(OUT / "window_three_quartet_support.tsv")}
    distance = {row["window_id"]: row for row in read(DIST)}
    topology = {row["window_id"]: row for row in read(TOPO)}
    assert set(output) == set(topology)
    for locus, source in distance.items():
        row = output[locus]
        assert row["start"] == source["start"]
        assert row["end"] == source["end"]
        assert float(row["A_distance"]) == float(source["A_distance"])
        assert float(row["A_relative"]) == float(source["A_relative"])


def test_new_figures_and_frozen_inversion_coordinates():
    for stem in ("anopheles_2la_three_quartet_support", "anopheles_2la_full_topology_distance", "anopheles_2la_supplementary_genealogy_controls"):
        assert (ROOT / f"empirical/anopheles_2la/figures/{stem}.pdf").exists()
        assert (ROOT / f"empirical/anopheles_2la/figures/{stem}.png").exists()
    audit = (OUT / "quartet_definition_audit.md").read_text()
    assert "q_species" in audit and "q_arrangement" in audit and "q_third" in audit
    report = (OUT / "validation_report.md").read_text()
    assert "20524058-42165532" in report
