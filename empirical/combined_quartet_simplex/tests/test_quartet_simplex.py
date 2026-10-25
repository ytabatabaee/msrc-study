import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from empirical.combined_quartet_simplex.build_quartet_simplex import interval_state, validate_qss, xy, SQRT3_2, TOPOLOGY_MAPPINGS, coordinate_groups, marker_area, MARKER_MAX_AREA, halo_alpha, halo_area, HALO_RADIUS_FACTOR


def test_simplex_vertices_and_barycenter():
    x, y = xy(__import__('numpy').array([1., 0., 0., 1/3]), __import__('numpy').array([0., 1., 0., 1/3]), __import__('numpy').array([0., 0., 1., 1/3]))
    assert (x[0], y[0]) == (0.5, SQRT3_2)
    assert (x[1], y[1]) == (0.0, 0.0)
    assert (x[2], y[2]) == (1.0, 0.0)
    assert math.isclose(x[3], .5) and math.isclose(y[3], SQRT3_2 / 3)


def test_interval_classification_and_boundary():
    intervals = [("R1", 100, 200)]
    assert interval_state(120, 150, intervals) == ("inside", "R1")
    assert interval_state(1, 50, intervals) == ("outside", "")
    assert interval_state(1, 120, intervals) == ("boundary", "R1")


def test_topology_order_is_explicit_and_q_values_are_not_reordered():
    assert TOPOLOGY_MAPPINGS["atlantic_cod"] == ("baseline/species", "arrangement", "third")
    row = {"dataset": "x", "chromosome": "c", "window_id": "w", "q1": .2, "q2": .7, "q3": .1, "n_informative_quartets": 10}
    kept = validate_qss([row], "x")[0]
    assert (kept["q1"], kept["q2"], kept["q3"]) == (.2, .7, .1)


def test_qss_validation_rejects_bad_sum_and_drops_zero_information():
    good = {"dataset": "x", "chromosome": "c", "window_id": "w1", "q1": .5, "q2": .25, "q3": .25, "n_informative_quartets": 1}
    zero = {**good, "window_id": "w2", "n_informative_quartets": 0}
    assert len(validate_qss([good, zero], "x")) == 1
    bad = {**good, "q3": .3}
    try:
        validate_qss([bad], "x")
    except ValueError:
        pass
    else:
        raise AssertionError("invalid QSS sum was accepted")


def test_coordinate_groups_preserve_coordinates_and_multiplicity():
    rows = [
        {"region_class": "inside", "q1": 0.0, "q2": 1.0, "q3": 0.0},
        {"region_class": "inside", "q1": 0.0, "q2": 1.0, "q3": 0.0},
        # Rounding only affects grouping; the first original coordinate is retained.
        {"region_class": "outside", "q1": 1 / 3, "q2": 1 / 3, "q3": 1 / 3},
        {"region_class": "outside", "q1": 1 / 3 + 1e-13, "q2": 1 / 3 - 1e-13, "q3": 1 / 3},
    ]
    groups = coordinate_groups(rows)
    assert sorted(g["n_overlap"] for g in groups) == [2, 2]
    inside = next(g for g in groups if g["region_class"] == "inside")
    assert (inside["q1"], inside["q2"], inside["q3"]) == (0.0, 1.0, 0.0)
    assert sum(g["n_overlap"] for g in groups) == len(rows)


def test_log_marker_scaling_distinguishes_counts_and_handles_vertices():
    assert marker_area(1) < marker_area(10) < marker_area(100)
    assert marker_area(1000000) == MARKER_MAX_AREA
    x, y = xy(__import__('numpy').array([0.0]), __import__('numpy').array([1.0]), __import__('numpy').array([0.0]))
    assert (x[0], y[0]) == (0.0, 0.0)


def test_halo_scaling_is_independent_and_multiplicity_aware():
    assert halo_alpha(1) < halo_alpha(10) < halo_alpha(100)
    assert halo_area(1, 1.0) > marker_area(1)
    assert math.isclose(math.sqrt(halo_area(1, 1.0) / marker_area(1)), HALO_RADIUS_FACTOR)
