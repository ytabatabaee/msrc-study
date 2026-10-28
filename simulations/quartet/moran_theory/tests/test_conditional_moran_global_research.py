import numpy as np

from simulations.quartet.moran_theory.analyze_conditional_moran_global_research import (
    radial_rows,
)


def test_radial_family_has_equal_nondominant_support():
    for shape in ("balanced", "pectinate"):
        rows = radial_rows(shape)
        for row in rows:
            assert row["pattern_probability"] > 0
            assert abs(row["normalization_error"]) < 1e-12
            assert abs(row["nonselected_difference"]) < 1e-12


def test_radial_family_connects_center_to_each_vertex_direction():
    for shape in ("balanced", "pectinate"):
        rows = radial_rows(shape)
        for pattern, selected in {"0011": 0, "0101": 1, "0110": 2}.items():
            subset = [r for r in rows if r["pattern"] == pattern]
            q_start = np.array([subset[0]["q1"], subset[0]["q2"], subset[0]["q3"]])
            q_end = np.array([subset[-1]["q1"], subset[-1]["q2"], subset[-1]["q3"]])
            assert np.allclose(q_start, 1 / 3, atol=1e-12)
            assert q_end[selected] > 0.99

