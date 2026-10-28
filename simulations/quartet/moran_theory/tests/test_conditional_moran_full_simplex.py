import numpy as np

from simulations.quartet.moran_theory.analyze_conditional_moran_full_simplex import (
    PATTERN_TO_TOPOLOGY,
    vertex_closure_rows,
)


def test_vertex_certificate_is_normalized_and_patterns_are_observable():
    rows = vertex_closure_rows(M=4, K=2, root_lengths=(0.0, 1.0, 8.0))
    assert len(rows) == 18
    for row in rows:
        q = np.array([row["q1"], row["q2"], row["q3"]])
        assert np.all(q >= -1e-12)
        assert abs(q.sum() - 1.0) < 1e-11
        assert row["pattern_probability"] > 0.0
        assert row["pattern"] in PATTERN_TO_TOPOLOGY


def test_selected_topology_becomes_dominant_in_long_root_limit():
    rows = vertex_closure_rows(M=4, K=2, root_lengths=(1.0, 8.0))
    for shape in ("balanced", "pectinate"):
        for pattern, selected in PATTERN_TO_TOPOLOGY.items():
            short = next(r for r in rows if r["shape"] == shape and r["pattern"] == pattern and r["root_length"] == 1.0)
            long = next(r for r in rows if r["shape"] == shape and r["pattern"] == pattern and r["root_length"] == 8.0)
            assert long["selected_q"] > short["selected_q"]
            assert long["selected_q"] > 0.99


def test_vertex_map_uses_the_three_paired_patterns():
    assert PATTERN_TO_TOPOLOGY == {"0011": 0, "0101": 1, "0110": 2}
