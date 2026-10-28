import numpy as np

from ..verify_simplex_coverage import crossed_pattern_quartet, theoretical_edge_vector


def test_three_paired_patterns_select_the_three_topologies():
    for pattern, favored in (("0011", 0), ("0101", 1), ("0110", 2)):
        q = crossed_pattern_quartet(pattern, 8.0)
        assert q[favored] > 1 - 1e-6
        assert np.isclose(q.sum(), 1.0, atol=1e-12)


def test_edge_formula_matches_exact_operator():
    for pattern, favored in (("0011", 0), ("0101", 1), ("0110", 2)):
        for T in (0.0, 1.0, 4.0):
            assert np.allclose(crossed_pattern_quartet(pattern, T), theoretical_edge_vector(favored, T), atol=2e-12)
