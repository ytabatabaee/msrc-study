import numpy as np

from ..analyze_conditional_moran_simplex import conditional_case


def test_exact_moran_paired_patterns_approach_distinct_vertices():
    for z, favored in (("0011", 0), ("0101", 1), ("0110", 2)):
        q, pz = conditional_case(z, 20.0, M=4, K=1)
        assert pz > 0
        assert q[favored] > 1 - 1e-6
        assert np.isclose(q.sum(), 1.0, atol=1e-12)


def test_conditional_pattern_is_not_an_unconditional_statement():
    q, pz = conditional_case("0101", 1.0, M=10, K=5)
    assert pz < 1
    assert q[1] > q[0]
