import numpy as np
from scipy.linalg import expm

try:
    from .quartet_theory import (
        PATTERNS,
        U,
        all_pattern_probabilities,
        conditional_quartet_probabilities,
        hm_kernel,
        msc_quartet_probabilities,
        pair_generator,
        pair_kernel_closed_form,
    )
except ImportError:
    from quartet_theory import (
        PATTERNS,
        U,
        all_pattern_probabilities,
        conditional_quartet_probabilities,
        hm_kernel,
        msc_quartet_probabilities,
        pair_generator,
        pair_kernel_closed_form,
    )


def test_one_lineage_transition_matrix():
    observed = U(0.7, 0.13)
    expected = 0.5 * np.array([[1 + np.exp(-0.182), 1 - np.exp(-0.182)], [1 - np.exp(-0.182), 1 + np.exp(-0.182)]])
    assert np.allclose(observed, expected)
    assert np.allclose(observed.sum(axis=1), 1.0)


def test_pair_closed_form_matches_expm():
    for m, lam, t in ((0.0, 0.7, 0.3), (0.17, 0.8, 0.63), (2.0, 0.2, 1.4)):
        assert np.allclose(pair_kernel_closed_form(t, m, lam), expm(pair_generator(m, lam) * t), atol=1e-12)


def test_pair_survival_and_coalescence_conservation():
    E = pair_kernel_closed_form(0.8, 0.2, 0.7)
    assert np.all(E >= -1e-14)
    assert np.all(E.sum(axis=1) <= 1.0 + 1e-12)
    assert np.all(E.sum(axis=1) > 0.0)


def test_all_sixteen_patterns_conserve_probability():
    qs = all_pattern_probabilities(terminal_lengths=(0.2,) * 4, sister_lengths=(0.3, 0.3), shared_length=0.5, m=0.17, lam=0.8)
    assert set(qs) == set(PATTERNS)
    assert all(np.isclose(q.sum(), 1.0) for q in qs.values())


def test_old_Hm_formula_is_recovered():
    m, lam, T = 0.17, 0.8, 0.63
    A = lam / (lam + 2 * m) * (1 - np.exp(-2 * (lam + 2 * m) * T))
    expected = np.array([(1 - A) / 3, (1 + 2 * A) / 3, (1 - A) / 3])
    observed = hm_kernel(T, m, lam)[int("1010", 2)]
    assert np.allclose(observed, expected, atol=1e-12)
    full = conditional_quartet_probabilities("1010", sister_lengths=(0.0, 0.0), shared_length=T, m=m, lam=lam)
    assert np.allclose(full, expected, atol=1e-12)


def test_species_tree_prior_coalescence_is_concordant():
    q = conditional_quartet_probabilities("0000", sister_lengths=(0.8, 0.8), shared_length=0.0, m=0.0, lam=0.7)
    assert q[0] > q[1] and q[0] > q[2]


def test_msc_limit():
    q = conditional_quartet_probabilities("1010", sister_lengths=(0.4, 0.4), shared_length=0.8, m=1e5, lam=1.0)
    expected = msc_quartet_probabilities((1.0 / 2.0) * (0.4 + 0.4))
    assert np.allclose(q, expected, atol=3e-4)
