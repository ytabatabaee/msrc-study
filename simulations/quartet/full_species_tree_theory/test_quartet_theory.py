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
        conditional_scalar,
        scalar_coefficients,
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
        conditional_scalar,
        scalar_coefficients,
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


def test_Hm_symmetry_classes_are_exact():
    m, lam, T_A = 0.17, 0.8, 0.63
    H = hm_kernel(T_A, m, lam)
    alpha = scalar_coefficients(0.0, m, lam, T_A)["alpha"]
    u = np.ones(3) / 3
    for z in PATTERNS:
        row = H[int(z, 2)]
        if z.count("1") in (0, 1, 3, 4):
            expected = u
        else:
            j = {"0011": 0, "1100": 0, "0101": 1, "1010": 1, "0110": 2, "1001": 2}[z]
            expected = u + alpha * (np.eye(3)[j] - u)
        assert np.allclose(row, expected, atol=1e-12)


def test_scalar_matches_matrix_for_all_patterns_and_grid():
    for m in (0.0, 1e-4, 1e-2, 0.1, 1.0, 10.0, 100.0):
        for lam, T, T_A in ((0.2, 0.1, 0.0), (0.8, 0.63, 0.4), (2.0, 1.2, 1.7)):
            for z in PATTERNS:
                matrix = conditional_quartet_probabilities(z, sister_lengths=(T, T), shared_length=T_A, m=m, lam=lam)
                scalar = conditional_scalar(z, T, m, lam, T_A)
                assert np.max(np.abs(matrix - scalar)) < 1e-12


def test_discordant_asymmetry_identity_and_reversal_condition():
    m, lam, T, T_A = 0.1, 0.8, 0.63, 0.5
    q = conditional_scalar("0101", T, m, lam, T_A)
    x = scalar_coefficients(T, m, lam, T_A)
    assert np.isclose(q[1] - conditional_scalar("0110", T, m, lam, T_A)[1], x["alpha"] * x["epsilon"] ** 2)
    assert (q[1] > q[0]) == (x["alpha"] * (x["d"] ** 2 + x["e"] ** 2 - 2 * x["c"] ** 2) > 1 - x["s_D"] ** 2)


def test_species_tree_prior_coalescence_is_concordant():
    q = conditional_quartet_probabilities("0000", sister_lengths=(0.8, 0.8), shared_length=0.0, m=0.0, lam=0.7)
    assert q[0] > q[1] and q[0] > q[2]


def test_msc_limit():
    q = conditional_quartet_probabilities("1010", sister_lengths=(0.4, 0.4), shared_length=0.8, m=1e5, lam=1.0)
    expected = msc_quartet_probabilities((1.0 / 2.0) * (0.4 + 0.4))
    assert np.allclose(q, expected, atol=3e-4)
