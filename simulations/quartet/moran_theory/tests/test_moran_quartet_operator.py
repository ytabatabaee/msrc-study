import numpy as np
from scipy.linalg import expm
import sys

sys.path.insert(0, "/Users/ytabatabaee/Desktop/msrc-sim")

from ..moran_quartet_operator import (
    BACKGROUND_STATES,
    G4,
    PATTERNS,
    fixed_frequency_result,
    integrated_quartet_probabilities,
    moran_generator,
    path_joint_probabilities,
    pattern_probabilities_from_moran,
    residual_vectors,
)


def _pi(M, k):
    out = np.zeros(M + 1)
    out[k] = 1.0
    return out


def test_generators_have_zero_rows_and_expected_moran_rates():
    M = 20
    Q = moran_generator(M).toarray()
    assert np.allclose(Q.sum(axis=1), 0.0)
    assert np.isclose(Q[5, 6], 5 * 15 / (2 * 20))
    for k in (0, M):
        assert np.allclose(G4(k, M, 0.1).toarray().sum(axis=1), 0.0)


def test_probability_conservation_and_pattern_identity():
    M, T, rho = 12, 4.0, 0.17
    result = integrated_quartet_probabilities(M, rho, T, _pi(M, 4))
    assert np.isclose(result.quartet_probabilities.sum(), 1.0)
    assert np.isclose(result.pattern_probabilities.sum(), 1.0)
    assert np.allclose(result.joint.sum(axis=1), result.pattern_probabilities)
    assert np.allclose(result.pattern_probabilities, pattern_probabilities_from_moran(M, T, _pi(M, 4)))


def test_t_zero_has_uniform_quartets_and_terminal_pattern_law():
    M = 20
    result = integrated_quartet_probabilities(M, 0.5, 0.0, _pi(M, 5))
    assert np.allclose(result.quartet_probabilities, 1 / 3)
    assert np.allclose(result.conditional_quartet_probabilities, 1 / 3)
    assert np.allclose(result.pattern_probabilities, pattern_probabilities_from_moran(M, 0.0, _pi(M, 5)))


def test_rho_zero_keeps_background_deterministic_along_a_fixed_path():
    M = 10
    result = integrated_quartet_probabilities(M, 0.0, 3.0, _pi(M, 5))
    # With p=1/2 and rho=0, complementary two-two patterns resolve in their
    # associated topology; one-count classes remain uniform after fallback.
    assert result.conditional_quartet_probabilities[PATTERNS.index("1010"), 1] > 1 / 3


def test_fixed_frequency_recovers_compute_Hm_and_A_m():
    M, rho, T = 20, 0.34, 0.63
    result = fixed_frequency_result(M, rho, T)
    m = rho / 2.0
    lam = 2.0 / M
    A = lam / (lam + 2 * m) * (1 - np.exp(-2 * (lam + 2 * m) * T))
    expected = np.array([(1 - A) / 3, (1 + 2 * A) / 3, (1 - A) / 3])
    assert np.allclose(result.conditional_quartet_probabilities[PATTERNS.index("1010")], expected, atol=1e-12)

    from msrcsim.analytic import compute_Hm

    H, states = compute_Hm(T, m, m, lam, lam)
    assert np.allclose(result.conditional_quartet_probabilities[PATTERNS.index("1010")], H[states.index((1, 0, 1, 0))], atol=1e-12)


def test_operator_ordering_matches_direct_backward_product():
    M, rho = 12, 0.21
    ks = (3, 8)
    durations = (0.17, 0.29)
    observed = path_joint_probabilities(ks, durations, M, rho)
    U_back = expm(G4(ks[1], M, rho).toarray() * durations[1]) @ expm(G4(ks[0], M, rho).toarray() * durations[0])
    direct = np.zeros((16, 3))
    final_k = ks[-1]
    for z_index, z in enumerate(BACKGROUND_STATES):
        b = np.zeros(19)
        b[z_index] = (final_k / M) ** sum(z) * (1 - final_k / M) ** (4 - sum(z))
        for r, h in enumerate(residual_vectors()):
            direct[z_index, r] = h @ U_back @ b
    assert np.allclose(observed, direct, atol=1e-12)
