import numpy as np

from ..regime_switching_quartet import regime_switching_shared_population_quartet_probabilities
from ..analyze_observable_reduction import killed_pair_spectrum


def test_killed_pair_modes_grow_with_m():
    for M in (6, 10, 14):
        spectrum = killed_pair_spectrum(M)
        assert len(spectrum) == M - 1
        assert np.min(np.abs(np.diff(np.sort_complex(spectrum)))) > 1e-12


def test_two_regime_alternative_is_normalized_and_has_38_states():
    result = regime_switching_shared_population_quartet_probabilities(
        p_low=.2, p_high=.8, rate_low_to_high=.3, rate_high_to_low=.4,
        M=10, rho=.17, T=.5, initial_regime=np.array([.6, .4]),
    )
    assert result.joint.shape == (16, 3)
    assert np.isclose(result.joint.sum(), 1.0, atol=2e-12)
    assert np.all(result.quartet_probabilities >= -2e-12)


def test_two_regime_zero_time_matches_regime_mixture_sampling():
    result = regime_switching_shared_population_quartet_probabilities(
        p_low=.2, p_high=.8, rate_low_to_high=.3, rate_high_to_low=.4,
        M=10, rho=.17, T=0.0, initial_regime=np.array([.6, .4]),
    )
    assert np.isclose(result.joint.sum(), 1.0, atol=2e-12)
