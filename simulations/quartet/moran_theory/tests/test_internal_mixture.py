import numpy as np

from ..verify_internal_mixture import external_three_component_mixture, target_weights


def test_target_weights_are_simplex_weights():
    target = np.array([.2, .5, .3])
    weights = target_weights(target, 4.0)
    assert np.all(weights >= 0)
    assert np.isclose(weights.sum(), 1)


def test_external_mixture_identity_is_exact():
    for target in (np.array([.96, .02, .02]), np.array([.2, .5, .3]), np.array([1 / 3] * 3)):
        assert np.allclose(external_three_component_mixture(target, 8.0), target, atol=2e-12)


def test_two_state_selector_is_not_three_state_selector():
    from ..verify_internal_mixture import two_state_selector_obstruction
    result = two_state_selector_obstruction()
    assert result["root_regime_states"] == 2
    assert result["independently_selectable_vertex_components"] == 2
