import numpy as np

from ..moran_quartet_operator import PATTERNS
from ..unbalanced_species_tree import (
    unbalanced_species_tree_conditional_quartet_probabilities,
    unbalanced_species_tree_quartet_probabilities,
)


def _kwargs(k=5):
    pi = np.zeros(11); pi[k] = 1.0
    return dict(
        M_root=10, M_intermediate=10, M_sister=10,
        M_terminal=(10, 10, 10, 10), rho=.17,
        root_length=.5, intermediate_length=.2, sister_length=.3,
        terminal_lengths=(0.0, 0.0, .3, .5), pi_old=pi,
    )


def test_joint_marginals_and_conditional_normalization():
    result = unbalanced_species_tree_conditional_quartet_probabilities(**_kwargs())
    unconditional = unbalanced_species_tree_quartet_probabilities(**_kwargs())
    assert result.joint.shape == (16, 3)
    assert result.pattern_probabilities.shape == (16,)
    assert np.isclose(result.joint.sum(), 1.0, atol=2e-12)
    assert np.allclose(result.joint.sum(axis=0), unconditional.quartet_probabilities, atol=2e-12)
    assert np.allclose(result.joint.sum(axis=1), result.pattern_probabilities, atol=2e-12)
    assert np.allclose(result.conditional_quartet_probabilities.sum(axis=1), 1.0, atol=2e-12)


def test_zero_probability_patterns_have_nan_conditionals():
    result = unbalanced_species_tree_conditional_quartet_probabilities(**_kwargs(k=0))
    zero = result.pattern_probabilities < 1e-14
    assert zero.any()
    assert np.all(np.isnan(result.conditional_quartet_probabilities[zero]))
    assert np.all(np.isfinite(result.conditional_quartet_probabilities[~zero]))


def test_unbalanced_conditional_effect_is_present_for_crossed_pattern():
    result = unbalanced_species_tree_conditional_quartet_probabilities(**_kwargs())
    i = PATTERNS.index("0101")
    assert result.conditional_quartet_probabilities[i, 1] > result.conditional_quartet_probabilities[i, 0]
    assert result.conditional_quartet_probabilities[i, 1] != result.conditional_quartet_probabilities[i, 2]


def test_complementary_patterns_swap_the_discordant_labels():
    result = unbalanced_species_tree_conditional_quartet_probabilities(**_kwargs())
    a = result.conditional_quartet_probabilities[PATTERNS.index("0101")]
    b = result.conditional_quartet_probabilities[PATTERNS.index("0110")]
    assert np.isclose(a[0], b[0], atol=2e-12)
    assert np.isclose(a[1], b[2], atol=2e-12)
    assert np.isclose(a[2], b[1], atol=2e-12)
