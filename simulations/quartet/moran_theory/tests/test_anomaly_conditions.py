import numpy as np

from ..analyze_anomaly_conditions import base_condition, conditional_q
from ..moran_quartet_operator import PATTERNS
from ..unbalanced_species_tree import unbalanced_species_tree_conditional_quartet_probabilities


def test_zero_switching_bound_is_sufficient_for_crossed_pattern():
    condition = base_condition(rho=0.0, intermediate_length=4.0)
    q = conditional_q(condition, "0101")
    assert q[1] > q[0]


def test_same_background_case_has_exchangeable_discordants_and_no_anomaly():
    condition = base_condition(rho=0.0)
    q = conditional_q(condition, "0000")
    assert np.isclose(q[1], q[2], atol=2e-12)
    assert q[0] >= q[1]


def test_exact_base_anomaly_is_positive_and_pattern_has_mass():
    condition = base_condition()
    result = unbalanced_species_tree_conditional_quartet_probabilities(
        M_root=condition["M_root"], M_intermediate=condition["M_intermediate"],
        M_sister=condition["M_sister"], M_terminal=condition["M_terminal"],
        rho=condition["rho"], root_length=condition["root_length"],
        intermediate_length=condition["intermediate_length"],
        sister_length=condition["sister_length"], terminal_lengths=(0., 0., .3, .5),
        pi_old=np.eye(11)[condition["k0"]], selection=condition["selection"],
    )
    i = PATTERNS.index("0101")
    q = result.conditional_quartet_probabilities[i]
    assert result.pattern_probabilities[i] > 0
    assert q[1] > q[0]
    assert q[1] != q[2]


def test_neutral_background_complement_preserves_delta():
    condition = base_condition(p0=0.3, rho=0.17)
    complement = base_condition(p0=0.7, rho=0.17)
    a = conditional_q(condition, "0101")
    b = conditional_q(complement, "1010")
    assert np.allclose(a, b, atol=2e-11)
    assert np.isclose(a[1] - a[2], b[1] - b[2], atol=2e-11)


def test_selected_background_complement_uses_fitness_transform():
    condition = base_condition(p0=0.3, selection=0.25)
    complement = base_condition(p0=0.7, selection=-0.25 / 1.25)
    a = conditional_q(condition, "0101")
    b = conditional_q(complement, "1010")
    assert np.allclose(a, b, atol=2e-10)
