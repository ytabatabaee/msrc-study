import numpy as np

from ..moran_quartet_operator import G3, joint_operator_3
from ..unbalanced_species_tree import unbalanced_species_tree_quartet_probabilities


def _pi(M: int, k: int) -> np.ndarray:
    out = np.zeros(M + 1)
    out[k] = 1.0
    return out


def _params(M: int = 12) -> dict:
    return dict(
        M_root=M,
        M_intermediate=M,
        M_sister=M,
        M_terminal=(M, M, M, M),
        rho=0.17,
        root_length=0.8,
        intermediate_length=0.4,
        sister_length=0.3,
        terminal_lengths=(0.2, 0.2, 0.2, 0.2),
        pi_old=_pi(M, M // 2),
    )


def test_g3_has_eight_transient_and_three_absorbing_states_at_boundaries():
    for k in (0, 5, 12):
        G = G3(k, 12, 0.17).toarray()
        assert G.shape == (11, 11)
        assert np.allclose(G.sum(axis=1), 0.0)
        assert np.allclose(G[8:, :], 0.0)
        assert np.all(G[:8, 8:] >= 0.0)
    assert joint_operator_3(12, 0.17).shape == (13 * 11, 13 * 11)


def test_pectinate_probability_conservation_and_root_decomposition():
    result = unbalanced_species_tree_quartet_probabilities(**_params())
    assert np.isclose(result.quartet_probabilities.sum(), 1.0, atol=2e-12)
    assert np.all(result.quartet_probabilities >= -2e-12)
    assert np.all(result.early_topology_contributions >= -2e-12)
    assert np.all(result.root_topology_contributions >= -2e-12)
    assert np.isclose(
        result.early_topology_contributions.sum() + result.root_topology_contributions.sum(),
        1.0,
        atol=2e-12,
    )
    assert np.isclose(result.diagnostics["root_entry_unresolved_mass"], result.root_topology_contributions.sum(), atol=2e-12)


def test_ordinary_msc_limit_has_the_expected_sister_branch_formula():
    # With all chromosomes in one background, rho=0, and no later branch
    # duration, the only topology bias is the 1--2 sister branch.
    M = 40
    t = 0.7
    result = unbalanced_species_tree_quartet_probabilities(
        M_root=M,
        M_intermediate=M,
        M_sister=M,
        M_terminal=(M, M, M, M),
        rho=0.0,
        root_length=0.0,
        intermediate_length=0.0,
        sister_length=M * t,
        terminal_lengths=(0.0, 0.0, 0.0, 0.0),
        pi_old=_pi(M, M),
    )
    expected = np.array([1.0 - (2.0 / 3.0) * np.exp(-t), np.exp(-t) / 3.0, np.exp(-t) / 3.0])
    assert np.allclose(result.quartet_probabilities, expected, atol=2e-10)


def test_zero_duration_branches_are_literal_and_conservative():
    result = unbalanced_species_tree_quartet_probabilities(
        M_root=10,
        M_intermediate=8,
        M_sister=6,
        M_terminal=(5, 7, 9, 11),
        rho=0.0,
        root_length=0.0,
        intermediate_length=0.0,
        sister_length=0.0,
        terminal_lengths=(0.0, 0.0, 0.0, 0.0),
        pi_old=_pi(10, 4),
    )
    assert np.isclose(result.quartet_probabilities.sum(), 1.0, atol=2e-12)


def test_swapping_taxa_one_and_two_preserves_pectinate_topology_labels():
    base = _params(10)
    base["M_terminal"] = (8, 12, 10, 10)
    base["terminal_lengths"] = (0.1, 0.4, 0.2, 0.3)
    swapped = dict(base)
    swapped["M_terminal"] = (12, 8, 10, 10)
    swapped["terminal_lengths"] = (0.4, 0.1, 0.2, 0.3)
    q = unbalanced_species_tree_quartet_probabilities(**base).quartet_probabilities
    q_swapped = unbalanced_species_tree_quartet_probabilities(**swapped).quartet_probabilities
    assert np.isclose(q[0], q_swapped[0], atol=2e-12)
    assert np.isclose(q[1], q_swapped[2], atol=2e-12)
    assert np.isclose(q[2], q_swapped[1], atol=2e-12)
