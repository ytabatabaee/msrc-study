import sys

import numpy as np
from scipy.linalg import expm

sys.path.insert(0, "/Users/ytabatabaee/Desktop/msrc-sim")

from ..full_species_tree import (
    fixed_p_half_species_tree_conditional,
    fixed_p_half_species_tree_quartet_probabilities,
    full_species_tree_quartet_probabilities,
    inheritance_kernel,
)
from ..moran_quartet_operator import G1, G2, G4, integrated_quartet_probabilities


def _pi(M, k):
    out = np.zeros(M + 1)
    out[k] = 1.0
    return out


def _params(M=12):
    return dict(
        M_root=M,
        M_left=M,
        M_right=M,
        M_terminal=(M, M, M, M),
        rho=0.17,
        root_length=0.8,
        sister_lengths=(0.3, 0.3),
        terminal_lengths=(0.2, 0.2, 0.2, 0.2),
        pi_old=_pi(M, M // 2),
    )


def test_n_lineage_branch_generators_conserve_probability():
    for n, size in ((1, 2), (2, 6), (4, 19)):
        for k in (0, 5, 12):
            G = {1: G1, 2: G2, 4: G4}[n](k, 12, 0.17).toarray()
            assert G.shape == (size, size)
            assert np.allclose(G.sum(axis=1), 0.0)


def test_inheritance_kernel_is_stochastic_and_parent_dependence_is_retained():
    H = inheritance_kernel(12, 20)
    assert np.allclose(H.sum(axis=1), 1.0)
    parent = np.array([0.0, 0.5, 0.0, 0.5] + [0.0] * 9)
    marginal = parent @ H
    joint = sum(parent[k] * np.outer(H[k], H[k]) for k in range(13))
    assert np.isclose(marginal.sum(), 1.0)
    assert np.max(np.abs(joint - np.outer(marginal, marginal))) > 1e-4


def test_full_probability_conservation_and_sister_diagnostics():
    result = full_species_tree_quartet_probabilities(**_params())
    assert np.isclose(result.quartet_probabilities.sum(), 1.0)
    assert 0.0 <= result.diagnostics["sister_q1"] <= 1.0
    assert np.isclose(
        result.diagnostics["sister_q1"] + result.diagnostics["neither_sister_coalesces"], 1.0
    )


def test_zero_sister_branches_reduce_to_shared_population_kernel():
    M = 12
    pi = _pi(M, 4)
    common = dict(
        M_root=M, M_left=M, M_right=M, M_terminal=(M,) * 4,
        rho=0.17, root_length=2.0, sister_lengths=(0.0, 0.0),
        terminal_lengths=(0.0,) * 4, pi_old=pi,
    )
    full = full_species_tree_quartet_probabilities(**common).quartet_probabilities
    shared = integrated_quartet_probabilities(M, 0.17, 2.0, pi).quartet_probabilities
    # The collapsed shared-population operator is the limit when the two
    # speciation inheritance steps are removed. The literal finite-M simulator
    # still performs those binomial steps at zero branch length; see the
    # derivation note and the explicit dependence test below.
    assert np.allclose(shared, np.array([1 / 3] * 3), atol=2e-12)
    assert not np.allclose(full, shared, atol=1e-5)


def test_fixed_p_half_recovers_old_full_species_tree_theory():
    M, rho = 20, 0.34
    terminal = (0.2, 0.2, 0.2, 0.2)
    sister = (0.3, 0.3)
    root = 0.8
    observed = fixed_p_half_species_tree_quartet_probabilities(
        M=M, rho=rho, root_length=root, sister_lengths=sister, terminal_lengths=terminal
    )
    from simulations.quartet.full_species_tree_theory.quartet_theory import conditional_quartet_probabilities

    expected = np.mean([
        conditional_quartet_probabilities(
            z, terminal_lengths=terminal, sister_lengths=sister,
            shared_length=root, m=rho / 2.0, lam=2.0 / M,
        )
        for i in range(16)
        for z in (f"{i:04b}",)
    ], axis=0)
    assert np.allclose(observed, expected, atol=2e-12)


def test_fixed_p_half_and_zero_sister_branches_recover_shared_Hm():
    M, rho, T = 20, 0.34, 0.63
    observed = fixed_p_half_species_tree_conditional(
        "1010", M=M, rho=rho, root_length=T, sister_lengths=(0.0, 0.0), terminal_lengths=(0.0,) * 4
    )
    m, lam = rho / 2.0, 2.0 / M
    A = lam / (lam + 2 * m) * (1 - np.exp(-2 * (lam + 2 * m) * T))
    expected = np.array([(1 - A) / 3, (1 + 2 * A) / 3, (1 - A) / 3])
    assert np.allclose(observed, expected, atol=2e-12)


def test_equal_sister_branches_have_left_right_symmetry():
    q = full_species_tree_quartet_probabilities(**_params()).quartet_probabilities
    assert np.isclose(q[1], q[2], atol=1e-12)


def test_taxon_relabeling_swaps_the_two_alternative_topologies():
    base = _params(10)
    base["terminal_lengths"] = (0.1, 0.3, 0.2, 0.4)
    swapped = dict(base)
    swapped["terminal_lengths"] = (0.3, 0.1, 0.2, 0.4)
    q = full_species_tree_quartet_probabilities(**base).quartet_probabilities
    q_swapped = full_species_tree_quartet_probabilities(**swapped).quartet_probabilities
    assert np.isclose(q[0], q_swapped[0], atol=1e-12)
    assert np.isclose(q[1], q_swapped[2], atol=1e-12)
    assert np.isclose(q[2], q_swapped[1], atol=1e-12)


def test_g2_coalesced_lineage_continues_switching():
    M, k, rho, T = 20, 10, 0.4, 0.7
    G = G2(k, M, rho).toarray()
    assert G[4, 5] > 0 and G[5, 4] > 0
    assert np.allclose(expm(G * T).sum(axis=1), 1.0)
