import numpy as np

from ..matrix_free_quartet import (
    matrix_free_dual_lift,
    matrix_free_operator,
    matrix_free_shared_population_quartet_probabilities,
)
from ..full_species_tree import full_species_tree_quartet_probabilities
from ..moran_quartet_operator import dual_lift, integrated_quartet_probabilities


def test_matrix_free_shared_population_matches_sparse_reference():
    for selection in (0.0, 0.25, -0.2):
        M = 8
        pi = np.zeros(M + 1)
        pi[[2, 5]] = [0.35, 0.65]
        sparse = integrated_quartet_probabilities(
            M, 0.17, 0.6, pi, selection=selection,
        )
        matrix_free = matrix_free_shared_population_quartet_probabilities(
            M, 0.17, 0.6, pi, selection=selection,
        )
        assert np.allclose(matrix_free.joint, sparse.joint, atol=2e-12)
        assert np.allclose(
            matrix_free.quartet_probabilities,
            sparse.quartet_probabilities,
            atol=2e-12,
        )
        assert np.allclose(
            matrix_free.pattern_probabilities,
            sparse.pattern_probabilities,
            atol=2e-12,
        )


def test_matrix_free_operator_keeps_count_dimension():
    op = matrix_free_operator(10, 0.17, 4)
    assert op.shape == (11 * 19, 11 * 19)


def test_matrix_free_dual_lift_matches_all_genealogy_state_sizes():
    rng = np.random.default_rng(14)
    for n, d in ((1, 2), (2, 6), (3, 11), (4, 19)):
        M = 6
        values = rng.normal(size=(M + 1) * d)
        expected = dual_lift(M, 0.11, 0.4, values, n, selection=0.15)
        observed = matrix_free_dual_lift(M, 0.11, 0.4, values, n, selection=0.15)
        assert np.allclose(observed, expected, atol=3e-12)


def test_equal_mean_count_distributions_need_not_give_equal_species_tree_laws():
    M = 10
    fixed = np.eye(M + 1)[5]
    mixture = np.zeros(M + 1)
    mixture[[2, 8]] = 0.5
    kwargs = dict(
        M_root=M, M_left=M, M_right=M, M_terminal=(M, M, M, M),
        rho=0.17, root_length=0.5, sister_lengths=(0.3, 0.1),
        terminal_lengths=(0.2, 0.4, 0.1, 0.3), selection=0.0,
    )
    fixed_q = full_species_tree_quartet_probabilities(pi_old=fixed, **kwargs).quartet_probabilities
    mixture_q = full_species_tree_quartet_probabilities(pi_old=mixture, **kwargs).quartet_probabilities
    assert np.max(np.abs(fixed_q - mixture_q)) > 1e-5
