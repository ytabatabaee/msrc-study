import numpy as np

from ..analyze_conditional_moran_image import PATTERNS, balanced_q, jacobian


def test_all_paired_patterns_have_positive_probability_and_normalize():
    from ..full_species_tree import full_species_tree_quartet_probabilities
    result = full_species_tree_quartet_probabilities(
        M_root=10, M_left=10, M_right=10, M_terminal=(10,) * 4, rho=.1,
        root_length=1, sister_lengths=(.2, .3), terminal_lengths=(.2, .1, .2, .3), pi_old=np.eye(11)[5])
    for z in PATTERNS:
        assert result.pattern_probabilities[int(z, 2)] > 0
        assert np.isclose(result.conditional_quartet_probabilities[int(z, 2)].sum(), 1)


def test_balanced_z0101_has_stable_rank_two_local_image():
    dets = [np.linalg.det(jacobian(lambda x: balanced_q("0101", x), h)) for h in (1e-3, 2e-4, 1e-4)]
    assert all(d < -1e-4 for d in dets)
