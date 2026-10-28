import numpy as np

from ..analyze_two_regime_reachability import jacobian_rank
from ..two_regime_species_tree import (
    RegimeBranch, TwoRegimeModel, balanced_two_regime_quartet_probabilities,
    unbalanced_two_regime_quartet_probabilities,
)


def test_symmetric_fixed_frequency_has_equal_discordant_topologies():
    b = lambda t: RegimeBranch(t, .2, .4, 10)
    result = balanced_two_regime_quartet_probabilities(
        model=TwoRegimeModel(.5, .5, .17), root_branch=b(.5),
        left_pair_branch=b(.3), right_pair_branch=b(.3),
        terminal_branches=tuple(b(.2) for _ in range(4)), initial_regime=np.array([.5, .5]))
    assert np.isclose(result.quartet_probabilities[1], result.quartet_probabilities[2], atol=2e-12)


def test_nonidentity_inheritance_can_create_unconditional_anomaly():
    branch = lambda t: RegimeBranch(t, 0.005686926929444283, 0.05337861031750239, 6)
    kernels = [
        np.array([[1., 0.], [0., 1.]]),
        np.array([[0.23390338334029215, 0.7660966166597079], [0.7481304407795402, 0.25186955922045984]]),
        np.eye(2), np.eye(2), np.eye(2),
        np.array([[0.887581609411756, 0.11241839058824399], [0.673264840197339, 0.326735159802661]]),
    ]
    result = unbalanced_two_regime_quartet_probabilities(
        model=TwoRegimeModel(.33275589777120923, .9865803248558471, .014644711631432385),
        root_branch=branch(.3681308221355468), intermediate_branch=branch(1.124333871562),
        sister_branch=branch(.7871955468288001), terminal_branches=tuple(branch(x) for x in (.0, .0, .0, .0)),
        initial_regime=np.array([.0770560885125247, .9229439114874753]),
        H_intermediate_pair=kernels[0], H_intermediate_3=kernels[1], H_root_intermediate=kernels[2],
        H_root_4=kernels[3], H_pair_1=kernels[4], H_pair_2=kernels[5])
    q = result.quartet_probabilities
    assert q[2] > q[0]
    assert abs(q[1] - q[2]) > .05
    assert np.isclose(q.sum(), 1.0, atol=2e-11)


def test_jacobian_has_two_local_simplex_directions():
    rank, singular = jacobian_rank()
    assert rank == 2
    assert singular[1] > 1e-7
