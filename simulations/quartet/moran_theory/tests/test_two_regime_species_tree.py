import numpy as np

from ..prescribed_frequency_quartet import (
    PrescribedTrajectory,
    balanced_prescribed_quartet_probabilities,
    unbalanced_prescribed_quartet_probabilities,
)
from ..two_regime_species_tree import (
    RegimeBranch,
    TwoRegimeModel,
    balanced_two_regime_quartet_probabilities,
    regime_operator,
    unbalanced_two_regime_quartet_probabilities,
)


def _b(duration):
    return RegimeBranch(duration, .2, .4, 10)


def test_operator_dimensions_are_4_12_22_38():
    model = TwoRegimeModel(.2, .8, .17)
    assert [regime_operator(model, _b(.1), n).shape[0] for n in (1, 2, 3, 4)] == [4, 12, 22, 38]


def test_balanced_p_equal_matches_prescribed_frequency():
    model = TwoRegimeModel(.5, .5, .17)
    observed = balanced_two_regime_quartet_probabilities(
        model=model, root_branch=_b(.5), left_pair_branch=_b(.3), right_pair_branch=_b(.1),
        terminal_branches=(_b(.2), _b(.4), _b(.1), _b(.3)), initial_regime=np.array([.3, .7]),
    ).quartet_probabilities
    t = PrescribedTrajectory.constant(.5, 1.0)
    expected = balanced_prescribed_quartet_probabilities(
        M_root=10, M_left=10, M_right=10, M_terminal=(10,) * 4, rho=.17,
        terminal_trajectories=(PrescribedTrajectory.constant(.5, .2), PrescribedTrajectory.constant(.5, .4), PrescribedTrajectory.constant(.5, .1), PrescribedTrajectory.constant(.5, .3)),
        sister_trajectories=(PrescribedTrajectory.constant(.5, .3), PrescribedTrajectory.constant(.5, .1)),
        root_trajectory=PrescribedTrajectory.constant(.5, .5),
    ).quartet_probabilities
    assert np.allclose(observed, expected, atol=2e-12)


def test_unbalanced_joint_law_conserves_and_marginalizes():
    result = unbalanced_two_regime_quartet_probabilities(
        model=TwoRegimeModel(.2, .8, .17), root_branch=_b(.5), intermediate_branch=_b(.2),
        sister_branch=_b(.3), terminal_branches=(_b(.2), _b(.4), _b(.1), _b(.3)),
        initial_regime=np.array([.6, .4]), with_joint=True,
        H_root_intermediate=np.array([[.8, .2], [.1, .9]]),
    )
    assert np.isclose(result.joint.sum(), 1.0, atol=2e-11)
    assert np.allclose(result.joint.sum(axis=0), result.quartet_probabilities, atol=2e-11)
    assert np.all(result.quartet_probabilities >= -2e-11)


def test_zero_regime_switching_is_a_weighted_fixed_regime_mixture():
    model = TwoRegimeModel(.2, .8, .17)
    observed = balanced_two_regime_quartet_probabilities(
        model=model, root_branch=RegimeBranch(.5, 0, 0, 10),
        left_pair_branch=RegimeBranch(.3, 0, 0, 10), right_pair_branch=RegimeBranch(.1, 0, 0, 10),
        terminal_branches=tuple(RegimeBranch(x, 0, 0, 10) for x in (.2, .4, .1, .3)),
        initial_regime=np.array([.6, .4]),
    ).quartet_probabilities
    assert np.isclose(observed.sum(), 1.0, atol=2e-11)
