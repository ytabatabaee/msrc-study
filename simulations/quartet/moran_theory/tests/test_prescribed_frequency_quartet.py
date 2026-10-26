import numpy as np
from scipy.integrate import solve_ivp

from ..full_species_tree import fixed_p_half_species_tree_quartet_probabilities
from ..moran_quartet_operator import G1, G2, G3, G4
from ..prescribed_frequency_quartet import (
    PrescribedTrajectory,
    balanced_prescribed_quartet_probabilities,
    prescribed_generator,
    unbalanced_prescribed_quartet_probabilities,
)
from ..prescribed_frequency_quartet import _propagate_row


def _t(*segments):
    return PrescribedTrajectory(tuple((float(p), float(dt)) for p, dt in segments))


def _balanced_kwargs(with_joint=False):
    return dict(
        M_root=10, M_left=10, M_right=10, M_terminal=(10, 10, 10, 10), rho=0.17,
        terminal_trajectories=(_t((.5, .2)), _t((.5, .4)), _t((.5, .1)), _t((.5, .3))),
        sister_trajectories=(_t((.5, .3)), _t((.5, .1))),
        root_trajectory=_t((.5, .5)), with_joint=with_joint,
    )


def test_prescribed_generators_match_fixed_count_generators():
    for n, old in ((1, G1), (2, G2), (3, G3), (4, G4)):
        observed = prescribed_generator(n, .5, 10, .17)
        expected = old(5, 10, .17).toarray()
        assert np.allclose(observed, expected, atol=1e-14)


def test_balanced_constant_frequency_matches_existing_fixed_frequency_result():
    observed = balanced_prescribed_quartet_probabilities(**_balanced_kwargs()).quartet_probabilities
    expected = fixed_p_half_species_tree_quartet_probabilities(
        M=10, rho=.17, root_length=.5, sister_lengths=(.3, .1),
        terminal_lengths=(.2, .4, .1, .3),
    )
    assert np.allclose(observed, expected, atol=2e-12)


def test_piecewise_propagation_matches_independent_ode_integration():
    trajectory = _t((.2, .13), (.8, .21), (.35, .17))
    row = np.array([.3, .7])
    observed = row.copy()
    for p, duration in trajectory.segments:
        if duration:
            solution = solve_ivp(
                lambda _t, y: prescribed_generator(1, p, 10, .23).T @ y,
                (0.0, duration), observed, rtol=2e-12, atol=2e-13,
            )
            observed = solution.y[:, -1]
    from_module = _propagate_row(row, trajectory, 1, 10, .23)
    assert np.allclose(observed, from_module, atol=2e-11)


def test_balanced_joint_marginal_and_unbalanced_conservation():
    balanced = balanced_prescribed_quartet_probabilities(**_balanced_kwargs(with_joint=True))
    assert np.isclose(balanced.joint.sum(), 1.0, atol=2e-12)
    assert np.allclose(balanced.joint.sum(axis=0), balanced.quartet_probabilities, atol=2e-12)
    assert np.allclose(balanced.joint.sum(axis=1), balanced.pattern_probabilities, atol=2e-12)

    trajectories = (_t((.5, .2)), _t((.5, .4)), _t((.5, .1)), _t((.5, .3)))
    unbalanced = unbalanced_prescribed_quartet_probabilities(
        M_root=10, M_intermediate=10, M_sister=10, M_terminal=(10, 10, 10, 10), rho=.17,
        terminal_trajectories=trajectories, sister_trajectory=_t((.5, .3)),
        intermediate_trajectory=_t((.5, .2)), root_trajectory=_t((.5, .5)), with_joint=True,
    )
    assert np.isclose(unbalanced.joint.sum(), 1.0, atol=2e-12)
    assert np.all(unbalanced.quartet_probabilities >= -2e-12)


def test_fixation_excludes_impossible_terminal_patterns():
    trajectories = (_t((0.0, .1)), _t((1.0, .1)), _t((.5, .1)), _t((.5, .1)))
    result = balanced_prescribed_quartet_probabilities(
        M_root=10, M_left=10, M_right=10, M_terminal=(10, 10, 10, 10), rho=.2,
        terminal_trajectories=trajectories,
        sister_trajectories=(_t((.5, .1)), _t((.5, .1))),
        root_trajectory=_t((.5, .1)), with_joint=True,
    )
    assert np.isclose(result.pattern_probabilities.sum(), 1.0, atol=2e-12)
    assert np.count_nonzero(result.pattern_probabilities) == 4
