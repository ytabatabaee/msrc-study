import numpy as np

from ..certify_conditional_moran_jacobian import (
    balanced_conditional_derivatives,
    pectinate_conditional_derivatives,
)


def test_analytic_jacobians_are_nonzero_and_normalized():
    for fun in (balanced_conditional_derivatives, pectinate_conditional_derivatives):
        q, jac, joint, pz = fun()
        assert np.isclose(q.sum(), 1.0)
        assert np.isclose(joint.sum(), pz)
        assert pz > 0
        assert abs(np.linalg.det(jac)) > 1e-5


def test_analytic_derivatives_match_exact_solver_finite_difference():
    from ..analyze_conditional_moran_image import balanced_q, pectinate_q

    for analytic, numeric in (
        (balanced_conditional_derivatives, lambda x: balanced_q("0101", x)),
        (pectinate_conditional_derivatives, lambda x: pectinate_q("0101", x, M=10)),
    ):
        _, jac, _, _ = analytic()
        x = np.array([1.0, .2])
        cols = []
        for j in range(2):
            e = np.zeros(2); e[j] = 1e-5
            cols.append((numeric(x + e)[:2] - numeric(x - e)[:2]) / (2e-5))
        assert np.max(np.abs(jac - np.column_stack(cols))) < 1e-7


def test_analytic_values_match_the_existing_conditional_recursions():
    from ..analyze_conditional_moran_image import balanced_q, pectinate_q

    q_bal, _, _, _ = balanced_conditional_derivatives()
    q_pec, _, _, _ = pectinate_conditional_derivatives()
    assert np.max(np.abs(q_bal - balanced_q("0101", (1.0, .2)))) < 1e-12
    assert np.max(np.abs(q_pec - pectinate_q("0101", (1.0, .2), M=10))) < 1e-12


def test_historical_pectinate_m4_value_is_reproduced_explicitly():
    from ..analyze_conditional_moran_image import pectinate_q

    q, jac, _, _ = pectinate_conditional_derivatives(M=4)
    assert np.max(np.abs(q - pectinate_q("0101", (1.0, .2), M=4))) < 1e-12
    assert np.isclose(np.linalg.det(jac), -6.494903570936153e-4, rtol=1e-8, atol=1e-12)
