import numpy as np

from ..certify_conditional_moran_jacobian import balanced_conditional_derivatives, pectinate_conditional_derivatives


def test_rigorous_audit_inputs_are_explicit_and_normalized():
    for function in (balanced_conditional_derivatives, pectinate_conditional_derivatives):
        q, jac, joint, probability = function(M=10)
        assert np.isclose(q.sum(), 1.0)
        assert np.isclose(joint.sum(), probability)
        assert probability > 0
        assert np.isfinite(jac).all()


def test_pectinate_m4_and_m10_are_distinct_parameterizations():
    q4, j4, _, _ = pectinate_conditional_derivatives(M=4)
    q10, j10, _, _ = pectinate_conditional_derivatives(M=10)
    assert not np.isclose(np.linalg.det(j4), np.linalg.det(j10))
    assert not np.allclose(q4, q10)
