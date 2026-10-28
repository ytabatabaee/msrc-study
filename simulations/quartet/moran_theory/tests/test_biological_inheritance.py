import numpy as np

from ..analyze_biological_inheritance import founding_classification_kernel, identity_kernel, symmetric_kernel


def test_inheritance_kernels_are_stochastic_and_bounded():
    for H in (identity_kernel(), symmetric_kernel(.1), founding_classification_kernel(.2, .8, 10)):
        assert H.shape == (2, 2)
        assert np.all(H >= 0)
        assert np.allclose(H.sum(axis=1), 1)


def test_foundation_kernel_is_derived_before_any_anomaly_search():
    H = founding_classification_kernel(.2, .8, 10)
    assert H[0, 0] > H[1, 0]
    assert H[0, 1] < H[1, 1]


def test_symmetric_inheritance_parameterization():
    H = symmetric_kernel(.02)
    assert np.allclose(H, H.T)
    assert np.allclose(H.sum(axis=1), 1)
