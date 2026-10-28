import numpy as np

from ..introgression_quartet_comparison import fit_introgression, introgression_quartet


def test_introgression_formula_normalizes():
    q = introgression_quartet(.4, .3, .2, 2)
    assert np.all(q >= 0)
    assert np.isclose(q.sum(), 1)


def test_nontrivial_ms_r_c_witness_has_an_introgression_match():
    target = np.array([.35417018362959773, .26868939582086815, .37714042054953417])
    params, fitted, error = fit_introgression(target, 2)
    assert 0 < params[0] < 1
    assert params[1] > 0 and params[2] > 0
    assert error < 1e-10
    assert np.allclose(fitted, target, atol=1e-10)
