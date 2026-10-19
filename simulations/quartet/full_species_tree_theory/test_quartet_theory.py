import numpy as np

try:
    from .quartet_theory import (
        PATTERNS,
        conditional_quartet_probabilities,
        msc_quartet_probabilities,
        shared_population_kernel,
    )
except ImportError:  # direct pytest invocation from this directory
    from quartet_theory import (
        PATTERNS,
        conditional_quartet_probabilities,
        msc_quartet_probabilities,
        shared_population_kernel,
    )


BASE = dict(sister_length=0.0, root_length_before_mrca=0.8, rho=4.0, p=0.5, ne=1.0)


def test_probability_conservation_all_patterns():
    for z in PATTERNS:
        q = conditional_quartet_probabilities(z, **BASE)
        assert np.all(q >= 0)
        assert np.isclose(q.sum(), 1.0)


def test_fast_switching_converges_to_msc():
    # The unrooted internal branch has length 2*T for two equal sister
    # branches, hence tau = 2*T/(2Ne) = T/Ne.
    tau = BASE["root_length_before_mrca"] / BASE["ne"]
    expected = msc_quartet_probabilities(tau)
    q = conditional_quartet_probabilities("0101", **{**BASE, "rho": 1e5})
    assert np.allclose(q, expected, atol=2e-4)


def test_no_structure_is_msc():
    tau = BASE["root_length_before_mrca"] / BASE["ne"]
    expected = msc_quartet_probabilities(tau)
    for z in PATTERNS:
        q = conditional_quartet_probabilities(z, **{**BASE, "rho": 1e5})
        assert np.allclose(q, expected, atol=2e-4)


def test_previous_shared_population_formula():
    m = 0.17
    lam = 0.8
    t = 0.63
    A = lam / (lam + 2 * m) * (1 - np.exp(-2 * (lam + 2 * m) * t))
    expected = np.array([(1 - A) / 3, (1 + 2 * A) / 3, (1 - A) / 3])
    observed = shared_population_kernel("1010", duration=t, m01=m, m10=m, lambda0=lam, lambda1=lam)
    assert np.allclose(observed, expected, atol=1e-12)


def test_species_boundary_can_create_concordant_bias():
    q = conditional_quartet_probabilities(
        "1100", sister_length=0.3, root_length_before_mrca=0.0,
        rho=0.0, p=0.5, ne=1.0,
    )
    assert q[0] > q[1]
    assert q[0] > q[2]
