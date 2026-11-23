import math
import numpy as np
from scipy.special import gammaln, xlogy, xlog1py
from ..full_species_tree import inheritance_kernel

def test_small_binomial_formula():
 for parent,daughter in [(4,6),(20,12),(40,40)]:
  h=inheritance_kernel(parent,daughter)
  expected=np.array([[math.comb(daughter,a)*(k/parent)**a*(1-k/parent)**(daughter-a) for a in range(daughter+1)] for k in range(parent+1)])
  np.testing.assert_allclose(h,expected,atol=2e-14,rtol=2e-13)

def test_large_normalization_endpoints_and_log_formula():
 m=2000;h=inheritance_kernel(m,m)
 assert np.isfinite(h).all() and (h>=0).all()
 np.testing.assert_allclose(h.sum(axis=1),1,atol=3e-13)
 assert h[0,0]==h[-1,-1]==1
 a=np.arange(m+1)
 for k in (1,100,1000,1999):
  p=k/m
  expected=np.exp(gammaln(m+1)-gammaln(a+1)-gammaln(m-a+1)+xlogy(a,p)+xlog1py(m-a,-p))
  np.testing.assert_allclose(h[k],expected,atol=2e-12,rtol=2e-11)
