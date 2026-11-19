import unittest
from empirical.scripts.joint_genomic_tracks.fdm import fdm_site,pool_fdm,window_fdm,jackknife_ratio,derived_frequency_from_alt

def dsuite_source_reference(p1,p2,p3,po=0.):
    if (p1==p2==p3==0.) or (p1==p2==p3==1.):return 0.,0.
    a=(1-p1)*p2*p3*(1-po);b=p1*(1-p2)*p3*(1-po)
    if p1<=p2:
        if p2>p3:d=((1-p1)*p2*p2*(1-po))-(p1*(1-p2)*p2*(1-po))
        else:d=((1-p1)*p3*p3*(1-po))-(p1*(1-p3)*p3*(1-po))
    elif p1>p3:d=-(((1-p1)*p2*p1*(1-po))-(p1*(1-p2)*p1*(1-po)))
    else:d=-(((1-p3)*p2*p3*(1-po))-(p3*(1-p2)*p3*(1-po)))
    return a-b,d

class FdMTests(unittest.TestCase):
    def test_piecewise_branches_match_dsuite_source(self):
        for p in ((.1,.7,.2),(.1,.3,.8),(.8,.3,.2),(.6,.3,.9)):
            got=fdm_site(*p,0.);ref=dsuite_source_reference(*p,0.)
            self.assertAlmostEqual(got.numerator,ref[0]);self.assertAlmostEqual(got.denominator,ref[1])
    def test_ref_alt_reversal_and_outgroup_polarization(self):
        # Original: REF ancestral (O=0), ALT frequencies .2,.7,.6.
        orig=[derived_frequency_from_alt(x,0.) for x in (.2,.7,.6)]
        # Recoded: alleles swap, so ALT frequencies complement and O is fixed ALT.
        flipped=[derived_frequency_from_alt(1-x,1.) for x in (.2,.7,.6)]
        for x,y in zip(orig,flipped):self.assertAlmostEqual(x,y)
        a,b=fdm_site(*orig),fdm_site(*flipped)
        self.assertAlmostEqual(a.numerator,b.numerator);self.assertAlmostEqual(a.denominator,b.denominator)
        self.assertIsNone(derived_frequency_from_alt(.4,.5))
    def test_positive_and_negative_signs(self):
        self.assertGreater(fdm_site(.1,.7,.6).numerator,0)
        self.assertLess(fdm_site(.8,.2,.6).numerator,0)
        self.assertLess(fdm_site(.8,.2,.6).numerator/fdm_site(.8,.2,.6).denominator,0)
    def test_missing_monomorphic_and_undefined(self):
        self.assertFalse(fdm_site(None,.2,.3).included);self.assertFalse(fdm_site(0,0,0).included)
        self.assertIsNone(pool_fdm([(0,0,0)]).fdm)
    def test_pool_before_division(self):
        rows=[(.1,.8,.7,0),(.7,.3,.2,0)];r=pool_fdm(rows);c=[fdm_site(*x) for x in rows]
        self.assertAlmostEqual(r.fdm,sum(x.numerator for x in c)/sum(x.denominator for x in c))
        self.assertNotAlmostEqual(r.fdm,sum(x.numerator/x.denominator for x in c)/2)
    def test_window_coordinate_and_block_aggregation(self):
        r,b=window_fdm([(1,.1,.8,.7,0),(100001,.2,.8,.6,0),(600001,.8,.2,.6,0)],100000,500000)
        self.assertEqual([x['start'] for x in r],[1,100001,600001]);self.assertEqual(len(b),2)
        self.assertIsNotNone(jackknife_ratio({1:(1,2),500001:(2,3)})['fdm'])
    def test_bad_frequency_rejected(self):
        with self.assertRaises(ValueError):fdm_site(-.1,.2,.3)
if __name__=='__main__':unittest.main()
