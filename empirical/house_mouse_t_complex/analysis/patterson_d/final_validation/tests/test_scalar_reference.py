#!/usr/bin/env python3
"""Independent tiny scalar reference for diploid population ABBA/BABA."""
import unittest

def scalar_site(groups, outgroup, ref='A', alt='G'):
    """Return one-site ABBA/BABA or None when site is unpolarizable.

    groups is [P1,P2,P3], each a list of diploid GT tuples. Missing alleles
    are None. Outgroup is list of diploid GT tuples.
    """
    called = [g for gt in outgroup if gt is not None and len(gt)==2 and None not in gt for g in gt]
    if len(called) != 2*len(outgroup):
        return None
    if all(a == 0 for a in called):
        flip = False
    elif all(a == 1 for a in called):
        flip = True
    else:
        return None
    ps=[]
    for group in groups:
        dosages=[]
        for gt in group:
            if gt is None or len(gt)!=2 or None in gt:
                continue
            dosages.extend(gt)
        if not dosages:
            return None
        p=sum(dosages)/(2*len([gt for gt in group if gt is not None and len(gt)==2 and None not in gt]))
        ps.append(1-p if flip else p)
    p1,p2,p3=ps
    return (1-p1)*p2*p3, p1*(1-p2)*p3

class ScalarReferenceTests(unittest.TestCase):
    def test_pure_abba(self):
        self.assertEqual(scalar_site([[(0,0)],[(1,1)],[(1,1)]],[(0,0)]),(1.0,0.0))
    def test_pure_baba(self):
        self.assertEqual(scalar_site([[(1,1)],[(0,0)],[(1,1)]],[(0,0)]),(0.0,1.0))
    def test_balanced_patterns_and_pooling(self):
        a=scalar_site([[(0,0)],[(1,1)],[(1,1)]],[(0,0)])
        b=scalar_site([[(1,1)],[(0,0)],[(1,1)]],[(0,0)])
        self.assertEqual((a[0]+b[0],a[1]+b[1]),(1.0,1.0))
        self.assertEqual((a[0]+b[0]-a[1]-b[1])/(sum(a)+sum(b)),0.0)
    def test_ref_alt_recode_invariance(self):
        ref=scalar_site([[(0,1)],[(0,0)],[(1,1)]],[(0,0)])
        recoded=scalar_site([[(0,1)],[(1,1)],[(0,0)]],[(1,1)])
        self.assertEqual(ref,recoded)
    def test_fixed_alt_outgroup_complements_ancestral_state(self):
        x=scalar_site([[(1,1)],[(0,0)],[(0,0)]],[(1,1)])
        self.assertEqual(x,(1.0,0.0))
    def test_polymorphic_outgroup_excluded(self):
        self.assertIsNone(scalar_site([[(0,0)],[(1,1)],[(1,1)]],[(0,1)]))
    def test_missing_and_heterozygote_population_frequency(self):
        x=scalar_site([[(0,1),(None,None)],[(1,1),(0,1)],[(1,1)]],[(0,0)])
        # p1=.5, p2=.75, p3=1 => ABBA=.375, BABA=.125
        self.assertAlmostEqual(x[0],.375)
        self.assertAlmostEqual(x[1],.125)
    def test_empty_population_excluded(self):
        self.assertIsNone(scalar_site([[(None,None)],[(1,1)],[(1,1)]],[(0,0)]))

if __name__=='__main__':
    unittest.main()
