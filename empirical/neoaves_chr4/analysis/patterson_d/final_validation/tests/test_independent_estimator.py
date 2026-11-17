"""Independent scalar ABBA/BABA and package-integration tests."""
import csv, importlib.util, math, sys, unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[6]
BASE=ROOT/'empirical/neoaves_chr4/analysis/patterson_d'
SPEC=importlib.util.spec_from_file_location('patterson_d',BASE/'scripts/patterson_d.py')
PD=importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(PD)

def scalar_reference(seqs, transversions_only=False):
    """Column-wise reference independent of the production branching logic."""
    order=('P1','P2','P3','O'); n=min(len(seqs[x]) for x in order)
    abba=baba=excluded=informative=0
    transitions={frozenset(('A','G')),frozenset(('C','T'))}
    for i in range(n):
        col=tuple(seqs[x][i].upper() for x in order)
        if not set(col)<=set('ACGT'):
            excluded+=1; continue
        states=set(col)
        if len(states)!=2: continue
        anc=col[3]
        derived=next(iter(states-{anc}),None)
        if derived is None: continue
        if transversions_only and frozenset((anc,derived)) in transitions: continue
        if col==(anc,derived,derived,anc): abba+=1; informative+=1
        elif col==(derived,anc,derived,anc): baba+=1; informative+=1
    return abba,baba,excluded,informative

def seqs_for(columns):
    names=('P1','P2','P3','O')
    return {n:''.join(c[i] for c in columns) for i,n in enumerate(names)}

class ScalarEstimatorTests(unittest.TestCase):
    def compare(self, cols, tx=False):
        seqs=seqs_for(cols)
        self.assertEqual(scalar_reference(seqs,tx),PD.count_sites(seqs,tx))
    def test_pure_abba(self): self.compare([('A','G','G','A')]); self.assertEqual(scalar_reference(seqs_for([('A','G','G','A')])),(1,0,0,1))
    def test_pure_baba(self): self.compare([('G','A','G','A')]); self.assertEqual(scalar_reference(seqs_for([('G','A','G','A')])),(0,1,0,1))
    def test_balanced(self): self.compare([('A','C','C','A'),('G','A','G','A')]); self.assertEqual(PD.patterson_d(1,1),0)
    def test_missing_outgroup_and_ambiguous_ingroups(self): self.compare([('A','C','C','N'),('N','A','G','A')]); self.assertEqual(scalar_reference(seqs_for([('A','C','C','N'),('N','A','G','A')])),(0,0,2,0))
    def test_multistate_is_uninformative(self): self.compare([('A','C','G','A')]); self.assertEqual(scalar_reference(seqs_for([('A','C','G','A')])),(0,0,0,0))
    def test_gap_is_excluded(self): self.compare([('A','C','C','A'),('A','-','C','A')]); self.assertEqual(scalar_reference(seqs_for([('A','C','C','A'),('A','-','C','A')])),(1,0,1,1))
    def test_zero_information_is_undefined(self): self.compare([('A','A','A','A')]); self.assertIsNone(PD.patterson_d(0,0))
    def test_allele_relabeling_preserves_pattern(self):
        before=scalar_reference(seqs_for([('A','C','C','A'),('G','A','G','A')]))
        after=scalar_reference(seqs_for([('T','G','G','T'),('C','T','C','T')]))
        self.assertEqual(before[:2],after[:2])
    def test_transversions_only_uses_ancestral_derived_pair(self):
        seqs=seqs_for([('A','G','G','A'),('A','C','C','A'),('C','T','T','C')])
        self.assertEqual(scalar_reference(seqs,True),(1,0,0,1))
        self.assertEqual(PD.count_sites(seqs,True),scalar_reference(seqs,True))
    def test_coordinate_boundary_inclusive(self):
        intervals=((25030000,32670000),(33510000,34470000),(44130000,56810000))
        self.assertEqual(PD.region_for(25030000),'PNAS_1')
        self.assertEqual(PD.region_for(32670000),'PNAS_1')
        self.assertEqual(PD.region_for(32670001),'outside_all_PNAS')
        self.assertTrue(all(any(lo<=x<=hi for lo,hi in intervals) for x in (33510000,34470000,44130000,56810000)))

class IntegrationTests(unittest.TestCase):
    def tsv(self,p):
        with open(p) as f:return list(csv.DictReader(f,delimiter='\t'))
    def test_pooling_identities_and_exact_quartet_outputs(self):
        checks=self.tsv(BASE/'final_validation/results/estimator_validation.tsv')
        self.assertTrue(checks); self.assertTrue(all(r['status']=='PASS' for r in checks))
        quart=self.tsv(BASE/'final_validation/results/quartet_topology_validation.tsv')
        self.assertEqual(len(quart),6); exact=[r for r in quart if r['source_tree']=='63K.tre']; self.assertEqual(len(exact),3); self.assertTrue(all(int(r['n_exact_quartet_trees'])>0 for r in exact))
    def test_transition_transversion_partition(self):
        rows=self.tsv(BASE/'final_validation/results/final_results.tsv')
        for d in ('Columbea_primary','Columbea_alt_outgroup','N61_primary','N61_alt_outgroup','N62_primary','N62_alt_outgroup'):
            x={r['filter']:r for r in rows if r['design']==d}
            for cls in ('inside','outside'):
                self.assertEqual(int(x['all'][f'ABBA_{cls}']),int(x['transitions_only'][f'ABBA_{cls}'])+int(x['transversions_only'][f'ABBA_{cls}']))
                self.assertEqual(int(x['all'][f'BABA_{cls}']),int(x['transitions_only'][f'BABA_{cls}'])+int(x['transversions_only'][f'BABA_{cls}']))
    def test_primary_results_keep_distinct_nulls(self):
        rows=self.tsv(BASE/'final_validation/results/final_results.tsv')
        for d in ('Columbea_primary','N61_primary','N62_primary'):
            r=next(x for x in rows if x['design']==d and x['filter']=='all')
            self.assertTrue(r['outside_CI95_low'] and r['delta_CI95_low'])
            a,b=int(r['ABBA_inside']),int(r['BABA_inside'])
            self.assertAlmostEqual(float(r['inside_D']),(a-b)/(a+b))

if __name__=='__main__': unittest.main()
