import csv
from pathlib import Path
import math
ROOT=Path(__file__).resolve().parents[6]
V=ROOT/'empirical/anopheles_2la/analysis/patterson_d/final_validation'
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def test_paired_mac_retains_informative_contributions():
 rows=read(V/'results/estimator_validation.tsv')
 assert rows
 for r in rows:
  assert math.isclose(float(r['primary_informative_contributions']),float(r['MAC2_informative_contributions']),abs_tol=1e-6)
  assert math.isclose(float(r['primary_ABBA']),float(r['MAC2_ABBA']),abs_tol=1e-8)
  assert math.isclose(float(r['primary_BABA']),float(r['MAC2_BABA']),abs_tol=1e-8)
def test_primary_inside_effect_signs_are_orientation_dependent():
 rows=read(V/'results/final_results.tsv')
 d={(r['test_id'],r['region']):float(r['D_pooled']) for r in rows}
 assert d['A_merus','inside'] < -.5
 assert d['B_merus','inside'] > .5
 assert abs(d['A_merus','outside']) < .05
 assert abs(d['B_merus','outside']) < .05
def test_sample_level_results_explicitly_not_run():
 p=(V/'results/sample_sensitivity_report.md').read_text()
 assert 'NOT_RUN_RAW_GENOTYPES_UNAVAILABLE_IN_VALIDATION_RUNTIME' in p
 assert 'not mathematically equivalent' in p
