import csv, importlib.util
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
def load():
 p=ROOT/'empirical/cross_system_boundary_detection/stage13c/scripts/13c_finalize_validation.py';s=importlib.util.spec_from_file_location('s13c',p);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def test_interval_set_iou_handles_disconnected_predictions():
 m=load();assert m.iou([(1,2),(9,10)],[(1,10)])[0]==4/10
def test_stage12b_regression_all_exact():
    rows=list(csv.DictReader((ROOT/'empirical/cross_system_boundary_detection/stage13c/results/stage12b_regression.tsv').open(),delimiter='\t'));assert len(rows)==4 and all(r['exact_within_tolerance']=='true' for r in rows)
def test_lg02_discrepancy_has_both_implementations():
 rows=list(csv.DictReader((ROOT/'empirical/cross_system_boundary_detection/stage13c/results/lg02_detector_discrepancy.tsv').open(),delimiter='\t'));assert {r['implementation'] for r in rows}=={'stage12b_B_robust','stage13b_cross_system'}
def test_nulls_have_explicit_replicates_and_hypothesis():
 rows=list(csv.DictReader((ROOT/'empirical/cross_system_boundary_detection/stage13c/results/change_point_null_results.tsv').open(),delimiter='\t'));assert rows and all(int(r['n_replicates'])>0 and r['null_hypothesis'] for r in rows)
