import csv
import importlib.util
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
SCRIPT=ROOT/'empirical/cross_system_boundary_detection/stage13b/scripts/13b_corrected_boundary_detection.py'
spec=importlib.util.spec_from_file_location('stage13b_detector',SCRIPT)
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)

def test_interval_union_merges_overlaps_and_adjacency():
    assert m.interval_set_length([(1,3),(3,5)]) == 5
    assert m.interval_set_length([(1,10),(3,4)]) == 10
    assert m.interval_set_length([(1,2),(5,6)]) == 4

def test_disconnected_iou_is_set_iou_not_hull_iou():
    got=m.union_intersection([(1,2),(9,10)],[(1,10)])
    assert got == 4

def test_empty_prediction_has_zero_metrics():
    x=m.interval_metrics([],[(1,10)])
    assert x['iou_union']==0 and x['f1']==0

def test_frozen_stage12b_regression_is_exact():
    p=ROOT/'empirical/cross_system_boundary_detection/stage13b/results/stage12b_regression.tsv'
    rows=list(csv.DictReader(p.open(),delimiter='\t'))
    assert len(rows)==4
    assert all(r['exact_within_tolerance']=='true' for r in rows)

def test_stage12b_result_is_preserved():
    p=ROOT/'empirical/atlantic_cod/results/boundary_detection_stage12b/model_comparison.tsv'
    assert p.exists()
