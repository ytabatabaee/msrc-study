from pathlib import Path
import csv

ROOT=Path(__file__).resolve().parents[3]
def test_mouse_distance_track_and_tree_audit():
    p=ROOT/'empirical/house_mouse_t_complex/results/stage4_distance/stage4_window_distance.tsv'
    rows=list(csv.DictReader(p.open(),delimiter='\t'))
    assert len(rows)==4046
    assert all(r['n_tips']=='55' for r in rows)
    assert 'A_relative' in rows[0]
def test_mouse_intervals_are_frozen_source():
    rows=list(csv.DictReader((ROOT/'empirical/house_mouse_t_complex/results/t_complex_inversion_intervals.tsv').open(),delimiter='\t'))
    assert [r['inversion_name'] for r in rows]==['In(17)1','In(17)2','In(17)3','In(17)4']
