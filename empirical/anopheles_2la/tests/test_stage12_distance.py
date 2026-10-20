from pathlib import Path
import csv

ROOT=Path(__file__).resolve().parents[3]
def test_distance_outputs_and_units():
    p=ROOT/'empirical/anopheles_2la/results/stage12_distance/stage12_window_distance.tsv'
    assert p.exists()
    rows=list(csv.DictReader(p.open(),delimiter='\t'))
    assert len(rows)>=900
    assert 'mean_pair_distance' in rows[0]
    assert 'A_distance' in rows[0]
def test_fixed_interval_and_gaps_are_recorded():
    rows=list(csv.DictReader((ROOT/'data/anopheles_2la/processed/stage2r_window_grid.tsv').open(),delimiter='\t'))
    assert any(r['region_class']=='inside' for r in rows)
    assert any(r['region_class']=='outside' for r in rows)
