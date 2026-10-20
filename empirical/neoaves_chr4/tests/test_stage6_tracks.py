import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / 'empirical/neoaves_chr4/results/stage6_genomic_tracks'

def read(name):
    with (RES / name).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def test_support_normalization_and_contrasts():
    rows = read('window_quartet_support.tsv')
    assert len(rows) == 18109
    for r in rows:
        q = [float(r[k]) for k in ('q1_species_history','q2_alternative_1','q3_alternative_2')]
        assert abs(sum(q)-float(r['support_sum'])) < 1e-9
        assert all(-1e-9 <= x <= 1+1e-9 for x in q)
        assert abs(float(r['D_q2_minus_q1']) - (q[1]-q[0])) < 1e-9
        assert abs(float(r['D_q3_minus_q1']) - (q[2]-q[0])) < 1e-9

def test_coordinates_and_fixed_distance():
    rows = read('window_distance_tracks.tsv')
    assert len(rows) > 6000
    mids = [float(r['midpoint']) for r in rows]
    assert mids == sorted(mids)
    for r in rows:
        if r['fixed_mapped_set_complete'] == 'True':
            assert int(r['n_fixed_taxa']) == 47
            assert math.isfinite(float(r['B_fixed_mapped']))

def test_no_unsupported_arrangement_a():
    text = (RES / 'data_availability_audit.md').read_text()
    assert 'arrangement-informed A_distance/A_relative are not constructed' in text
