import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / 'empirical/neoaves_chr4/results/stage7_arrangement_divergence'
STAGE6 = ROOT / 'empirical/neoaves_chr4/results/stage6_genomic_tracks'

def read(name):
    with (RES / name).open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def test_state_audit_and_eligibility():
    rows = read('region_arrangement_eligibility.tsv')
    assert len(rows) == 11
    assert {r['canonical_region_id'] for r in rows if r['eligible'] == 'True'} == {'CR_01','CR_10','CR_11'}
    assert all(int(r['n_state_taxa']) == 6 for r in rows)
    assert all(int(r['n_same_pairs']) >= 0 for r in rows)

def test_distance_contrast_and_denominator():
    rows = read('window_arrangement_distance.tsv')
    assert len(rows) > 17000
    for r in rows[::101]:
        a = float(r['A_distance'])
        same = float(r['mean_same_state'])
        diff = float(r['mean_different_state'])
        b = float(r['B_region'])
        assert abs(a - (diff - same)) < 1e-9
        assert abs(float(r['A_relative']) - a / b) < 1e-9
        assert int(r['n_eligible_pairs']) == int(r['n_same_pairs']) + int(r['n_different_state_pairs'])
        assert b > 0

def test_no_unsupported_binary_relabeling():
    text = (RES / 'structural_state_provenance.md').read_text()
    assert 'multistate labels' in text
    assert 'not a universal binary arrangement assignment' in text
    assert 'independent of the local ML trees' in text

def test_stage6_quartet_support_is_reused_without_change():
    source = list(csv.DictReader((STAGE6 / 'window_quartet_support.tsv').open(), delimiter='\t'))
    assert len(source) > 18000
    script = (ROOT / 'empirical/neoaves_chr4/scripts/07_arrangement_divergence.py').read_text()
    assert 'support=read(SUPPORT)' in script
    for r in source[::137]:
        for field in ('q1_species_history', 'q2_alternative_1', 'q3_alternative_2', 'D_q2_minus_q1', 'D_q3_minus_q1'):
            assert math.isfinite(float(r[field]))
        assert abs(float(r['support_sum']) - 1.0) < 1e-9

def test_stage7_manifest_records_provenance():
    manifest = (RES / 'manifest.json').read_text()
    assert 'reproduction_command' in manifest
    assert 'arrangement_A_primary_supported' in manifest
