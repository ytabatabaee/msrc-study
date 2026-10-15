from __future__ import annotations
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "empirical/house_mouse_t_complex"
AUDIT = BASE / "results/recombination_suppression_source_audit.tsv"
REPORT = BASE / "results/recombination_suppression_report.md"
FIG_PNG = BASE / "figures/house_mouse_t_complex_suppression.png"
FIG_PDF = BASE / "figures/house_mouse_t_complex_suppression.pdf"
Q_SCAN = BASE / "results/stage2_fixed_quartet_scan.tsv"

def test_audit_contains_direct_quantitative_and_negative_evidence():
    text = AUDIT.read_text()
    assert "silver_artzt_1981" in text and "0.1-0.5 cM" in text
    assert "hammer_bliss_silver_1991" in text and "100-200-fold" in text
    assert "kelemen_vicoso_2018" in text and "No direct crossover frequency" in text

def test_coordinate_and_rate_claims_are_explicitly_gated():
    text = REPORT.read_text()
    assert "rather than coordinate-resolved" in text
    assert "coordinate-resolved recombination profile" in text
    assert "No lift-over is required" in text
    assert "not direct suppression rates" in text

def test_frozen_quartet_scan_preserved():
    assert hashlib.sha256(Q_SCAN.read_bytes()).hexdigest() == "9ad62b579fae816d04b05b234f075e2b6f143198bdb5e05df0df592bddfbcbcd"

def test_suppression_figure_exists_and_is_nonempty():
    assert FIG_PNG.stat().st_size > 20_000
    assert FIG_PDF.stat().st_size > 10_000
