#!/usr/bin/env python3
"""Regression checks for the published/approximate inversion overlay."""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
INTERVALS = ROOT / "empirical/house_mouse_t_complex/results/t_complex_inversion_intervals.tsv"
PLOT = ROOT / "empirical/house_mouse_t_complex/scripts/02m_make_visualizations.py"
Q_SCAN = ROOT / "empirical/house_mouse_t_complex/results/stage2_fixed_quartet_scan.tsv"
RECOMB = ROOT / "empirical/house_mouse_t_complex/results/stage2_recombination_state_500kb.tsv"


def rows():
    lines = INTERVALS.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line.strip()]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_interval_order_and_coordinate_gate():
    intervals = rows()
    assert [row["inversion_name"] for row in intervals] == ["In(17)1", "In(17)2", "In(17)3", "In(17)4"]
    assert all(int(row["start_bp"]) < int(row["end_bp"]) for row in intervals)
    assert all(row["coordinate_status"] == "approximate_schematic" for row in intervals)
    assert all(row["comparison_eligible"] == "false" for row in intervals)


def test_plot_uses_shared_axis_and_preserves_tracks():
    source = PLOT.read_text()
    assert "def plot_inversion_panel" in source
    assert "sharex=axq" in source
    assert "plot_inversion_panel(axinv" in source
    assert 'axq.tick_params(axis="x", labelbottom=False)' in source
    assert 'axd.tick_params(axis="x", labelbottom=False)' in source
    assert "plot_q_axis(axq, qrows, \"ALL_TIPS\"" in source
    assert "plot_recomb_axes(axes, recomb)" in source
    assert digest(Q_SCAN) == "9ad62b579fae816d04b05b234f075e2b6f143198bdb5e05df0df592bddfbcbcd"
    assert digest(RECOMB) == "8410d60b20aec421c65b55c1150df4c36e32f395be84cf7291fa70c4e471c6d4"


def test_plot_module_interval_loader():
    spec = importlib.util.spec_from_file_location("house_mouse_visualizations", PLOT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    loaded = module.load_inversion_intervals()
    assert len(loaded) == 4
    assert int(loaded[0]["start_bp"]) == 5_000_000
    assert int(loaded[-1]["end_bp"]) == 40_000_000
