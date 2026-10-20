"""Focused tests for the fire-ant patristic distance extension."""
from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "empirical/fire_ants_chr16/scripts/09_patristic_distance_track.py"
spec = importlib.util.spec_from_file_location("stage9_distance", SCRIPT)
module = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(module)


def test_group_pair_partition_is_six_pairs() -> None:
    assert len(module.ALL_PAIRS) == 6
    assert len(module.OPPOSITE) == 4
    assert len(module.SAME) == 2
    assert set(module.OPPOSITE).isdisjoint(module.SAME)


def test_patristic_distance_and_balanced_contrast() -> None:
    tree = module.parse_tree("((a:1,b:1):2,(c:1,d:1):2);")
    paths = module.ancestor_depths(tree)
    assert module.distance_for_pair(paths, "a", "b") == 2.0
    assert module.distance_for_pair(paths, "a", "c") == 6.0
    assert module.finite_mean([1.0, 3.0]) == 2.0


def test_zero_denominator_is_missing() -> None:
    assert module.finite_mean([]) is None


def test_gap_plot_helper_is_callable() -> None:
    assert callable(module.plot_series)
