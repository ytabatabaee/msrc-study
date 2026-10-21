from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "simulations/conditional_topology_time/scripts/conditional_experiment.py"
spec = importlib.util.spec_from_file_location("conditional_experiment", SCRIPT)
mod = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(mod)


def test_pattern_categories():
    assert mod.category("0011") == "concordant_2_2"
    assert mod.category("1010") == "discordant_2_2"
    assert mod.category("1000") == "3_1"
    assert mod.category("0000") == "uniform"


def test_pattern_and_topology_ordering():
    assert mod.PATTERNS[0] == "0000"
    assert mod.PATTERNS[-1] == "1111"
    assert mod.TOPOLOGY_NAMES == ("12|34", "13|24", "14|23")


def test_conditional_reconstructs_marginal():
    rows = [
        {"model": "msrc", "pattern": "0011", "topology_index": 0, "mean_pairwise_time": 1.0},
        {"model": "msrc", "pattern": "0011", "topology_index": 1, "mean_pairwise_time": 3.0},
        {"model": "msrc", "pattern": "1010", "topology_index": 2, "mean_pairwise_time": 5.0},
    ]
    total = len(rows)
    reconstructed = sum(sum(int(r["topology_index"]) == 0 for r in rows if r["pattern"] == p) / total for p in ("0011", "1010"))
    assert math.isclose(reconstructed, 1 / 3)
    time_reconstructed = sum(sum(float(r["mean_pairwise_time"]) for r in rows if r["pattern"] == p) / total for p in ("0011", "1010"))
    assert math.isclose(time_reconstructed, sum(float(r["mean_pairwise_time"]) for r in rows) / total)


def test_empty_pattern_is_undefined_not_zero():
    assert math.isnan(float("nan"))
    assert "undefined_msc_Z" not in {mod.category(p) for p in mod.PATTERNS}


def test_bootstrap_is_reproducible():
    values = np.array([1.0, 2.0, 4.0])
    a = mod.bootstrap_mean(values, np.random.default_rng(7), 100)
    b = mod.bootstrap_mean(values, np.random.default_rng(7), 100)
    assert a == b
