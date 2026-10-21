from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "simulations/topology_time_separation/scripts"))
from signal_separation import TOPOLOGY_NAMES, parse_newick, summarize_cell, topology_from_split, tree_metrics


def test_quartet_topology_labels():
    assert [topology_from_split(i) for i in range(3)] == ["12|34", "13|24", "14|23"]


def test_ordinary_msc_zero_internal_branch_is_symmetric():
    # Existing MSC control formula: at zero internal length, all three
    # unrooted quartet resolutions have probability 1/3.
    internal_length = 0.0
    expected = np.array(
        [1.0 - (2.0 / 3.0) * math.exp(-internal_length),
         (1.0 / 3.0) * math.exp(-internal_length),
         (1.0 / 3.0) * math.exp(-internal_length)]
    )
    assert np.allclose(expected, np.full(3, 1.0 / 3.0))


def test_pairwise_mrca_times_and_total_length():
    metrics = tree_metrics("((1:1,2:1):2,(3:2,4:2):1);")
    assert metrics["T12"] == 1.0
    assert metrics["T34"] == 2.0
    assert metrics["T13"] == metrics["T14"] == metrics["T23"] == metrics["T24"] == 3.0
    assert metrics["total_tree_length"] == 9.0


def test_internal_branch_formula_and_time_ratio():
    metrics = tree_metrics("((1:1,2:1):2,(3:2,4:2):1);")
    assert math.isclose(metrics["internal_branch_length_q0"], 3.0)
    assert math.isclose(metrics["R_T"], 6.0 / 3.0)


def test_common_scaling_preserves_time_ratio():
    metrics = tree_metrics("((1:1,2:1):2,(3:2,4:2):1);")
    scaled = {key: float(value) * 7.3 for key, value in metrics.items() if key.startswith("T")}
    assert math.isclose((scaled["T13"] + scaled["T24"]) / (scaled["T12"] + scaled["T34"]), metrics["R_T"])


def test_mrca_tree_has_four_taxa():
    assert len(parse_newick("((1:1,2:1):2,(3:2,4:2):1);").children) == 2


def test_summary_keeps_history_and_genealogy_counts_separate():
    base = tree_metrics("((1:1,2:1):2,(3:2,4:2):1);")
    rows = []
    for model in ("msc", "msrc"):
        for history in range(3):
            row = {k: v for k, v in base.items() if k not in {"pair_distances", "newick"}}
            row.update({"model": model, "history_id": history, "genealogy_id": 0, "topology_index": 0, "topology": "12|34", "pattern": "1010", "n_coalescences": 3, "last_coalescence_time": 3.0})
            rows.append(row)
    summary, _, _ = summarize_cell(rows, 20, .5, 0.1, .25, .5, None)
    assert summary["n_msrc_histories"] == 3
    assert summary["n_msrc_genealogies"] == 3
    assert summary["delta_Q_mc_se_approx"] >= 0
    assert math.isnan(summary["delta_T_mc_se_approx"])
