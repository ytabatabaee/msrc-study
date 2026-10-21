from __future__ import annotations

import importlib.util
from pathlib import Path

import math

ROOT = Path(__file__).resolve().parents[3]
SCRIPT = ROOT / "empirical/atlantic_cod/scripts/08c_tcmm_baseline_projection.py"
spec = importlib.util.spec_from_file_location("stage8c", SCRIPT)
stage8c = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(stage8c)


def test_baseline_has_exactly_12_labels_and_fixed_topology():
    tree = stage8c.load_baseline()
    assert stage8c.terminals(tree) == stage8c.TAXA
    assert len(stage8c.topology_keys(tree)) == 11


def test_descendant_sets_are_node_identities():
    tree = stage8c.parse_tree("((a:1,b:1):1,(c:1,d:1):1);")
    keys = stage8c.topology_keys(tree)
    assert ("a", "b") in keys
    assert ("c", "d") in keys
    assert ("a", "b", "c", "d") in keys


def test_ultrametric_and_metric_audits_are_distinct():
    ultrametric = stage8c.parse_tree("((a:1,b:1):2,(c:1,d:1):2);")
    non_ultrametric = stage8c.parse_tree("((a:1,b:1):2,(c:1,d:2):2);")
    assert stage8c.audit_tree(ultrametric)["root_to_tip_deviation"] == 0
    assert stage8c.audit_tree(non_ultrametric)["root_to_tip_deviation"] > 0


def test_ultrametric_mrca_time_is_half_patristic_distance():
    tree = stage8c.parse_tree("((a:1,b:1):2,(c:1,d:1):2);")
    mrca = stage8c.pairwise_values(tree)
    patristic = stage8c.patristic_values(tree)
    for pair in mrca:
        assert math.isclose(2 * mrca[pair], patristic[pair])


def test_stage8b_outputs_are_preserved():
    assert stage8c.STAGE8B_NODE.exists()
    assert stage8c.STAGE8B_FIT.exists()


def test_tcmm_outputs_have_expected_lambda_sensitivity_files():
    run_dir = ROOT / "empirical/atlantic_cod/results/stage8c_tcmm/raw/tcmm_runs"
    for treatment in ("outside", "all"):
        for lam in ("0p0", "0p001"):
            assert (run_dir / f"{treatment}_lambda_{lam}.trees").exists()
