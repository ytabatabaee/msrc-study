from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "analysis/suppression_effect"


def read(name):
    with (OUT / name).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_within_msrc_reference_is_kappa_one():
    rows = read("suppression_effects.tsv")
    assert len(rows) == 12
    for rep in {r["replicate"] for r in rows}:
        ref = next(r for r in rows if r["replicate"] == rep and r["kappa"] == "1.0")
        assert float(ref["delta_astral_rf_vs_kappa1"]) == 0.0
        assert float(ref["delta_gene_tree_rf_vs_kappa1"]) == 0.0


def test_bootstrap_preserves_history_blocks():
    rows = read("astral_bootstrap.tsv")
    assert len(rows) == 80
    assert {r["n_blocks"] for r in rows} == {"25"}
    assert "independent forward histories" in {r["block_type"] for r in rows}
    assert "ordinary MSC pseudo-blocks" in {r["block_type"] for r in rows}


def test_split_and_branch_outputs_are_nonempty():
    assert len(read("split_recovery.tsv")) > 100
    rows = read("branch_length_comparison.tsv")
    assert rows
    assert all(float(r["absolute_error"]) >= 0 for r in rows)
