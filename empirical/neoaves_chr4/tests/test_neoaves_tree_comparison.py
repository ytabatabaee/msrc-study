"""Lightweight regression checks for the Neoaves tree comparison outputs."""
import csv
from pathlib import Path

from empirical.neoaves_chr4.scripts.make_neoaves_tree_comparison_figures import (
    ALL_TREATMENTS, CLUSTER, RESULTS, DISPLAY_TAXA, load_tree, prune_tree, tree_splits,
)


def test_treatment_rf_patterns_from_completed_trees():
    expected = {"48": [0, 6, 0, 6, 6], "363": [0, 2, 0, 2, 2]}
    for scale, treatments in ALL_TREATMENTS.items():
        base = tree_splits(load_tree(treatments[0], "cu"))
        observed = []
        for treatment in treatments:
            splits = tree_splits(load_tree(treatment, "cu"))
            observed.append(len(set(base) ^ set(splits)))
        assert observed == expected[scale]


def test_pruning_preserves_induced_topology_and_taxa():
    tree = load_tree("J48_T0", "cu")
    display = prune_tree(tree, DISPLAY_TAXA)
    assert {leaf.taxon.label for leaf in display.leaf_node_iter()} == set(DISPLAY_TAXA)
    assert len(tree_splits(display)) == len(DISPLAY_TAXA) - 2


def test_display_mapping_marks_terminal_cu_as_unavailable():
    path = RESULTS / "neoaves_display_tree_edge_mapping.tsv"
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    terminals = [r for r in rows if r["source_edge_kind"] == "terminal_cu_unavailable"]
    assert terminals
    assert all(r["display_CU_length"] == "NA" for r in terminals)
    assert any(r["source_edge_kind"] == "exact_internal" for r in rows)
    assert any(r["source_edge_kind"] == "composite_internal_path" for r in rows)


def test_branch_comparisons_do_not_encode_absent_splits_as_zero():
    path = RESULTS / "neoaves_cross_taxon_branch_lengths.tsv"
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows
    assert all(r["shared_internal_split"] == "1" for r in rows)
    assert all(r["split_key"] for r in rows)
