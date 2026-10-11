import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from rf_metrics import TreeValidationError, parse_newick, rf_distance


def test_identical_trees_have_zero_rf():
    result = rf_distance(parse_newick("((A,B),(C,D));"), parse_newick("((A,B),(C,D));"))
    assert result.distance == 0 and result.normalized == 0


def test_rooting_and_branch_lengths_are_ignored():
    a = parse_newick("((A:1,B:2):3,(C:4,D:5):6):7;")
    b = parse_newick("(A,(B,(C,D)));" )
    result = rf_distance(a, b)
    assert result.distance == 0 and result.normalized == 0


def test_one_quartet_nni_is_maximal():
    result = rf_distance(parse_newick("((A,B),(C,D));"), parse_newick("((A,C),(B,D));"))
    assert result.distance == 2 and result.normalized == 1


def test_unresolved_tree_uses_actual_resolution_denominator():
    result = rf_distance(parse_newick("((A,B,C),D);"), parse_newick("((A,B),(C,D));"))
    assert result.distance == 1
    assert result.denominator == 1
    assert result.normalized == 1


def test_duplicate_taxa_are_rejected():
    with pytest.raises(TreeValidationError, match="duplicate"):
        parse_newick("((A,A),(C,D));")


def test_mismatched_taxa_are_rejected():
    with pytest.raises(TreeValidationError, match="taxon sets differ"):
        rf_distance(parse_newick("((A,B),(C,D));"), parse_newick("((A,B),(C,E));"))


def test_dendropy_independent_validation():
    dendropy = pytest.importorskip("dendropy")
    from dendropy.calculate import treecompare
    namespace = dendropy.TaxonNamespace(["A", "B", "C", "D"])
    a = dendropy.Tree.get(data="((A,B),(C,D));", schema="newick", taxon_namespace=namespace)
    b = dendropy.Tree.get(data="((A,C),(B,D));", schema="newick", taxon_namespace=namespace)
    a.is_rooted = False
    b.is_rooted = False
    a.encode_bipartitions()
    b.encode_bipartitions()
    expected = treecompare.symmetric_difference(a, b)
    observed = rf_distance(parse_newick("((A,B),(C,D));"), parse_newick("((A,C),(B,D));")).distance
    assert observed == expected
