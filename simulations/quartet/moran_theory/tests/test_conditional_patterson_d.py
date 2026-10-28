from types import SimpleNamespace

import numpy as np

from simulations.quartet.moran_theory.conditional_patterson_d import (
    classify_genealogy,
    ratio_d,
    summarize_genealogies,
)


def _result(pair, t1=1.0, t2=3.0):
    other = tuple(x for x in ("P1", "P2", "P3", "O") if x not in pair)
    coals = (
        SimpleNamespace(descendant_set_1=pair[0], descendant_set_2=pair[1]),
        SimpleNamespace(descendant_set_1=",".join(pair), descendant_set_2=",".join(other)),
    )
    return SimpleNamespace(coalescences=coals, coalescence_times=(t1, t2))


def test_explicit_branch_mapping_for_abba_and_baba():
    assert classify_genealogy(_result(("P1", "P3"))).topology == 1
    assert classify_genealogy(_result(("P2", "P3"))).topology == 2
    assert classify_genealogy(_result(("P1", "P2"))).topology == 0


def test_branch_length_d_is_distinct_from_topology_d():
    summary = summarize_genealogies([_result(("P1", "P3"), 1, 2), _result(("P2", "P3"), 1, 5)])
    assert np.isclose(summary.d_topology, 0.0)
    assert summary.d_branch < 0.0
    assert summary.d_patterson > 0.0


def test_zero_denominator_is_nan():
    assert np.isnan(ratio_d(0.0, 0.0))
