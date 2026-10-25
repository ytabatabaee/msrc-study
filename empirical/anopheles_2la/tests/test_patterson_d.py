#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np
import importlib.util

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from patterson_d_core import calculate_window, jackknife, jackknife_difference  # noqa: E402

spec = importlib.util.spec_from_file_location("run_patterson_d", Path(__file__).resolve().parents[1] / "scripts" / "run_patterson_d.py")
run_patterson_d = importlib.util.module_from_spec(spec); spec.loader.exec_module(run_patterson_d)

def diploid(*rows):
    return np.asarray(rows, dtype=np.int8)

def test_limits_and_symmetry():
    groups = {"P1": [0], "P2": [1], "P3": [2]}
    d_plus = calculate_window(diploid([[0,0],[1,1],[1,1],[0,0]]), groups, [3])
    assert d_plus.d == 1.0
    d_minus = calculate_window(diploid([[1,1],[0,0],[1,1],[0,0]]), groups, [3])
    assert d_minus.d == -1.0
    equal = calculate_window(diploid([[0,0],[1,1],[1,1],[0,0]], [[1,1],[0,0],[1,1],[0,0]]), groups, [3])
    assert abs(equal.d) < 1e-12
    swapped = calculate_window(diploid([[1,1],[0,0],[1,1],[0,0]]), {"P1":[1],"P2":[0],"P3":[2]}, [3])
    assert swapped.d == 1.0

def test_ref_alt_recode_and_missing_heterozygote():
    gt = diploid([[0,1],[1,1],[1,1],[0,0]], [[-1,-1],[0,1],[1,1],[0,0]])
    recoded = 1 - gt
    recoded[gt < 0] = -1
    groups = {"P1":[0],"P2":[1],"P3":[2]}
    a = calculate_window(gt, groups, [3]); b = calculate_window(recoded, groups, [3])
    assert abs(a.d - b.d) < 1e-12 and a.n_pass == 1

def test_multiple_diploid_samples_are_counted_by_complete_genotype():
    gt = diploid([[0, 0], [1, 1], [1, 1], [0, 0], [0, 1], [1, 1]])
    result = calculate_window(gt, {"P1": [0, 4], "P2": [1], "P3": [2]}, [3])
    assert result.n_pass == 1 and result.n_informative == 1

def test_exclusions_and_zero_denominator():
    gt = diploid([[0,0],[0,0],[0,0],[0,0]], [[0,2],[0,0],[0,0],[0,0]], [[0,0],[-1,-1],[0,0],[0,0]], [[0,0],[0,0],[0,0],[0,0]])
    result = calculate_window(gt, {"P1":[0],"P2":[1],"P3":[2]}, [3])
    assert result.d is None and result.exclusions["multiallelic"] == 1 and result.exclusions["insufficient_group_calls"] >= 1
    polymorphic = calculate_window(diploid([[0,0],[0,0],[1,1],[0,0],[1,1]]), {"P1":[0],"P2":[1],"P3":[2]}, [3,4])
    assert polymorphic.d is None and polymorphic.exclusions["outgroup_polymorphic"] == 1
    partial_outgroup = calculate_window(diploid([[0,0],[1,1],[1,1],[-1,0]]), {"P1":[0],"P2":[1],"P3":[2]}, [3])
    assert partial_outgroup.exclusions["outgroup_missing"] == 1

def test_jackknife():
    result = jackknife([(10.0, 0.0, 10), (0.0, 10.0, 10), (5.0, 5.0, 10)])
    assert result["n_blocks"] == 3 and result["se"] is not None and result["pooled_d"] == 0.0
    filtered = jackknife([(10.0, 0.0, 10), (0.0, 0.0, 0), (0.0, 10.0, 10)])
    assert filtered["n_blocks"] == 2 and filtered["pooled_d"] == 0.0
    difference = jackknife_difference([(10.0, 0.0, 0.0, 10.0, 10, 10), (0.0, 10.0, 10.0, 0.0, 10, 10), (0.0, 0.0, 0.0, 0.0, 0, 0)])
    assert difference["n_blocks"] == 2 and difference["delta_d"] == 0.0

def test_pilot_queries_are_exactly_eight_independent_50kb_windows():
    frozen = run_patterson_d.windows(True)
    regions = run_patterson_d.query_regions(frozen, True, 1_000_000)
    assert len(regions) == 8 and all(len(group) == 1 for group in regions)
    assert sum(int(group[0]["end"]) - int(group[0]["start"]) + 1 for group in regions) == 400_000
    assert [group[0]["window_id"] for group in regions] == [w["window_id"] for w in frozen]
