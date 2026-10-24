#!/usr/bin/env python3
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from patterson_d_core import calculate_window, jackknife  # noqa: E402

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

def test_exclusions_and_zero_denominator():
    gt = diploid([[0,0],[0,0],[0,0],[0,0]], [[0,2],[0,0],[0,0],[0,0]], [[0,0],[-1,-1],[0,0],[0,0]], [[0,0],[0,0],[0,0],[0,0]])
    result = calculate_window(gt, {"P1":[0],"P2":[1],"P3":[2]}, [3])
    assert result.d is None and result.exclusions["multiallelic"] == 1 and result.exclusions["insufficient_group_calls"] >= 1
    polymorphic = calculate_window(diploid([[0,0],[0,0],[1,1],[0,0],[1,1]]), {"P1":[0],"P2":[1],"P3":[2]}, [3,4])
    assert polymorphic.d is None and polymorphic.exclusions["outgroup_polymorphic"] == 1

def test_jackknife():
    result = jackknife([(10.0, 0.0, 10), (0.0, 10.0, 10), (5.0, 5.0, 10)])
    assert result["n_blocks"] == 3 and result["se"] is not None and result["pooled_d"] == 0.0
