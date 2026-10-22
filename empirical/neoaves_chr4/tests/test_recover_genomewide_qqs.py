import importlib.util
import math
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "recover_genomewide_qqs.py"
spec = importlib.util.spec_from_file_location("recover_genomewide_qqs", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_published_qqs_formula_and_normalization():
    result = module.qqs([8.0, 5.0, 3.0, 1.0])
    assert result is not None
    assert math.isclose(result["q1"] + result["q2"] + result["q3"], 1.0)
    assert result["d1"] == 7.0
    assert result["d2"] == 4.0
    assert result["d3"] == 2.0


def test_zero_denominator_is_invalid():
    assert module.qqs([1.0, 1.0, 1.0, 1.0]) is None


def test_mean_and_aggregate_estimands_can_differ():
    rows = [
        {"x1": 2.0, "x2": 1.0, "x3": 1.0, "x4": 0.0, "q1": 0.5, "q2": 0.25, "q3": 0.25},
        {"x1": 2.0, "x2": 2.0, "x3": 0.0, "x4": 0.0, "q1": 0.5, "q2": 0.5, "q3": 0.0},
    ]
    mean, aggregate = module.aggregate_q(rows)
    assert mean != aggregate
    assert math.isclose(sum(mean), 1.0)
    assert math.isclose(sum(aggregate), 1.0)
