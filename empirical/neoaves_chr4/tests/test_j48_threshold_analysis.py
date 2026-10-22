import importlib.util
import math
import sys
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "j48_threshold_analysis.py"
spec = importlib.util.spec_from_file_location("j48_threshold_analysis", SCRIPT)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)


def test_canonical_unrooted_splits_ignore_child_order():
    left = module.parse_newick("((a,b),(c,d));")
    right = module.parse_newick("((d,c),(b,a));")
    assert module.split_keys(left, {"a", "b", "c", "d"}) == module.split_keys(right, {"a", "b", "c", "d"})


def test_all_three_quartet_topologies_are_distinct():
    trees = [module.parse_newick(x) for x in ("((a,b),(c,d));", "((a,c),(b,d));", "((a,d),(b,c));")]
    keys = [module.split_keys(tree, {"a", "b", "c", "d"}) for tree in trees]
    assert len({frozenset(key) for key in keys}) == 3


def test_pairwise_crossover_and_invalid_cases():
    result = module.crossover([0.8, 0.1, 0.1], [0.1, 0.8, 0.1], 0, 1)
    assert result["status"] == "valid_interior"
    assert math.isclose(result["epsilon"], 0.5)
    assert module.crossover([0.5, 0.5, 0.0], [0.5, 0.5, 0.0], 0, 1)["status"] == "denominator_zero"
    assert module.crossover([0.8, 0.1, 0.1], [0.7, 0.2, 0.1], 0, 1)["status"] == "outside_[0,1]"
