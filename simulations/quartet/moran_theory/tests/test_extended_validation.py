import numpy as np

from ..validate_unbalanced_extended import _theory, conditions


def test_extended_grid_covers_requested_complexity_without_cartesian_product():
    grid = conditions()
    assert len(grid) == 23
    assert {row["p0"] for row in grid} >= {0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0}
    assert {row["rho"] for row in grid} >= {0.0, 0.05, 0.2, 0.5, 1.0}
    assert any(len(set(row["M_terminal"])) > 1 for row in grid)
    assert any(row["intermediate_length"] > row["sister_length"] for row in grid)
    assert any(row["sister_length"] > row["intermediate_length"] for row in grid)


def test_asymmetric_theory_cell_is_conservative():
    row = next(row for row in conditions() if row["condition"] == "M_asym_3")
    result = _theory(row)
    assert np.all(result.quartet_probabilities >= -2e-12)
    assert np.isclose(result.quartet_probabilities.sum(), 1.0, atol=2e-12)
    assert np.isclose(result.early_topology_contributions.sum() + result.root_topology_contributions.sum(), 1.0, atol=2e-12)


def test_terminal_duration_mapping_is_explicit_for_pectinate_tree():
    row = next(row for row in conditions() if row["condition"] == "dur_I_long")
    assert (0.0, 0.0, row["sister_length"], row["sister_length"] + row["intermediate_length"]) == (
        0.0, 0.0, 0.1, 0.9
    )
