import numpy as np

from simulations.quartet.moran_theory.analyze_conditional_moran_asymmetry import (
    certificate_rows,
)


def test_asymmetry_certificate_is_normalized_and_conditioning_is_positive():
    rows = certificate_rows()
    assert len(rows) == 4
    for row in rows:
        assert row["P_Z"] > 0
        assert abs(row["normalization_error"]) < 1e-12
        assert min(row["q1"], row["q2"], row["q3"]) >= -1e-12


def test_exact_solver_can_move_r_to_both_sides_of_one_half():
    rows = {row["case"]: row for row in certificate_rows()}
    assert rows["r_below_half"]["r_q1_over_q1q3"] < 0.5
    assert rows["r_above_half"]["r_q1_over_q1q3"] > 0.5


def test_long_root_controls_q2_dominance_in_the_certificate_family():
    rows = {row["case"]: row for row in certificate_rows()}
    assert rows["dominant_probability_long_root"]["d_q2"] > rows["symmetric_zero_lower"]["d_q2"]
