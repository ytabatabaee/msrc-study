"""Audit the proposed three-component internal mixture construction."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .verify_simplex_coverage import crossed_pattern_quartet


def target_weights(target: np.ndarray, T: float) -> np.ndarray:
    target = np.asarray(target, dtype=float)
    floor = np.exp(-2 * T) / 3
    weights = (target - floor) / (1 - 3 * floor)
    if np.any(target < 0) or not np.isclose(target.sum(), 1):
        raise ValueError("target must be a simplex point")
    if np.any(weights < -1e-12):
        raise ValueError("T is too short for this target")
    return weights


def external_three_component_mixture(target: np.ndarray, T: float) -> np.ndarray:
    """Algebraic mixture of three separately parameterized MSRC models.

    This is intentionally not an internal two-regime model; it is a reference
    calculation showing why the target-weight formula is algebraically valid.
    """
    w = target_weights(target, T)
    components = np.array([crossed_pattern_quartet(z, T) for z in ("0011", "0101", "0110")])
    return w @ components


def two_state_selector_obstruction() -> dict[str, int | str]:
    return {
        "root_regime_states": 2,
        "independently_selectable_vertex_components": 2,
        "conclusion": "A three-component selector is not encoded by the two-state regime variable.",
    }


def run(output: Path):
    output.mkdir(parents=True, exist_ok=True)
    targets = {
        "near_Q1": np.array([.96, .02, .02]),
        "near_Q2": np.array([.02, .96, .02]),
        "near_Q3": np.array([.02, .02, .96]),
        "center": np.array([1 / 3, 1 / 3, 1 / 3]),
        "interior": np.array([.2, .5, .3]),
        "near_edge": np.array([.499999, .499999, .000002]),
    }
    T = 4.0
    with (output / "target_certificates.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t")
        w.writerow(("target", "q1", "q2", "q3", "T", "w1", "w2", "w3", "max_external_mixture_error"))
        for name, target in targets.items():
            if np.min(target) < np.exp(-2 * T) / 3:
                T_use = 8.0
            else:
                T_use = T
            fitted = external_three_component_mixture(target, T_use)
            w.writerow((name, *target, T_use, *target_weights(target, T_use), float(np.max(np.abs(fitted - target)))))
    with (output / "selector_obstruction.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t"); w.writerow(("quantity", "value"))
        for key, value in two_state_selector_obstruction().items(): w.writerow((key, value))


if __name__ == "__main__":
    run(Path(__file__).parent / "results" / "internal_mixture")
