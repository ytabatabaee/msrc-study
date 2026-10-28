"""Four-taxon MSC-with-introgression quartet comparison.

The model is a two-component mixture: with probability 1-gamma the locus
follows Q1 and with probability gamma it follows an alternative tree Q2 or Q3.
Each component uses the standard MSC quartet probabilities
(1-2 exp(-t)/3, exp(-t)/3, exp(-t)/3).
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


def msc_quartet(t: float, favored: int = 0) -> np.ndarray:
    e = np.exp(-float(t))
    q = np.full(3, e / 3.0)
    q[favored] = 1.0 - 2.0 * e / 3.0
    return q


def introgression_quartet(gamma: float, t_base: float, t_intro: float, intro_topology: int = 1) -> np.ndarray:
    if intro_topology not in (1, 2) or not 0 <= gamma <= 1 or t_base < 0 or t_intro < 0:
        raise ValueError("invalid introgression parameters")
    return (1 - gamma) * msc_quartet(t_base, 0) + gamma * msc_quartet(t_intro, intro_topology)


def fit_introgression(target: np.ndarray, intro_topology: int = 2):
    target = np.asarray(target, dtype=float)
    fit = least_squares(lambda x: introgression_quartet(*x, intro_topology) - target,
                        x0=(.5, .3, .3), bounds=([0, 0, 0], [1, 30, 30]),
                        xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=10000)
    return fit.x, introgression_quartet(*fit.x, intro_topology), float(np.max(np.abs(fit.fun)))


def write_witness(path: Path, target: np.ndarray, intro_topology: int):
    params, fitted, error = fit_introgression(target, intro_topology)
    with path.open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(("target_q1", "target_q2", "target_q3", "gamma", "t_base", "t_intro", "intro_topology", "max_abs_error"))
        w.writerow((*target, *params, intro_topology + 1, error))
    return params, fitted, error


if __name__ == "__main__":
    out = Path(__file__).parent / "results" / "biological_inheritance"
    out.mkdir(parents=True, exist_ok=True)
    targets = {
        "legacy_arbitrary_inheritance": (np.array([.35417018362959773, .26868939582086815, .37714042054953417]), 2),
        "identity_inheritance": (np.array([.4329387873725959, .2836441797083963, .2834170329190079]), 1),
    }
    for name, (target, topology) in targets.items():
        print(name, write_witness(out / f"introgression_{name}.tsv", target, topology))
