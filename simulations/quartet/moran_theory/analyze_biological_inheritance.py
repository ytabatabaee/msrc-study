"""Targeted search under biologically constrained two-regime inheritance."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.stats import binom

from .two_regime_species_tree import (
    RegimeBranch, TwoRegimeModel, balanced_two_regime_quartet_probabilities,
    unbalanced_two_regime_quartet_probabilities,
)


def identity_kernel():
    return np.eye(2)


def symmetric_kernel(epsilon: float):
    if not 0 <= epsilon <= 1:
        raise ValueError("epsilon must be in [0,1]")
    return np.array([[1 - epsilon, epsilon], [epsilon, 1 - epsilon]], dtype=float)


def founding_classification_kernel(p_low: float, p_high: float, M: int, threshold: float | None = None):
    """Kernel from binomial founding counts followed by low/high classification.

    A daughter has X~Binomial(M,p_parent); it is classified low when X/M is
    at or below the midpoint of the two regime frequencies. This is an
    explicit approximation: the daughter branch is then represented by the
    fixed regime frequency associated with its class.
    """
    threshold = (p_low + p_high) / 2 if threshold is None else threshold
    cutoff = int(np.floor(M * threshold + 1e-12))
    out = np.zeros((2, 2), dtype=float)
    for parent, p in enumerate((p_low, p_high)):
        out[parent, 0] = binom.cdf(cutoff, M, p)
        out[parent, 1] = 1.0 - out[parent, 0]
    return out


def _result(shape, inheritance, p_low, p_high, rho, a, b, M, pi_low, durations):
    model = TwoRegimeModel(p_low, p_high, rho)
    branch = lambda t: RegimeBranch(t, a, b, M)
    root, b1, b2, terminal = durations
    pi = np.array([pi_low, 1 - pi_low])
    if shape == "balanced":
        return balanced_two_regime_quartet_probabilities(
            model=model, root_branch=branch(root), left_pair_branch=branch(b1), right_pair_branch=branch(b2),
            terminal_branches=tuple(branch(x) for x in terminal), initial_regime=pi,
            H_root_left=inheritance, H_root_right=inheritance, H_left_1=inheritance,
            H_left_2=inheritance, H_right_3=inheritance, H_right_4=inheritance)
    return unbalanced_two_regime_quartet_probabilities(
        model=model, root_branch=branch(root), intermediate_branch=branch(b1), sister_branch=branch(b2),
        terminal_branches=tuple(branch(x) for x in terminal), initial_regime=pi,
        H_intermediate_pair=inheritance, H_intermediate_3=inheritance,
        H_root_intermediate=inheritance, H_root_4=inheritance,
        H_pair_1=inheritance, H_pair_2=inheritance)


def targeted_grid():
    for shape in ("balanced", "unbalanced"):
        for p_low, p_high in ((.1, .9), (.2, .8), (.35, .65)):
            for rho in (0.0, .17, 1.0):
                for a, b in ((0.0, 0.0), (.2, .4), (.5, .5)):
                    for pi_low in (.1, .5):
                        durations = (.5, .8, .2, (.1, .4, .2, .3)) if shape == "unbalanced" else (.5, .2, .8, (.1, .4, .2, .3))
                        for M in (6, 10):
                            yield shape, p_low, p_high, rho, a, b, M, pi_low, durations


def run(output: Path):
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    for index, (shape, pl, ph, rho, a, b, M, pi, durations) in enumerate(targeted_grid()):
        models = [("identity", identity_kernel()), ("symmetric_0.02", symmetric_kernel(.02)),
                  ("symmetric_0.10", symmetric_kernel(.10)),
                  ("founding", founding_classification_kernel(pl, ph, M))]
        for name, H in models:
            q = _result(shape, H, pl, ph, rho, a, b, M, pi, durations).quartet_probabilities
            rows.append(dict(index=index, shape=shape, inheritance_model=name, p_low=pl, p_high=ph,
                             rho=rho, a=a, b=b, M=M, pi_low=pi, q1=q[0], q2=q[1], q3=q[2],
                             anomaly=max(q[1], q[2]) - q[0], asymmetry=abs(q[1] - q[2]),
                             normalization_error=q.sum() - 1))
    with (output / "search.tsv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    with (output / "summary.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(("inheritance_model", "shape", "max_anomaly", "max_asymmetry", "min_anomaly"))
        for model in ("identity", "symmetric_0.02", "symmetric_0.10", "founding"):
            subset = [r for r in rows if r["inheritance_model"] == model]
            for shape in ("balanced", "unbalanced"):
                s = [r for r in subset if r["shape"] == shape]
                w.writerow((model, shape, max(r["anomaly"] for r in s), max(r["asymmetry"] for r in s), min(r["anomaly"] for r in s)))
    _plot(output, rows)
    return rows


def _xy(q):
    return .5 * (2 * q[0] + q[2]), np.sqrt(3) * q[0] / 2


def _plot(output, rows):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.4), constrained_layout=True)
    colors = {"identity": "#0072B2", "symmetric_0.02": "#56B4E9", "symmetric_0.10": "#009E73", "founding": "#D55E00"}
    tri = np.array([_xy((1, 0, 0)), _xy((0, 1, 0)), _xy((0, 0, 1)), _xy((1, 0, 0))])
    for ax, shape in zip(axes, ("balanced", "unbalanced")):
        ax.plot(tri[:, 0], tri[:, 1], color="#555", lw=.8)
        for model, color in colors.items():
            s = [r for r in rows if r["shape"] == shape and r["inheritance_model"] == model]
            ax.scatter([_xy((r["q1"], r["q2"], r["q3"]))[0] for r in s],
                       [_xy((r["q1"], r["q2"], r["q3"]))[1] for r in s],
                       s=9, alpha=.45, color=color, edgecolors="none", label=model)
        ax.plot([1, .25], [np.sqrt(3)/2, 0], "--", color="#999", lw=.7)
        ax.set_aspect("equal"); ax.axis("off"); ax.set_title(shape)
        ax.text(.5, .91, "q1", ha="center", fontsize=8); ax.text(-.04, -.04, "q2", ha="right", fontsize=8); ax.text(1.04, -.04, "q3", ha="left", fontsize=8)
    axes[1].legend(frameon=False, fontsize=6, loc="lower left")
    fig.savefig(output / "biological_inheritance_simplex.pdf")
    fig.savefig(output / "biological_inheritance_simplex.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    run(Path(__file__).parent / "results" / "biological_inheritance")
