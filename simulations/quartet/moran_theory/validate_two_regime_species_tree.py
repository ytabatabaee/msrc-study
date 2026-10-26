"""Deterministic validation and compact diagnostics for the two-regime model."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from .prescribed_frequency_quartet import PrescribedTrajectory, balanced_prescribed_quartet_probabilities
from .two_regime_species_tree import (
    RegimeBranch, TwoRegimeModel, balanced_two_regime_quartet_probabilities,
    unbalanced_two_regime_quartet_probabilities,
)


def fixed_balanced(p: float) -> np.ndarray:
    b = lambda t: PrescribedTrajectory.constant(p, t)
    return balanced_prescribed_quartet_probabilities(
        M_root=10, M_left=10, M_right=10, M_terminal=(10,) * 4, rho=.17,
        terminal_trajectories=(b(.2), b(.4), b(.1), b(.3)),
        sister_trajectories=(b(.3), b(.1)), root_trajectory=b(.5),
    ).quartet_probabilities


def run(outdir: Path) -> list[dict[str, object]]:
    outdir.mkdir(parents=True, exist_ok=True)
    model = TwoRegimeModel(.2, .8, .17)
    active = lambda t: RegimeBranch(t, .2, .4, 10)
    result = balanced_two_regime_quartet_probabilities(
        model=model, root_branch=active(.5), left_pair_branch=active(.3), right_pair_branch=active(.1),
        terminal_branches=(active(.2), active(.4), active(.1), active(.3)),
        initial_regime=np.array([.6, .4]), with_joint=True,
    )
    rows = [{"check": "balanced_probability_sum", "value": float(result.quartet_probabilities.sum()), "error": abs(float(result.quartet_probabilities.sum() - 1))}]
    rows.append({"check": "balanced_joint_sum", "value": float(result.joint.sum()), "error": abs(float(result.joint.sum() - 1))})
    for z, probs in zip(result.pattern_probabilities, result.conditional_quartet_probabilities):
        rows.append({"check": "conditional_pattern", "value": float(z), "error": float(np.nanmax(np.abs(probs.sum() - 1)))})

    no_switch = lambda t: RegimeBranch(t, 0, 0, 10)
    mixture = balanced_two_regime_quartet_probabilities(
        model=model, root_branch=no_switch(.5), left_pair_branch=no_switch(.3), right_pair_branch=no_switch(.1),
        terminal_branches=(no_switch(.2), no_switch(.4), no_switch(.1), no_switch(.3)),
        initial_regime=np.array([.6, .4]),
    ).quartet_probabilities
    fixed_mix = .6 * fixed_balanced(.2) + .4 * fixed_balanced(.8)
    rows.extend([
        {"check": "p_equal_weighted_fixed_regime", "value": float(np.max(np.abs(mixture - fixed_mix))), "error": float(np.max(np.abs(mixture - fixed_mix)))},
        {"check": "conditional_delta_max", "value": float(np.nanmax(np.abs(result.conditional_quartet_probabilities[:, 1] - result.conditional_quartet_probabilities[:, 2]))), "error": 0.0},
        {"check": "conditional_anomaly_max", "value": float(np.nanmax(np.maximum(result.conditional_quartet_probabilities[:, 1], result.conditional_quartet_probabilities[:, 2]) - result.conditional_quartet_probabilities[:, 0])), "error": 0.0},
    ])
    with (outdir / "validation.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("check", "value", "error"), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    _plot(outdir, result.quartet_probabilities, mixture, fixed_mix)
    return rows


def _plot(outdir: Path, regime: np.ndarray, no_switch: np.ndarray, fixed_mix: np.ndarray) -> None:
    import matplotlib.pyplot as plt
    labels = ("q1", "q2", "q3")
    x = np.arange(3); width = .25
    fig, ax = plt.subplots(figsize=(5.4, 3.1), constrained_layout=True)
    ax.bar(x - width, regime, width, label="switching regimes", color="#0072B2")
    ax.bar(x, no_switch, width, label="no regime switching", color="#D55E00")
    ax.bar(x + width, fixed_mix, width, label="fixed-regime mixture", color="#009E73", alpha=.8)
    ax.set_xticks(x, labels); ax.set_ylim(0, 1); ax.set_ylabel("quartet probability")
    ax.legend(frameon=False, fontsize=8); ax.spines[["top", "right"]].set_visible(False)
    fig.savefig(outdir / "two_regime_quartet_comparison.pdf")
    fig.savefig(outdir / "two_regime_quartet_comparison.png", dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results" / "two_regime")
    args = parser.parse_args(); run(args.output); print(args.output)


if __name__ == "__main__":
    main()
