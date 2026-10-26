"""Plot the compact prescribed-frequency validation diagnostics."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from simulations.quartet.moran_theory.full_species_tree import fixed_p_half_species_tree_quartet_probabilities
from simulations.quartet.moran_theory.prescribed_frequency_quartet import PrescribedTrajectory, balanced_prescribed_quartet_probabilities


def main() -> None:
    out = Path(__file__).parent / "results"
    out.mkdir(exist_ok=True)
    constant = PrescribedTrajectory.constant(0.5, 0.5)
    prescribed = balanced_prescribed_quartet_probabilities(
        M_root=10, M_left=10, M_right=10, M_terminal=(10,) * 4, rho=.17,
        terminal_trajectories=(constant,) * 4,
        sister_trajectories=(PrescribedTrajectory.constant(.5, .3), PrescribedTrajectory.constant(.5, .1)),
        root_trajectory=constant,
    ).quartet_probabilities
    fixed = fixed_p_half_species_tree_quartet_probabilities(
        M=10, rho=.17, root_length=.5, sister_lengths=(.3, .1), terminal_lengths=(.5,) * 4,
    )
    varying = balanced_prescribed_quartet_probabilities(
        M_root=10, M_left=10, M_right=10, M_terminal=(10,) * 4, rho=.17,
        terminal_trajectories=tuple(PrescribedTrajectory(((.5, .1), (p, d))) for p, d in ((.2, .2), (.8, .15), (.35, .25), (.65, .1))),
        sister_trajectories=(PrescribedTrajectory(((.45, .2), (.55, .1))), PrescribedTrajectory(((.6, .2),))),
        root_trajectory=PrescribedTrajectory(((.5, .3), (.4, .2))),
    ).quartet_probabilities
    labels = ["q1", "q2", "q3"]
    fig, axes = plt.subplots(1, 3, figsize=(8.0, 2.7), constrained_layout=True)
    colors = ["#0072B2", "#D55E00", "#009E73"]
    axes[0].plot(labels, fixed, "o-", color="#555555", label="fixed-frequency")
    axes[0].plot(labels, prescribed, "s", color="#0072B2", label="prescribed")
    axes[0].set_title("Constant frequency")
    axes[0].set_ylabel("quartet probability")
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].bar(labels, np.abs(prescribed - fixed), color=colors)
    axes[1].set_title("Absolute error")
    axes[1].set_ylabel("absolute difference")
    axes[1].set_yscale("log")
    axes[2].bar(labels, varying, color=colors)
    axes[2].set_title("Piecewise history")
    axes[2].set_ylabel("quartet probability")
    for ax in axes:
        if ax is not axes[1]:
            ax.set_ylim(bottom=0)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=8)
    fig.savefig(out / "prescribed_frequency_validation.pdf")
    fig.savefig(out / "prescribed_frequency_validation.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
