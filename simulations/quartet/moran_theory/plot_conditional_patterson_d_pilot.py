"""Plot the small conditional Patterson's D pilot without smoothing."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_rows(path: Path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def plot(input_path: Path, output_stem: Path):
    rows = read_rows(input_path)
    labels = ["H1", "H2", "H3", "MSC"]
    selected = rows[:-1] + rows[-1:]
    d_branch = [float(row["D_branch_q2_minus_q3"]) for row in selected]
    d_sites = [float(row["D_sites"]) for row in selected]
    q2 = [float(row["q2"]) for row in selected]
    q3 = [float(row["q3"]) for row in selected]
    x = np.arange(len(labels))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), constrained_layout=True)
    width = 0.24
    axes[0].bar(x - width, d_branch, width, label="branch-length D")
    axes[0].bar(x, d_sites, width, label="sequence D")
    axes[0].bar(x + width, [float(row["D_topology_q2_minus_q3"]) for row in selected], width,
                label="topology diagnostic")
    axes[0].axhline(0, color="#444", lw=.8)
    axes[0].set_xticks(x, labels)
    axes[0].set_ylabel("signed statistic")
    axes[0].set_title("Conditional D diagnostics")
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].plot(x, q2, "o-", label="q2 / BABA topology")
    axes[1].plot(x, q3, "o-", label="q3 / ABBA topology")
    axes[1].set_xticks(x, labels)
    axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("topology frequency")
    axes[1].set_title("Topology component")
    axes[1].legend(frameon=False, fontsize=7)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem.with_suffix(".pdf"))
    fig.savefig(output_stem.with_suffix(".png"), dpi=240)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=Path(__file__).parent / "results" / "conditional_patterson_d" / "pilot_results.tsv")
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "results" / "conditional_patterson_d" / "conditional_patterson_d_pilot")
    args = parser.parse_args()
    plot(args.input, args.output)


if __name__ == "__main__":  # pragma: no cover
    main()
