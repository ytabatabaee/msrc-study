"""Summarize completed conditional-MS​​RC Patterson's D history blocks.

This postprocessor only combines completed TSV rows; it never simulates or
selects histories.  The pilot directory is accepted so the same code can
produce a reproducible small-run diagnostic before the Bridges-2 array runs.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


NUMERIC = ("q1", "q2", "q3", "D_patterson_branch", "D_topology_q2_minus_q3",
           "D_branch_q2_minus_q3", "D_patterson_se", "ABBA_branch_sum",
           "BABA_branch_sum", "history_acceptance_fraction")


def read_rows(paths):
    rows = []
    for path in sorted(paths):
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("condition") == "conditional_msrc":
                    for key in NUMERIC:
                        if row.get(key, "") not in ("", "nan", "NaN"):
                            row[key] = float(row[key])
                    rows.append(row)
    if not rows:
        raise ValueError("No conditional_msrc rows found in the supplied result files")
    return rows


def summarize(rows):
    def stats(key):
        values = np.asarray([row[key] for row in rows], dtype=float)
        return {"mean": float(np.mean(values)), "median": float(np.median(values)),
                "sd": float(np.std(values, ddof=1)) if len(values) > 1 else float("nan"),
                "min": float(np.min(values)), "max": float(np.max(values))}
    d = np.asarray([row["D_patterson_branch"] for row in rows], dtype=float)
    se = np.asarray([row["D_patterson_se"] for row in rows], dtype=float)
    return {
        "n_histories": len(rows),
        "D_patterson_branch": stats("D_patterson_branch"),
        "D_topology_q2_minus_q3": stats("D_topology_q2_minus_q3"),
        "D_branch_q2_minus_q3": stats("D_branch_q2_minus_q3"),
        "mean_ABBA_branch_sum": float(np.mean([row["ABBA_branch_sum"] for row in rows])),
        "mean_BABA_branch_sum": float(np.mean([row["BABA_branch_sum"] for row in rows])),
        "fraction_genealogy_CI_excludes_zero": float(np.mean(np.abs(d) > 1.96 * se)),
        "fraction_negative_D": float(np.mean(d < 0)),
    }


def write_summary(rows, output: Path):
    output.mkdir(parents=True, exist_ok=True)
    summary = summarize(rows)
    with (output / "history_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["metric", "value"])
        writer.writerow(["n_histories", summary["n_histories"]])
        for key, value in summary.items():
            if isinstance(value, dict):
                for stat, stat_value in value.items():
                    writer.writerow([f"{key}_{stat}", stat_value])
            elif key != "n_histories":
                writer.writerow([key, value])
    return summary


def plot(rows, output: Path):
    import matplotlib.pyplot as plt
    output.mkdir(parents=True, exist_ok=True)
    d = np.asarray([r["D_patterson_branch"] for r in rows])
    qdiff = np.asarray([r["D_topology_q2_minus_q3"] for r in rows])
    abba = np.asarray([r["ABBA_branch_sum"] for r in rows])
    baba = np.asarray([r["BABA_branch_sum"] for r in rows])
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.2), constrained_layout=True)
    axes[0].hist(d, bins=min(12, max(4, len(d) // 3)), color="#386cb0", alpha=.75)
    axes[0].axvline(0, color="black", lw=.8); axes[0].set(xlabel="Patterson's D (ABBA−BABA)", ylabel="Histories")
    axes[1].scatter(qdiff, d, s=22, alpha=.75, color="#7b3294", edgecolors="none")
    axes[1].axhline(0, color="black", lw=.8); axes[1].axvline(0, color="black", lw=.8)
    axes[1].set(xlabel="q2−q3 diagnostic", ylabel="Patterson's D")
    axes[2].scatter(abba, baba, s=22, alpha=.75, color="#fdae61", edgecolors="none")
    lim = max(float(np.max(abba)), float(np.max(baba))) * 1.05
    axes[2].plot([0, lim], [0, lim], "k--", lw=.8)
    axes[2].set(xlabel="ABBA branch sum", ylabel="BABA branch sum")
    fig.savefig(output / "conditional_patterson_d_history_summary.pdf")
    fig.savefig(output / "conditional_patterson_d_history_summary.png", dpi=220)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    paths = []
    for item in args.inputs:
        paths.extend(item.rglob("pilot_results.tsv") if item.is_dir() else [item])
    rows = read_rows(paths)
    summary = write_summary(rows, args.output)
    plot(rows, args.output)
    print(summary)


if __name__ == "__main__":
    main()
