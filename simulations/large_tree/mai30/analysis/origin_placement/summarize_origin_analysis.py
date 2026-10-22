#!/usr/bin/env python3
"""Summarize the executed forward audit and create honest status figures."""

from __future__ import annotations

import csv
import statistics
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(name):
    with (HERE / name).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(name, rows):
    with (HERE / name).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def main():
    persistence = read("persistence_results.tsv")
    exposure = read("quartet_exposure.tsv")
    summary = []
    for key in sorted({(r["condition"], r["replicate"], r["scenario"]) for r in persistence}):
        rows = [r for r in persistence if (r["condition"], r["replicate"], r["scenario"]) == key and r.get("status") == "executed_forward_history"]
        if not rows:
            continue
        summary.append({"condition": key[0], "replicate": key[1], "scenario": key[2], "n_histories": len(rows), "mean_origin_polymorphic_fraction": statistics.mean(float(r["origin_branch_polymorphic_fraction"]) for r in rows), "fraction_origin_end_polymorphic": statistics.mean(float(r["origin_branch_end_polymorphic"]) for r in rows), "fraction_origin_ever_fixed": statistics.mean(float(r["origin_branch_ever_fixed"]) for r in rows), "fraction_origin_ever_lost": statistics.mean(float(r["origin_branch_ever_lost"]) for r in rows), "fraction_root_terminal_polymorphic": statistics.mean(float(r["root_terminal_polymorphic"]) for r in rows)})
    write("persistence_summary.tsv", summary)
    write("quartet_margins.tsv", [{"status": "not_executed", "reason": "30-taxon backward genealogy discovery was stopped after local runtime stress; no population-level quartet margins are reported"}])
    write("quartet_reversal_results.tsv", [{"status": "not_executed", "reason": "no validated genealogy-level quartet estimates"}])
    write("astral_results.tsv", [{"status": "not_executed", "reason": "no origin-placement genealogy collection met the discovery gate"}])

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from collections import OrderedDict
    figdir = HERE / "figures"; figdir.mkdir(exist_ok=True)
    colors = {"A_existing_internal": "#377eb8", "B_deep_internal": "#4daf4a", "C_root_adjacent": "#ff7f00", "D_older_root_extended": "#984ea3"}
    labels = {"A_existing_internal": "A existing", "B_deep_internal": "B deep", "C_root_adjacent": "C root-adjacent", "D_older_root_extended": "D older/extended"}

    agg = defaultdict(list)
    for row in exposure:
        agg[row["scenario"]].append((float(row["descendant_taxa"]), float(row["fraction_structurally_exposed"])))
    scenarios = list(labels)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for i, s in enumerate(scenarios):
        vals = agg[s]; ax.scatter([i] * len(vals), [v[0] for v in vals], color=colors[s], s=35)
    ax.set_xticks(range(4), [labels[s] for s in scenarios], rotation=25, ha="right"); ax.set_ylabel("descendant taxa"); ax.set_title("Origin placement on the four Mai pilot trees")
    fig.tight_layout(); fig.savefig(figdir / "origin_placement_schematic.pdf"); fig.savefig(figdir / "origin_placement_schematic.png", dpi=220); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for row in persistence:
        if row.get("status") != "executed_forward_history": continue
        s = row["scenario"]; ax.scatter(float(row["origin_age"]), float(row["origin_branch_polymorphic_fraction"]), color=colors[s], label=labels[s], alpha=.8)
    handles, labels_seen = ax.get_legend_handles_labels(); unique = dict(zip(labels_seen, handles)); ax.legend(unique.values(), unique.keys(), frameon=False, fontsize=8)
    ax.set_xlabel("origin age (backward time units)"); ax.set_ylabel("fraction of origin-branch records polymorphic"); ax.set_title("Forward persistence audit; one history per scenario/tree")
    fig.tight_layout(); fig.savefig(figdir / "persistence_origin_age.pdf"); fig.savefig(figdir / "persistence_origin_age.png", dpi=220); plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    for i, s in enumerate(scenarios):
        vals = [v[1] for v in agg[s]]; ax.scatter([i] * len(vals), vals, color=colors[s], s=35)
    ax.set_xticks(range(4), [labels[s] for s in scenarios], rotation=25, ha="right"); ax.set_ylabel("fraction of all 27,405 quartets with ≥2 descendants"); ax.set_ylim(0, 1.05); ax.set_title("Structural quartet exposure")
    fig.tight_layout(); fig.savefig(figdir / "quartet_exposure.pdf"); fig.savefig(figdir / "quartet_exposure.png", dpi=220); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5)); ax.text(.5, .55, "No validated genealogy-level\nquartet margins executed", ha="center", va="center", fontsize=15); ax.text(.5, .35, "rho=.1 and 1.0 remain planned stress points", ha="center", va="center", fontsize=10); ax.set_axis_off(); fig.tight_layout(); fig.savefig(figdir / "quartet_reversal_phase_diagram.pdf"); fig.savefig(figdir / "quartet_reversal_phase_diagram.png", dpi=220); plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4.5)); ax.text(.5, .55, "Origin-placement ASTRAL validation\nnot executed", ha="center", va="center", fontsize=15); ax.text(.5, .35, "No topology claim is made from the forward audit", ha="center", va="center", fontsize=10); ax.set_axis_off(); fig.tight_layout(); fig.savefig(figdir / "astral_accuracy.pdf"); fig.savefig(figdir / "astral_accuracy.png", dpi=220); plt.close(fig)


if __name__ == "__main__":
    main()
