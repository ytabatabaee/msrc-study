#!/usr/bin/env python3
"""Summarize retained mechanism diagnostics and make audit figures."""

from __future__ import annotations

import csv
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
OUT = ROOT / "simulations/large_tree/mai30/analysis/mechanism_diagnostics"

TREATMENTS = {
    "no_suppression": "msrc_no_suppression",
    "moderate_suppression": "msrc_moderate_suppression",
    "strong_suppression": "msrc_strong_suppression",
}


def treatment_from_group(group):
    for part in group.split("__"):
        if part.startswith("msrc_"):
            return part
    return ""


def read(path):
    with (OUT / path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows):
    with (OUT / path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)


def main():
    switches = read("switching_diagnostics.tsv")
    persistence = read("persistence_probabilities.tsv")
    patterns = read("arrangement_pattern_probabilities.tsv")
    conditionals = read("conditional_quartet_probabilities.tsv")
    astral = read("astral_results.tsv")
    species_rows = {r["group"]: r for r in read("../results/mai30_msrc_vs_species.tsv")} if (OUT / "../results/mai30_msrc_vs_species.tsv").exists() else {}

    by_treatment = defaultdict(list)
    for row in switches:
        by_treatment[row["treatment"]].append(row)
    switching_summary = []
    for treatment, rows in sorted(by_treatment.items()):
        values = [float(r["mean_switches_per_genealogy"]) for r in rows]
        zeros = [float(r["zero_switch_fraction"]) for r in rows]
        switching_summary.append({
            "treatment": treatment, "rho": rows[0]["rho"], "n_histories": len(rows),
            "mean_switches_per_genealogy": statistics.mean(values),
            "sd_history_mean_switches": statistics.stdev(values),
            "mean_zero_switch_fraction": statistics.mean(zeros),
            "mean_max_switches": statistics.mean(float(r["max_switches"]) for r in rows),
        })
    write("switching_summary.tsv", switching_summary)

    persistence_summary = []
    for treatment, branch in sorted({(r["treatment"], r["branch_id"]) for r in persistence}):
        subset = [r for r in persistence if r["treatment"] == treatment and r["branch_id"] == branch]
        persistence_summary.append({
            "treatment": treatment, "rho": subset[0]["rho"], "branch_id": branch,
            "n_histories": len(subset),
            "p_polymorphic_at_branch_start": statistics.mean(float(r["polymorphic_at_branch_start"]) for r in subset),
            "p_polymorphic_at_branch_end": statistics.mean(float(r["polymorphic_at_branch_end"]) for r in subset),
            "mean_fraction_recorded_generations_polymorphic": statistics.mean(float(r["fraction_recorded_generations_polymorphic"]) for r in subset),
            "mean_integrated_hazard_0_to_1": statistics.mean(float(r["integrated_hazard_0_to_1"]) for r in subset),
            "mean_integrated_hazard_1_to_0": statistics.mean(float(r["integrated_hazard_1_to_0"]) for r in subset),
        })
    write("persistence_summary.tsv", persistence_summary)

    pattern_summary = []
    for treatment, category in sorted({(treatment_from_group(r["group"]), r["category"]) for r in patterns}):
        subset = [r for r in patterns if treatment_from_group(r["group"]) == treatment and r["category"] == category]
        pattern_summary.append({"treatment": treatment, "category": category, "n_quartet_group_pattern_rows": len(subset), "mean_pattern_probability": statistics.mean(float(r["probability"]) for r in subset)})
    write("pattern_summary.tsv", pattern_summary)

    # Observed pilot points only; no unexecuted rate-grid values are inserted.
    rows = []
    for r in astral:
        treatment = next((part for part in r["group"].split("__") if part.startswith("msrc_")), "")
        if treatment in {"msrc_no_suppression", "msrc_moderate_suppression", "msrc_strong_suppression"}:
            condition = r["group"].split("__")[0]
            replicate = next(part.removeprefix("replicate_") for part in r["group"].split("__") if part.startswith("replicate_"))
            rows.append({"group": r["group"], "condition": condition, "replicate": replicate, "treatment": treatment, "rho": {"msrc_no_suppression": 0.02, "msrc_moderate_suppression": 0.01, "msrc_strong_suppression": 0.002}[treatment], "gene_tree_normalized_rf": species_rows.get(r["group"], {}).get("mean_normalized_rf", ""), "astral_normalized_rf": r["normalized_rf"], "analysis_status": "completed_pilot"})
    write("gene_tree_vs_astral_summary.tsv", rows)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = {"msrc_no_suppression": "#377eb8", "msrc_moderate_suppression": "#ff8c42", "msrc_strong_suppression": "#d73027"}
    labels = {"msrc_no_suppression": "κ=1 (ρ=.02)", "msrc_moderate_suppression": "κ=.5 (ρ=.01)", "msrc_strong_suppression": "κ=.1 (ρ=.002)"}
    OUT.joinpath("figures").mkdir(exist_ok=True)

    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for treatment, rs in by_treatment.items():
        ax.scatter([float(r["rho"]) for r in rs], [float(r["mean_switches_per_genealogy"]) for r in rs], s=18, alpha=.45, color=colors[treatment], label=labels[treatment])
    ax.set_xscale("log"); ax.set_xlabel("effective switching intensity ρ"); ax.set_ylabel("switches per genealogy"); ax.set_title("Mai30 retained pilot: realized switching")
    ax.legend(frameon=False, fontsize=8); fig.tight_layout()
    fig.savefig(OUT / "figures/switching_behavior.pdf"); fig.savefig(OUT / "figures/switching_behavior.png", dpi=220); plt.close(fig)

    # A compact regime map for the four retained species-tree replicates.  It
    # displays observed pilot switching only; the unexecuted rate-grid rows
    # are intentionally absent.
    replicate_order = ["007", "042", "074", "005"]
    treatment_order = ["msrc_no_suppression", "msrc_moderate_suppression", "msrc_strong_suppression"]
    matrix = []
    for replicate in replicate_order:
        row = []
        for treatment in treatment_order:
            values = [float(r["mean_switches_per_genealogy"]) for r in switches if r["replicate"] == replicate and r["treatment"] == treatment]
            row.append(statistics.mean(values) if values else float("nan"))
        matrix.append(row)
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    image = ax.imshow(matrix, aspect="auto", cmap="magma")
    ax.set_xticks(range(3), [labels[t] for t in treatment_order], rotation=25, ha="right")
    ax.set_yticks(range(4), [f"replicate {r}" for r in replicate_order])
    ax.set_xlabel("retained pilot treatment")
    ax.set_title("Observed switching regimes; retained pilot only")
    fig.colorbar(image, ax=ax, label="mean switches per genealogy")
    fig.tight_layout(); fig.savefig(OUT / "figures/observed_switching_regimes.pdf"); fig.savefig(OUT / "figures/observed_switching_regimes.png", dpi=220); plt.close(fig)

    ps = [r for r in persistence_summary if r["branch_id"] == "N8"]
    cats = [r for r in pattern_summary if r["category"] == "discordant_2_2"]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 4.2))
    axes[0].bar([labels[r["treatment"]] for r in ps], [float(r["p_polymorphic_at_branch_end"]) for r in ps], color=[colors[r["treatment"]] for r in ps]); axes[0].set_ylabel("P(polymorphic at N8 younger end)"); axes[0].tick_params(axis="x", rotation=35)
    axes[1].bar([labels[r["treatment"]] for r in cats], [float(r["mean_pattern_probability"]) for r in cats], color=[colors[r["treatment"]] for r in cats]); axes[1].set_ylabel("mean P(discordant 2:2 Z)"); axes[1].tick_params(axis="x", rotation=35)
    fig.suptitle("Persistence and terminal arrangement sorting"); fig.tight_layout(); fig.savefig(OUT / "figures/persistence_sorting.pdf"); fig.savefig(OUT / "figures/persistence_sorting.png", dpi=220); plt.close(fig)

    # Use the species-tree topology row for an interpretable conditional-vs-
    # marginal panel.  Values are aggregated over panel quartets and groups.
    cond_plot = []
    for treatment in colors:
        for category in ["uniform_4_0", "3_1", "concordant_2_2", "discordant_2_2"]:
            rs = [r for r in conditionals if treatment_from_group(r["group"]) == treatment and r["category"] == category and r["topology"] == "species_tree"]
            if rs:
                cond_plot.append((treatment, category, statistics.mean(float(r["conditional_probability"]) for r in rs)))
    fig, ax = plt.subplots(figsize=(8.5, 4.5))
    categories = ["uniform_4_0", "3_1", "concordant_2_2", "discordant_2_2"]
    for i, category in enumerate(categories):
        for treatment in colors:
            value = next((v for t, c, v in cond_plot if t == treatment and c == category), None)
            if value is not None: ax.scatter(i + list(colors).index(treatment)*.18-.18, value, color=colors[treatment], label=labels[treatment] if i == 0 else None)
    ax.set_xticks(range(4), categories, rotation=25); ax.set_ylabel("P(species-tree quartet | Z category)"); ax.set_title("Conditional topology composition; pilot panel average"); ax.legend(frameon=False, fontsize=8); fig.tight_layout(); fig.savefig(OUT / "figures/conditional_vs_marginal.pdf"); fig.savefig(OUT / "figures/conditional_vs_marginal.png", dpi=220); plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.8, 4.5))
    for treatment in colors:
        rs = [r for r in rows if r["treatment"] == treatment]
        ax.scatter([float(r["gene_tree_normalized_rf"]) for r in rs], [float(r["astral_normalized_rf"]) for r in rs], color=colors[treatment], label=labels[treatment], alpha=.8)
    ax.set_xlabel("gene-tree normalized RF to species tree"); ax.set_ylabel("ASTRAL normalized RF"); ax.set_title("Gene-tree distortion versus ASTRAL accuracy"); ax.legend(frameon=False, fontsize=8); fig.tight_layout(); fig.savefig(OUT / "figures/gene_tree_vs_astral.pdf"); fig.savefig(OUT / "figures/gene_tree_vs_astral.png", dpi=220); plt.close(fig)


if __name__ == "__main__":
    main()
