#!/usr/bin/env python3
"""Summarize checkpointed history-level quartet validation."""
from __future__ import annotations

import csv
import gzip
import json
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
RHOS = (0.0, 0.0001, 0.001, 0.002, 0.01, 0.02)
TOPOLOGIES = ("Q1", "Q2", "Q3")
MULTIPLE_COMPARISON_Z = 4.5


def read_json(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def bootstrap(values: np.ndarray, rng: np.random.Generator, reps: int = 100) -> tuple[float, float]:
    if len(values) < 2:
        return float("nan"), float("nan")
    # History-level normal bootstrap equivalent.  The individual history is
    # the resampling unit; this avoids treating the many quartet rows within a
    # history as independent loci.
    se = float(values.std(ddof=1) / math.sqrt(len(values)))
    return float(values.mean() - 1.96 * se), float(values.mean() + 1.96 * se)


def main() -> None:
    parser = __import__("argparse").ArgumentParser()
    parser.add_argument("--run-dir", type=Path, default=HERE / "validated_100")
    parser.add_argument("--out-dir", type=Path, default=HERE)
    args = parser.parse_args()
    files = sorted(args.run_dir.glob("*/history_*.json.gz"))
    histories = [read_json(path) for path in files]
    if not histories:
        raise SystemExit(f"no history checkpoints under {args.run_dir}")
    rng = np.random.default_rng(20261023)

    groups = defaultdict(list)
    for h in histories:
        groups[(h["condition"], h["replicate"], h["scenario"])].append(h)
    persistence = []
    for key, rows in sorted(groups.items()):
        n = len(rows)
        fields = ("origin_polymorphic_fraction", "origin_end_polymorphic", "ever_fixed", "ever_lost", "terminal_polymorphic")
        item = {"condition": key[0], "replicate": key[1], "scenario": key[2], "n_histories": n}
        for field in fields:
            values = np.array([float(r[field]) for r in rows])
            item[field + "_mean"] = float(values.mean())
            item[field + "_se"] = float(values.std(ddof=1) / math.sqrt(n)) if n > 1 else float("nan")
        item["status"] = "complete" if n == 100 else "partial"
        persistence.append(item)
    write(args.out_dir / "persistence_validation.tsv", persistence)

    switch_groups = defaultdict(list)
    for h in histories:
        for row in h["genealogies"]:
            switch_groups[(h["condition"], h["replicate"], h["scenario"], float(row["rho"]))].append(row)
    switching = []
    for key, rows in sorted(switch_groups.items()):
        switches = np.array([float(r["n_switches"]) for r in rows])
        coals = np.array([float(r["n_coalescences"]) for r in rows])
        switching.append({"condition": key[0], "replicate": key[1], "scenario": key[2], "rho": key[3],
                          "n_history_genealogies": len(rows), "mean_switches": switches.mean(),
                          "switch_se": switches.std(ddof=1) / math.sqrt(len(rows)) if len(rows) > 1 else float("nan"),
                          "fraction_zero_switches": float(np.mean(switches == 0)), "mean_coalescences": coals.mean()})
    write(args.out_dir / "switching_diagnostics.tsv", switching)

    qgroups = defaultdict(list)
    conditional_groups = defaultdict(list)
    for h in histories:
        for row in h["genealogies"]:
            key = (h["condition"], h["replicate"], h["scenario"], float(row["rho"]), row["quartet"])
            qgroups[key].append(row)
            conditional_groups[key + (row["pattern"],)].append(row)

    qrows = []
    ci_rows = []
    reversal_rows = []
    for key, rows in sorted(qgroups.items()):
        condition, replicate, scenario, rho, quartet = key
        species = rows[0]["species_topology"]
        alternatives = [x for x in TOPOLOGIES if x != species]
        indicators = {x: np.array([r["topology"] == x for r in rows], dtype=float) for x in TOPOLOGIES}
        probs = {x: float(indicators[x].mean()) for x in TOPOLOGIES}
        cis = {x: bootstrap(indicators[x], rng) for x in TOPOLOGIES}
        deltas = {x: indicators[species] - indicators[x] for x in alternatives}
        delta_ci = {x: bootstrap(deltas[x], rng) for x in alternatives}
        dominant = max(TOPOLOGIES, key=lambda x: probs[x])
        qrows.append({"condition": condition, "replicate": replicate, "scenario": scenario, "rho": rho, "quartet": quartet,
                      "species_topology": species, "p_Q1": probs["Q1"], "p_Q2": probs["Q2"], "p_Q3": probs["Q3"],
                      "quartet_margin": probs[species] - max(probs[x] for x in alternatives), "dominant_topology": dominant,
                      "n_histories": len(rows)})
        ci_rows.append({"condition": condition, "replicate": replicate, "scenario": scenario, "rho": rho, "quartet": quartet,
                        "species_topology": species, "alternative_1": alternatives[0], "alternative_2": alternatives[1],
                        "delta_species_alt1": probs[species] - probs[alternatives[0]],
                        "delta_species_alt1_low": delta_ci[alternatives[0]][0], "delta_species_alt1_high": delta_ci[alternatives[0]][1],
                        "delta_species_alt2": probs[species] - probs[alternatives[1]],
                        "delta_species_alt2_low": delta_ci[alternatives[1]][0], "delta_species_alt2_high": delta_ci[alternatives[1]][1],
                        "negative_point_margin": int(qrows[-1]["quartet_margin"] < 0), "n_histories": len(rows)})
        supported = []
        corrected_supported = []
        for alt in alternatives:
            other = next(x for x in alternatives if x != alt)
            alt_other = indicators[alt] - indicators[other]
            alt_other_ci = bootstrap(alt_other, rng)
            if delta_ci[alt][1] < 0 and alt_other_ci[0] > 0:
                supported.append(alt)
            # The prespecified many-quartet gate uses a conservative
            # Bonferroni-style z threshold rather than treating one negative
            # row among thousands as validation.
            alt_delta = float(deltas[alt].mean())
            alt_se = float(deltas[alt].std(ddof=1) / math.sqrt(len(rows))) if len(rows) > 1 else float("nan")
            other_delta = float(alt_other.mean())
            other_se = float(alt_other.std(ddof=1) / math.sqrt(len(rows))) if len(rows) > 1 else float("nan")
            if alt_delta + MULTIPLE_COMPARISON_Z * alt_se < 0 and other_delta - MULTIPLE_COMPARISON_Z * other_se > 0:
                corrected_supported.append(alt)
        reversal_rows.append({"condition": condition, "replicate": replicate, "scenario": scenario, "rho": rho,
                              "quartet": quartet, "n_histories": len(rows), "point_margin": qrows[-1]["quartet_margin"],
                              "point_dominant": dominant, "supported_discordant_dominant": int(bool(supported)),
                              "supported_topology": supported[0] if supported else "",
                              "multiple_comparison_supported": int(bool(corrected_supported)),
                              "multiple_comparison_topology": corrected_supported[0] if corrected_supported else ""})
    write(args.out_dir / "quartet_probabilities.tsv", qrows)
    write(args.out_dir / "quartet_margin_confidence_intervals.tsv", ci_rows)
    write(args.out_dir / "quartet_reversal_validation.tsv", reversal_rows)

    conditional = []
    for key, rows in sorted(conditional_groups.items()):
        condition, replicate, scenario, rho, quartet, pattern = key
        species = rows[0]["species_topology"]
        counts = {x: sum(r["topology"] == x for r in rows) for x in TOPOLOGIES}
        total = len(rows)
        conditional.append({"condition": condition, "replicate": replicate, "scenario": scenario, "rho": rho, "quartet": quartet,
                            "pattern": pattern, "species_topology": species, "p_Q1": counts["Q1"] / total,
                            "p_Q2": counts["Q2"] / total, "p_Q3": counts["Q3"] / total, "n_histories": total,
                            "natural_pattern_frequency": total / len(qgroups[(condition, replicate, scenario, rho, quartet)])})
    write(args.out_dir / "conditional_quartet_probabilities.tsv", conditional)

    figures = args.out_dir / "figures"; figures.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for key, rows in sorted(groups.items()):
        label = f"{key[1]} {key[2][0]}"
        x = [float(r["rho"]) for r in switching if (r["condition"], r["replicate"], r["scenario"]) == key]
        y = [float(r["mean_switches"]) for r in switching if (r["condition"], r["replicate"], r["scenario"]) == key]
        if x: axes[0].plot(x, y, marker="o", label=label)
        p = next(r for r in persistence if (r["condition"], r["replicate"], r["scenario"]) == key)
        axes[1].errorbar([label], [p["origin_polymorphic_fraction_mean"]], yerr=[p["origin_polymorphic_fraction_se"]], fmt="o")
    axes[0].set_xscale("symlog", linthresh=1e-4); axes[0].set_xlabel("effective switching intensity rho"); axes[0].set_ylabel("mean switches per genealogy"); axes[0].legend(fontsize=6)
    axes[1].set_ylabel("origin polymorphic fraction"); axes[1].tick_params(axis="x", rotation=90)
    for suffix in ("pdf", "png"):
        fig.savefig(figures / f"validated_persistence_switching.{suffix}", dpi=220)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 5))
    for key in sorted({(r["condition"], r["replicate"], r["scenario"]) for r in qrows}):
        rows = [r for r in qrows if (r["condition"], r["replicate"], r["scenario"]) == key]
        means = []
        for rho in RHOS:
            vals = [float(r["quartet_margin"]) for r in rows if float(r["rho"]) == rho]
            means.append(float(np.mean(vals)) if vals else float("nan"))
        ax.plot(RHOS, means, marker="o", label=f"{key[1]} {key[2][0]}")
    ax.set_xscale("symlog", linthresh=1e-4); ax.axhline(0, color="black", linewidth=.8); ax.set_xlabel("effective switching intensity rho"); ax.set_ylabel("mean quartet margin"); ax.legend(fontsize=6)
    for suffix in ("pdf", "png"):
        fig.savefig(figures / f"validated_quartet_margins.{suffix}", dpi=220)
    plt.close(fig)

    manifest = json.loads((args.run_dir / "run_manifest.json").read_text())
    supported = sum(int(r["multiple_comparison_supported"]) for r in reversal_rows)
    raw_supported = sum(int(r["supported_discordant_dominant"]) for r in reversal_rows)
    report = ["# Validated Mai30 Quartet-Reversal Experiment", "",
              f"The summarized run used simulator commit `{manifest['simulator_commit']}`, {len(histories)} completed history checkpoints, six rho values, and one genealogy per history/rho.",
              "", "No full 30-taxon ASTRAL validation was run: the prespecified gate requires an independently supported discordant dominant quartet.", "",
              f"Uncorrected negative-dominant rows passing the nominal history-level interval rule: **{raw_supported}**; rows passing the conservative multiple-comparison gate (z={MULTIPLE_COMPARISON_Z}): **{supported}**.",
              "", "Intervals use the forward-history unit; the multiple-comparison gate is required before calling a quartet reversal validated. Conditional rows are descriptive and are not used to infer natural pattern frequencies.", ""]
    (args.out_dir / "results_text.md").write_text("\n".join(report))
    (args.out_dir / "report.md").write_text("\n".join(report + ["", "The 1,000-history expansion remains gated on these initial results and computational cost. Raw compressed forward histories are retained in the checkpoint directory and are not included in the committed summary outputs.", ""]))


if __name__ == "__main__":
    main()
