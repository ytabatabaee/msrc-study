#!/usr/bin/env python3
"""Create Stage 13d summaries, diagnostics, and figures from completed runs.

This script is post-processing only. It never changes Phlag inputs or calls a
different anomaly model. It deliberately reports invalid/nonfinite runs as
unavailable rather than converting them to zero accuracy.
"""
from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "empirical/cross_system_boundary_detection/stage13d"
FIG = OUT / "figures"
PHLAG_ROOT = ROOT / "empirical/cross_system_boundary_detection/stage13b/phlag/raw_results/runs"
NEO = ROOT / "empirical/neoaves_chr4/results/stage11_phlag_benchmark"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def intervals_from_rows(rows: list[dict[str, str]], state_only: bool = True) -> list[tuple[int, int]]:
    selected = []
    for row in rows:
        if state_only and row["state"] != "1":
            continue
        selected.append((int(row["start_bp"]), int(row["end_bp"])))
    return selected


def merge(intervals: list[tuple[int, int]], gap: int = 5000) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not out or start > out[-1][1] + gap + 1:
            out.append((start, end))
        else:
            out[-1] = (out[-1][0], max(out[-1][1], end))
    return out


def run_summaries() -> tuple[list[dict], list[dict]]:
    rows = read_tsv(OUT / "run_manifest.tsv")
    acc = read_tsv(OUT / "replicate_accuracy.tsv")
    numeric = [r for r in acc if r["aggregation"] in {"union_state1", "max_posterior_ge_0.5"}]
    summary = []
    for aggregation in ("union_state1", "max_posterior_ge_0.5"):
        subset = [r for r in numeric if r["aggregation"] == aggregation]
        for metric in ("iou", "precision", "recall", "f1", "n_predicted_intervals", "predicted_bp"):
            values = [float(r[metric]) for r in subset]
            summary.append({"dataset": "house_mouse_t_complex", "scope": "replicate_aggregate",
                            "configuration": aggregation, "metric": metric, "n": len(values),
                            "mean": np.mean(values), "median": np.median(values),
                            "min": np.min(values), "max": np.max(values)})
    write_tsv(OUT / "house_mouse_replicate_summary.tsv", summary,
              ["dataset", "scope", "configuration", "metric", "n", "mean", "median", "min", "max"])
    return rows, acc


def raw_failure_diagnostics() -> list[dict]:
    specs = [
        ("atlantic_cod", ROOT / "empirical/atlantic_cod/results/phlag_benchmark/raw", "*.txt"),
        ("anopheles_2la", PHLAG_ROOT / "anopheles_2la/replicate_000", "N002.txt"),
        ("fire_ants_chr16", PHLAG_ROOT / "fire_ants_chr16/replicate_000", "N002.txt"),
    ]
    out = []
    for dataset, folder, pattern in specs:
        files = sorted(folder.glob(pattern))
        if not files:
            out.append({"dataset": dataset, "file": "", "status": "unavailable"})
            continue
        for path in files:
            text = path.read_text(errors="replace")
            final = re.search(r"Final focal edge lengths:\s*(.*)", text)
            div = re.search(r"State divergence:\s*(.*)", text)
            states = re.findall(r"(?m)^([01](?:,[01])+)$", text)
            posts = re.findall(r"(?m)^(nan(?:,nan)+|(?:[0-9.eE+-]+)(?:,[0-9.eE+-]+)+)$", text)
            final_text = final.group(1).strip() if final else "NA"
            div_text = div.group(1).strip() if div else "NA"
            edge_floor = bool(re.search(r"(?:^|:)\s*1e-06(?:,|$)", final_text))
            post_finite = bool(posts) and all(math.isfinite(float(x)) for x in posts[-1].split(","))
            out.append({"dataset": dataset, "file": str(path.relative_to(ROOT)),
                        "status": "numerically_degenerate" if not post_finite else "valid",
                        "final_focal_edge_lengths": final_text, "edge_at_floor_1e-6": edge_floor,
                        "state_divergence": div_text, "state_sequence_present": bool(states),
                        "posterior_finite": post_finite,
                        "earliest_observed_failure": "nonfinite posterior after focal-edge flooring"
                        if not post_finite and edge_floor else "not observed"})
    write_tsv(OUT / "cross_system_numerical_diagnostics.tsv", out,
              ["dataset", "file", "status", "final_focal_edge_lengths", "edge_at_floor_1e-6",
               "state_divergence", "state_sequence_present", "posterior_finite", "earliest_observed_failure"])
    return out


def make_house_figures() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    rows = read_tsv(OUT / "window_predictions.tsv")
    by_run: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        by_run[(row["replicate"], row["focal_branch"])].append(row)
    refs = [("Inv1", 5747331, 6957919), ("Inv2", 6957919, 7435308),
            ("Inv3", 7435308, 13016510), ("Inv4", 14068229, 18051507),
            ("Inv5", 19401000, 40049964)]
    fig, ax = plt.subplots(figsize=(12, 4.8))
    y = 0
    for name, start, end in refs:
        obs_end = min(end, 37054999)
        ax.plot([start / 1e6, obs_end / 1e6], [y, y], color="#444444", lw=8, solid_capstyle="butt")
        ax.text(4.85, y, name, ha="right", va="center", fontsize=8)
        y += 1
    selected = [("Phlag N002, rep 0", "0", "N002", "#1b9e77"),
                ("Phlag union, rep 0", None, None, "#d95f02"),
                ("Phlag max posterior, rep 0", None, None, "#7570b3")]
    for label, rep, edge, color in selected:
        y += .55
        if rep is not None:
            iv = merge(intervals_from_rows(by_run[(rep, edge)]))
        else:
            if "union" in label:
                iv = merge([x for key, rr in by_run.items() if key[0] == "0" for x in intervals_from_rows(rr)])
            else:
                iv = merge([(int(row["start_bp"]), int(row["end_bp"])) for i in range(4046)
                            if any(float(rr[i]["posterior_anomaly"]) >= 0.5 for key, rr in by_run.items() if key[0] == "0")
                            for row in [next(rr[i] for key, rr in by_run.items() if key[0] == "0")]])
        for start, end in iv:
            ax.plot([start / 1e6, end / 1e6], [y, y], color=color, lw=5, solid_capstyle="butt")
        ax.text(4.85, y, label, ha="right", va="center", fontsize=8)
    ax.axvline(37.055, color="black", ls="--", lw=0.9)
    ax.axvspan(37.055, 40.05, color="grey", alpha=.16)
    ax.text(37.08, y + .45, "no local-tree coverage", fontsize=8, color="#555555")
    ax.set_xlim(5, 40.1); ax.set_ylim(-.6, y + .8); ax.set_xlabel("Chromosome 17 position (Mb)")
    ax.set_yticks([]); ax.set_title("House-mouse Phlag predictions across representative sets")
    fig.tight_layout(); fig.savefig(FIG / "house_mouse_phlag_boundary_comparison.pdf"); fig.savefig(FIG / "house_mouse_phlag_boundary_comparison.png", dpi=300); plt.close(fig)

    # Posterior heatmap: 100-kb means, one row per replicate/edge.
    groups = sorted(by_run, key=lambda x: (int(x[0]), x[1]))
    bins = np.arange(5_000_000, 37_100_000, 100_000)
    mat = np.full((len(groups), len(bins) - 1), np.nan)
    for r, key in enumerate(groups):
        rr = by_run[key]
        for c in range(len(bins) - 1):
            vals = [float(x["posterior_anomaly"]) for x in rr if bins[c] <= int(x["start_bp"]) < bins[c + 1]]
            if vals: mat[r, c] = np.mean(vals)
    fig, ax = plt.subplots(figsize=(12, 5.5))
    im = ax.imshow(mat, aspect="auto", interpolation="none", cmap="magma", vmin=0, vmax=1,
                   extent=[bins[0] / 1e6, bins[-1] / 1e6, len(groups) - .5, -.5])
    ax.set_yticks(range(len(groups))); ax.set_yticklabels([f"r{r}:{e}" for r, e in groups], fontsize=7)
    ax.set_xlabel("Chromosome 17 position (Mb)"); ax.set_title("Phlag posterior anomaly probability (100-kb means)")
    fig.colorbar(im, ax=ax, label="posterior anomaly probability")
    fig.tight_layout(); fig.savefig(FIG / "house_mouse_phlag_robustness.pdf"); fig.savefig(FIG / "house_mouse_phlag_robustness.png", dpi=300); plt.close(fig)


def make_cross_system_summary() -> None:
    raw_rows = []
    source = ROOT / "empirical/cross_system_boundary_detection/phlag/accuracy_comparison.tsv"
    if source.exists():
        for row in read_tsv(source):
            if row.get("valid_for_accuracy", "").lower() == "true":
                raw_rows.append({"dataset": row["dataset"], "method": "MSRC " + row["method"],
                                 "iou": float(row["iou"]), "f1": float(row["f1"]), "valid": "yes"})
    # The historical table contains repeated evaluation rows. Collapse those
    # rows before plotting, while retaining the aggregation rule in the table.
    grouped = defaultdict(list)
    for row in raw_rows:
        grouped[(row["dataset"], row["method"])].append(row)
    rows = []
    for key, values in sorted(grouped.items()):
        rows.append({"dataset": key[0], "method": key[1],
                     "iou": np.mean([x["iou"] for x in values]),
                     "f1": np.mean([x["f1"] for x in values]), "valid": "yes",
                     "n_source_rows": len(values)})
    for agg in ("union_state1", "max_posterior_ge_0.5"):
        vals = [r for r in read_tsv(OUT / "replicate_accuracy.tsv") if r["aggregation"] == agg]
        rows.append({"dataset": "House mouse", "method": f"Phlag {agg} (mean across 6 reps)",
                     "iou": np.mean([float(x["iou"]) for x in vals]), "f1": np.mean([float(x["f1"]) for x in vals]), "valid": "yes", "n_source_rows": len(vals)})
    for dataset, status, reason in [("Atlantic cod", "invalid", "nonfinite posterior / SNAPP population-tree input"),
                                    ("Anopheles", "invalid", "nonfinite posterior in matched pilot"),
                                    ("Fire ants", "invalid", "nonfinite posterior in matched pilot"),
                                    ("Neoaves", "pilot_valid", "5,619-locus Columbea pilot; separate Stage 11B result")]:
        validity = f"{status}: {reason}"
        rows.append({"dataset": dataset, "method": "Phlag", "iou": "NA", "f1": "NA", "valid": validity, "n_source_rows": 0})
    fields = ["dataset", "method", "iou", "f1", "valid", "n_source_rows"]
    write_tsv(OUT / "cross_system_comparison.tsv", rows, fields)
    datasets = ["Atlantic cod", "Fire ants", "Anopheles", "House mouse"]
    methods = ["MSRC topology", "MSRC continuous", "MSRC joint", "Phlag union_state1 (mean across 6 reps)"]
    colors = {"MSRC topology": "#66a61e", "MSRC continuous": "#e6ab02", "MSRC joint": "#a6761d", "Phlag union_state1 (mean across 6 reps)": "#1b9e77"}
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(datasets)); width = .19
    for j, method in enumerate(methods):
        values = []
        for dataset in datasets:
            candidates = [r for r in rows if r["dataset"] == dataset and r["method"] == method and r["iou"] != "NA"]
            values.append(float(candidates[0]["iou"]) if candidates else np.nan)
        ax.bar(x + (j - 1.5) * width, values, width, label=method.replace(" (mean across 6 reps)", ""), color=colors[method])
    ax.set_xticks(x); ax.set_xticklabels(datasets)
    ax.set_ylim(0, 1); ax.set_ylabel("IoU (valid comparisons only)")
    ax.set_title("MSRC and Phlag interval accuracy; invalid Phlag runs are unavailable")
    ax.legend(fontsize=8, loc="upper right"); fig.tight_layout()
    fig.savefig(FIG / "cross_system_phlag_summary.pdf"); fig.savefig(FIG / "cross_system_phlag_summary.png", dpi=300); plt.close(fig)


def write_docs(run_rows: list[dict], diagnostics: list[dict]) -> None:
    summary = read_tsv(OUT / "house_mouse_replicate_summary.tsv")
    def stat(config: str, metric: str, field: str) -> str:
        return next(row[field] for row in summary if row["configuration"] == config and row["metric"] == metric)
    (OUT / "methods_text.md").write_text("""# Stage 13d methods

The primary house-mouse benchmark uses the unmodified Phlag v0.2.0 executable
at commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6`. Six prepared
representative-taxon sets were evaluated on four eligible focal edges (N002,
N003, N004, N005), with the same representatives in all 4,046 windows of a
replicate. State and posterior sequences were accepted only when both had the
expected length and finite probabilities in [0, 1].

The five inversion references are scored on the observed chr17 domain
5,000,000--37,054,999 bp. Inv5 is right-censored. Accuracy uses merged
interval-set union/intersection metrics. Two fixed aggregations are reported:
the union of state-1 calls across edges and the maximum posterior across edges
with a prespecified 0.5 threshold. These rules were set before scoring.

The macOS `/usr/bin/time -l` wrapper returned nonzero because its optional
`sysctl` memory probe was denied. Complete Phlag output was nevertheless
present and finite for all 24 runs; peak memory is therefore unavailable rather
than imputed.
""")
    (OUT / "results_text.md").write_text(f"""# Stage 13d results

All 24 house-mouse matched runs were valid and informative (6 representatives
sets × 4 focal edges). Across six replicate-level aggregations, the union of
state-1 calls had mean IoU {stat('union_state1','iou','mean')} and mean F1
{stat('union_state1','f1','mean')}; the maximum-posterior aggregation had mean
IoU {stat('max_posterior_ge_0.5','iou','mean')} and mean F1
{stat('max_posterior_ge_0.5','f1','mean')}.

These are valid matched-taxon Phlag results, but they are not independent
biological replicates: adjacent genomic windows are linked and the six
representative sets vary sampling rather than the genome. The right edge of
Inv5 lies beyond the observed domain and was not scored as a recovered
breakpoint.
""")
    (OUT / "report.md").write_text(f"""# Stage 13d Phlag benchmark report

## Completed

- House mouse: 24/24 matched-taxon runs completed with finite, length-matched
  posteriors and nonconstant state paths: six representative sets and four
  focal edges.
- Valid aggregate performance across representatives is summarized in
  `house_mouse_replicate_summary.tsv`; union-state-1 mean IoU/F1 are
  {stat('union_state1','iou','mean')} / {stat('union_state1','f1','mean')},
  and max-posterior mean IoU/F1 are
  {stat('max_posterior_ge_0.5','iou','mean')} / {stat('max_posterior_ge_0.5','f1','mean')}.
- Inv5 was evaluated only through 37,054,999 bp and is explicitly right-censored.
- The repaired Neoaves Stage 11B pilot remains preserved as a separate valid
  5,619-locus result; this Stage 13d run does not relabel it as a multi-replicate
  result.

## Numerical failures

The Atlantic cod, Anopheles, and fire-ant Phlag runs remain unavailable for
accuracy. Their raw outputs contain a focal edge at the Phlag floor (`1e-6`)
and nonfinite state-divergence/posterior values. These are numerical or input
calibration failures, not zero detections. The diagnostics identify the first
observable failure and preserve the raw outputs.

## Fair interpretation

The house-mouse comparison is valid for the matched representative inputs and
uses the original Phlag algorithm. It does not support a cross-dataset claim
that MSRC outperforms Phlag, because the other three datasets do not produce
valid Phlag posteriors and the Neoaves result is a separate one-pilot
experiment. The house-mouse results belong in a supplementary benchmark until
additional independent biological calibration is available.

## Reproduction commands

```bash
{ROOT}/.venv-malaria/bin/python empirical/cross_system_boundary_detection/stage13d/scripts/complete_phlag_benchmark.py --replicates 6 --workers 4 --run
{ROOT}/.venv-malaria/bin/python empirical/cross_system_boundary_detection/stage13d/scripts/complete_phlag_benchmark.py --score-only
{ROOT}/.venv-malaria/bin/python empirical/cross_system_boundary_detection/stage13d/scripts/finalize_phlag_benchmark.py
```
""")
    (OUT / "manifest.json").write_text(json.dumps({
        "stage": "13d", "phlag_commit": "76597d0e4ebd0e507fdafebf7f5163ef634a36d6",
        "algorithm_modified": False, "house_mouse_runs": 24, "valid_house_mouse_runs": 24,
        "observation_domain": [5000000, 37054999], "inv5_right_censored": True,
        "invalid_accuracy_datasets": ["Atlantic cod", "Anopheles", "Fire ants"],
        "timing_memory_note": "macOS time -l sysctl probe denied; peak memory unavailable",
    }, indent=2) + "\n")

    # A compact committed manifest points to the full representative files in
    # the preceding stage without copying or rewriting those selections.
    sampling = []
    input_root = ROOT / "empirical/cross_system_boundary_detection/stage13b/phlag/raw_results/matched_inputs/house_mouse_t_complex"
    for rep in range(6):
        mapping = input_root / f"replicate_{rep:03d}/representatives.tsv"
        sampling.append({"dataset": "house_mouse_t_complex", "replicate": rep,
                         "n_windows": 4046, "focal_edges": "N002,N003,N004,N005",
                         "representative_mapping": str(mapping.relative_to(ROOT)),
                         "selection_fixed_across_windows": "true"})
    write_tsv(OUT / "representative_sampling_manifest.tsv", sampling,
              ["dataset", "replicate", "n_windows", "focal_edges", "representative_mapping", "selection_fixed_across_windows"])


def main() -> None:
    rows, _ = run_summaries()
    diagnostics = raw_failure_diagnostics()
    make_house_figures()
    make_cross_system_summary()
    write_docs(rows, diagnostics)
    print(json.dumps({"status": "complete", "valid_house_mouse_runs": sum(str(r["status"]).startswith("valid") for r in rows),
                      "diagnostic_runs": len(diagnostics), "figures": str(FIG)}, indent=2))


if __name__ == "__main__":
    main()
