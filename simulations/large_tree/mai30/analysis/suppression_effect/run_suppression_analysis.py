#!/usr/bin/env python3
"""Audit and quantify the Mai30 MSRC suppression effect without changing inputs."""

from __future__ import annotations

import argparse
import csv
import json
import random
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[5]
MAI = ROOT / "simulations/large_tree/mai30"
RESULTS = MAI / "analysis/results"
OUT = MAI / "analysis/suppression_effect"
ASTRAL_DEFAULT = Path("/Users/ytabatabaee/Desktop/ASTER/bin/astral")
TREATMENT_KAPPA = {"msrc_no_suppression": 1.0, "msrc_moderate_suppression": 0.5, "msrc_strong_suppression": 0.1}

import sys
sys.path.insert(0, str(MAI / "analysis/scripts"))
from compute_original_baseline import parse_newick  # noqa: E402
from run_astral_scaled_pilot import normalized_rf, parse_lengths_correct  # noqa: E402


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def group_parts(group: str) -> tuple[str, str]:
    parts = group.split("__")
    condition = parts[0]
    rep = next(p.removeprefix("replicate_") for p in parts if p.startswith("replicate_"))
    return condition, rep


def true_species(condition: str, rep: str):
    return MAI / "datasets/raw" / condition / rep / "s_tree.trees"


def scaled_ne(condition: str, rep: str) -> int:
    config = next((MAI / "datasets/processed/mai30_scaled_wf_runs").glob(f"{condition}__replicate_{rep}__*__freq_0.5/history_000/simulation_config.json"))
    return int(json.loads(config.read_text())["species_tree"]["default_effective_population_size"])


def split_label(split: frozenset[str]) -> str:
    return ",".join(sorted(split))


def audit_rows() -> tuple[list[dict], dict[str, dict]]:
    species_rows = {r["group"]: r for r in read_tsv(RESULTS / "mai30_msrc_vs_species.tsv")}
    astral_rows = {r["group"]: r for r in read_tsv(RESULTS / "mai30_scaled_wf_astral.tsv")}
    groups = sorted(species_rows)
    by_cell = {}
    for group in groups:
        sr = species_rows[group]
        ar = astral_rows[group]
        condition, rep = group_parts(group)
        kappa = TREATMENT_KAPPA[sr["treatment"]]
        by_cell[group] = {"species": sr, "astral": ar, "condition": condition, "replicate": rep, "kappa": kappa}
    output = []
    for group, data in by_cell.items():
        sr, ar = data["species"], data["astral"]
        condition, rep, kappa = data["condition"], data["replicate"], data["kappa"]
        reference = next(v for v in by_cell.values() if v["condition"] == condition and v["replicate"] == rep and v["kappa"] == 1.0)
        row = {
            "condition": condition, "replicate": rep, "treatment": sr["treatment"], "kappa": kappa,
            "n_gene_trees": sr["n_gene_trees"], "n_independent_histories": sr["n_independent_histories"],
            "gene_tree_normalized_rf": float(sr["mean_normalized_rf"]),
            "quartet_concordance": float(sr["mean_sampled_quartet_concordance"]),
            "species_split_recovery": float(sr["mean_species_split_concordance"]),
            "astral_normalized_rf": float(ar["normalized_rf"]),
            "astral_exact_recovery": int(ar["exact_topology_recovery"]),
            "astral_common_split_abs_error": float(ar["mean_absolute_cu_error_common_splits"]),
            "astral_common_split_signed_error": float(ar["mean_cu_bias_common_splits"]),
            "astral_common_split_count": int(ar["common_internal_splits"]),
            "delta_gene_tree_rf_vs_kappa1": float(sr["mean_normalized_rf"]) - float(reference["species"]["mean_normalized_rf"]),
            "delta_quartet_concordance_vs_kappa1": float(sr["mean_sampled_quartet_concordance"]) - float(reference["species"]["mean_sampled_quartet_concordance"]),
            "delta_astral_rf_vs_kappa1": float(ar["normalized_rf"]) - float(reference["astral"]["normalized_rf"]),
            "delta_astral_abs_branch_error_vs_kappa1": float(ar["mean_absolute_cu_error_common_splits"]) - float(reference["astral"]["mean_absolute_cu_error_common_splits"]),
            "delta_astral_signed_branch_error_vs_kappa1": float(ar["mean_cu_bias_common_splits"]) - float(reference["astral"]["mean_cu_bias_common_splits"]),
            "delta_exact_recovery_vs_kappa1": int(ar["exact_topology_recovery"]) - int(reference["astral"]["exact_topology_recovery"]),
            "reference_ordinary_msc": "separate table; not the within-MSRC contrast",
        }
        output.append(row)
    return output, by_cell


def split_rows(by_cell: dict[str, dict]) -> list[dict]:
    rows = []
    groups = list(by_cell)
    groups += [f"{condition}__replicate_{rep}__msc_control" for condition, rep in sorted({(d["condition"], d["replicate"]) for d in by_cell.values()})]
    for group in sorted(groups):
        if "msc_control" in group:
            condition, rep = group_parts(group)
            method = "MSC"
            astral_path = MAI / "datasets/processed/mai30_scaled_msc_controls" / condition / rep / "astral_species_tree.nwk"
        else:
            data = by_cell[group]
            condition, rep = data["condition"], data["replicate"]
            method = data["species"]["treatment"]
            astral_path = Path(data["astral"]["astral_output"])
        true_taxa, true_splits, true_lengths = parse_lengths_correct(true_species(condition, rep).read_text())
        inferred_taxa, inferred_splits, inferred_lengths = parse_lengths_correct(astral_path.read_text())
        for split in sorted(true_splits, key=lambda x: tuple(sorted(x))):
            rows.append({"condition": condition, "replicate": rep, "method": method, "split_type": "true_split", "split": split_label(split), "recovered": int(split in inferred_splits), "estimated_culength": inferred_lengths.get(split, ""), "true_culength": true_lengths[split] / (2 * scaled_ne(condition, rep)), "shared_inferred_split": int(split in inferred_splits)})
        for split in sorted(inferred_splits - true_splits, key=lambda x: tuple(sorted(x))):
            rows.append({"condition": condition, "replicate": rep, "method": method, "split_type": "incorrect_inferred_split", "split": split_label(split), "recovered": 0, "estimated_culength": inferred_lengths.get(split, ""), "true_culength": "", "shared_inferred_split": 0})
    return rows


def branch_rows(by_cell: dict[str, dict]) -> list[dict]:
    out = []
    groups = list(by_cell)
    groups += [f"{condition}__replicate_{rep}__msc_control" for condition, rep in sorted({(d["condition"], d["replicate"]) for d in by_cell.values()})]
    for group in sorted(groups):
        if "msc_control" in group:
            condition, rep = group_parts(group)
            method = "MSC"
            astral_path = MAI / "datasets/processed/mai30_scaled_msc_controls" / condition / rep / "astral_species_tree.nwk"
        else:
            data = by_cell[group]
            condition, rep = data["condition"], data["replicate"]
            method = data["species"]["treatment"]
            astral_path = Path(data["astral"]["astral_output"])
        true_taxa, true_splits, true_lengths = parse_lengths_correct(true_species(condition, rep).read_text())
        inferred_taxa, inferred_splits, inferred_lengths = parse_lengths_correct(astral_path.read_text())
        for split in sorted(true_splits & inferred_splits, key=lambda x: tuple(sorted(x))):
            true_cu = true_lengths[split] / (2 * scaled_ne(condition, rep))
            estimate = inferred_lengths[split]
            out.append({"condition": condition, "replicate": rep, "method": method, "split": split_label(split), "true_culength": true_cu, "estimated_culength": estimate, "signed_error": estimate - true_cu, "absolute_error": abs(estimate - true_cu), "topology_correct": 1})
    return out


def tree_blocks(group: str) -> tuple[list[list[str]], str, str, Path]:
    condition, rep = group_parts(group)
    if "msc_control" in group:
        source = MAI / "datasets/processed/mai30_scaled_msc_controls" / condition / rep / "msc_gene_trees.nwk"
        lines = [x for x in source.read_text().splitlines() if x.strip()]
        return [lines[i:i + 20] for i in range(0, len(lines), 20)], "ordinary MSC pseudo-blocks", condition, true_species(condition, rep)
    source = MAI / "datasets/processed/mai30_scaled_wf_pilot" / group / "gene_trees.nwk"
    mapping = read_tsv(source.parent / "locus_history_mapping.tsv")
    lines = [x for x in source.read_text().splitlines() if x.strip()]
    grouped = defaultdict(list)
    for row, line in zip(mapping, lines):
        grouped[row["history_id"]].append(line)
    return list(grouped.values()), "independent forward histories", condition, true_species(condition, rep)


def bootstrap_astral(by_cell: dict[str, dict], output: Path, astral: Path, reps: int, seed: int) -> None:
    rows = []
    for group in sorted(list(by_cell) + [
        f"{condition}__replicate_{rep}__msc_control" for condition, rep in sorted({(d["condition"], d["replicate"]) for d in by_cell.values()})
    ]):
        if "msc_control" in group:
            condition, rep = group_parts(group.replace("__msc_control", "__msc_no_suppression"))
            method = "MSC"
        else:
            if group not in by_cell:
                continue
            condition, rep = group_parts(group)
            method = by_cell[group]["species"]["treatment"]
        blocks, block_type, condition, truth = tree_blocks(group)
        if len(blocks) < 2:
            continue
        true_taxa, true_splits, true_lengths = parse_lengths_correct(truth.read_text())
        n_scaled = int(json.loads((MAI / "datasets/processed/mai30_scaled_wf_runs" / (group if "msc_control" not in group else next(g for g,d in by_cell.items() if d["condition"] == condition and d["replicate"] == rep)) / "history_000/simulation_config.json").read_text())["species_tree"]["default_effective_population_size"])
        rng = random.Random(seed + len(rows))
        with tempfile.TemporaryDirectory(prefix="mai30_astral_bootstrap_") as tmp:
            tmp = Path(tmp)
            for b in range(reps):
                selected = [rng.choice(blocks) for _ in range(25)]
                gene_file = tmp / f"{b}.nwk"
                output_tree = tmp / f"{b}.out.nwk"
                gene_file.write_text("\n".join(line for block in selected for line in block) + "\n")
                cmd = [str(astral), "-i", str(gene_file), "-o", str(output_tree), "-t", "2", "-u", "1", "--length", "CULength"]
                result = subprocess.run(cmd, capture_output=True, text=True, check=False)
                if result.returncode != 0:
                    raise RuntimeError(f"ASTRAL bootstrap failed for {group} replicate {b}: {result.stderr[-500:]}")
                taxa, splits, lengths = parse_lengths_correct(output_tree.read_text())
                common = true_splits & splits
                errors = [lengths[s] - true_lengths[s] / (2 * n_scaled) for s in common]
                rows.append({"condition": condition, "replicate": rep, "method": method, "bootstrap": b, "block_type": block_type, "n_blocks": len(blocks), "normalized_rf": normalized_rf(splits, true_splits, len(true_taxa)), "exact_recovery": int(splits == true_splits), "common_splits": len(common), "mean_abs_branch_error": float(np.mean(np.abs(errors))) if errors else np.nan, "mean_signed_branch_error": float(np.mean(errors)) if errors else np.nan})
    write_tsv(output / "astral_bootstrap.tsv", rows)


def make_figures(out: Path, effects: list[dict], splits: list[dict], branches: list[dict], bootstrap: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(parents=True, exist_ok=True)
    methods = ["MSC", "msrc_no_suppression", "msrc_moderate_suppression", "msrc_strong_suppression"]
    colors = {"MSC": "#444444", "msrc_no_suppression": "#377eb8", "msrc_moderate_suppression": "#ff8c42", "msrc_strong_suppression": "#d73027"}
    for ext in ("pdf", "png"):
        fig, ax = plt.subplots(figsize=(9, 5))
        for i, rep in enumerate(["042", "005", "074", "007"]):
            condition = next(r["condition"] for r in effects if r["replicate"] == rep)
            for j, method in enumerate(methods):
                if method == "MSC":
                    candidates = [r for r in bootstrap if r["replicate"] == rep and r["method"] == method]
                    value = np.mean([float(r["normalized_rf"]) for r in candidates]) if candidates else np.nan
                else:
                    candidates = [r for r in effects if r["replicate"] == rep and r["treatment"] == method]
                    value = candidates[0]["astral_normalized_rf"] if candidates else np.nan
                ax.scatter(i + (j - 1.5) * .16, value, color=colors[method], s=38)
        ax.set_xticks(range(4), ["042", "005", "074", "007"]); ax.set_ylabel("ASTRAL normalized RF"); ax.set_title("ASTRAL topology accuracy"); ax.grid(axis="y", alpha=.25)
        for method in methods: ax.scatter([], [], color=colors[method], label=method)
        ax.legend(frameon=False, fontsize=8); fig.tight_layout(); fig.savefig(out / f"astral_topology_accuracy.{ext}", dpi=220, bbox_inches="tight"); plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 4.5))
        for rep in ["042", "005", "074", "007"]:
            xs = {r["treatment"]: float(r["delta_astral_rf_vs_kappa1"]) for r in effects if r["replicate"] == rep}
            ax.scatter([rep, rep], [xs.get("msrc_moderate_suppression", np.nan), xs.get("msrc_strong_suppression", np.nan)], color=["#ff8c42", "#d73027"], s=38)
        ax.axhline(0, color="black", lw=.7); ax.set_ylabel("ASTRAL RF difference from κ=1"); ax.set_xlabel("Mai replicate"); ax.set_title("Within-MSRC suppression effect"); fig.tight_layout(); fig.savefig(out / f"suppression_effect.{ext}", dpi=220, bbox_inches="tight"); plt.close(fig)

        rep_splits = [r for r in splits if r["replicate"] == "007"]
        methods_split = methods
        split_names = sorted({r["split"] for r in rep_splits})
        matrix = np.full((len(methods_split), len(split_names)), np.nan)
        for i, method in enumerate(methods_split):
            for j, split in enumerate(split_names):
                rows = [r for r in rep_splits if r["method"] == method and r["split"] == split and r["split_type"] == "true_split"]
                if rows: matrix[i, j] = int(rows[0]["recovered"])
        fig, ax = plt.subplots(figsize=(14, 3.8)); ax.imshow(matrix, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1); ax.set_yticks(range(len(methods_split)), methods_split); ax.set_xticks(range(len(split_names)), [f"S{i+1}" for i in range(len(split_names))], rotation=90); ax.set_title("Replicate 007: true species splits recovered by ASTRAL"); fig.tight_layout(); fig.savefig(out / f"replicate_007_split_recovery.{ext}", dpi=220, bbox_inches="tight"); plt.close(fig)

        fig, ax = plt.subplots(figsize=(9, 5))
        branch_methods = ["MSC", "msrc_no_suppression", "msrc_moderate_suppression", "msrc_strong_suppression"]
        branch_labels = ["MSC", "κ=1", "κ=0.5", "κ=0.1"]
        branch_colors = ["#444444", colors["msrc_no_suppression"], colors["msrc_moderate_suppression"], colors["msrc_strong_suppression"]]
        for j, method in enumerate(branch_methods):
            vals = [float(r["absolute_error"]) for r in branches if r["method"] == method]
            ax.boxplot(vals, positions=[j], widths=.5, patch_artist=True, boxprops={"facecolor": branch_colors[j]})
        ax.set_xticks(range(4), branch_labels); ax.set_ylabel("absolute CULength error on matching splits"); ax.set_title("ASTRAL branch-length error"); ax.grid(axis="y", alpha=.25); fig.tight_layout(); fig.savefig(out / f"branch_length_error.{ext}", dpi=220, bbox_inches="tight"); plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--bootstrap-reps", type=int, default=5)
    p.add_argument("--bootstrap-seed", type=int, default=20261130)
    p.add_argument("--astral", type=Path, default=ASTRAL_DEFAULT)
    p.add_argument("--skip-bootstrap", action="store_true")
    args = p.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    effects, by_cell = audit_rows()
    write_tsv(OUT / "suppression_effects.tsv", effects)
    write_tsv(OUT / "split_recovery.tsv", split_rows(by_cell))
    write_tsv(OUT / "branch_length_comparison.tsv", branch_rows(by_cell))
    if not args.skip_bootstrap:
        bootstrap_astral(by_cell, OUT, args.astral, args.bootstrap_reps, args.bootstrap_seed)
    bootstrap = read_tsv(OUT / "astral_bootstrap.tsv")
    make_figures(OUT / "figures", effects, read_tsv(OUT / "split_recovery.tsv"), read_tsv(OUT / "branch_length_comparison.tsv"), bootstrap)
    (OUT / "analysis_manifest.json").write_text(json.dumps({"bootstrap_replicates": args.bootstrap_reps, "bootstrap_seed": args.bootstrap_seed, "independent_unit": "MSRC forward frequency history; MSC sequential pseudo-blocks", "astral_settings": "-t 2 -u 1 --length CULength"}, indent=2) + "\n")


if __name__ == "__main__":
    main()
