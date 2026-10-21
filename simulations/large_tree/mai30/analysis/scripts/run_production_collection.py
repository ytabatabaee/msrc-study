#!/usr/bin/env python3
"""Run one resumable Mai30 production collection from the manifest."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import yaml

SCRIPT = Path(__file__).resolve()
STUDY_ROOT = SCRIPT.parents[2]
sys.path.insert(0, str(SCRIPT.parent))
from compute_original_baseline import parse_newick  # noqa: E402


def load_row(manifest: Path, task_id: int) -> dict[str, str]:
    with manifest.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    matches = [row for row in rows if int(row["task_id"]) == task_id]
    if len(matches) != 1:
        raise ValueError(f"task_id {task_id} is not unique in {manifest}")
    return matches[0]


def config_for(row: dict[str, str], output: Path, seed: int) -> dict:
    newick = (STUDY_ROOT / row["species_tree_path"]).read_text().strip()
    n = int(row["scaled_effective_population_size"])
    return {
        "mode": "arbitrary_tree",
        "seed": seed,
        "num_loci": int(row["loci_per_history"]),
        "species_tree": {"newick": newick, "root_extension": int(row["scaled_root_extension_generations"]), "default_effective_population_size": n, "ultrametric_tolerance": 1e-6},
        "rearrangement": {"id": "mai_pilot_inv", "type": "inversion", "origin_branch": row["origin_branch"], "origin_time_from_branch_start": int(row["origin_time_from_branch_start"]), "initial_copy_count": int(round(2 * n * float(row["initial_frequency"]))), "selection": {"model": "genic", "coefficient": 0.0}},
        "recombination": {"baseline_rate": 0.02, "effective_cross_arrangement_fraction": float(row["cross_arrangement_fraction"])},
        "genome": {"length_bp": int(row["loci_per_history"]), "rearrangement_interval": {"start": 0, "end": int(row["loci_per_history"]), "id": "mai_pilot_inv"}, "breakpoint_rate_per_bp": 0.0},
        "population_process": {"model": "wright_fisher"},
        "ordinary_msc_outside": True,
        "output": {"directory": str(output), "plots": False},
    }


def validate_history(directory: Path, expected_loci: int) -> dict:
    summary = json.loads((directory / "simulation_summary.json").read_text())
    lines = [line.strip() for line in (directory / "local_gene_trees.nwk").read_text().splitlines() if line.strip()]
    if len(lines) != expected_loci or summary.get("num_gene_trees") != expected_loci or summary.get("num_genealogy_blocks") != expected_loci:
        raise ValueError(f"invalid tree count in {directory}")
    taxa = frozenset(summary["taxa"])
    for line in lines:
        tree = parse_newick(line)
        if tree.taxa != taxa or tree.n_internal != len(taxa) - 1 or tree.root_height < 0:
            raise ValueError(f"invalid tree in {directory}")
    if summary.get("num_coalescent_events") != expected_loci * (len(taxa) - 1):
        raise ValueError(f"incomplete coalescence in {directory}")
    return summary


def run(args: argparse.Namespace) -> None:
    row = load_row(args.manifest, args.task_id)
    output_root = Path(args.output_root or STUDY_ROOT / row["output_root"]).resolve()
    collection = output_root / row["collection_id"]
    run_root = collection / "histories"
    run_root.mkdir(parents=True, exist_ok=True)
    histories = int(args.histories or row["n_histories"])
    loci = int(row["loci_per_history"])
    summaries = []
    for history in range(histories):
        out = run_root / f"history_{history:03d}"
        marker = out / "history_complete.json"
        if marker.exists():
            summaries.append(validate_history(out, loci))
            continue
        out.mkdir(parents=True, exist_ok=True)
        seed = int(row["seed_base"]) + history * 1000
        config = config_for(row, out, seed)
        config_path = out / "resolved_config.yaml"
        config_path.write_text(yaml.safe_dump(config, sort_keys=False))
        log_path = out / "simulation.log"
        command = [sys.executable, "-m", "msrcsim.arbitrary_cli", "--config", str(config_path)]
        completed = subprocess.run(command, cwd=STUDY_ROOT, stdout=log_path.open("w"), stderr=subprocess.STDOUT, check=False)
        if completed.returncode != 0:
            raise RuntimeError(f"simulation failed for task {args.task_id}, history {history}; see {log_path}")
        summary = validate_history(out, loci)
        marker.write_text(json.dumps({"status": "PASS", "task_id": args.task_id, "history_id": history, "seed": seed, "num_gene_trees": loci}, indent=2) + "\n")
        summaries.append(summary)

    combined = collection / "gene_trees.nwk"
    mapping = collection / "locus_history_mapping.tsv"
    with combined.open("w") as gene_handle, mapping.open("w", newline="") as map_handle:
        writer = csv.DictWriter(map_handle, fieldnames=["collection_id", "history_id", "locus_id", "source_tree"], delimiter="\t")
        writer.writeheader()
        locus = 0
        for history in range(histories):
            source = run_root / f"history_{history:03d}" / "local_gene_trees.nwk"
            for line in source.read_text().splitlines():
                if line.strip():
                    gene_handle.write(line.strip() + "\n")
                    writer.writerow({"collection_id": row["collection_id"], "history_id": history, "locus_id": locus, "source_tree": str(source.relative_to(STUDY_ROOT))})
                    locus += 1
    if locus != histories * loci:
        raise ValueError(f"combined collection has {locus} trees, expected {histories * loci}")
    summary = {"status": "PASS", "task_id": args.task_id, "collection_id": row["collection_id"], "n_independent_histories": histories, "n_gene_trees": locus, "n_taxa": summaries[0]["n_taxa"], "total_coalescent_events": sum(s["num_coalescent_events"] for s in summaries), "gene_tree_path": str(combined.relative_to(STUDY_ROOT)), "mapping_path": str(mapping.relative_to(STUDY_ROOT)), "scale_factor": float(row["scale_factor"]), "msrc_diploid_ne": int(row["scaled_effective_population_size"])}
    (collection / "collection_complete.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, sort_keys=True))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=STUDY_ROOT / "datasets/processed/mai30_production_manifest.tsv")
    parser.add_argument("--task-id", type=int, required=True)
    parser.add_argument("--histories", type=int)
    parser.add_argument("--output-root", type=Path)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
