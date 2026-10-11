#!/usr/bin/env python3
"""Run and combine the resumable scaled-WF MSRC pilot histories.

Each independent history produces 20 loci.  Twenty-five histories are combined
to 500 loci, while the mapping remains explicit so those loci are not treated
as 500 independent rearrangement histories.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import yaml

SCRIPT = Path(__file__).resolve()
STUDY_ROOT = SCRIPT.parents[2]
CONFIG_ROOT = STUDY_ROOT / "datasets/configs/mai30_msrc_pilot"
RUN_ROOT = STUDY_ROOT / "datasets/processed/mai30_scaled_wf_runs"
COMBINED_ROOT = STUDY_ROOT / "datasets/processed/mai30_scaled_wf_pilot"
sys.path.insert(0, str(SCRIPT.parent))
from compute_original_baseline import parse_newick  # noqa: E402


def config_jobs(frequency: str, histories: int) -> list[dict]:
    jobs = []
    for template in sorted(CONFIG_ROOT.glob(f"*__{frequency}.yaml")):
        with template.open() as handle:
            config = yaml.safe_load(handle)
        group = template.stem
        for history_index in range(histories):
            jobs.append({
                "template": template,
                "config": config,
                "group": group,
                "history_index": history_index,
                "seed": int(config["seed"]) + history_index * 1000,
            })
    return jobs


def run_one(job: dict) -> dict:
    group = job["group"]
    history_index = job["history_index"]
    run_dir = RUN_ROOT / group / f"history_{history_index:03d}"
    tree_file = run_dir / "local_gene_trees.nwk"
    summary_file = run_dir / "simulation_summary.json"
    run_dir.mkdir(parents=True, exist_ok=True)
    if tree_file.exists() and summary_file.exists():
        summary = json.loads(summary_file.read_text())
        if summary.get("num_gene_trees") == 20:
            return {"group": group, "history_index": history_index, "status": "reused", **summary}

    config = json.loads(json.dumps(job["config"]))
    config["seed"] = job["seed"]
    config["output"]["directory"] = str(run_dir)
    config["output"]["plots"] = False
    config_path = run_dir / "resolved_config.yaml"
    config_path.write_text(yaml.safe_dump(config, sort_keys=False))
    log_path = run_dir / "simulation.log"
    command = [sys.executable, "-m", "msrcsim.arbitrary_cli", "--config", str(config_path)]
    with log_path.open("w") as log:
        completed = subprocess.run(command, cwd=STUDY_ROOT, stdout=log, stderr=subprocess.STDOUT, check=False)
    if completed.returncode != 0:
        raise RuntimeError(f"failed {group}/history_{history_index:03d}; see {log_path}")
    summary = json.loads(summary_file.read_text())
    return {"group": group, "history_index": history_index, "status": "ran", **summary}


def validate_history(group: str, history_index: int) -> dict:
    run_dir = RUN_ROOT / group / f"history_{history_index:03d}"
    summary = json.loads((run_dir / "simulation_summary.json").read_text())
    gene_path = run_dir / "local_gene_trees.nwk"
    lines = [line.strip() for line in gene_path.read_text().splitlines() if line.strip()]
    if len(lines) != 20 or summary.get("num_gene_trees") != 20 or summary.get("num_genealogy_blocks") != 20:
        raise ValueError(f"wrong block count in {run_dir}")
    expected_taxa = frozenset(summary["taxa"])
    for index, line in enumerate(lines, start=1):
        tree = parse_newick(line)
        if tree.taxa != expected_taxa:
            raise ValueError(f"taxon set mismatch in {gene_path}:{index}")
        if tree.n_internal != len(expected_taxa) - 1:
            raise ValueError(f"non-bifurcating tree in {gene_path}:{index}")
        if tree.root_height < 0:
            raise ValueError(f"invalid branch length in {gene_path}:{index}")
    if summary.get("num_coalescent_events") != 20 * (len(expected_taxa) - 1):
        raise ValueError(f"incomplete coalescence count in {run_dir}")
    return summary


def combine_group(group: str, histories: int) -> dict:
    summaries = [validate_history(group, i) for i in range(histories)]
    out = COMBINED_ROOT / group
    out.mkdir(parents=True, exist_ok=True)
    combined = out / "gene_trees.nwk"
    mapping_path = out / "locus_history_mapping.tsv"
    with combined.open("w") as gene_handle, mapping_path.open("w", newline="") as map_handle:
        writer = csv.DictWriter(map_handle, fieldnames=["group", "history_id", "locus_id", "source_tree"], delimiter="\t")
        writer.writeheader()
        locus_id = 0
        for history_index in range(histories):
            source = RUN_ROOT / group / f"history_{history_index:03d}" / "local_gene_trees.nwk"
            for within_history, line in enumerate(source.read_text().splitlines()):
                if not line.strip():
                    continue
                gene_handle.write(line.strip() + "\n")
                writer.writerow({
                    "group": group,
                    "history_id": f"{group}__history_{history_index:03d}",
                    "locus_id": locus_id,
                    "source_tree": str(source),
                })
                locus_id += 1
    if locus_id != histories * 20:
        raise ValueError(f"combined collection has {locus_id} trees, expected {histories * 20}: {group}")
    summary = {
        "group": group,
        "n_independent_histories": histories,
        "n_gene_trees": locus_id,
        "n_taxa": summaries[0]["n_taxa"],
        "total_coalescent_events": sum(s["num_coalescent_events"] for s in summaries),
        "history_runtime_seconds": sum(s["runtime_seconds"] for s in summaries),
        "history_peak_rss_bytes_max": max(s["peak_rss_bytes"] for s in summaries),
        "output_contract": "PASS",
        "gene_tree_path": str(combined),
        "mapping_path": str(mapping_path),
    }
    (out / "collection_summary.json").write_text(json.dumps(summary, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frequency", default="freq_0.5", choices=["freq_0.5", "freq_0.2"])
    parser.add_argument("--histories", type=int, default=25)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    jobs = config_jobs(args.frequency, args.histories)
    if not jobs:
        raise SystemExit(f"no configs found for {args.frequency}")
    RUN_ROOT.mkdir(parents=True, exist_ok=True)
    results = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(run_one, job) for job in jobs]
        for future in as_completed(futures):
            result = future.result()
            results.append(result)
            print(f"{result['group']} history {result['history_index']:03d}: {result['status']}", flush=True)
    groups = sorted({result["group"] for result in results})
    collections = [combine_group(group, args.histories) for group in groups]
    manifest = COMBINED_ROOT / f"pilot_manifest_{args.frequency}.tsv"
    with manifest.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(collections[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(collections)
    print(f"wrote {manifest}")


if __name__ == "__main__":
    main()
