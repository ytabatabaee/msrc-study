#!/usr/bin/env python3
"""Build epsilon mixtures from matched complete MSC/MSRC source pools."""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
POOL_ROOT = ROOT / "datasets/processed/mai30_scaled_wf_pilot"
MSC_ROOT = ROOT / "datasets/processed/mai30_scaled_msc_controls"
OUT_ROOT = ROOT / "datasets/processed/mai30_epsilon_mixtures"
EPSILONS = (0.0, 0.1, 0.25, 0.5, 0.75, 1.0)
TREATMENTS = ("msrc_no_suppression", "msrc_moderate_suppression", "msrc_strong_suppression")


def selected_complete_groups() -> list[tuple[str, str]]:
    out = []
    for rep in ("007", "042"):
        condition = "outgroup.0.species.1.5.genes.1.5"
        for treatment in TREATMENTS:
            group = f"{condition}__replicate_{rep}__{treatment}__freq_0.5"
            summary = POOL_ROOT / group / "collection_summary.json"
            if summary.exists() and json.loads(summary.read_text())["n_gene_trees"] == 500:
                out.append((rep, group))
    return out


def source_history_map(msrc_group: str) -> list[str]:
    path = POOL_ROOT / msrc_group / "locus_history_mapping.tsv"
    with path.open() as handle:
        return [row["history_id"] for row in csv.DictReader(handle, delimiter="\t")]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=3, help="mixture resamples per epsilon")
    parser.add_argument("--seed", type=int, default=20261011)
    args = parser.parse_args()
    manifest_rows = []
    for rep, msrc_group in selected_complete_groups():
        condition = msrc_group.split("__replicate_")[0]
        msc_path = MSC_ROOT / condition / rep / "msc_gene_trees.nwk"
        msrc_path = POOL_ROOT / msrc_group / "gene_trees.nwk"
        msc_trees = [x for x in msc_path.read_text().splitlines() if x.strip()]
        msrc_trees = [x for x in msrc_path.read_text().splitlines() if x.strip()]
        histories = source_history_map(msrc_group)
        if len(msc_trees) != 500 or len(msrc_trees) != 500:
            raise ValueError(f"source pool is not complete: {msc_group}")
        for epsilon in EPSILONS:
            n_msrc = int(math.floor(500 * epsilon + 0.5))
            n_msc = 500 - n_msrc
            for mixture_replicate in range(args.replicates):
                seed = args.seed + int(round(epsilon * 1000)) * 10000 + mixture_replicate + (0 if rep == "007" else 1000000) + TREATMENTS.index(msrc_group.split("__")[2]) * 100000
                rng = random.Random(seed)
                msrc_indices = rng.sample(range(500), n_msrc)
                msc_indices = rng.sample(range(500), n_msc)
                records = [("MSRC", i, msrc_trees[i], histories[i]) for i in msrc_indices]
                records += [("MSC", i, msc_trees[i], "MSC_direct") for i in msc_indices]
                rng.shuffle(records)
                name = f"replicate_{rep}__{msrc_group.split('__')[2]}__epsilon_{epsilon:g}__mixture_{mixture_replicate:02d}"
                out = OUT_ROOT / name
                out.mkdir(parents=True, exist_ok=True)
                with (out / "gene_trees.nwk").open("w") as tree_handle, (out / "manifest.tsv").open("w", newline="") as manifest_handle:
                    writer = csv.DictWriter(manifest_handle, fieldnames=["destination_locus", "source_type", "source_index", "source_pool", "history_id", "epsilon", "seed"], delimiter="\t")
                    writer.writeheader()
                    for locus, (source_type, source_index, tree, history_id) in enumerate(records):
                        tree_handle.write(tree + "\n")
                        source_pool = str(msrc_path if source_type == "MSRC" else msc_path)
                        writer.writerow({"destination_locus": locus, "source_type": source_type, "source_index": source_index, "source_pool": source_pool, "history_id": history_id, "epsilon": epsilon, "seed": seed})
                        manifest_rows.append({"mixture": name, "destination_locus": locus, "source_type": source_type, "source_index": source_index, "source_pool": source_pool, "history_id": history_id, "epsilon": epsilon, "seed": seed})
                (out / "mixture_config.json").write_text(json.dumps({"replicate": rep, "condition": condition, "msrc_group": msrc_group, "epsilon": epsilon, "n_msrc": n_msrc, "n_msc": n_msc, "mixture_replicate": mixture_replicate, "sampling": "without replacement within each realized 500-tree source pool; source-pool conditional", "seed": seed}, indent=2))
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    with (OUT_ROOT / "mixture_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"built {len(manifest_rows) // 500} mixture datasets and {len(manifest_rows)} manifest rows")


if __name__ == "__main__":
    main()
