#!/usr/bin/env python3
"""Build one deterministic J48 mixture on a Bridges-2 compute node.

This command is intentionally a data-preparation step.  It streams the
verified named-recollapsed source, prunes to the frozen 48-taxon set, writes
temporary background/chr4 strata, then samples with replacement to realize the
requested record weight.  It never launches ASTRAL.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "empirical" / "neoaves_chr4"
SOURCE = BASE / "data" / "stage4d" / "named_recollapsed.tre"
LOCUS_MANIFEST = BASE / "results" / "stage4d_locus_manifest.tsv"
TAXA = BASE / "data" / "stage4d" / "jarvis48_taxa.tsv"
OUT = BASE / "results" / "stage4e" / "threshold_mixtures"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage4d_core import emit_newick, parse_newick, prune_tree  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def boolish(value: str) -> bool:
    return value.lower() in {"true", "1", "yes"}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_mixture(mixture_id: str, epsilon: float, seed: int, output: Path) -> dict[str, object]:
    if not 0.0 <= epsilon <= 1.0:
        raise ValueError("epsilon must be in [0, 1]")
    manifest = read_tsv(LOCUS_MANIFEST)
    if len(manifest) != 63430:
        raise ValueError("unexpected locus manifest count")
    taxa = {row["stiller2024_taxon"] for row in read_tsv(TAXA)}
    if len(taxa) != 48:
        raise ValueError("unexpected 48-taxon mapping")
    strata = {"background": OUT / f".{mixture_id}.background.tre", "affected_chr4": OUT / f".{mixture_id}.affected.tre"}
    OUT.mkdir(parents=True, exist_ok=True)
    counts = {key: 0 for key in strata}
    with SOURCE.open() as source, strata["background"].open("w") as background, strata["affected_chr4"].open("w") as affected:
        for row, line in zip(manifest, source, strict=True):
            tree = prune_tree(parse_newick(line), taxa)
            if tree is None:
                continue
            children = tree[0] if tree[0] else []
            if len(children) < 2:
                continue
            rendered = emit_newick(tree) + "\n"
            stream = affected if boolish(row["is_chr4"]) else background
            key = "affected_chr4" if boolish(row["is_chr4"]) else "background"
            stream.write(rendered)
            counts[key] += 1
    total = counts["background"] + counts["affected_chr4"]
    if total == 0 or not all(counts.values()):
        raise ValueError(f"empty stratum after 48-taxon pruning: {counts}")
    target = 62945
    n_affected = round(target * epsilon)
    n_background = target - n_affected
    rng = random.Random(seed)
    selected = [("background", rng.randrange(counts["background"])) for _ in range(n_background)]
    selected += [("affected_chr4", rng.randrange(counts["affected_chr4"])) for _ in range(n_affected)]
    rng.shuffle(selected)
    offsets: dict[str, list[int]] = {}
    for key, path in strata.items():
        positions = []
        with path.open("rb") as handle:
            while True:
                position = handle.tell(); line = handle.readline()
                if not line: break
                positions.append(position)
        offsets[key] = positions
    with output.open("w") as sink:
        handles = {key: path.open() for key, path in strata.items()}
        try:
            for key, index in selected:
                handles[key].seek(offsets[key][index])
                sink.write(handles[key].readline())
        finally:
            for handle in handles.values(): handle.close()
    for path in strata.values(): path.unlink(missing_ok=True)
    return {"mixture_id": mixture_id, "epsilon_requested": epsilon, "seed": seed,
            "background_distinct_loci": counts["background"], "affected_distinct_loci": counts["affected_chr4"],
            "total_records": target, "affected_record_fraction": n_affected / target,
            "output": str(output), "output_sha256": sha256(output), "source_sha256": sha256(SOURCE)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mixture-id", required=True)
    parser.add_argument("--epsilon", type=float, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or OUT / f"{args.mixture_id}.tre"
    result = write_mixture(args.mixture_id, args.epsilon, args.seed, output)
    print(result)


if __name__ == "__main__":
    main()
