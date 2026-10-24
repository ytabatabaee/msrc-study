#!/usr/bin/env python3
"""Create the frozen 25-job J48 replication manifest.

With ``--realize``, this script is also the cluster-side, deterministic input
builder.  It never launches ASTRAL and it refuses to replace an existing
mixture input with different content.
"""
from __future__ import annotations

import argparse, csv, hashlib, json, random, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / "empirical" / "neoaves_chr4"
SOURCE = BASE / "data" / "stage4d" / "named_recollapsed.tre"
LOCUS = BASE / "results" / "stage4d_locus_manifest.tsv"
TAXA = BASE / "data" / "stage4d" / "jarvis48_taxa.tsv"
RESULTS = BASE / "results" / "stage4e"
MANIFEST = RESULTS / "j48_mixture_replication_manifest.tsv"
sys.path.insert(0, str(ROOT / "empirical" / "neoaves_chr4" / "scripts"))
from stage4d_core import emit_newick, parse_newick, prune_tree  # noqa: E402

TOTAL = 62945
EPSILONS = (0.02, 0.05, 0.075, 0.10, 0.15)

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def tag(epsilon: float) -> str:
    return f"E{epsilon:.3f}".replace(".", "p")

def planned_rows(source: Path = SOURCE, locus: Path = LOCUS, taxa: Path = TAXA) -> list[dict[str, object]]:
    if not source.exists() or not locus.exists() or not taxa.exists():
        raise FileNotFoundError("source tree, locus manifest, and taxon mapping are required")
    source_hash, locus_hash, taxa_hash = sha256(source), sha256(locus), sha256(taxa)
    for index, epsilon in enumerate(EPSILONS):
        n_chr4 = round(TOTAL * epsilon)
        for replicate in range(1, 6):
            seed = 76000 + index * 100 + replicate
            mixture_id = f"J48_MIX_REP_{tag(epsilon)}_R{replicate:02d}"
            yield {
                "mixture_id": mixture_id, "epsilon_requested": f"{epsilon:.3f}",
                "replicate": replicate, "seed": seed, "total_records": TOTAL,
                "n_chr4_records": n_chr4, "n_background_records": TOTAL - n_chr4,
                "realized_epsilon": f"{n_chr4 / TOTAL:.12f}",
                "source_tree_sha256": source_hash, "locus_manifest_sha256": locus_hash,
                "taxon_mapping_sha256": taxa_hash, "status": "prepared_not_run",
                "input_relative_path": f"inputs/{mixture_id}.tre",
                "result_relative_path": f"results/{mixture_id}",
            }

def write_manifest(path: Path = MANIFEST) -> None:
    fields = list(next(iter(planned_rows())).keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(planned_rows())

def read_manifest(path: Path) -> list[dict[str, str]]:
    result = rows(path)
    if len(result) != 25 or len({r["mixture_id"] for r in result}) != 25:
        raise ValueError("replication manifest must contain 25 unique jobs")
    return result

def realize(row: dict[str, str], input_root: Path, source: Path, locus: Path, taxa_path: Path) -> dict[str, object]:
    if sha256(source) != row["source_tree_sha256"] or sha256(locus) != row["locus_manifest_sha256"] or sha256(taxa_path) != row["taxon_mapping_sha256"]:
        raise ValueError("source checksum differs from frozen manifest")
    output = input_root / Path(row["input_relative_path"]).name
    sidecar = output.with_suffix(".json")
    if output.exists() or sidecar.exists():
        if sidecar.exists() and output.exists() and json.loads(sidecar.read_text()).get("input_sha256") == sha256(output):
            return json.loads(sidecar.read_text()) | {"status": "already_validated"}
        raise FileExistsError(f"refusing to replace existing mixture input: {output}")
    locus_rows, frozen_taxa = rows(locus), {r["stiller2024_taxon"] for r in rows(taxa_path)}
    if len(locus_rows) != 63430 or len(frozen_taxa) != 48:
        raise ValueError("unexpected source dimensions")
    strata = {"background": [], "chr4": []}
    with source.open() as f:
        for locus_row, line in zip(locus_rows, f, strict=True):
            tree = prune_tree(parse_newick(line), frozen_taxa)
            if tree is None or len(tree[0] or []) < 2:
                continue
            strata["chr4" if locus_row["is_chr4"].lower() in {"true", "1", "yes"} else "background"].append(emit_newick(tree) + "\n")
    target = int(row["total_records"]); n_chr4 = int(row["n_chr4_records"])
    rng = random.Random(int(row["seed"]))
    selected = [strata["background"][rng.randrange(len(strata["background"]))] for _ in range(target - n_chr4)]
    selected += [strata["chr4"][rng.randrange(len(strata["chr4"]))] for _ in range(n_chr4)]
    rng.shuffle(selected)
    input_root.mkdir(parents=True, exist_ok=True)
    partial = output.with_suffix(".tre.partial")
    with partial.open("x") as f:
        f.writelines(selected)
    partial.replace(output)
    record = {**row, "input_sha256": sha256(output), "n_background_distinct": len(strata["background"]), "n_chr4_distinct": len(strata["chr4"]), "status": "prepared"}
    sidecar.write_text(json.dumps(record, indent=2) + "\n")
    return record

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write-manifest", action="store_true")
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--realize", action="store_true")
    parser.add_argument("--job-index", type=int)
    parser.add_argument("--input-root", type=Path)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--locus-manifest", type=Path, default=LOCUS)
    parser.add_argument("--taxa", type=Path, default=TAXA)
    args = parser.parse_args()
    if args.write_manifest:
        write_manifest(args.manifest)
    if args.realize:
        if args.job_index is None or args.input_root is None:
            raise SystemExit("--realize requires --job-index and --input-root")
        job_rows = read_manifest(args.manifest)
        if not 0 <= args.job_index < len(job_rows): raise SystemExit("job index outside manifest")
        print(json.dumps(realize(job_rows[args.job_index], args.input_root, args.source, args.locus_manifest, args.taxa), indent=2))
    if not args.write_manifest and not args.realize:
        write_manifest(args.manifest)

if __name__ == "__main__": main()
