#!/usr/bin/env python3
"""Collect and validate returned Stage 4E ASTER/ASTRAL4 cluster outputs."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from stage4d_core import ROOT, RESULTS, named_splits, sha256  # noqa: E402

BASE = ROOT / "empirical/neoaves_chr4"
TREE_DIR = BASE / "results/stage4e/trees"


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def classify(output_dir: Path, treatment: str, expected_taxa: int, mode: str) -> dict[str, str]:
    exit_code = output_dir / "exit_code.txt"
    tree = output_dir / ("tree.cu.nwk" if mode == "topology_cu" else "tree.su.nwk")
    # The CU/SU runner writes mode-specific logs; accept the historical
    # single-log names as well.
    stdout_paths = [output_dir / f"{mode}.stdout.log", output_dir / "stdout.log"]
    stderr_paths = [output_dir / f"{mode}.stderr.log", output_dir / "stderr.log"]
    row = {"treatment": treatment, "mode": mode, "status": "missing", "output_dir": str(output_dir), "tree_sha256": "", "error": "output_not_present"}
    if not output_dir.exists():
        return row
    if not exit_code.exists():
        row["status"] = "incomplete"
        return row
    code = exit_code.read_text().strip()
    if code != "0":
        row["status"] = "failed" if code else "incomplete"
        return row
    if not tree.exists() or tree.stat().st_size == 0:
        row["status"] = "incomplete"
        return row
    try:
        taxa, _ = named_splits(tree.read_text())
    except Exception as exc:
        row["status"] = "failed"
        row["error"] = f"Newick parse failed: {exc}"
        return row
    if len(taxa) != expected_taxa:
        row["status"] = "failed"
        row["error"] = f"Expected {expected_taxa} taxa, observed {len(taxa)}"
        return row
    logs = "\n".join(path.read_text() for path in stdout_paths + stderr_paths if path.exists())
    if logs and any(marker in logs.lower() for marker in ("error", "exception", "killed")):
        row["status"] = "failed"
        row["error"] = "Log contains error/exception/killed marker"
        return row
    row["status"] = "complete"
    row["tree_sha256"] = sha256(tree)
    return row


def copy_validated_tree(source: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and sha256(dest) != sha256(source):
        raise RuntimeError(f"Refusing to overwrite differing tree: {dest}")
    shutil.copyfile(source, dest)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cluster-root", type=Path, default=BASE / "results/stage4e/cluster")
    parser.add_argument("--manifest", type=Path, default=RESULTS / "stage4e_input_manifest.tsv")
    parser.add_argument("--status-output", type=Path, default=RESULTS / "stage4e_cu_su_collect_status.tsv")
    parser.add_argument("--legacy-status-output", type=Path,
                        help="Explicitly replace the historical single-unit status table.")
    args = parser.parse_args()
    rows = []
    for item in read_manifest(args.manifest):
        expected_taxa = 48 if item["scope"] == "jarvis48" else 363
        output_dir = args.cluster_root / item["treatment"]
        for mode in ("topology_cu", "topology_su"):
            row = classify(output_dir, item["treatment"], expected_taxa, mode)
            rows.append(row)
            if row["status"] == "complete":
                suffix = "cu" if mode == "topology_cu" else "su"
                copy_validated_tree(output_dir / f"tree.{suffix}.nwk", TREE_DIR / f"{item['treatment']}.{suffix}.nwk")
    fields = ["treatment", "mode", "status", "output_dir", "tree_sha256", "error"]
    out = args.status_output
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    if args.legacy_status_output:
        shutil.copyfile(out, args.legacy_status_output)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
