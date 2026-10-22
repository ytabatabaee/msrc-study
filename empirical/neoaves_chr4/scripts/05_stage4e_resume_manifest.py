#!/usr/bin/env python3
"""Classify Stage 4E outputs without submitting or overwriting any job."""
from __future__ import annotations
import argparse
import csv
from pathlib import Path

def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    with args.manifest.open() as handle:
        for item in csv.DictReader(handle, delimiter="\t"):
            out = args.output_root / item["treatment"]
            annotated = out / "tree.annotated.nwk"
            derived = all((out / name).exists() for name in ("tree.cu.nwk", "tree.su.nwk", "branch_lengths.tsv"))
            if annotated.exists() and derived and (out / "exit_code.txt").read_text().strip() == "0":
                state = "valid_annotated_and_derived"
            elif annotated.exists():
                state = "annotated_exists_extract_only"
            elif (out / "tree.su.nwk").exists():
                state = "legacy_su_extract_only"
            elif (out / "tree.cu.nwk").exists():
                state = "cu_only_incomplete"
            elif out.exists():
                state = "corrupted_or_incomplete"
            else:
                state = "no_output_run_inference"
            rows.append({"treatment": item["treatment"], "scope": item["scope"],
                         "output_dir": str(out), "state": state})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["treatment", "scope", "output_dir", "state"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    print(f"Wrote {args.output}")

if __name__ == "__main__":
    main()
