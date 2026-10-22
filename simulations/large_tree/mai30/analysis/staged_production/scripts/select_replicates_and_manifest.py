#!/usr/bin/env python3
"""Select the prespecified Mai30 staged sample and write its 90-task manifest.

Selection uses only the original Mai gene-tree discordance summary.  The
existing 900-task production manifest supplies the already-audited population
mapping, origin choice, and deterministic seeds; no MSRC result enters the
selection.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


CONDITIONS = (
    "outgroup.0.species.0.15.genes.0.15",
    "outgroup.0.species.1.5.genes.1.5",
    "outgroup.0.species.5.genes.5",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"no rows for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def select(rows: list[dict[str, str]], condition: str) -> list[dict[str, str]]:
    candidates = sorted(
        (row for row in rows if row["condition"] == condition),
        key=lambda row: (float(row["mean_normalized_rf"]), row["replicate"]),
    )
    if len(candidates) != 100:
        raise ValueError(f"{condition}: expected 100 original Mai replicates, found {len(candidates)}")
    selected: list[dict[str, str]] = []
    for stratum in range(5):
        block = candidates[stratum * 20 : (stratum + 1) * 20]
        for position, row in ((0, block[0]), (19, block[-1])):
            out = dict(row)
            out["selection_stratum"] = stratum + 1
            out["selection_position_in_stratum"] = position + 1
            out["selection_rule"] = "sort original mean_normalized_rf; select both stratum endpoints in five equal strata"
            selected.append(out)
    return selected


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--study-root", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--selection-output", type=Path)
    parser.add_argument("--manifest-output", type=Path)
    args = parser.parse_args()
    root = args.study_root.resolve()
    processed = root / "datasets/processed"
    selection_output = args.selection_output or root / "analysis/staged_production/selected_replicates.tsv"
    manifest_output = args.manifest_output or root / "analysis/staged_production/staged_manifest.tsv"

    baseline_path = root / "analysis/staged_production/original_baseline_all.tsv"
    if not baseline_path.exists():
        raise FileNotFoundError(
            f"{baseline_path} is missing; run build_original_selection_summary.py first. "
            "The four-replicate validated baseline is intentionally not overwritten."
        )
    baseline = read_tsv(baseline_path)
    selected = [row for condition in CONDITIONS for row in select(baseline, condition)]
    if len(selected) != 30 or {r["condition"] for r in selected} != set(CONDITIONS):
        raise AssertionError("selection did not produce 10 replicates per condition")
    write_tsv(selection_output, selected)

    production = read_tsv(processed / "mai30_production_manifest.tsv")
    by_key = {(row["condition"], row["replicate"], row["treatment"]): row for row in production}
    staged = []
    task_id = 0
    for selected_row in selected:
        key = (selected_row["condition"], selected_row["replicate"])
        for treatment in ("msrc_no_suppression", "msrc_moderate_suppression", "msrc_strong_suppression"):
            if (*key, treatment) not in by_key:
                raise KeyError(f"missing production metadata for {key} {treatment}")
            row = dict(by_key[(*key, treatment)])
            row["task_id"] = str(task_id)
            row["output_root"] = "datasets/processed/mai30_staged_runs"
            row["status"] = "planned"
            row["selection_stratum"] = selected_row["selection_stratum"]
            row["selection_rule"] = selected_row["selection_rule"]
            staged.append(row)
            task_id += 1
    if len(staged) != 90:
        raise AssertionError(f"expected 90 staged tasks, found {len(staged)}")
    write_tsv(manifest_output, staged)
    print(f"wrote {selection_output} ({len(selected)} replicates)")
    print(f"wrote {manifest_output} ({len(staged)} collections)")


if __name__ == "__main__":
    main()
