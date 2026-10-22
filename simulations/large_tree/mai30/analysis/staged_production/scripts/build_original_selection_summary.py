#!/usr/bin/env python3
"""Summarize all original no-outgroup Mai gene-tree collections.

This is an audit-only pass over existing ``truegenetrees`` files.  It does
not simulate or alter the validated four-replicate baseline summary.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "analysis/scripts"))
from compute_original_baseline import summarize_replicate  # noqa: E402

CONDITIONS = (
    "outgroup.0.species.0.15.genes.0.15",
    "outgroup.0.species.1.5.genes.1.5",
    "outgroup.0.species.5.genes.5",
)


def main():
    output = ROOT / "analysis/staged_production/original_baseline_all.tsv"
    rows = []
    for condition in CONDITIONS:
        for rep_dir in sorted((ROOT / "datasets/raw" / condition).glob("[0-9][0-9][0-9]")):
            if not (rep_dir / "truegenetrees").exists():
                continue
            rows.append(summarize_replicate(condition, rep_dir.name))
    if len(rows) != 300:
        raise ValueError(f"expected 300 original no-outgroup summaries, found {len(rows)}")
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {output} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
