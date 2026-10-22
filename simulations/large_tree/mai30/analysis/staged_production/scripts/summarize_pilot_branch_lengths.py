#!/usr/bin/env python3
"""Summarize the corrected four-replicate branch table.

This does not touch the frozen pilot tables; it writes a staged audit copy
with explicit summary statistics and per-split values.
"""

from __future__ import annotations

import csv
import math
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "analysis/suppression_effect/branch_length_comparison.tsv"
OUT = ROOT / "analysis/staged_production/results"


def read_rows():
    with SOURCE.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    rows = read_rows()
    for row in rows:
        if "relative_error" not in row:
            true = float(row["true_culength"])
            row["relative_error"] = float(row["signed_error"]) / true if true else ""
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "pilot_branch_lengths.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    summary = []
    for method in sorted({r["method"] for r in rows}):
        subset = [r for r in rows if r["method"] == method]
        signed = [float(r["signed_error"]) for r in subset]
        absolute = [float(r["absolute_error"]) for r in subset]
        relative = [float(r["relative_error"]) for r in subset if r["relative_error"] != ""]
        summary.append({
            "method": method,
            "matching_internal_splits": len(subset),
            "mean_signed_cu_error": statistics.mean(signed),
            "mean_absolute_cu_error": statistics.mean(absolute),
            "median_absolute_cu_error": statistics.median(absolute),
            "mean_relative_error": statistics.mean(relative) if relative else "",
            "median_relative_error": statistics.median(relative) if relative else "",
            "max_true_culength": max(float(r["true_culength"]) for r in subset),
        })
    with (OUT / "pilot_branch_length_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(summary)
    print(f"wrote {OUT / 'pilot_branch_lengths.tsv'} ({len(rows)} rows)")
    print(f"wrote {OUT / 'pilot_branch_length_summary.tsv'} ({len(summary)} methods)")


if __name__ == "__main__":
    main()
