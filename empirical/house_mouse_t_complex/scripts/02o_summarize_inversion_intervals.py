#!/usr/bin/env python3
"""Summarize quartet support and recombination classes inside/outside mm10 inversions."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
RESULTS = ROOT / "empirical/house_mouse_t_complex/results"
INTERVALS = RESULTS / "t_complex_inversion_intervals.tsv"
Q_SCAN = RESULTS / "stage2_fixed_quartet_scan.tsv"
RECOMB = RESULTS / "stage2_recombination_state_5kb.tsv"
OUTPUT = RESULTS / "t_complex_inversion_interval_summary.tsv"
SPECIES = ("domesticus", "musculus", "castaneus")
CLASSES = ("NO_RECENT_RECOMBINATION", "RECENT_OR_OLDER", "VERY_RECENT_OR_EXTENSIVE")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def side(row: dict[str, str], start: int, end: int) -> str:
    ws, we = int(row["start_bp"]), int(row["end_bp"])
    if ws >= start and we <= end:
        return "inside"
    if we < start or ws > end:
        return "outside"
    return "boundary_intersecting"


def main() -> None:
    intervals = read_tsv(INTERVALS)
    qrows = [row for row in read_tsv(Q_SCAN) if row["treatment"] == "ALL_TIPS"]
    q_by_id = {row["locus_id"]: row for row in qrows if int(row["n_resolved_quartets"]) > 0}
    recomb_by_id: dict[str, list[dict[str, str]]] = {}
    for row in read_tsv(RECOMB):
        recomb_by_id.setdefault(row["locus_id"], []).append(row)

    output: list[dict[str, str]] = []
    for interval in intervals:
        start, end = int(interval["start_bp"]), int(interval["end_bp"])
        assignments = {locus_id: side(row, start, end) for locus_id, row in q_by_id.items()}
        boundary_count = sum(value == "boundary_intersecting" for value in assignments.values())
        for group in ("inside", "outside"):
            locus_ids = [locus_id for locus_id, value in assignments.items() if value == group]
            qgroup = [q_by_id[locus_id] for locus_id in locus_ids]
            means = {
                key: sum(float(row[key]) for row in qgroup) / len(qgroup) if qgroup else float("nan")
                for key in ("q_species", "q_t_alt", "q_other")
            }
            means["mean_delta"] = means["q_t_alt"] - means["q_species"] if qgroup else float("nan")
            fields = {
                "inversion_name": interval["inversion_name"],
                "published_structural_id": interval["published_structural_id"],
                "comparison": group,
                "n_usable_5kb_windows": str(len(qgroup)),
                "n_boundary_intersecting_windows_excluded": str(boundary_count),
                "mean_q_species": f"{means['q_species']:.10f}",
                "mean_q_t_alt": f"{means['q_t_alt']:.10f}",
                "mean_q_other": f"{means['q_other']:.10f}",
                "mean_D_q_t_alt_minus_q_species": f"{means['mean_delta']:.10f}",
            }
            for species in SPECIES:
                species_rows = [row for locus_id in locus_ids for row in recomb_by_id.get(locus_id, []) if row["subspecies"] == species]
                counts = Counter(row["recombination_class"] for row in species_rows)
                for classification in CLASSES:
                    fields[f"{species}_{classification}_n"] = str(counts[classification])
                    fields[f"{species}_{classification}_fraction"] = f"{counts[classification] / len(species_rows):.10f}" if species_rows else "nan"
            output.append(fields)

    with OUTPUT.open("w", newline="") as handle:
        fieldnames = list(output[0])
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)


if __name__ == "__main__":
    main()
