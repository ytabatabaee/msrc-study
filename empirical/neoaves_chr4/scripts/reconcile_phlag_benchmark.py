#!/usr/bin/env python3
"""Rebuild the reconciled Neoaves Phlag comparison from retained outputs.

This script never invokes Phlag or fits a detector. It reads the validated
Stage 11B predictions, the matched four-taxon detector output, the preserved
original multi-taxon detector output, and the canonical structural candidate
regions used only for post-fit evaluation.
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical/neoaves_chr4"
RES = EMP / "results/stage11_phlag_benchmark"
MULTI = EMP / "results/stage11_original_multitaxon_detector"
FIG = EMP / "figures/neoaves_chr4_phlag_comparison"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    fields = list(rows[0])
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def union(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], end))
        else:
            result.append((start, end))
    return result


def metrics(predicted: list[tuple[int, int]], reference: list[tuple[int, int]]) -> dict[str, float | int]:
    predicted = union(predicted)
    reference = union(reference)
    intersection = sum(
        max(0, min(p_end, r_end) - max(p_start, r_start))
        for p_start, p_end in predicted
        for r_start, r_end in reference
    )
    predicted_bp = sum(end - start for start, end in predicted)
    reference_bp = sum(end - start for start, end in reference)
    union_bp = predicted_bp + reference_bp - intersection
    precision = intersection / predicted_bp if predicted_bp else 0.0
    recall = intersection / reference_bp if reference_bp else 0.0
    return {
        "n_predicted_intervals": len(predicted),
        "predicted_bp": predicted_bp,
        "reference_bp": reference_bp,
        "intersection_bp": intersection,
        "union_bp": union_bp,
        "iou": intersection / union_bp if union_bp else 0.0,
        "precision": precision,
        "recall": recall,
        "f1": 2 * precision * recall / (precision + recall) if precision + recall else 0.0,
        "false_positive_bp": predicted_bp - intersection,
        "false_negative_bp": reference_bp - intersection,
    }


def parse_intervals(path: Path, start_field: str = "start", end_field: str = "end") -> list[tuple[int, int]]:
    return [(int(row[start_field]), int(row[end_field])) for row in read_tsv(path)]


def parse_prediction_string(value: str) -> list[tuple[int, int]]:
    return [
        (int(part.split("-", 1)[0]), int(part.split("-", 1)[1]))
        for part in value.split(";")
        if part.strip()
    ]


def detector_intervals() -> dict[str, list[tuple[int, int]]]:
    matched = read_tsv(RES / "matched_detector_predictions.tsv")[0]
    multi = read_tsv(MULTI / "our_detector_predictions.tsv")
    return {
        "official_phlag_matched_four_taxon": parse_intervals(RES / "phlag_predicted_intervals.tsv"),
        "our_matched_four_taxon_topology": [(int(matched["start"]), int(matched["end"]))],
        "our_original_multitaxon_topology": parse_prediction_string(next(row["predicted_intervals"] for row in multi if row["method"] == "topology")),
        "our_original_multitaxon_topology_plus_B": parse_prediction_string(next(row["predicted_intervals"] for row in multi if row["method"] == "topology_plus_B")),
    }


def rebuild_tables(intervals: dict[str, list[tuple[int, int]]], references: list[dict[str, str]]) -> None:
    reference_union = [(int(row["canonical_start"]), int(row["canonical_end"])) for row in references]
    rows = []
    for method, predicted in intervals.items():
        matched = "matched_four_taxon" in method
        evaluation = "matched_input_post_fit" if matched else "unmatched_input_post_fit"
        source = "matched four-taxon manifest (5,619 loci)" if matched else "original multi-taxon observations (6,096 loci)"
        targets = [("all_canonical_regions_union", reference_union)]
        targets.extend((row["canonical_region_id"], [(int(row["canonical_start"]), int(row["canonical_end"]))]) for row in references)
        for target, reference in targets:
            row = metrics(predicted, reference)
            rows.append({
                "method": method,
                "evaluation": evaluation,
                "source_observations": source,
                "target": target,
                **row,
            })
    write_tsv(RES / "interval_comparison.tsv", rows)
    write_tsv(RES / "boundary_accuracy.tsv", rows)


def rebuild_figure(intervals: dict[str, list[tuple[int, int]]], references: list[dict[str, str]]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    tracks = [
        ("Canonical structural candidates", "#777777", [(int(r["canonical_start"]), int(r["canonical_end"])) for r in references]),
        ("Official Phlag · matched 4-taxon", "#4c78a8", intervals["official_phlag_matched_four_taxon"]),
        ("Our topology · matched 4-taxon", "#f58518", intervals["our_matched_four_taxon_topology"]),
        ("Our topology · original multi-taxon", "#54a24b", intervals["our_original_multitaxon_topology"]),
        ("Our topology+B · original multi-taxon", "#e45756", intervals["our_original_multitaxon_topology_plus_B"]),
    ]
    xmin, xmax = 12.0, 60.0
    fig, axes = plt.subplots(len(tracks), 1, figsize=(13, 5.8), sharex=True, gridspec_kw={"hspace": 0.12})
    for axis, (label, color, ranges) in zip(axes, tracks):
        axis.set_xlim(xmin, xmax)
        axis.set_ylim(-0.5, 0.5)
        axis.set_yticks([])
        axis.grid(axis="x", color="#eeeeee", linewidth=0.6)
        for start, end in ranges:
            axis.plot([start / 1e6, end / 1e6], [0, 0], color=color, linewidth=9, solid_capstyle="butt")
        axis.text(0.005, 0.5, label, transform=axis.transAxes, va="center", ha="left", fontsize=9,
                  bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 1.5})
        for spine in ("top", "right", "left"):
            axis.spines[spine].set_visible(False)
    axes[-1].set_xlabel("GalGal6 chromosome 4 position (Mb)")
    fig.suptitle("Neoaves chromosome 4: reconciled Phlag and detector intervals", fontsize=13)
    fig.text(0.01, 0.01, "Structural regions are synteny-based candidate annotations and are used only for post-fit evaluation.\n"
             "Matched and original multi-taxon tracks use different observation domains.", fontsize=8)
    fig.savefig(FIG.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(FIG.with_suffix(".png"), dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    references = read_tsv(EMP / "results/canonical_structural_regions.tsv")
    intervals = detector_intervals()
    rebuild_tables(intervals, references)
    rebuild_figure(intervals, references)
    print("reconciled", ", ".join(f"{name}={len(value)}" for name, value in intervals.items()))


if __name__ == "__main__":
    main()
