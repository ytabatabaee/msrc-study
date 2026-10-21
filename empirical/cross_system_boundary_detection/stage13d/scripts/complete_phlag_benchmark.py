#!/usr/bin/env python3
"""Complete the prepared house-mouse matched-taxon Phlag benchmark.

This script reuses representative-selected inputs from Stage 13b, runs the
unmodified Phlag executable, validates state/posterior lengths, and scores
only finite runs.  Representative selection and Phlag execution never read
the inversion coordinates; those are loaded only by the scoring/plotting
phase.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
SOURCE = ROOT / "empirical/cross_system_boundary_detection/stage13b/phlag"
OUT = ROOT / "empirical/cross_system_boundary_detection/stage13d"
DATASET = "house_mouse_t_complex"
MOUSE_INPUTS = SOURCE / "raw_results/matched_inputs" / DATASET
MOUSE_RAW = OUT / "raw_results"
PHLAG_DEFAULT = Path("/private/tmp/phlag-venv/bin/phlag")
PHLAG_COMMIT = "76597d0e4ebd0e507fdafebf7f5163ef634a36d6"
OBS_START, OBS_END = 5_000_000, 37_054_999


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def parse_phlag_output(path: Path, expected: int) -> tuple[list[str], list[float], str]:
    lines = [line.strip() for line in path.read_text(errors="replace").splitlines()
             if line.strip() and not line.lstrip().startswith("#") and "," in line]
    state_lines = [line for line in lines if all(token in {"0", "1"} for token in line.split(","))]
    if len(state_lines) != 1:
        raise ValueError(f"expected exactly one binary state sequence, found {len(state_lines)}")
    state = state_lines[0].split(",")
    state_index = lines.index(state_lines[0])
    if state_index + 1 >= len(lines):
        raise ValueError("posterior sequence is missing")
    posterior = [float(token) for token in lines[state_index + 1].split(",")]
    if len(state) != expected or len(posterior) != expected:
        raise ValueError(f"Phlag length mismatch: states={len(state)}, posterior={len(posterior)}, expected={expected}")
    if not all(math.isfinite(value) for value in posterior):
        raise ValueError("nonfinite posterior probability")
    if not all(0.0 <= value <= 1.0 for value in posterior):
        raise ValueError("posterior probability outside [0,1]")
    status = "valid_no_detection" if len(set(state)) == 1 else "valid_informative"
    return state, posterior, status


def inversion_intervals() -> list[dict[str, int | str]]:
    # The fifth interval is documented in provenance but is not in the older
    # four-row TSV.  It is added here explicitly for this versioned benchmark.
    return [
        {"name": "Inv1", "start": 5_747_331, "end": 6_957_919},
        {"name": "Inv2", "start": 6_957_919, "end": 7_435_308},
        {"name": "Inv3", "start": 7_435_308, "end": 13_016_510},
        {"name": "Inv4", "start": 14_068_229, "end": 18_051_507},
        {"name": "Inv5", "start": 19_401_000, "end": 40_049_964},
    ]


def observable_reference(ref: dict[str, int | str]) -> tuple[int, int]:
    return max(int(ref["start"]), OBS_START), min(int(ref["end"]), OBS_END)


def merge_intervals(intervals: list[tuple[int, int]], max_gap: int = 5_000) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1] + max_gap + 1:
            merged.append((start, end))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
    return merged


def interval_length(interval: tuple[int, int]) -> int:
    return max(0, interval[1] - interval[0] + 1)


def union_length(intervals: list[tuple[int, int]]) -> int:
    return sum(interval_length(interval) for interval in merge_intervals(intervals, max_gap=0))


def intersection_length(left: list[tuple[int, int]], right: list[tuple[int, int]]) -> int:
    return sum(max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
               for a_start, a_end in left for b_start, b_end in right)


def score(predicted: list[tuple[int, int]], references: list[tuple[int, int]]) -> dict[str, object]:
    # Accuracy is set-based. Merge only overlapping or directly adjacent
    # windows; physical-gap tolerances belong to the separately reported
    # interval-construction step, not to the primary IoU calculation.
    predicted = merge_intervals(predicted, max_gap=0)
    references = merge_intervals(references, max_gap=0)
    pred_bp, ref_bp = union_length(predicted), union_length(references)
    overlap = intersection_length(predicted, references)
    union = pred_bp + ref_bp - overlap
    precision = overlap / pred_bp if pred_bp else 0.0
    recall = overlap / ref_bp if ref_bp else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {"predicted_intervals": ";".join(f"{a}-{b}" for a, b in predicted),
            "n_predicted_intervals": len(predicted), "predicted_bp": pred_bp,
            "reference_bp": ref_bp, "overlap_bp": overlap,
            "false_positive_bp": pred_bp - overlap, "false_negative_bp": ref_bp - overlap,
            "iou": overlap / union if union else 0.0, "precision": precision,
            "recall": recall, "f1": f1}


def coordinate_rows() -> list[dict[str, str]]:
    rows = read_tsv(ROOT / "empirical/house_mouse_t_complex/results/stage4_distance/stage4_window_distance.tsv")
    if len(rows) != 4046:
        raise RuntimeError(f"expected 4046 frozen mouse coordinate rows, found {len(rows)}")
    return rows


def selected_jobs(replicates: int) -> list[dict[str, object]]:
    jobs = []
    for rep in range(replicates):
        base = MOUSE_INPUTS / f"replicate_{rep:03d}"
        gene = base / "gene_trees.nwk"
        species = base / "species_tree.nwk"
        mapping = base / "window_mapping.tsv"
        if not all(path.exists() for path in (gene, species, mapping)):
            raise RuntimeError(f"missing prepared matched input for replicate {rep}: {base}")
        mappings = read_tsv(mapping)
        if len(mappings) != 4046:
            raise RuntimeError(f"replicate {rep} has {len(mappings)} loci, expected 4046")
        for edge in ("N002", "N003", "N004", "N005"):
            jobs.append({"replicate": rep, "focal_branch": edge, "gene": gene,
                         "species": species, "mapping": mapping, "n_windows": len(mappings)})
    return jobs


def run_job(job: dict[str, object], phlag: Path) -> dict[str, object]:
    rep = int(job["replicate"]); edge = str(job["focal_branch"])
    run_dir = MOUSE_RAW / f"replicate_{rep:03d}"; run_dir.mkdir(parents=True, exist_ok=True)
    output = run_dir / f"{edge}.txt"
    stdout = run_dir / f"{edge}.stdout"; stderr = run_dir / f"{edge}.stderr"
    command = ["/usr/bin/time", "-l", str(phlag), "-s", str(job["species"]),
               "-g", str(job["gene"]), "-e", edge, "-o", str(output)]
    started = time.time()
    with stdout.open("w") as out_handle, stderr.open("w") as err_handle:
        result = subprocess.run(command, stdout=out_handle, stderr=err_handle, text=True, check=False)
    runtime = time.time() - started
    text = stderr.read_text(errors="replace")
    match = re.search(r"maximum resident set size:\s*(\d+)", text)
    rss = int(match.group(1)) if match else "NA"
    status, error = "execution_failure", ""
    n_state1 = "NA"; posterior_max = "NA"
    # On macOS, `/usr/bin/time -l` can return nonzero after a successful child
    # because its optional sysctl memory probe is denied.  Parse the Phlag
    # output first; a complete finite state/posterior pair is authoritative for
    # scientific validity, while the wrapper status is retained in `error`.
    if output.exists():
        try:
            states, posterior, status = parse_phlag_output(output, int(job["n_windows"]))
            n_state1, posterior_max = states.count("1"), max(posterior)
            if result.returncode != 0:
                error = f"resource_wrapper_returncode={result.returncode}"
        except Exception as exc:
            status, error = "numerically_degenerate", str(exc)
    elif result.returncode != 0:
        error = f"returncode={result.returncode}"
    return {"dataset": DATASET, "replicate": rep, "focal_branch": edge,
            "status": status, "returncode": result.returncode,
            "phlag_exit_status": "not_separately_observable",
            "timing_wrapper_returncode": result.returncode,
            "n_windows": job["n_windows"], "n_state1": n_state1,
            "posterior_finite": status.startswith("valid"), "posterior_max": posterior_max,
            "runtime_seconds": runtime, "peak_rss_bytes": rss,
            "output": str(output), "error": error, "command": " ".join(command)}


def load_run_rows() -> list[dict[str, str]]:
    path = OUT / "run_manifest.tsv"
    return read_tsv(path) if path.exists() else []


def reclassify_existing(rows: list[dict[str, str]]) -> list[dict[str, object]]:
    """Reclassify completed output files independently of wrapper exit codes."""
    updated = []
    for row in rows:
        item = dict(row)
        output = Path(item.get("output", ""))
        if not output.is_absolute():
            output = OUT / output
        try:
            states, posterior, status = parse_phlag_output(output, 4046)
            item.update(status=status, n_state1=states.count("1"), posterior_finite=True,
                        posterior_max=max(posterior), error="")
        except Exception as exc:
            item.update(status="numerically_degenerate", posterior_finite=False,
                        posterior_max="NA", error=str(exc))
        updated.append(item)
    return updated


def predictions_and_accuracy(run_rows: list[dict[str, str]], coords: list[dict[str, str]]) -> None:
    prediction_rows = []
    branch_accuracy = []
    by_rep: dict[str, list[dict]] = {}
    for run in run_rows:
        if not run["status"].startswith("valid"):
            continue
        path = Path(run["output"])
        if not path.is_absolute():
            path = OUT / path
        states, posterior, _ = parse_phlag_output(path, 4046)
        rep = run["replicate"]; edge = run["focal_branch"]
        intervals = []
        for index, (state, prob, coord) in enumerate(zip(states, posterior, coords)):
            start, end = int(coord["start_bp"]), int(coord["end_bp"])
            prediction_rows.append({"replicate": rep, "focal_branch": edge, "window_index": index,
                                    "locus_id": coord["locus_id"], "start_bp": start, "end_bp": end,
                                    "state": state, "posterior_anomaly": prob})
            if state == "1": intervals.append((start, end))
        by_rep.setdefault(rep, []).append({"run": run, "states": states, "posterior": posterior, "intervals": intervals})
        for ref in inversion_intervals():
            ref_obs = observable_reference(ref)
            metrics = score(intervals, [ref_obs]) if ref_obs[0] <= ref_obs[1] else score([], [])
            branch_accuracy.append({"replicate": rep, "focal_branch": edge, "reference": ref["name"],
                                    "status": run["status"], "observable_reference_start": ref_obs[0],
                                    "observable_reference_end": ref_obs[1], **metrics})
    write_tsv(OUT / "window_predictions.tsv", prediction_rows,
              ["replicate", "focal_branch", "window_index", "locus_id", "start_bp", "end_bp", "state", "posterior_anomaly"])
    write_tsv(OUT / "branch_accuracy.tsv", branch_accuracy,
              ["replicate", "focal_branch", "reference", "status", "observable_reference_start", "observable_reference_end",
               "predicted_intervals", "n_predicted_intervals", "predicted_bp", "reference_bp", "overlap_bp",
               "false_positive_bp", "false_negative_bp", "iou", "precision", "recall", "f1"])
    aggregate_rows = []
    references = [(int(x["start"]), min(int(x["end"]), OBS_END)) for x in inversion_intervals()]
    references = [(max(a, OBS_START), b) for a, b in references if max(a, OBS_START) <= b]
    for rep, entries in sorted(by_rep.items(), key=lambda item: int(item[0])):
        union_intervals = [interval for entry in entries for interval in entry["intervals"]]
        max_intervals = [(int(coords[i]["start_bp"]), int(coords[i]["end_bp"]))
                         for i in range(4046) if any(entry["posterior"][i] >= 0.5 for entry in entries)]
        for aggregation, intervals in (("union_state1", union_intervals), ("max_posterior_ge_0.5", max_intervals)):
            aggregate_rows.append({"replicate": rep, "aggregation": aggregation, **score(intervals, references)})
    write_tsv(OUT / "replicate_accuracy.tsv", aggregate_rows,
              ["replicate", "aggregation", "predicted_intervals", "n_predicted_intervals", "predicted_bp", "reference_bp", "overlap_bp",
               "false_positive_bp", "false_negative_bp", "iou", "precision", "recall", "f1"])


def write_status_and_docs(run_rows: list[dict[str, object]], args: argparse.Namespace) -> None:
    fields = ["dataset", "replicate", "focal_branch", "status", "returncode", "phlag_exit_status",
              "timing_wrapper_returncode", "n_windows", "n_state1", "posterior_finite", "posterior_max",
              "runtime_seconds", "peak_rss_bytes", "output", "error", "command"]
    manifest_rows = []
    for source_row in run_rows:
        row = dict(source_row)
        output = Path(str(row.get("output", "")))
        try:
            row["output"] = str(output.relative_to(OUT)) if output.is_absolute() else str(output)
        except ValueError:
            row["output"] = str(output)
        row["command"] = "phlag -s <species_tree.nwk> -g <gene_trees.nwk> -e " + str(row.get("focal_branch", ""))
        row.setdefault("phlag_exit_status", "not_separately_observable")
        row.setdefault("timing_wrapper_returncode", row.get("returncode", "NA"))
        manifest_rows.append(row)
    write_tsv(OUT / "run_manifest.tsv", manifest_rows, fields)
    write_tsv(OUT / "numerical_diagnostics.tsv", manifest_rows, fields)
    (OUT / "manifest.json").write_text(json.dumps({
        "stage": "13d", "dataset": DATASET, "phlag_commit": PHLAG_COMMIT,
        "phlag_algorithm_modified": False,
        "source_inputs": "empirical/cross_system_boundary_detection/stage13b/phlag/raw_results/matched_inputs/house_mouse_t_complex",
        "replicates": args.replicates, "focal_edges": ["N002", "N003", "N004", "N005"],
        "observation_domain": [OBS_START, OBS_END], "five_inversion_intervals": inversion_intervals(),
        "max_gap_bp": 5000, "execution_workers": args.workers,
    }, indent=2) + "\n")
    (OUT / "report.md").write_text("""# Stage 13d house-mouse Phlag completion

This analysis reuses the prepared Stage 13b matched-taxon inputs and runs the
unmodified upstream Phlag executable at commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6`.
Representative choices are fixed across all 4,046 windows within each
replicate. Accuracy is calculated only for finite, length-matched posterior
outputs. The five structural intervals are scored on the observable domain
chr17:5,000,000–37,054,999; Inv5 is therefore right-censored.

The run manifest separates valid informative, valid no-detection, numerical
degeneracy, and execution failures. Per-edge results and two fixed aggregation
rules are reported separately. Structural coordinates are read only during
post-inference scoring.
""")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phlag", type=Path, default=PHLAG_DEFAULT)
    parser.add_argument("--replicates", type=int, default=6)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--score-only", action="store_true")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    coords = coordinate_rows()
    if not args.score_only:
        if not args.phlag.exists():
            raise SystemExit(f"Phlag executable not found: {args.phlag}")
        jobs = selected_jobs(args.replicates)
        if not args.run:
            print(json.dumps({"status": "preflight_only", "jobs": len(jobs), "workers": args.workers,
                              "phlag": str(args.phlag), "input_windows": 4046}, indent=2))
            return
        rows: list[dict[str, object]] = []
        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(run_job, job, args.phlag) for job in jobs]
            for future in as_completed(futures):
                row = future.result(); rows.append(row)
                print(f"replicate={row['replicate']} edge={row['focal_branch']} status={row['status']}", flush=True)
        rows.sort(key=lambda row: (int(row["replicate"]), str(row["focal_branch"])))
        write_status_and_docs(rows, args)
    else:
        rows = load_run_rows()
        if not rows:
            raise SystemExit("No Stage 13d run_manifest.tsv available")
        rows = reclassify_existing(rows)
        write_status_and_docs(rows, args)
    predictions_and_accuracy([dict(row) for row in rows], coords)
    print(json.dumps({"status": "complete", "runs": len(rows), "valid_runs": sum(str(row["status"]).startswith("valid") for row in rows),
                      "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
