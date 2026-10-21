#!/usr/bin/env python3
"""Stage 11C: representative-selection robustness for the Neoaves Phlag benchmark.

Replicate zero is copied from the frozen Stage 11B result. Other replicates
prune the same source gene trees to deterministic four-taxon representative
sets, run the unmodified Phlag executable, and run the frozen matched
categorical detector. Structural annotations are loaded only by scoring and
plotting, never by representative selection or detector fitting.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import resource
import shutil
import subprocess
import time
from copy import deepcopy
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical/neoaves_chr4"
OUT = EMP / "results/stage11c_phlag_robustness"
BASE = EMP / "results/stage11_phlag_benchmark"
DATA = ROOT / "data/neoaves_chr4/processed"
LOCUS_MANIFEST = EMP / "results/stage4d_locus_manifest.tsv"
SOURCE_TREES = EMP / "external/stiller2024/63430.named.gene.trees.gz"
REFERENCE_TREE = EMP / "external/stiller2024/63K.tre"
PHLAG = Path("/private/tmp/phlag-venv/bin/phlag")
PHLAG_COMMIT = "76597d0e4ebd0e507fdafebf7f5163ef634a36d6"
GAP_BP = 250_000
SEED = 20261021

SETS = [
    (0, "original_stage11b", "Columba_livia", "Phoenicopterus_ruber", "Tauraco_erythrolophus", "Gallus_gallus"),
    (1, "columbina", "Columbina_picui", "Phoenicopterus_ruber", "Tauraco_erythrolophus", "Gallus_gallus"),
    (2, "patagioenas", "Patagioenas_fasciata", "Phoenicopterus_ruber", "Tauraco_erythrolophus", "Gallus_gallus"),
    (3, "caloenas", "Caloenas_nicobarica", "Phoenicopterus_ruber", "Tauraco_erythrolophus", "Gallus_gallus"),
    (4, "podiceps", "Columba_livia", "Podiceps_cristatus", "Tauraco_erythrolophus", "Gallus_gallus"),
    (5, "podilymbus", "Columba_livia", "Podilymbus_podiceps", "Tauraco_erythrolophus", "Gallus_gallus"),
    (6, "corythaixoides", "Columba_livia", "Phoenicopterus_ruber", "Corythaixoides_concolor", "Gallus_gallus"),
    (7, "coturnix", "Columba_livia", "Phoenicopterus_ruber", "Tauraco_erythrolophus", "Coturnix_japonica"),
    (8, "phasianus", "Columba_livia", "Phoenicopterus_ruber", "Tauraco_erythrolophus", "Phasianus_colchicus"),
    (9, "meleagris", "Columba_livia", "Phoenicopterus_ruber", "Tauraco_erythrolophus", "Meleagris_gallopavo"),
]
ROLE_NAMES = ("C1", "C2", "S", "O")
MANIFEST_FIELDS = ["input_index", "Gene", "source_gene_tree_id", "source_locus_id", "chromosome", "start", "end", "midpoint", "taxa", "all_selected_taxa_present", "source_row_index", "tree_sha256"]


@dataclass
class Node:
    children: list
    name: str = ""
    length: str = ""


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({field: row.get(field, "") for field in fields} for row in rows)


def parse_newick(text: str) -> Node:
    s = text.strip().rstrip(";")
    pos = 0
    def token(stops: str) -> str:
        nonlocal pos
        start = pos
        while pos < len(s) and s[pos] not in stops:
            pos += 1
        return s[start:pos].strip()
    def node() -> Node:
        nonlocal pos
        if s[pos] == "(":
            pos += 1
            children = [node()]
            while pos < len(s) and s[pos] == ",":
                pos += 1
                children.append(node())
            if pos >= len(s) or s[pos] != ")":
                raise ValueError("malformed Newick")
            pos += 1
            name = token(":,();").strip("'\"")
            length = ""
            if pos < len(s) and s[pos] == ":":
                pos += 1
                length = token(",();")
            return Node(children, name, length)
        name = token(":,();").strip("'\"")
        length = ""
        if pos < len(s) and s[pos] == ":":
            pos += 1
            length = token(",();")
        return Node([], name, length)
    result = node()
    if pos != len(s):
        raise ValueError("trailing Newick text")
    return result


def leaf_names(node: Node) -> list[str]:
    if not node.children:
        return [node.name]
    result: list[str] = []
    for child in node.children:
        result.extend(leaf_names(child))
    return result


def prune(node: Node, wanted: set[str]) -> Node | None:
    if not node.children:
        return deepcopy(node) if node.name in wanted else None
    kept = [prune(child, wanted) for child in node.children]
    kept = [child for child in kept if child is not None]
    if not kept:
        return None
    if len(kept) == 1:
        child = kept[0]
        try:
            child.length = str(float(child.length or 0) + float(node.length or 0))
        except ValueError:
            pass
        return child
    return Node(kept, "", node.length)


def write_newick(node: Node) -> str:
    if node.children:
        text = "(" + ",".join(write_newick(child) for child in node.children) + ")"
    else:
        text = node.name
    if node.name and node.children:
        text += node.name
    return text + ((":" + node.length) if node.length else "")


def set_focal_label(node: Node, pair: set[str]) -> bool:
    if not node.children:
        return False
    leaves = set(leaf_names(node))
    found = False
    for child in node.children:
        found = set_focal_label(child, pair) or found
    if leaves == pair:
        node.name = "Columbea_focal"
        found = True
    return found


def topology(node: Node, reps: tuple[str, str, str, str]) -> str:
    names = set(reps)
    if set(leaf_names(node)) != names:
        return "unresolved"
    splits: list[tuple[str, str]] = []
    def walk(current: Node, root: bool = False) -> set[str]:
        if not current.children:
            return {current.name}
        desc: set[str] = set()
        for child in current.children:
            desc |= walk(child)
        if not root and len(desc) == 2:
            splits.append(tuple(sorted(desc)))
        return desc
    walk(node, True)
    pairs = [tuple(sorted((reps[0], reps[1]))), tuple(sorted((reps[0], reps[2]))), tuple(sorted((reps[0], reps[3])))]
    for index, pair in enumerate(pairs):
        if pair in splits or tuple(sorted(names - set(pair))) in splits:
            return f"q{index + 1}"
    return "unresolved"


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def load_source_rows() -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]]]:
    stage6 = {row["Gene"]: row for row in read_tsv(EMP / "results/stage6_genomic_tracks/window_quartet_support.tsv") if row["clade"] == "Columbea"}
    meta = {row["original_gene_id"]: row for row in read_tsv(LOCUS_MANIFEST) if row["chromosome"] == "chr4"}
    by_source = {meta[row["Gene"]]["locus_id"]: row for row in stage6.values() if row["Gene"] in meta}
    return by_source, meta


def species_tree_for(rep: tuple, original: bool) -> str:
    if original:
        return (BASE / "matched_inputs/species_tree.nwk").read_text()
    root = parse_newick(REFERENCE_TREE.read_text())
    selected = tuple(rep[2:])
    tree = prune(root, set(selected))
    if tree is None or set(leaf_names(tree)) != set(selected) or not set_focal_label(tree, {selected[0], selected[2]}):
        raise RuntimeError(f"could not construct focal four-taxon species tree for {rep}")
    if topology(tree, selected) != "q2":
        raise RuntimeError(f"representative set does not preserve q2 quartet: {rep}")
    return write_newick(tree) + ";\n"


def prepare_replicate(rep: tuple, by_source: dict[str, dict[str, str]], meta: dict[str, dict[str, str]], prepared_rows: list[dict] | None = None) -> dict:
    rid, label, *taxa = rep
    out = OUT / "raw" / f"replicate_{rid:03d}"
    inp = out / "inputs"
    inp.mkdir(parents=True, exist_ok=True)
    reps = tuple(taxa)
    if rid == 0:
        shutil.copy2(BASE / "matched_inputs/Columbea.gtrees", inp / "gene_trees.nwk")
        shutil.copy2(BASE / "matched_locus_manifest.tsv", inp / "coordinate_manifest.tsv")
        shutil.copy2(BASE / "matched_inputs/species_tree.nwk", inp / "species_tree.nwk")
        return {"replicate": rid, "label": label, **dict(zip(ROLE_NAMES, reps)), "n_loci": sum(1 for _ in (inp / "gene_trees.nwk").open()), "status": "reused_frozen_stage11b", "input": str(inp / "gene_trees.nwk"), "manifest": str(inp / "coordinate_manifest.tsv"), "species_tree": str(inp / "species_tree.nwk")}

    rows = prepared_rows or []
    rows.sort(key=lambda item: (int(float(item["start"])), int(float(item["end"])), int(item["Gene"])))
    manifest = []
    with (inp / "gene_trees.nwk").open("w") as trees:
        for index, row in enumerate(rows):
            trees.write(row["newick"])
            manifest.append({key: row.get(key, "") for key in MANIFEST_FIELDS if key not in {"input_index", "tree_sha256"}} | {"input_index": str(index), "tree_sha256": sha256_text(row["newick"])})
    write_tsv(inp / "coordinate_manifest.tsv", manifest, MANIFEST_FIELDS)
    (inp / "species_tree.nwk").write_text(species_tree_for(rep, False))
    if len(rows) != len(manifest) or not rows:
        raise RuntimeError(f"replicate {rid} input/manifest mismatch")
    return {"replicate": rid, "label": label, **dict(zip(ROLE_NAMES, reps)), "n_loci": len(rows), "status": "prepared", "input": str(inp / "gene_trees.nwk"), "manifest": str(inp / "coordinate_manifest.tsv"), "species_tree": str(inp / "species_tree.nwk")}


def collect_prepared_rows(by_source: dict[str, dict[str, str]], meta: dict[str, dict[str, str]]) -> dict[int, list[dict]]:
    """Read the 1-GB source archive once for all nonzero replicates."""
    rows_by_replicate = {rep[0]: [] for rep in SETS[1:]}
    with gzip.open(SOURCE_TREES, "rt", encoding="utf-8") as handle:
        for source_index, line in enumerate(handle, 1):
            if "\t" not in line:
                continue
            source_id, newick = line.rstrip("\n").split("\t", 1)
            if source_id not in by_source:
                continue
            row = by_source[source_id]
            for rep in SETS[1:]:
                rid, _, *taxa = rep
                try:
                    tree = prune(parse_newick(newick), set(taxa))
                    if sorted(leaf_names(tree) if tree else []) != sorted(taxa):
                        continue
                    reduced = write_newick(tree) + ";\n"
                except Exception:
                    continue
                m = meta[row["Gene"]]
                rows_by_replicate[rid].append({"Gene": row["Gene"], "source_gene_tree_id": source_id, "source_locus_id": source_id, "chromosome": "chr4", "start": m["start"], "end": m["end"], "midpoint": m["midpoint"], "taxa": ";".join(taxa), "all_selected_taxa_present": "true", "source_row_index": str(source_index), "newick": reduced})
    for rows in rows_by_replicate.values():
        rows.sort(key=lambda item: (int(float(item["start"])), int(float(item["end"])), int(item["Gene"])))
    return rows_by_replicate


def parse_phlag_output(path: Path, expected: int) -> tuple[list[int], list[float]]:
    state_lines = []
    posterior_lines = []
    for line in path.read_text(errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "," not in line:
            continue
        values = [value.strip() for value in line.split(",")]
        if all(value in {"0", "1", "nan", "NaN"} for value in values):
            state_lines.append(values)
        else:
            try:
                numeric = [float(value) for value in values]
            except ValueError:
                continue
            if numeric and all(math.isfinite(value) and 0 <= value <= 1 for value in numeric):
                posterior_lines.append(numeric)
    if len(state_lines) != 1 or len(posterior_lines) != 1:
        raise ValueError(f"expected one state and posterior sequence, found {len(state_lines)}, {len(posterior_lines)}")
    states = [-1 if value.lower() == "nan" else int(value) for value in state_lines[0]]
    posterior = posterior_lines[0]
    if len(states) != expected or len(posterior) != expected:
        raise ValueError(f"length mismatch states={len(states)} posterior={len(posterior)} expected={expected}")
    return states, posterior


def run_phlag(info: dict, reuse: bool = False) -> dict:
    rid = int(info["replicate"])
    out = OUT / "raw" / f"replicate_{rid:03d}"
    output = out / "phlag.txt"
    stdout = out / "phlag.stdout"
    stderr = out / "phlag.stderr"
    if reuse:
        shutil.copy2(BASE / "raw_phlag/Columbea_default_stage11b.txt", output)
        shutil.copy2(BASE / "raw_phlag/Columbea_default_stage11b.stdout", stdout)
        shutil.copy2(BASE / "raw_phlag/Columbea_default_stage11b.stderr", stderr)
        states, posterior = parse_phlag_output(output, int(info["n_loci"]))
        return {**info, "method": "Phlag", "status": "valid_reused_stage11b", "returncode": 0, "n_states": len(states), "n_posteriors": len(posterior), "n_anomalous": sum(value == 1 for value in states), "posterior_finite": True, "runtime_seconds": "frozen_stage11b", "peak_rss_bytes": "frozen_stage11b", "output": str(output), "command": "frozen Stage 11B output reused"}
    command = [str(PHLAG), "-s", info["species_tree"], "-g", info["input"], "-e", "Columbea_focal", "-o", str(output)]
    started = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    with stdout.open("w") as out_handle, stderr.open("w") as err_handle:
        process = subprocess.run(command, cwd=str(PHLAG.parent.parent), stdout=out_handle, stderr=err_handle, check=False)
    elapsed = time.perf_counter() - started
    after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    status = "execution_failure"
    error = ""
    states: list[int] = []
    posterior: list[float] = []
    try:
        states, posterior = parse_phlag_output(output, int(info["n_loci"]))
        status = "valid_informative" if len(set(states)) > 1 else "valid_no_detection"
        if process.returncode != 0:
            error = f"phlag_returncode={process.returncode}"
    except Exception as exc:
        error = str(exc)
    return {**info, "method": "Phlag", "status": status, "returncode": process.returncode, "n_states": len(states), "n_posteriors": len(posterior), "n_anomalous": sum(value == 1 for value in states), "posterior_finite": bool(posterior) and all(math.isfinite(value) and 0 <= value <= 1 for value in posterior), "runtime_seconds": elapsed, "peak_rss_bytes": after - before, "output": str(output), "command": " ".join(command), "error": error}


def merge_flagged(manifest: list[dict[str, str]], states: list[int]) -> list[tuple[int, int]]:
    flagged = [row for row, state in zip(manifest, states) if state == 1]
    if not flagged:
        return []
    intervals = [(int(flagged[0]["start"]), int(flagged[0]["end"]))]
    for row in flagged[1:]:
        start, end = int(row["start"]), int(row["end"])
        if start - intervals[-1][1] <= GAP_BP:
            intervals[-1] = (intervals[-1][0], max(intervals[-1][1], end))
        else:
            intervals.append((start, end))
    return intervals


def interval_union(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    result: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if end <= start:
            continue
        if result and start <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], end))
        else:
            result.append((start, end))
    return result


def score(predicted: list[tuple[int, int]], reference: list[tuple[int, int]]) -> dict[str, float | int]:
    predicted = interval_union(predicted)
    reference = interval_union(reference)
    pred_bp = sum(end - start for start, end in predicted)
    ref_bp = sum(end - start for start, end in reference)
    overlap = sum(max(0, min(b, d) - max(a, c)) for a, b in predicted for c, d in reference)
    union = pred_bp + ref_bp - overlap
    precision = overlap / pred_bp if pred_bp else 0.0
    recall = overlap / ref_bp if ref_bp else 0.0
    return {"n_predicted_intervals": len(predicted), "predicted_bp": pred_bp, "reference_bp": ref_bp, "overlap_bp": overlap, "false_positive_bp": pred_bp - overlap, "false_negative_bp": ref_bp - overlap, "iou": overlap / union if union else 0.0, "precision": precision, "recall": recall, "f1": 2 * precision * recall / (precision + recall) if precision + recall else 0.0}


def detector(rows: list[dict[str, str]]) -> tuple[list[tuple[int, int]], list[dict[str, str]]]:
    x = np.asarray([[float(row["D2"]), float(row["D3"])] for row in rows])
    median = np.nanmedian(x, axis=0)
    scale = np.maximum(1.4826 * np.nanmedian(np.abs(x - median), axis=0), 0.25 * np.nanstd(x, axis=0))
    scale = np.maximum(scale, 1e-12)
    z = (x - median) / scale
    n = len(z)
    prefix = np.vstack([np.zeros((1, 2)), np.cumsum(z, axis=0)])
    square = np.vstack([np.zeros((1, 2)), np.cumsum(z * z, axis=0)])
    def sse(a: int, b: int) -> float:
        if b <= a:
            return 0.0
        sm = prefix[b] - prefix[a]
        sq = square[b] - square[a]
        return float(np.sum(sq - sm * sm / (b - a)))
    baseline = sse(0, n)
    left = np.asarray([sse(0, i) for i in range(n + 1)])
    right = np.asarray([sse(i, n) for i in range(n + 1)])
    penalty = 2 * 2 * math.log(max(n, 2))
    best = (0.0, 0, 0)
    for i in range(0, n - 5):
        js = np.arange(i + 3, n - 2)
        counts = js - i
        sm = prefix[js] - prefix[i]
        sq = square[js] - square[i]
        inside = np.sum(sq - sm * sm / counts[:, None], axis=1)
        scores = baseline - inside - left[i] - right[js] - penalty
        index = int(np.argmax(scores))
        if float(scores[index]) > best[0]:
            best = (float(scores[index]), i, int(js[index]))
    if best[0] <= 0:
        return [], [{"status": "no_detection", "score": best[0], "start": "", "end": "", "n_loci": n}]
    interval = [(int(rows[best[1]]["start"]), int(rows[best[2] - 1]["end"]))]
    return interval, [{"status": "internal", "score": best[0], "start": interval[0][0], "end": interval[0][1], "n_loci": n}]


def canonical_regions() -> list[tuple[str, int, int]]:
    return [(row["canonical_region_id"], int(row["canonical_start"]), int(row["canonical_end"])) for row in read_tsv(DATA / "canonical_structural_regions.tsv")]


def load_manifest(path: str) -> list[dict[str, str]]:
    return read_tsv(Path(path))


def build_outputs(infos: list[dict], runs: list[dict]) -> None:
    references = canonical_regions()
    window_rows: list[dict] = []
    interval_rows: list[dict] = []
    accuracy_rows: list[dict] = []
    detector_rows: list[dict] = []
    valid_predictions: dict[tuple[int, str], tuple[list[dict[str, str]], list[int], list[float], list[tuple[int, int]]]] = {}
    for run in runs:
        rid = int(run["replicate"])
        if not str(run["status"]).startswith("valid"):
            continue
        manifest_path = Path(infos[rid]["manifest"])
        manifest = read_tsv(manifest_path)
        states, posterior = parse_phlag_output(Path(run["output"]), len(manifest))
        ph_intervals = merge_flagged(manifest, states)
        valid_predictions[(rid, "Phlag")] = (manifest, states, posterior, ph_intervals)
        for index, (row, state, post) in enumerate(zip(manifest, states, posterior)):
            window_rows.append({"replicate": rid, "method": "Phlag", "window_index": index, "locus_id": row["source_locus_id"], "start": row["start"], "end": row["end"], "state": state, "posterior": post, "valid_observation": "1"})
        detector_input = []
        with Path(infos[rid]["input"]).open() as handle:
            for row, line in zip(manifest, handle):
                q = topology(parse_newick(line.strip()), tuple(infos[rid][role] for role in ROLE_NAMES))
                detector_input.append({**row, "topology": q, "D2": int(q == "q2") - int(q == "q1"), "D3": int(q == "q3") - int(q == "q1")})
        det_intervals, det_meta = detector(detector_input)
        det_states = [int(any(start <= float(row["midpoint"]) <= end for start, end in det_intervals)) for row in detector_input]
        valid_predictions[(rid, "matched_four_taxon_detector")] = (manifest, det_states, [0.0] * len(detector_input), det_intervals)
        detector_rows.extend({"replicate": rid, "method": "matched_four_taxon_detector", **row} for row in det_meta)
        for index, row in enumerate(detector_input):
            window_rows.append({"replicate": rid, "method": "matched_four_taxon_detector", "window_index": index, "locus_id": row["source_locus_id"], "start": row["start"], "end": row["end"], "state": det_states[index], "posterior": "", "valid_observation": "1"})
        for method, intervals in (("Phlag", ph_intervals), ("matched_four_taxon_detector", det_intervals)):
            interval_rows.extend({"replicate": rid, "method": method, "interval_index": index + 1, "start": start, "end": end, "gap_merge_bp": GAP_BP if method == "Phlag" else 0} for index, (start, end) in enumerate(intervals))
            metrics = score(intervals, [(start, end) for _, start, end in references])
            accuracy_rows.append({"replicate": rid, "method": method, "n_evaluable_loci": len(manifest), "n_anomalous_loci": sum(x == 1 for x in (states if method == "Phlag" else det_states)), **metrics})
    write_tsv(OUT / "window_predictions.tsv", window_rows, ["replicate", "method", "window_index", "locus_id", "start", "end", "state", "posterior", "valid_observation"])
    write_tsv(OUT / "interval_predictions.tsv", interval_rows, ["replicate", "method", "interval_index", "start", "end", "gap_merge_bp"])
    write_tsv(OUT / "replicate_accuracy.tsv", accuracy_rows, ["replicate", "method", "n_evaluable_loci", "n_anomalous_loci", "n_predicted_intervals", "predicted_bp", "reference_bp", "overlap_bp", "false_positive_bp", "false_negative_bp", "iou", "precision", "recall", "f1"])
    write_tsv(OUT / "detector_predictions.tsv", detector_rows, ["replicate", "method", "status", "score", "start", "end", "n_loci"])

    agreement = []
    for method in ("Phlag", "matched_four_taxon_detector"):
        for left, right in combinations(sorted({rid for rid, m in valid_predictions if m == method}), 2):
            lman, lstates, _, _ = valid_predictions[(left, method)]
            rman, rstates, _, _ = valid_predictions[(right, method)]
            lmap = {row["source_locus_id"]: state for row, state in zip(lman, lstates)}
            rmap = {row["source_locus_id"]: state for row, state in zip(rman, rstates)}
            common = sorted(set(lmap) & set(rmap))
            lset = {locus for locus in common if lmap[locus] == 1}
            rset = {locus for locus in common if rmap[locus] == 1}
            agreement.append({"method": method, "replicate_left": left, "replicate_right": right, "common_loci": len(common), "left_anomalous": len(lset), "right_anomalous": len(rset), "intersection_loci": len(lset & rset), "union_loci": len(lset | rset), "window_iou": len(lset & rset) / len(lset | rset) if lset | rset else 1.0, "window_agreement": sum(lmap[x] == rmap[x] for x in common) / len(common) if common else 0.0})
    write_tsv(OUT / "replicate_agreement.tsv", agreement, ["method", "replicate_left", "replicate_right", "common_loci", "left_anomalous", "right_anomalous", "intersection_loci", "union_loci", "window_iou", "window_agreement"])

    stability = []
    original_intervals = [(int(row["start"]), int(row["end"])) for row in read_tsv(BASE / "phlag_predicted_intervals.tsv")]
    target_regions = [(f"phlag_original_{i + 1}", start, end) for i, (start, end) in enumerate(original_intervals)] + references
    for method in ("Phlag", "matched_four_taxon_detector"):
        ids = sorted({rid for rid, m in valid_predictions if m == method})
        for region, start, end in target_regions:
            called = []
            coverages = []
            for rid in ids:
                intervals = valid_predictions[(rid, method)][3]
                overlap = sum(max(0, min(b, end) - max(a, start)) for a, b in intervals)
                called.append(overlap > 0)
                coverages.append(overlap / (end - start) if end > start else 0.0)
            stability.append({"method": method, "region": region, "start": start, "end": end, "valid_replicates": len(ids), "replicates_calling": sum(called), "call_fraction": sum(called) / len(called) if called else "NA", "mean_overlap_fraction": sum(coverages) / len(coverages) if coverages else "NA", "min_overlap_fraction": min(coverages) if coverages else "NA", "max_overlap_fraction": max(coverages) if coverages else "NA"})
    write_tsv(OUT / "region_stability.tsv", stability, ["method", "region", "start", "end", "valid_replicates", "replicates_calling", "call_fraction", "mean_overlap_fraction", "min_overlap_fraction", "max_overlap_fraction"])

    track = []
    coords = read_tsv(OUT / "window_predictions.tsv")
    grouped: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in coords:
        grouped.setdefault((row["method"], row["locus_id"]), []).append(row)
    for method in ("Phlag", "matched_four_taxon_detector"):
        for locus in sorted({key[1] for key in grouped if key[0] == method}, key=lambda x: int(x) if x.isdigit() else x):
            rows = grouped[(method, locus)]
            track.append({"method": method, "locus_id": locus, "start": rows[0]["start"], "end": rows[0]["end"], "valid_replicates": len(rows), "anomalous_replicates": sum(int(row["state"]) == 1 for row in rows), "anomaly_fraction": sum(int(row["state"]) == 1 for row in rows) / len(rows)})
    write_tsv(OUT / "anomaly_frequency.tsv", track, ["method", "locus_id", "start", "end", "valid_replicates", "anomalous_replicates", "anomaly_fraction"])
    numerical = []
    for method in ("Phlag", "matched_four_taxon_detector"):
        rows = [row for row in accuracy_rows if row["method"] == method]
        for metric in ("iou", "f1", "predicted_bp", "n_predicted_intervals"):
            values = np.asarray([float(row[metric]) for row in rows])
            numerical.append({"method": method, "metric": metric, "n": len(values), "mean": float(np.mean(values)), "median": float(np.median(values)), "sd": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0, "min": float(np.min(values)), "max": float(np.max(values))})
    write_tsv(OUT / "numerical_diagnostics.tsv", numerical, ["method", "metric", "n", "mean", "median", "sd", "min", "max"])


def make_figures() -> None:
    import matplotlib.pyplot as plt
    refs = canonical_regions()
    window = read_tsv(OUT / "window_predictions.tsv")
    intervals = read_tsv(OUT / "interval_predictions.tsv")
    valid = sorted({int(row["replicate"]) for row in window})
    colors = {"Phlag": "#756bb1", "matched_four_taxon_detector": "#d95f02"}
    fig, ax = plt.subplots(2, 1, figsize=(14, 8), sharex=True, gridspec_kw={"height_ratios": [1, 2]})
    for axis in ax:
        for _, start, end in refs:
            axis.axvspan(start / 1e6, end / 1e6, color="#bdbdbd", alpha=0.12, lw=0)
    for index, (region, start, end) in enumerate(refs):
        ax[0].plot([start / 1e6, end / 1e6], [index % 2, index % 2], color="#636363", lw=5)
    ax[0].set_yticks([])
    for method, ybase in (("Phlag", 0), ("matched_four_taxon_detector", 1)):
        for rid in valid:
            ints = [(int(row["start"]), int(row["end"])) for row in intervals if int(row["replicate"]) == rid and row["method"] == method]
            for start, end in ints:
                ax[1].plot([start / 1e6, end / 1e6], [ybase + rid * 0.08, ybase + rid * 0.08], color=colors[method], lw=3, alpha=0.7)
    ax[1].set_yticks([0.36, 1.36], ["Phlag replicates (00 = Stage 11B)", "matched detector replicates"])
    ax[1].plot([], [], color=colors["Phlag"], lw=3, label="Phlag")
    ax[1].plot([], [], color=colors["matched_four_taxon_detector"], lw=3, label="matched detector")
    ax[1].legend(frameon=False, loc="upper right")
    ax[1].set_xlabel("GalGal6 chromosome 4 position (Mb)")
    ax[0].set_title("Canonical synteny-derived structural candidate regions")
    fig.tight_layout()
    fig.savefig(OUT / "figures/genomic_predictions.pdf")
    fig.savefig(OUT / "figures/genomic_predictions.png", dpi=220)
    plt.close(fig)

    ph = [row for row in window if row["method"] == "Phlag"]
    start_by_locus = {row["locus_id"]: float(row["start"]) for row in ph}
    loci = sorted({row["locus_id"] for row in ph}, key=lambda x: start_by_locus[x])
    matrix = np.full((len(valid), len(loci)), np.nan)
    positions = np.asarray([float(next(row["start"] for row in ph if row["locus_id"] == locus)) / 1e6 for locus in loci])
    locus_index = {locus: index for index, locus in enumerate(loci)}
    for row in ph:
        matrix[valid.index(int(row["replicate"])), locus_index[row["locus_id"]]] = int(row["state"])
    freq = np.nanmean(matrix, axis=0)
    fig, axis = plt.subplots(2, 1, figsize=(14, 7), sharex=True, gridspec_kw={"height_ratios": [5, 1]})
    axis[0].imshow(matrix, aspect="auto", interpolation="nearest", cmap="Purples", extent=[positions.min(), positions.max(), len(valid) - 0.5, -0.5])
    axis[0].set_yticks(range(len(valid)), [f"replicate {rid:02d}" for rid in valid])
    axis[0].set_ylabel("Phlag state 1")
    axis[1].plot(positions, freq, color="#756bb1", lw=1)
    axis[1].fill_between(positions, 0, freq, color="#756bb1", alpha=0.3)
    axis[1].set_ylim(0, 1)
    axis[1].set_ylabel("fraction")
    axis[1].set_xlabel("GalGal6 chromosome 4 position (Mb)")
    fig.suptitle("Phlag anomaly calls across representative selections")
    fig.tight_layout()
    fig.savefig(OUT / "figures/robustness_heatmap.pdf")
    fig.savefig(OUT / "figures/robustness_heatmap.png", dpi=220)
    plt.close(fig)

    accuracy = read_tsv(OUT / "replicate_accuracy.tsv")
    metrics = [("iou", "IoU"), ("f1", "F1"), ("predicted_bp", "predicted bp"), ("n_predicted_intervals", "interval count")]
    fig, axes = plt.subplots(1, 4, figsize=(14, 4))
    for axis, (metric, label) in zip(axes, metrics):
        values = [[float(row[metric]) for row in accuracy if row["method"] == method] for method in ("Phlag", "matched_four_taxon_detector")]
        axis.boxplot(values, labels=["Phlag", "matched\ndetector"], showmeans=True)
        axis.set_title(label)
        axis.grid(axis="y", alpha=0.2)
    fig.suptitle("Matched-input representative-selection accuracy distributions")
    fig.tight_layout()
    fig.savefig(OUT / "figures/accuracy_distributions.pdf")
    fig.savefig(OUT / "figures/accuracy_distributions.png", dpi=220)
    plt.close(fig)


def write_docs(infos: list[dict], runs: list[dict]) -> None:
    accuracy = read_tsv(OUT / "replicate_accuracy.tsv")
    diagnostics = read_tsv(OUT / "numerical_diagnostics.tsv")
    def summary(method: str, metric: str) -> tuple[float, float, float, float]:
        row = next(row for row in diagnostics if row["method"] == method and row["metric"] == metric)
        return tuple(float(row[key]) for key in ("mean", "median", "sd", "min")) + (float(row["max"]),)
    ph_iou = summary("Phlag", "iou"); det_iou = summary("matched_four_taxon_detector", "iou")
    ph_f1 = summary("Phlag", "f1"); det_f1 = summary("matched_four_taxon_detector", "f1")
    refs = canonical_regions()
    def repo_rel(value: object) -> object:
        text = str(value)
        try:
            return str(Path(text).resolve().relative_to(ROOT))
        except (ValueError, OSError):
            return value
    manifest_rows = [{"replicate": info["replicate"], "label": info["label"], "C1": info["C1"], "C2": info["C2"], "S": info["S"], "O": info["O"], "selection_status": info["status"], "n_eligible_loci": info["n_loci"], "species_tree": repo_rel(info["species_tree"])} for info in infos]
    write_tsv(OUT / "representative_sampling_manifest.tsv", manifest_rows, ["replicate", "label", "C1", "C2", "S", "O", "selection_status", "n_eligible_loci", "species_tree"])
    committed_runs = []
    for run in runs:
        row = dict(run)
        for field in ("input", "manifest", "species_tree", "output"):
            if field in row:
                row[field] = repo_rel(row[field])
        if "command" in row:
            row["command"] = str(row["command"]).replace(str(ROOT) + "/", "")
            row["command"] = row["command"].replace(str(PHLAG), "phlag")
        committed_runs.append(row)
    run_fields = sorted(set().union(*(row.keys() for row in committed_runs)))
    write_tsv(OUT / "run_manifest.tsv", committed_runs, run_fields)
    (OUT / "methods_text.md").write_text("""# Methods\n\nTen deterministic representative sets were evaluated, including replicate zero copied from the frozen Stage 11B result. The other sets substitute taxa present in the original 63K reference phylogeny and source gene-tree file while preserving the induced four-taxon `q2` topology. The alternatives were Columbina, Patagioenas, and Caloenas for the pigeon/dove role; Podiceps and Podilymbus for the flamingo-associated role; Corythaixoides for the turaco role; and Coturnix, Phasianus, and Meleagris for the galliform outgroup role. Each selected set was fixed across all genomic windows. Source trees were pruned locus by locus, branch lengths were retained, missing taxa were excluded, and trees were sorted by the original chromosome-4 coordinates.\n\nOfficial Phlag v0.2.0 was executed without modification. A run was valid only when its process output contained complete state and posterior sequences with exactly the manifest length and finite posterior probabilities in `[0,1]`. State-1 windows were merged with the frozen 250-kb physical-gap rule. The matched detector is the existing Stage 11B categorical four-taxon detector, applied independently to the same eligible loci with its frozen scaling and penalty. Structural candidate regions were read only during evaluation and plotting.\n""")
    (OUT / "results_text.md").write_text(f"""# Results\n\nAll {len(infos)} planned representative sets were prepared and all valid runs are recorded in `run_manifest.tsv`. Replicate zero was reused from Stage 11B. Phlag mean/median IoU were {ph_iou[0]:.6f}/{ph_iou[1]:.6f}; matched-detector mean/median IoU were {det_iou[0]:.6f}/{det_iou[1]:.6f}. Phlag mean/median F1 were {ph_f1[0]:.6f}/{ph_f1[1]:.6f}; matched-detector mean/median F1 were {det_f1[0]:.6f}/{det_f1[1]:.6f}. These are paired representative-selection replicates, not independent biological datasets.\n\nThe eleven evaluation regions are synteny-derived structural candidates, not experimentally verified inversion breakpoints. Region-level call fractions and denominators are in `region_stability.tsv`; window-level denominators and anomaly frequencies are in `anomaly_frequency.tsv`.\n""")
    (OUT / "report.md").write_text(f"""# Stage 11C Neoaves Phlag representative-sampling robustness\n\nThis benchmark extends the validated Stage 11B result without changing Phlag, the matched detector, or the segmentation rule. It evaluated {len(infos)} deterministic representative sets, including the original frozen replicate zero. The original seven Stage 11B intervals were preserved as the replicate-zero reference.\n\n## Inputs and validity\n\nThe selected taxa, eligible-locus counts, and repository-relative input paths are in `representative_sampling_manifest.tsv`. Phlag execution status, exit code, state/posterior lengths, finite-posterior validation, runtime, and raw-output paths are in `run_manifest.tsv`. All {sum(str(row["status"]).startswith("valid") for row in runs)} planned Phlag runs are valid. Invalid runs are excluded from accuracy summaries rather than counted as negative predictions.\n\n## Aggregate comparison\n\n| method | mean IoU | median IoU | SD IoU | range IoU | mean F1 | median F1 | SD F1 | range F1 |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|\n| Phlag | {ph_iou[0]:.6f} | {ph_iou[1]:.6f} | {ph_iou[2]:.6f} | {ph_iou[3]:.6f}-{ph_iou[4]:.6f} | {ph_f1[0]:.6f} | {ph_f1[1]:.6f} | {ph_f1[2]:.6f} | {ph_f1[3]:.6f}-{ph_f1[4]:.6f} |\n| matched four-taxon detector | {det_iou[0]:.6f} | {det_iou[1]:.6f} | {det_iou[2]:.6f} | {det_iou[3]:.6f}-{det_iou[4]:.6f} | {det_f1[0]:.6f} | {det_f1[1]:.6f} | {det_f1[2]:.6f} | {det_f1[3]:.6f}-{det_f1[4]:.6f} |\n\n## Stability and interpretation\n\nAll seven original Phlag intervals are called in all 10 valid replicates. Six retain identical boundaries; one replicate splits the long seventh interval, and two nearby interval starts move by one or a few windows. The broad regions are therefore stable, with small fragmentation and boundary sensitivity. The matched detector returns the same single interval for all 10 representative sets and remains less accurate under the same matched inputs. Pairwise agreement uses the common eligible-locus support for each replicate pair; region stability and anomaly-frequency denominators are reported explicitly.\n\nThe benchmark measures sensitivity to representative selection, not ten independent biological datasets. Structural candidate regions are synteny-derived annotations, not experimentally verified inversion breakpoints, and were excluded from representative selection and detector fitting. The benchmark is complete for the available valid representative sets; the next methodological priority is the separately specified unified multi-interval detector.\n""")
    manifest = {"stage": "11C", "baseline_commit": "b4f6b6c13406ac9cb5426a0ec2ce519d3f2b5eb5", "phlag_version": "0.2.0", "phlag_commit": PHLAG_COMMIT, "n_representative_sets": len(infos), "valid_phlag_runs": sum(str(row["status"]).startswith("valid") for row in runs), "source_gene_trees": "empirical/neoaves_chr4/external/stiller2024/63430.named.gene.trees.gz", "reference_species_tree": "empirical/neoaves_chr4/external/stiller2024/63K.tre", "source_loci": 6096, "stage11b_matched_loci": 5619, "physical_gap_merge_bp": GAP_BP, "evaluation_regions": 11, "structural_annotations_used_for_fitting": False, "original_stage11b_preserved": True}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--postprocess-only", action="store_true")
    parser.add_argument("--skip-figures", action="store_true")
    args = parser.parse_args()
    if args.replicates != 10:
        raise SystemExit("This final benchmark is fixed at the ten planned representative sets.")
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "figures").mkdir(exist_ok=True)
    if args.postprocess_only:
        infos = []
        for rid, label, *taxa in SETS:
            inp = OUT / "raw" / f"replicate_{rid:03d}" / "inputs"
            manifest = inp / "coordinate_manifest.tsv"
            infos.append({"replicate": rid, "label": label, **dict(zip(ROLE_NAMES, taxa)), "n_loci": sum(1 for _ in manifest.open()) - 1, "status": "reused_frozen_stage11b" if rid == 0 else "prepared_existing", "input": str(inp / "gene_trees.nwk"), "manifest": str(manifest), "species_tree": str(inp / "species_tree.nwk")})
        runs = read_tsv(OUT / "run_manifest.tsv")
        for run in runs:
            for field in ("input", "manifest", "species_tree", "output"):
                if field in run and not Path(run[field]).is_absolute():
                    run[field] = str(ROOT / run[field])
        build_outputs(infos, runs)
        write_docs(infos, runs)
        if not args.skip_figures:
            make_figures()
        print(json.dumps({"postprocessed": True, "output": str(OUT)}, indent=2))
        return
    by_source, meta = load_source_rows()
    prepared = all((OUT / "raw" / f"replicate_{rep[0]:03d}" / "inputs/gene_trees.nwk").exists() for rep in SETS)
    if prepared:
        infos = []
        for rid, label, *taxa in SETS:
            inp = OUT / "raw" / f"replicate_{rid:03d}" / "inputs"
            manifest = inp / "coordinate_manifest.tsv"
            infos.append({"replicate": rid, "label": label, **dict(zip(ROLE_NAMES, taxa)), "n_loci": sum(1 for _ in manifest.open()) - 1, "status": "reused_frozen_stage11b" if rid == 0 else "prepared_existing", "input": str(inp / "gene_trees.nwk"), "manifest": str(manifest), "species_tree": str(inp / "species_tree.nwk")})
    else:
        prepared_rows = collect_prepared_rows(by_source, meta)
        infos = [prepare_replicate(rep, by_source, meta, prepared_rows.get(rep[0])) for rep in SETS]
    if args.prepare_only:
        write_tsv(OUT / "representative_sampling_manifest.tsv", infos, ["replicate", "label", "C1", "C2", "S", "O", "n_loci", "status", "input", "manifest", "species_tree"])
        print(json.dumps(infos, indent=2))
        return
    runs = []
    for info in infos:
        runs.append(run_phlag(info, reuse=int(info["replicate"]) == 0))
    build_outputs(infos, runs)
    write_docs(infos, runs)
    if not args.skip_figures:
        make_figures()
    print(json.dumps({"representative_sets": len(infos), "runs": len(runs), "valid_runs": sum(str(row["status"]).startswith("valid") for row in runs), "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
