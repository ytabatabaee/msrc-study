#!/usr/bin/env python3
"""Run ASTER/ASTRAL-IV and calculate topology and CULength summaries."""

from __future__ import annotations

import csv
import json
import re
import subprocess
from pathlib import Path

import numpy as np

SCRIPT = Path(__file__).resolve()
STUDY_ROOT = SCRIPT.parents[2]
PILOT_ROOT = STUDY_ROOT / "datasets/processed/mai30_scaled_wf_pilot"
RESULTS = STUDY_ROOT / "analysis/results/mai30_scaled_wf_astral.tsv"
ASTRAL = Path("/Users/ytabatabaee/Desktop/ASTER/bin/astral")
sys_taxa = None


def canonical(side: frozenset[str], taxa: frozenset[str]) -> frozenset[str] | None:
    other = taxa - side
    if len(side) < 2 or len(other) < 2:
        return None
    return side if tuple(sorted(side)) <= tuple(sorted(other)) else other


def parse_with_lengths(text: str) -> tuple[frozenset[str], frozenset[frozenset[str]], dict[frozenset[str], float]]:
    tokens = re.findall(r"[(),:;]|[^(),:;\s]+", text.strip())
    pos = 0

    def subtree():
        nonlocal pos
        if tokens[pos] == "(":
            pos += 1
            children = [subtree()]
            while tokens[pos] == ",":
                pos += 1
                children.append(subtree())
            if tokens[pos] != ")":
                raise ValueError("malformed Newick")
            pos += 1
            if pos < len(tokens) and tokens[pos] not in ",);:":
                pos += 1
            internal = True
        else:
            name = tokens[pos]
            if name in ",);:":
                raise ValueError("missing taxon")
            pos += 1
            children = [(frozenset([name]), [], {}, False)]
            internal = False
        length = 0.0
        if pos < len(tokens) and tokens[pos] == ":":
            pos += 1
            length = float(tokens[pos])
            if not np.isfinite(length) or length < 0:
                raise ValueError("invalid branch length")
            pos += 1
        taxa = frozenset().union(*(x[0] for x in children))
        splits = frozenset().union(*(x[1] for x in children))
        lengths = {}
        for child_taxa, child_splits, child_lengths, _ in children:
            splits = splits | child_splits
            lengths.update(child_lengths)
            split = canonical(child_taxa, taxa | frozenset())
            if split is not None:
                lengths[split] = length
        if internal:
            splits = splits | frozenset({taxa})
        return taxa, splits, lengths, internal

    taxa, _splits, lengths, _ = subtree()
    if pos < len(tokens) and tokens[pos] == ";":
        pos += 1
    if pos != len(tokens):
        raise ValueError(f"unparsed Newick tokens: {tokens[pos:]}")
    # Rebuild canonical splits with the complete taxon set; the recursive
    # helper above stores child sides, and this second pass handles the root.
    return taxa, frozenset(s for s in lengths if s is not None), lengths


def parse_lengths_correct(text: str):
    """Parse clade edge lengths and canonicalize them against all taxa."""
    tokens = re.findall(r"[(),:;]|[^(),:;\s]+", text.strip())
    pos = 0

    def rec():
        nonlocal pos
        if tokens[pos] == "(":
            pos += 1
            children = [rec()]
            while tokens[pos] == ",":
                pos += 1
                children.append(rec())
            if tokens[pos] != ")":
                raise ValueError("malformed Newick")
            pos += 1
            if pos < len(tokens) and tokens[pos] not in ",);:":
                pos += 1
            internal = True
        else:
            name = tokens[pos]
            pos += 1
            children = []
            internal = False
            taxa = frozenset([name])
        length = 0.0
        if pos < len(tokens) and tokens[pos] == ":":
            pos += 1
            length = float(tokens[pos])
            pos += 1
        if internal:
            taxa = frozenset().union(*(child[0] for child in children))
        return taxa, children, length

    root = rec()
    if pos < len(tokens) and tokens[pos] == ";":
        pos += 1
    taxa = root[0]
    lengths = {}

    def collect(node):
        node_taxa, children, length = node
        if children:
            split = canonical(node_taxa, taxa)
            if split is not None:
                lengths[split] = length
            for child in children:
                collect(child)

    collect(root)
    splits = frozenset(lengths)
    return taxa, splits, lengths


def normalized_rf(inferred: frozenset[frozenset[str]], truth: frozenset[frozenset[str]], n: int) -> float:
    return len(inferred ^ truth) / (2 * (n - 3))


def main() -> None:
    rows = []
    processed_root = STUDY_ROOT / "datasets/processed"
    collections = []
    for collection in sorted(PILOT_ROOT.glob("*/collection_summary.json")):
        collections.append((json.loads(collection.read_text()), collection.parent, collection.parent / "gene_trees.nwk", collection.parent / "astral_species_tree.nwk"))
    controls_root = STUDY_ROOT / "datasets/processed/mai30_scaled_msc_controls"
    for control_file in sorted(controls_root.glob("*/*/msc_gene_trees.nwk")):
        condition = control_file.parent.parent.name
        rep = control_file.parent.name
        pilot_group = next((p for p in PILOT_ROOT.glob(f"{condition}__replicate_{rep}__*__freq_0.5") if p.is_dir()), None)
        if pilot_group is None:
            continue
        source_summary = json.loads((pilot_group / "collection_summary.json").read_text())
        summary = {"group": f"{condition}__replicate_{rep}__msc_control", "n_independent_histories": "direct_msc", "n_gene_trees": 500}
        collections.append((summary, control_file.parent, control_file, control_file.parent / "astral_species_tree.nwk"))
    for summary, group_dir, gene_file, output in collections:
        command = [str(ASTRAL), "-i", str(gene_file), "-o", str(output), "-t", "2", "-u", "1", "--length", "CULength"]
        completed = subprocess.run(command, cwd=STUDY_ROOT, capture_output=True, text=True, check=False)
        (group_dir / "astral.log").write_text(completed.stdout + completed.stderr)
        if completed.returncode != 0:
            raise RuntimeError(f"ASTRAL failed for {group_dir}: {completed.stderr[-1000:]}")
        inferred_taxa, inferred_splits, inferred_lengths = parse_lengths_correct(output.read_text())
        if summary["group"].endswith("__msc_control"):
            condition, rep = summary["group"].rsplit("__replicate_", 1)
            rep = rep.replace("__msc_control", "")
            pilot_group = next(p for p in PILOT_ROOT.glob(f"{condition}__replicate_{rep}__*__freq_0.5") if p.is_dir())
            source = processed_root / "mai30_scaled_wf_runs" / pilot_group.name / "history_000" / "species_tree.nwk"
        else:
            source = processed_root / "mai30_scaled_wf_runs" / summary["group"] / "history_000" / "species_tree.nwk"
        true_taxa, true_splits, true_lengths = parse_lengths_correct(source.read_text())
        if inferred_taxa != true_taxa:
            raise ValueError(f"ASTRAL taxon set mismatch: {group_dir}")
        n = len(true_taxa)
        # Recover N_scaled from the resolved history configuration.
        cfg_group = summary["group"]
        if cfg_group.endswith("__msc_control"):
            cfg_group = pilot_group.name
        cfg = json.loads((processed_root / "mai30_scaled_wf_runs" / cfg_group / "history_000" / "simulation_config.json").read_text())
        n_scaled = int(cfg["species_tree"]["default_effective_population_size"])
        true_cu = {split: length / (2.0 * n_scaled) for split, length in true_lengths.items() if split in true_splits}
        common = sorted(set(inferred_lengths) & set(true_cu), key=lambda x: tuple(sorted(x)))
        errors = [inferred_lengths[s] - true_cu[s] for s in common]
        abs_errors = [abs(x) for x in errors]
        rel_errors = [abs(e) / abs(true_cu[s]) for e, s in zip(errors, common) if abs(true_cu[s]) > 1e-12]
        rows.append({
            "group": summary["group"],
            "n_independent_histories": summary["n_independent_histories"],
            "n_gene_trees": summary["n_gene_trees"],
            "n_taxa": n,
            "normalized_rf": normalized_rf(inferred_splits, true_splits, n),
            "exact_topology_recovery": int(inferred_splits == true_splits),
            "common_internal_splits": len(common),
            "mean_absolute_cu_error_common_splits": float(np.mean(abs_errors)) if abs_errors else "",
            "mean_cu_bias_common_splits": float(np.mean(errors)) if errors else "",
            "mean_relative_cu_error_common_splits": float(np.mean(rel_errors)) if rel_errors else "",
            "astral_output": str(output),
            "astral_command": " ".join(command),
        })
    if rows:
        RESULTS.parent.mkdir(parents=True, exist_ok=True)
        with RESULTS.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        print(f"wrote {RESULTS}")


if __name__ == "__main__":
    main()
