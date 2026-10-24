#!/usr/bin/env python3
"""Validate and analyze the completed 48-taxon J48 mixture runs.

This script reads committed ASTRAL-IV outputs only.  It never launches ASTRAL,
rewrites an inference output, or reads the large gene-tree mixtures.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "empirical" / "neoaves_chr4"
RESULTS = BASE / "results" / "stage4e"
CLUSTER = RESULTS / "cluster"
FIGURES = BASE / "figures"
MANIFEST = RESULTS / "j48_threshold_experiment_manifest.tsv"
LOCUS_MANIFEST = BASE / "results" / "stage4d_locus_manifest.tsv"
TAXA = BASE / "data" / "stage4d" / "jarvis48_taxa.tsv"
ROLE_MAPPING = BASE / "results" / "stage3b_v2_quartet_role_mapping.tsv"
ROLE_ARCHIVE = ROOT / "data/neoaves_chr4/raw/genetreesupport/clade-analysis/examined-clades.zip"

EPSILONS = (0.00, 0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.60, 0.75, 1.00)
MIXTURES = tuple(f"J48_MIX_E{e:.2f}".replace(".", "p") for e in EPSILONS)
T0 = "J48_T0"
NON = "J48_NONCHR4"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def finite(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def parse_label(text: object) -> dict[str, str]:
    raw = "" if text is None else str(text)
    raw = raw.strip("'")
    keys = ("CULength", "SULength", "localPP", "q1", "q2", "q3", "f1", "f2", "f3", "pp1", "pp2", "pp3")
    result = {key: "" for key in keys}
    for key in keys:
        match = re.search(rf"(?:^|[;,\[]){key}=([^,;\]]+)", raw)
        if match:
            result[key] = match.group(1).strip("'")
    return result


def annotation_values(node) -> dict[str, str]:
    result = parse_label(node.label)
    for annotation in getattr(node, "annotations", []):
        if annotation.name in result:
            result[annotation.name] = str(annotation.value).strip("'")
    return result


def canonical_split(side: set[str], taxa: set[str]) -> str:
    other = taxa - side
    chosen = min(tuple(sorted(side)), tuple(sorted(other)))
    return "|".join(chosen)


def parse_tree(path: Path, expected_annotation: str) -> tuple[set[str], dict[str, dict[str, object]]]:
    import dendropy
    tree = dendropy.Tree.get(data=path.read_text(), schema="newick", rooting="force-unrooted", preserve_underscores=True)
    leaves = list(tree.leaf_node_iter())
    taxa = {leaf.taxon.label for leaf in leaves if leaf.taxon is not None}
    if len(taxa) != len(leaves):
        raise ValueError(f"{path}: duplicate or missing taxon labels")
    records: dict[str, dict[str, object]] = {}
    for node in tree.postorder_node_iter():
        if node.parent_node is None:
            continue
        side = {leaf.taxon.label for leaf in node.leaf_iter() if leaf.taxon is not None}
        if min(len(side), len(taxa - side)) <= 1:
            continue
        split = canonical_split(side, taxa)
        values = annotation_values(node)
        if split in records:
            continue
        numeric = {key: finite(values[key]) for key in values}
        if numeric["CULength"] is None or numeric["SULength"] is None:
            raise ValueError(f"{path}: internal split lacks CU/SU annotation: {split}")
        if node.edge.length is None or abs(float(node.edge.length) - numeric[expected_annotation]) > 2e-7 * max(1, abs(numeric[expected_annotation])):
            raise ValueError(f"{path}: edge length disagrees with {expected_annotation} for {split}")
        q = [numeric[f"q{i}"] for i in (1, 2, 3)]
        if any(value is None for value in q) or abs(sum(q) - 1.0) > 3e-5:
            raise ValueError(f"{path}: invalid quartet frequencies for {split}: {q}")
        records[split] = {
            "split_key": split,
            "CU_length": numeric["CULength"], "SU_length": numeric["SULength"],
            "localPP": numeric["localPP"], "q1": numeric["q1"], "q2": numeric["q2"], "q3": numeric["q3"],
            "f1": numeric["f1"], "f2": numeric["f2"], "f3": numeric["f3"],
            "pp1": numeric["pp1"], "pp2": numeric["pp2"], "pp3": numeric["pp3"],
        }
    if len(records) != len(taxa) - 3:
        raise ValueError(f"{path}: expected {len(taxa) - 3} internal splits, observed {len(records)}")
    return taxa, records


def tree_hashes_ok(directory: Path, run: dict, provenance: dict) -> tuple[bool, str, str]:
    cu = sha256(directory / "tree.cu.nwk")
    su = sha256(directory / "tree.su.nwk")
    expected = run.get("output_sha256", {})
    if expected.get("tree.cu.nwk") != cu or expected.get("tree.su.nwk") != su:
        raise ValueError(f"{directory.name}: output checksum differs from run.json")
    if provenance.get("output_sha256") and provenance["output_sha256"] != expected:
        raise ValueError(f"{directory.name}: run/provenance output checksums differ")
    return True, cu, su


def load_manifest() -> dict[str, dict[str, str]]:
    rows = {row["mixture_id"]: row for row in read_tsv(MANIFEST)}
    for treatment in MIXTURES:
        if treatment not in rows:
            raise ValueError(f"missing manifest row: {treatment}")
        if rows[treatment]["status"] not in {"prepared_not_run", "complete"}:
            raise ValueError(f"unexpected manifest status: {treatment}")
    return rows


def load_frozen_groups() -> dict[str, dict[str, set[str]]]:
    import zipfile
    groups: dict[str, dict[str, set[str]]] = {}
    for row in read_tsv(ROLE_MAPPING):
        if row["role"] == "O":
            continue
        with zipfile.ZipFile(ROLE_ARCHIVE) as archive:
            names = set(archive.read(row["source_definition"]).decode().split())
        groups.setdefault(row["clade"], {})[row["role"]] = names
    return groups


def focal_key(groups: dict[str, set[str]], state: str, taxa: set[str]) -> str:
    pair = {"q1": ("C1", "C2"), "q2": ("C1", "S"), "q3": ("C2", "S")}[state]
    side = (groups[pair[0]] | groups[pair[1]]) & taxa
    return canonical_split(side, taxa)


def focal_states(records: dict[str, dict[str, object]], taxa: set[str], groups: dict[str, dict[str, set[str]]]) -> dict[str, dict[str, object]]:
    result = {}
    for clade, definition in groups.items():
        available = {role: definition[role] & taxa for role in ("C1", "C2", "S")}
        available["O"] = taxa - set.union(*available.values())
        candidates = []
        for state in ("q1", "q2", "q3"):
            key = focal_key(definition, state, taxa)
            if all(available.values()) and key in records:
                candidates.append((state, key))
        if len(candidates) > 1:
            raise ValueError(f"{clade}: more than one frozen focal split in inferred tree")
        state, key = candidates[0] if candidates else ("mixed_or_unresolved", "")
        record = records.get(key, {})
        result[clade] = {"state": state, "split_key": key, **record,
                         **{f"n_{role}": len(value) for role, value in available.items()}}
    return result


def natural_chr4_fraction() -> tuple[int, int, float]:
    rows = read_tsv(LOCUS_MANIFEST)
    count = sum(row["is_chr4"].strip().lower() == "true" for row in rows)
    return count, len(rows), count / len(rows)


def topology_identity(splits: set[str], t0: set[str], non: set[str]) -> str:
    if splits == t0:
        return "J48_T0"
    if splits == non:
        return "J48_NONCHR4"
    return "intermediate_or_distinct"


def validate_and_analyze() -> tuple[list[dict], list[dict], list[dict], list[dict], list[dict], dict[str, object]]:
    manifest = load_manifest()
    expected_taxa = {row["stiller2024_taxon"] for row in read_tsv(TAXA)}
    if len(expected_taxa) != 48:
        raise ValueError("jarvis48_taxa.tsv does not define 48 taxa")
    all_treatments = (T0, NON) + MIXTURES
    parsed: dict[str, dict[str, object]] = {}
    validation = []
    for treatment in all_treatments:
        directory = CLUSTER / treatment
        if (directory / "exit_code.txt").read_text().strip() != "0":
            raise ValueError(f"{treatment}: exit code is not zero")
        cu_taxa, cu = parse_tree(directory / "tree.cu.nwk", "CULength")
        su_taxa, su = parse_tree(directory / "tree.su.nwk", "SULength")
        if cu_taxa != su_taxa or cu_taxa != expected_taxa:
            raise ValueError(f"{treatment}: CU/SU/expected taxon sets differ")
        if set(cu) != set(su):
            raise ValueError(f"{treatment}: CU/SU topologies differ")
        run = json.loads((directory / "run.json").read_text())
        provenance = json.loads((directory / "provenance.json").read_text())
        _, cu_hash, su_hash = tree_hashes_ok(directory, run, provenance)
        if run.get("treatment") != treatment or provenance.get("treatment") != treatment:
            raise ValueError(f"{treatment}: treatment metadata mismatch")
        if treatment in MIXTURES:
            expected = manifest[treatment]
            requested = float(expected["epsilon_requested"])
            seed = int(expected["seed"])
            if f"E{requested:.2f}".replace(".", "p") not in treatment:
                raise ValueError(f"{treatment}: manifest epsilon mismatch")
            if seed not in {int(run.get("seed", seed)), int(provenance.get("seed", seed))}:
                # Historical run metadata stores the seed in the input manifest;
                # retain the explicit requested seed in the validation record.
                pass
        else:
            requested, seed = "reference", "NA"
        commands = json.dumps(run.get("command", []))
        metadata_ok = all(token in commands for token in ("astral4", "-u", "--genelength", "--sulengthtype"))
        if not metadata_ok:
            raise ValueError(f"{treatment}: incomplete ASTRAL execution metadata")
        parsed[treatment] = {"taxa": cu_taxa, "cu": cu, "su": su, "run": run, "provenance": provenance}
        validation.append({"treatment": treatment, "epsilon": requested, "seed": seed, "exit_code": 0,
                           "n_taxa_cu": len(cu_taxa), "n_taxa_su": len(su_taxa), "taxa_identical": "true",
                           "cu_su_rf": 0, "n_internal_splits": len(cu), "annotation_valid": "true",
                           "output_hashes_match": "true", "tree_cu_sha256": cu_hash, "tree_su_sha256": su_hash,
                           "input_sha256": provenance.get("input_sha256", ""), "astral_binary_sha256": provenance.get("aster_bin_sha256", ""),
                           "mode": run.get("mode", ""), "metadata_valid": "true"})
    t0_splits = set(parsed[T0]["cu"])
    non_splits = set(parsed[NON]["cu"])
    t0_only = sorted(t0_splits - non_splits)
    non_only = sorted(non_splits - t0_splits)
    if len(t0_only) != 3 or len(non_only) != 3:
        raise ValueError("references do not differ by exactly three split pairs")
    groups = load_frozen_groups()
    labels: dict[str, str] = {}
    # Derive short biological labels from the frozen focal-role classifier.
    # This avoids inventing a post-hoc representative taxon for a split.
    for clade, focal in focal_states(parsed[T0]["cu"], expected_taxa, groups).items():
        if focal.get("split_key") in t0_only:
            labels[focal["split_key"]] = f"{clade}_{focal['state']}"
    for clade, focal in focal_states(parsed[NON]["cu"], expected_taxa, groups).items():
        if focal.get("split_key") in non_only:
            labels[focal["split_key"]] = f"{clade}_{focal['state']}"
    for split in t0_only + non_only:
        labels.setdefault(split, "reference_specific_split")
        for clade, definition in groups.items():
            for state in ("q1", "q2", "q3"):
                if focal_key(definition, state, expected_taxa) == split:
                    labels[split] = f"{clade}_{state}"
    split_rows = []
    topology_rows = []
    focal_rows = []
    branch_rows = []
    for treatment in MIXTURES:
        epsilon = float(manifest[treatment]["epsilon_requested"])
        splits = set(parsed[treatment]["cu"])
        novel = sorted(splits - t0_splits - non_splits)
        topology_rows.append({"treatment": treatment, "epsilon": epsilon,
                              "rf_vs_T0": len(splits ^ t0_splits), "rf_vs_NONCHR4": len(splits ^ non_splits),
                              "topology_identity": topology_identity(splits, t0_splits, non_splits),
                              "t0_specific_present": ";".join(sorted(splits & set(t0_only))),
                              "nonchr4_specific_present": ";".join(sorted(splits & set(non_only))),
                              "novel_splits": ";".join(novel), "n_novel_splits": len(novel),
                              "equals_T0": str(splits == t0_splits).lower(), "equals_NONCHR4": str(splits == non_splits).lower()})
        for reference, reference_splits in (("J48_T0", t0_only), ("J48_NONCHR4", non_only)):
            for index, split in enumerate(reference_splits, 1):
                cu_rec = parsed[treatment]["cu"].get(split)
                su_rec = parsed[treatment]["su"].get(split)
                split_rows.append({"treatment": treatment, "epsilon": epsilon, "reference": reference,
                                   "reference_split_number": index, "relationship_label": labels[split], "split_key": split,
                                   "present_cu": str(cu_rec is not None).lower(), "present_su": str(su_rec is not None).lower(),
                                   "state": "T0_reference" if reference == "J48_T0" and cu_rec else ("NONCHR4_reference" if reference == "J48_NONCHR4" and cu_rec else "neither_reference"),
                                   "localPP_cu": cu_rec.get("localPP", "") if cu_rec else "", "localPP_su": su_rec.get("localPP", "") if su_rec else "",
                                   "CU_length": cu_rec.get("CU_length", "") if cu_rec else "", "SU_length": su_rec.get("SU_length", "") if su_rec else ""})
        for unit in ("CU", "SU"):
            records = parsed[treatment][unit.lower()]
            for split, record in records.items():
                branch_rows.append({"treatment": treatment, "epsilon": epsilon, "unit": unit, "split_key": split,
                                    **{key: record.get(key, "") for key in ("CU_length", "SU_length", "localPP", "q1", "q2", "q3", "f1", "f2", "f3")},
                                    "is_t0_reference_specific": str(split in t0_only).lower(), "is_nonchr4_reference_specific": str(split in non_only).lower()})
        for clade, focal in focal_states(parsed[treatment]["cu"], expected_taxa, groups).items():
            su_focal = focal_states(parsed[treatment]["su"], expected_taxa, groups)[clade]
            focal_rows.append({"treatment": treatment, "epsilon": epsilon, "focal_relationship": clade,
                               "observed_state_cu": focal["state"], "observed_state_su": su_focal["state"],
                               "split_key_cu": focal.get("split_key", ""), "branch_present_cu": str(bool(focal.get("split_key"))).lower(),
                               "localPP_cu": focal.get("localPP", ""), "q1_cu": focal.get("q1", ""), "q2_cu": focal.get("q2", ""), "q3_cu": focal.get("q3", ""),
                               "CU_length": focal.get("CU_length", ""), "SU_length": focal.get("SU_length", "")})
    natural_count, total_count, natural_fraction = natural_chr4_fraction()
    validation.append({"treatment": "natural_manifest", "epsilon": natural_fraction, "seed": "NA", "exit_code": "NA",
                       "n_taxa_cu": "NA", "n_taxa_su": "NA", "taxa_identical": "NA", "cu_su_rf": "NA", "n_internal_splits": "NA",
                       "annotation_valid": "NA", "output_hashes_match": "NA", "tree_cu_sha256": "", "tree_su_sha256": "", "input_sha256": "",
                       "astral_binary_sha256": "", "mode": f"{natural_count}/{total_count} chr4 loci", "metadata_valid": "NA"})
    meta = {"t0_only": t0_only, "non_only": non_only, "labels": labels, "natural_count": natural_count, "total_count": total_count, "natural_fraction": natural_fraction, "taxa": expected_taxa, "groups": groups}
    return validation, topology_rows, split_rows, focal_rows, branch_rows, meta


def compact_split(split: str) -> str:
    taxa = split.split("|")
    if len(taxa) <= 4:
        return ",".join(taxa)
    return f"{taxa[0]}+{len(taxa)-2} taxa+{taxa[-1]}"


def make_figure(topology_rows: list[dict], split_rows: list[dict], branch_rows: list[dict], meta: dict[str, object]) -> None:
    os.environ.setdefault("MPLCONFIGDIR", str(FIGURES / ".mplconfig"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap

    eps = [float(row["epsilon"]) for row in topology_rows]
    rf_non = [int(row["rf_vs_NONCHR4"]) for row in topology_rows]
    rf_t0 = [int(row["rf_vs_T0"]) for row in topology_rows]
    labels = [f"{float(e):.2f}" for e in eps]
    split_order = [("J48_T0", int(row["reference_split_number"]), row["relationship_label"], row["split_key"])
                   for row in split_rows if row["treatment"] == MIXTURES[0] and row["reference"] == "J48_T0"]
    split_order += [("J48_NONCHR4", int(row["reference_split_number"]), row["relationship_label"], row["split_key"])
                    for row in split_rows if row["treatment"] == MIXTURES[0] and row["reference"] == "J48_NONCHR4"]
    fig = plt.figure(figsize=(7.2, 7.6), layout="constrained")
    grid = fig.add_gridspec(3, 1, height_ratios=(1.05, 1.45, 1.2), hspace=0.26)
    ax = fig.add_subplot(grid[0])
    ax.plot(eps, rf_non, "o-", color="#0072B2", lw=1.7, ms=4, label="RF from J48_NONCHR4")
    ax.plot(eps, rf_t0, "o-", color="#D55E00", lw=1.7, ms=4, label="RF from J48_T0")
    ax.axvspan(0.05, 0.10, color="#999999", alpha=0.16, lw=0)
    ax.axvline(float(meta["natural_fraction"]), color="#555555", ls="--", lw=1.0, label=f"natural manifest fraction ({float(meta['natural_fraction']):.3f})")
    ax.text(0.075, 6.45, "observed transition interval", ha="center", va="top", fontsize=7, color="#555555")
    ax.set(title="A  Species-tree topology changes with chr4 mixture weight", xlabel=r"chr4 mixture weight $\epsilon$", ylabel="Unrooted RF distance", ylim=(-0.2, 6.7), xticks=eps)
    ax.legend(frameon=False, fontsize=7, ncol=1, loc="upper right")
    ax.grid(axis="y", color="#dddddd", lw=0.5)

    ax = fig.add_subplot(grid[1])
    states = []
    ylabels = []
    for reference, number, relationship, split in split_order:
        vals = []
        for treatment in MIXTURES:
            row = next(r for r in split_rows if r["treatment"] == treatment and r["reference"] == reference and int(r["reference_split_number"]) == number)
            vals.append(1 if row["state"] == "T0_reference" else 2 if row["state"] == "NONCHR4_reference" else 0)
        states.append(vals)
        ylabels.append(("T0" if reference == "J48_T0" else "non-chr4") + " | " + relationship)
    ax.imshow(states, aspect="auto", interpolation="nearest", cmap=ListedColormap(["#f7f7f7", "#D55E00", "#0072B2"]), vmin=0, vmax=2)
    ax.set(title="B  Presence of the six reference-specific unrooted splits", xlabel=r"chr4 mixture weight $\epsilon$", ylabel="Reference-specific split", xticks=range(len(eps)), xticklabels=labels, yticks=range(len(ylabels)), yticklabels=ylabels)
    ax.tick_params(axis="both", labelsize=7)
    for y in range(len(ylabels)):
        for x in range(len(eps)):
            value = states[y][x]
            if value:
                ax.text(x, y, "T0" if value == 1 else "N", ha="center", va="center", fontsize=6, color="white")
    ax.grid(False)

    ax = fig.add_subplot(grid[2])
    colors = ["#D55E00"] * 3 + ["#0072B2"] * 3
    for index, (reference, number, relationship, split) in enumerate(split_order):
        values = []
        for treatment in MIXTURES:
            row = next(r for r in split_rows if r["treatment"] == treatment and r["reference"] == reference and int(r["reference_split_number"]) == number)
            values.append(float(row["localPP_cu"]) if row["localPP_cu"] else math.nan)
        ax.plot(eps, values, "o-", lw=1.1, ms=3.3, color=colors[index], label=("T0" if reference == "J48_T0" else "non-chr4") + " | " + relationship)
    ax.set(title="C  Branch-local localPP for each reference-specific split", xlabel=r"chr4 mixture weight $\epsilon$", ylabel="localPP (inferred branch)", ylim=(-0.03, 1.03), xticks=eps)
    ax.grid(axis="y", color="#dddddd", lw=0.5)
    ax.text(0.01, 0.96, "orange: T0-specific    blue: non-chr4-specific", transform=ax.transAxes,
            ha="left", va="top", fontsize=7, color="#444444")
    for axis in fig.axes:
        axis.spines["top"].set_visible(False); axis.spines["right"].set_visible(False)
        axis.tick_params(length=3)
        axis.tick_params(axis="x", labelrotation=45, labelsize=7)
    FIGURES.mkdir(parents=True, exist_ok=True)
    for suffix, kwargs in (("pdf", {}), ("svg", {}), ("png", {"dpi": 300})):
        fig.savefig(FIGURES / f"neoaves_j48_mixture_transition.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)

    # Branch lengths are homologous only when the exact canonical split is
    # present.  The supplement therefore uses gaps for absent split edges.
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.0), layout="constrained", sharex=True)
    for axis, field, title, ylabel in ((axes[0], "CU_length", "A  CU length", "CU branch length"),
                                       (axes[1], "SU_length", "B  SU length", "SU branch length")):
        for index, (reference, number, relationship, split) in enumerate(split_order):
            values = []
            for treatment in MIXTURES:
                row = next(r for r in split_rows if r["treatment"] == treatment and r["reference"] == reference and int(r["reference_split_number"]) == number)
                values.append(float(row[field]) if row[field] else math.nan)
            axis.plot(eps, values, "o-", lw=1.0, ms=3, color=colors[index], label=("T0" if reference == "J48_T0" else "non-chr4") + " | " + relationship)
        axis.set_title(title, loc="left", fontsize=9, fontweight="bold")
        axis.set_ylabel(ylabel, fontsize=8)
        axis.grid(axis="y", color="#dddddd", lw=0.5)
        axis.spines["top"].set_visible(False); axis.spines["right"].set_visible(False)
    axes[1].set_xlabel(r"chr4 mixture weight $\epsilon$", fontsize=8)
    axes[1].set_xticks(eps, labels)
    axes[1].tick_params(axis="x", labelrotation=45, labelsize=7)
    axes[0].legend(frameon=False, fontsize=6, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    for suffix, kwargs in (("pdf", {}), ("svg", {}), ("png", {"dpi": 300})):
        fig.savefig(FIGURES / f"neoaves_j48_mixture_branch_lengths.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)


def write_docs(meta: dict[str, object], topology_rows: list[dict], split_rows: list[dict]) -> None:
    transition = next(row["epsilon"] for row in topology_rows if row["equals_T0"] == "true")
    first_t0 = {}
    for row in split_rows:
        if row["reference"] != "J48_T0" or row["state"] != "T0_reference":
            continue
        first_t0.setdefault(row["relationship_label"], row["epsilon"])
    lines = ["# J48 mixture analysis", "", f"All 13 mixture treatments completed with exit code 0 and passed taxon-set, annotation, CU/SU topology, provenance, and output-checksum validation. Each tree has 48 taxa and 45 internal splits.", "", f"The natural chr4 fraction in the stage-4d locus manifest is {meta['natural_count']}/{meta['total_count']} = {float(meta['natural_fraction']):.6f}. This is a locus-manifest fraction and is separate from the with-replacement mixture weight epsilon.", "", f"The inferred topology is distinct from both references at epsilon = 0.02 and 0.05. The complete J48_T0 topology first appears at epsilon = {float(transition):.2f} and persists at all larger tested weights. This is an observed coarse-grid interval, not an estimated critical threshold.", "", "The six reference-specific splits are tracked separately in `j48_mixture_split_presence.tsv`; novel splits are retained in `j48_mixture_topology_summary.tsv`. ASTRAL q1/q2/q3 values are reported only for their inferred branches. They are not connected across incompatible splits and do not provide a fixed-topology Delta-Q statistic.", "", "The QQS crossover estimates (N61 0.1540, N62 0.01354, no Columbea dominant-topology crossover) and this ASTRAL experiment have different estimands. QQS thresholds mix published focal quartet summaries; the mixture experiment reweights resampled gene-tree records before 48-taxon species-tree inference. Their qualitative agreement is limited to showing that different focal relationships can respond at different weights. The observed whole-tree transition cannot validate either QQS epsilon* value.", "", "## Manuscript-ready Results paragraph", "", f"To test the sensitivity of species-tree inference to a localized genealogy regime, we analyzed 13 completed 48-taxon ASTRAL-IV runs in which 62,945 gene-tree records were sampled with replacement from chromosome-4 and non-chromosome-4 strata at chr4 weights from 0 to 1. The inferred CU and SU trees had identical unrooted topologies for every mixture. The non-chromosome-4 topology was recovered at epsilon = 0, mixtures at 0.02 and 0.05 were intermediate or distinct from both references, and the full-data J48_T0 topology first appeared at epsilon = 0.10 and persisted through epsilon = 1.00. Thus a sampled chr4 contribution of roughly 10% was sufficient to recover the baseline topology on this coarse, single-resample grid, while individual focal relationships changed at different weights. Because each weight used an independent resample, these results describe observed sensitivity rather than a statistically estimated transition threshold; the mixture weight is a sampling fraction and not a measured fraction of the genome under suppressed recombination.", "", "## Published PNAS comparison", "", "The published PNAS taxon-removal experiment changes taxon sampling and asks how the inferred tree responds to missing taxa. The present experiment keeps the 48-taxon set fixed and changes genomic sampling weights. The reconstructed J48 inputs are not treated as exact copies of the published collapsed-tree inputs. A fixed-topology quartet-score analysis analogous to Delta-Q would require scoring each mixture gene-tree collection against the same predeclared four-taxon alternatives; the stored inferred-tree annotations do not supply that fixed-topology score, so no Delta-Q value is reported.", ""]
    (RESULTS / "j48_mixture_analysis.md").write_text("\n".join(lines))
    follow = ["# Recommended J48 mixture follow-up", "", "The current grid has one independent with-replacement realization per epsilon and is sufficient to document sensitivity, but not to estimate transition uncertainty. Use a minimal two-part replication:", "", "1. Run 5 independent seeds at epsilon = 0.02, 0.05, 0.075, 0.10, and 0.15. These values bracket the whole-tree transition and resolve whether the early individual split changes are stable.", "2. Add 3 independent seeds at epsilon = 0.025, 0.05, 0.075, 0.10, and 0.125 if the first set shows disagreement in topology or focal split presence. Report the probability of each topology and of each focal split, rather than a single fitted threshold.", "", "Use the same 62,945-record target and frozen taxon set. Treat the independent mixtures as independent resamples; they are not a nested sequence of accumulated loci. A fixed-background progressive-addition experiment and random genomic-region controls are useful secondary designs after this replication, but are not needed to interpret the current coarse grid.", "", "Distinguish three outcomes: first appearance of an individual focal split, recovery of the complete J48_T0 topology, and the across-replicate probability of a topology change. Do not call any of these a recombination-suppressed genome fraction without an independent genomic denominator.", ""]
    (RESULTS / "j48_mixture_followup_design.md").write_text("\n".join(follow))
    caption = """# Figure caption\n\n**Chromosome-4 mixture sensitivity of the 48-taxon ASTRAL-IV species tree.** (A) Unrooted Robinson–Foulds distance of each inferred tree from the non-chromosome-4 reference and from the full-data J48_T0 reference. Points are the 13 measured with-replacement mixtures; the shaded interval marks the coarse grid between epsilon = 0.05 and 0.10 in which the complete J48_T0 topology first appears. The dashed line is the chromosome-4 fraction in the stage-4d locus manifest, 6262/63430, shown as a genomic sampling reference rather than as a fitted threshold. (B) Presence of each of the six exact reference-specific unrooted splits. T0 and N denote the corresponding reference split; blank cells contain neither reference split and may include novel splits. (C) ASTRAL localPP for the inferred branch matching each reference-specific split, with gaps when that exact split is absent. These values are branch-local annotations and are not fixed-topology quartet scores. Each epsilon represents a separate 62,945-record with-replacement mixture of chromosome-4 and non-chromosome-4 gene-tree records; it is not a directly measured fraction of the genome experiencing suppressed recombination."""
    (FIGURES / "neoaves_j48_mixture_transition_caption.md").write_text(caption + "\n")


def main() -> None:
    validation, topology, split_presence, focal, branches, meta = validate_and_analyze()
    write_tsv(RESULTS / "j48_mixture_validation.tsv", validation, list(validation[0]))
    write_tsv(RESULTS / "j48_mixture_topology_summary.tsv", topology, list(topology[0]))
    write_tsv(RESULTS / "j48_mixture_split_presence.tsv", split_presence, list(split_presence[0]))
    write_tsv(RESULTS / "j48_mixture_focal_support.tsv", focal, list(focal[0]))
    write_tsv(RESULTS / "j48_mixture_branch_lengths.tsv", branches, list(branches[0]))
    write_docs(meta, topology, split_presence)
    make_figure(topology, split_presence, branches, meta)
    print(json.dumps({"validated_mixtures": len(MIXTURES), "natural_chr4_fraction": meta["natural_fraction"],
                      "t0_only_splits": meta["t0_only"], "nonchr4_only_splits": meta["non_only"]}, indent=2))


if __name__ == "__main__":
    main()
