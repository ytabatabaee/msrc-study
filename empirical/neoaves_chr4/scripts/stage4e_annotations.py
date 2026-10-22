#!/usr/bin/env python3
"""Validated extraction of CU/SU annotations from one ASTRAL-IV tree."""
from __future__ import annotations

import hashlib
import math
import os
import re
import tempfile
import csv
from pathlib import Path

FIELDS = ["treatment", "branch_id", "branch_type", "split_key", "taxon_label",
          "CU_length", "SU_length", "localPP", "q1", "q2", "q3",
          "CU_estimated", "SU_estimated"]
ANNOTATION_KEYS = ("CULength", "SULength", "localPP", "q1", "q2", "q3")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(text)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(name, path)
    finally:
        if os.path.exists(name):
            os.unlink(name)


def finite_number(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _parse_annotation_text(text: object) -> dict[str, str]:
    """Parse ASTRAL's quoted ``[key=value;...]`` node labels.

    DendroPy may expose the surrounding quotes and brackets in ``node.label``;
    the parser deliberately accepts both that form and annotation objects.
    """
    raw = "" if text is None else str(text).strip()
    raw = raw.strip("'")
    if raw.startswith("[") and raw.endswith("]"):
        raw = raw[1:-1]
    result: dict[str, str] = {}
    for part in re.split(r"[;,]", raw):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        key, value = key.strip(), value.strip().strip("'")
        if key in ANNOTATION_KEYS:
            result[key] = value
    return result


def _canonical_split(left: set[str], all_taxa: set[str]) -> str:
    right = all_taxa - left
    if len(left) == 1:
        return "|".join(sorted(left))
    if len(right) == 1:
        return "|".join(sorted(right))
    a = tuple(sorted(left))
    b = tuple(sorted(right))
    chosen = min(a, b)
    return "|".join(chosen)


def _dendropy_tree(text: str):
    import dendropy
    return dendropy.Tree.get(data=text, schema="newick", rooting="force-unrooted",
                             preserve_underscores=True)


def extract_records(text: str, tolerance: float = 1e-8,
                    edge_annotation: str = "SULength") -> tuple[object, set[str], list[dict[str, object]]]:
    if edge_annotation not in {"CULength", "SULength"}:
        raise ValueError("edge_annotation must be CULength or SULength")
    tree = _dendropy_tree(text)
    leaves = list(tree.leaf_node_iter())
    taxa = {node.taxon.label for node in leaves if node.taxon is not None}
    if len(taxa) != len(leaves):
        raise ValueError("tree has missing or duplicate taxon labels")
    records: list[dict[str, object]] = []
    internal_index = 0
    seen_internal: set[str] = set()
    for node in tree.postorder_node_iter():
        if node.parent_node is None:
            continue
        side = {leaf.taxon.label for leaf in node.leaf_iter() if leaf.taxon is not None}
        if not side:
            raise ValueError("branch has no descendant taxa")
        terminal = len(side) == 1 or len(taxa - side) == 1
        branch_type = "terminal" if terminal else "internal"
        split_key = _canonical_split(side, taxa)
        if branch_type == "internal" and split_key in seen_internal:
            # A rooted Newick representation of an unrooted tree may encode
            # the central edge twice when its root has degree two.
            continue
        if branch_type == "internal":
            seen_internal.add(split_key)
        annotation = _parse_annotation_text(node.label)
        for item in getattr(node, "annotations", []):
            if item.name in ANNOTATION_KEYS:
                annotation[item.name] = str(item.value)
        cu = finite_number(annotation.get("CULength"))
        su = finite_number(annotation.get("SULength"))
        edge_length = finite_number(node.edge.length)
        if branch_type == "terminal" and edge_length is not None:
            su = edge_length
        if branch_type == "internal":
            if cu is None or su is None:
                raise ValueError(f"internal branch {split_key} lacks CU/SU annotation")
            expected_edge = cu if edge_annotation == "CULength" else su
            if edge_length is None or abs(edge_length - expected_edge) > tolerance * max(1.0, abs(expected_edge)):
                raise ValueError(f"internal {edge_annotation} disagrees with Newick edge for {split_key}")
            internal_index += 1
            branch_id = f"internal_{internal_index:04d}"
        else:
            taxon = next(iter(side)) if len(side) == 1 else next(iter(taxa - side))
            branch_id = f"terminal_{taxon}"
        records.append({"branch_id": branch_id, "branch_type": branch_type,
                        "split_key": split_key, "taxon_label": branch_id.removeprefix("terminal_"),
                        "CU_length": cu, "SU_length": su,
                        "localPP": finite_number(annotation.get("localPP")),
                        "q1": finite_number(annotation.get("q1")),
                        "q2": finite_number(annotation.get("q2")),
                        "q3": finite_number(annotation.get("q3")),
                        "edge_length": edge_length, "node": node})
    if not records:
        raise ValueError("tree has no branches")
    return tree, taxa, records


def _format(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return format(value, ".12g")
    return str(value)


def write_branch_table(path: Path, treatment: str, records: list[dict[str, object]]) -> None:
    rows = []
    for record in records:
        row = {field: "NA" for field in FIELDS}
        row.update({"treatment": treatment, **{key: _format(record.get(key)) for key in
                  ("branch_id", "branch_type", "split_key", "taxon_label", "CU_length",
                   "SU_length", "localPP", "q1", "q2", "q3")}})
        row["CU_estimated"] = "1" if record["branch_type"] == "internal" and record["CU_length"] is not None else "0"
        row["SU_estimated"] = "1" if record["SU_length"] is not None else "0"
        rows.append(row)
    output = []
    from io import StringIO
    buffer = StringIO()
    writer = csv.DictWriter(buffer, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    atomic_write(path, buffer.getvalue())


def write_derived_trees(path: Path, tree: object, records: list[dict[str, object]]) -> None:
    cu_tree = tree.clone(depth=2)
    su_tree = tree.clone(depth=2)
    # Cloned DendroPy nodes do not preserve Python identity, so match branches
    # in traversal order, which is preserved by clone(depth=2).
    original = [node for node in tree.postorder_node_iter() if node.parent_node is not None]
    record_map = {(record["branch_type"], record["split_key"]): record for record in records}
    all_taxa = {leaf.taxon.label for leaf in tree.leaf_node_iter()}
    for derived in (cu_tree, su_tree):
        derived_nodes = [node for node in derived.postorder_node_iter() if node.parent_node is not None]
        for source_node, derived_node in zip(original, derived_nodes):
            side = {leaf.taxon.label for leaf in source_node.leaf_iter()}
            branch_type = "terminal" if len(side) == 1 or len(all_taxa - side) == 1 else "internal"
            split_key = _canonical_split(side, all_taxa)
            record = record_map[(branch_type, split_key)]
            if derived is cu_tree:
                derived_node.edge.length = record["CU_length"] if record["branch_type"] == "internal" else None
            else:
                derived_node.edge.length = record["SU_length"] if record["branch_type"] == "internal" else record["edge_length"]
        if derived is cu_tree:
            atomic_write(path.with_name("tree.cu.nwk"), derived.as_string(schema="newick", suppress_rooting=True))
        else:
            atomic_write(path.with_name("tree.su.nwk"), derived.as_string(schema="newick", suppress_rooting=True))


def extract_and_write(input_path: Path, output_dir: Path, treatment: str,
                      source_kind: str = "single_run", tolerance: float = 1e-8) -> dict[str, object]:
    source_text = input_path.read_text()
    tree, taxa, records = extract_records(source_text, tolerance=tolerance)
    output_dir.mkdir(parents=True, exist_ok=True)
    annotated = output_dir / "tree.annotated.nwk"
    if annotated.exists() and annotated.read_text() != source_text:
        raise FileExistsError(f"refusing to overwrite differing annotated tree: {annotated}")
    if not annotated.exists():
        atomic_write(annotated, source_text)
    write_derived_trees(annotated, tree, records)
    write_branch_table(output_dir / "branch_lengths.tsv", treatment, records)
    return {"taxa": len(taxa), "branches": len(records), "source_sha256": sha256(input_path),
            "annotated_sha256": sha256(annotated), "source_kind": source_kind,
            "internal_branches": sum(r["branch_type"] == "internal" for r in records),
            "terminal_branches": sum(r["branch_type"] == "terminal" for r in records)}
