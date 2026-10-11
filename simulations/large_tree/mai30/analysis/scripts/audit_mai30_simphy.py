#!/usr/bin/env python3
"""Audit the locally recovered Mai D1 SimPhy outputs.

The script uses the files produced by SimPhy itself.  It does not infer
demographic values from condition names.  In particular, parameter.log.info
is cross-checked against the detailed log because the published zero-outgroup
directories omit replicate 001 from that table.
"""

from __future__ import annotations

import csv
import hashlib
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "datasets" / "raw"
OUT = ROOT / "datasets" / "processed"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_newick(text: str) -> tuple[list[str], float, float, float]:
    """Return taxa, root-to-tip minimum, maximum, and root height.

    SimPhy D1 trees have no internal labels.  This small parser intentionally
    accepts only the Newick features needed here and preserves taxon labels.
    """
    tokens = re.findall(r"[(),:;]|[^(),:;\s]+", text.strip())
    pos = 0

    def subtree() -> tuple[list[str], float, float]:
        nonlocal pos
        if tokens[pos] == "(":
            pos += 1
            children = [subtree()]
            while tokens[pos] == ",":
                pos += 1
                children.append(subtree())
            if tokens[pos] != ")":
                raise ValueError("Malformed Newick: expected ')' ")
            pos += 1
            if pos < len(tokens) and tokens[pos] not in ",);:":
                pos += 1  # optional internal label
        else:
            children = [([tokens[pos]], 0.0, 0.0)]
            pos += 1
        length = 0.0
        if pos < len(tokens) and tokens[pos] == ":":
            pos += 1
            length = float(tokens[pos])
            pos += 1
        taxa: list[str] = []
        heights: list[float] = []
        for child_taxa, child_min, child_max in children:
            taxa.extend(child_taxa)
            heights.extend([child_min + length, child_max + length])
        return taxa, min(heights), max(heights)

    taxa, minimum, maximum = subtree()
    if pos < len(tokens) and tokens[pos] == ";":
        pos += 1
    if pos != len(tokens):
        raise ValueError(f"Unparsed Newick tokens: {tokens[pos:]}")
    return taxa, minimum, maximum, maximum


def parse_parameter_table(path: Path) -> dict[str, dict[str, str]]:
    rows: dict[str, dict[str, str]] = {}
    with path.open() as handle:
        next(handle)
        for line in handle:
            fields = line.split()
            if len(fields) != 11:
                continue
            rows[fields[0]] = {
                "speciation_rate": fields[1],
                "extinction_rate": fields[2],
                "outgroup_ratio": fields[3],
                "number_of_leaves": fields[4],
                "generations": fields[5],
                "haploid_effective_population_size": fields[6],
                "global_substitution_rate": fields[7],
                "lineage_gamma_shape": fields[8],
                "gene_family_gamma_shape": fields[9],
                "individuals_per_species": fields[10],
                "parameter_source": str(path),
            }
    return rows


def parse_detailed_log(path: Path) -> dict[str, dict[str, str]]:
    text = path.read_text(errors="replace")
    rows: dict[str, dict[str, str]] = {}
    blocks = re.split(r"(?=^Replicate \d+ parameters:$)", text, flags=re.MULTILINE)
    patterns = {
        "speciation_rate": r"-Speciation rate: ([^\s]+)",
        "extinction_rate": r"-Extinction rate: ([^\s]+)",
        "generations": r"-Generations: ([^\s]+)",
        "haploid_effective_population_size": r"-Haploid efective population size: ([^\s]+)",
        "global_substitution_rate": r"-Global substitution rate: ([^\s]+)",
        "individuals_per_species": r"-Individuals per species: ([^\s]+)",
        "number_of_leaves": r"-Number of leaves: ([^\s]+)",
    }
    for block in blocks:
        match = re.search(r"^Replicate (\d+) parameters:", block, flags=re.MULTILINE)
        if not match:
            continue
        replicate = f"{int(match.group(1)):03d}"
        row = {"parameter_source": str(path)}
        for key, pattern in patterns.items():
            found = re.search(pattern, block)
            if found:
                row[key] = found.group(1)
        outgroup = re.search(r"-Outgroup addition: ([^\n]+)", block)
        row["outgroup_description"] = outgroup.group(1).strip() if outgroup else ""
        ratio = re.search(r"-Ingroup divergence/divergence to the ingroup ratio: ([^\s]+)", block)
        row["outgroup_ratio"] = ratio.group(1) if ratio else ("0" if "No addition" in row["outgroup_description"] else "")
        rows[replicate] = row
    return rows


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    replicate_path = OUT / "mai30_replicate_parameters.tsv"
    condition_path = OUT / "mai30_condition_summary.tsv"
    replicate_fields = [
        "condition", "replicate", "species_tree", "species_tree_sha256",
        "taxon_count", "root_height_generations", "root_to_tip_min_generations",
        "root_to_tip_max_generations", "branch_length_units", "command",
        "command_sha256", "log", "log_sha256", "outgroup_ratio",
        "generations", "haploid_effective_population_size",
        "global_substitution_rate", "lineage_gamma_shape", "gene_family_gamma_shape",
        "individuals_per_species", "parameter_source", "parameter_recovery_note",
    ]
    condition_fields = [
        "condition", "replicate_count", "taxon_count_min", "taxon_count_max",
        "root_height_min_generations", "root_height_max_generations",
        "haploid_ne_min", "haploid_ne_max", "outgroup_ratio", "command",
    ]
    replicate_rows: list[dict[str, str]] = []
    condition_rows: list[dict[str, str]] = []
    for condition_dir in sorted(p for p in RAW.glob("outgroup.*") if p.is_dir()):
        condition = condition_dir.name
        command_path = condition_dir / f"{condition}.command"
        command = command_path.read_text().strip()
        command_digest = sha256(command_path)
        log_path = condition_dir / "log.txt"
        log_digest = sha256(log_path)
        table = parse_parameter_table(condition_dir / "parameter.log.info")
        detailed = parse_detailed_log(condition_dir / "log.txt")
        condition_reps: list[dict[str, str]] = []
        for tree_path in sorted(condition_dir.glob("[0-9][0-9][0-9]/s_tree.trees")):
            replicate = tree_path.parent.name
            tree_text = tree_path.read_text().strip()
            taxa, minimum, maximum, root_height = parse_newick(tree_text)
            row = table.get(replicate, {}).copy()
            recovery_note = "parameter.log.info"
            malformed_table = (
                replicate in table
                and (float(table[replicate].get("number_of_leaves", "0")) <= 0
                     or float(table[replicate].get("generations", "0")) <= 0)
            )
            if (replicate not in table or malformed_table) and replicate in detailed:
                row.update(detailed[replicate])
                recovery_note = "parameter.log.info missing or malformed; recovered from detailed log.txt"
            if "generations" not in row:
                raise ValueError(f"No parameter record for {condition} replicate {replicate}")
            out = {
                "condition": condition,
                "replicate": replicate,
                "species_tree": str(tree_path),
                "species_tree_sha256": sha256(tree_path),
                "taxon_count": str(len(taxa)),
                "root_height_generations": f"{root_height:.10g}",
                "root_to_tip_min_generations": f"{minimum:.10g}",
                "root_to_tip_max_generations": f"{maximum:.10g}",
                "branch_length_units": "generations (SimPhy log)",
                "command": command,
                "command_sha256": command_digest,
                "log": str(log_path),
                "log_sha256": log_digest,
                "outgroup_ratio": row.get("outgroup_ratio", ""),
                "generations": row.get("generations", ""),
                "haploid_effective_population_size": row.get("haploid_effective_population_size", ""),
                "global_substitution_rate": row.get("global_substitution_rate", ""),
                "lineage_gamma_shape": row.get("lineage_gamma_shape", ""),
                "gene_family_gamma_shape": row.get("gene_family_gamma_shape", ""),
                "individuals_per_species": row.get("individuals_per_species", ""),
                "parameter_source": row.get("parameter_source", ""),
                "parameter_recovery_note": recovery_note,
            }
            replicate_rows.append(out)
            condition_reps.append(out)
        condition_rows.append({
            "condition": condition,
            "replicate_count": str(len(condition_reps)),
            "taxon_count_min": str(min(int(r["taxon_count"]) for r in condition_reps)),
            "taxon_count_max": str(max(int(r["taxon_count"]) for r in condition_reps)),
            "root_height_min_generations": f"{min(float(r['root_height_generations']) for r in condition_reps):.10g}",
            "root_height_max_generations": f"{max(float(r['root_height_generations']) for r in condition_reps):.10g}",
            "haploid_ne_min": f"{min(float(r['haploid_effective_population_size']) for r in condition_reps):.10g}",
            "haploid_ne_max": f"{max(float(r['haploid_effective_population_size']) for r in condition_reps):.10g}",
            "outgroup_ratio": sorted({r["outgroup_ratio"] for r in condition_reps})[0],
            "command": command,
        })
    for path, fields, rows in [
        (replicate_path, replicate_fields, replicate_rows),
        (condition_path, condition_fields, condition_rows),
    ]:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    print(f"wrote {replicate_path} ({len(replicate_rows)} replicates)")
    print(f"wrote {condition_path} ({len(condition_rows)} conditions)")


if __name__ == "__main__":
    main()
