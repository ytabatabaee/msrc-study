#!/usr/bin/env python3
"""Stage 0 audit for the deposited Microbotryum gene-tree resources.

The script is deliberately dependency-free. It audits files when present and
emits explicit unavailable statuses when Dryad files could not be retrieved;
it never fabricates tree, coordinate, or stratum results.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import re
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data/microbotryum/raw"
OUT = ROOT / "empirical/microbotryum/results"
META = ROOT / "data/microbotryum/metadata"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def split_newicks(text: str) -> list[str]:
    """Return semicolon-terminated Newick records, respecting quoted text."""
    records, start, quote = [], None, None
    for i, ch in enumerate(text):
        if ch in "'\"":
            quote = None if quote == ch else (ch if quote is None else quote)
        if ch == ";" and quote is None:
            if start is not None:
                records.append(text[start:i + 1].strip())
            start = None
        elif start is None and ch in "(>":
            start = i
    return [r for r in records if "(" in r and ")" in r]


def parse_newick(text: str) -> tuple[set[str], bool, bool]:
    """Validate parentheses and return tips, branch-length, resolution flags."""
    depth = 0
    tips: list[str] = []
    has_lengths = bool(re.search(r"\)[^,();:]*:[-+0-9.eE]+|[),][^,();:]*:[-+0-9.eE]+", text))
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth < 0:
                raise ValueError("unbalanced closing parenthesis")
    if depth != 0 or not text.endswith(";"):
        raise ValueError("unbalanced or unterminated Newick")
    body = text[:-1]
    for token in re.split(r"[(),]", body):
        token = token.strip()
        if not token:
            continue
        token = token.split(":", 1)[0].strip()
        if token and not re.fullmatch(r"[-+0-9.eE]+", token):
            tips.append(token.strip("'\""))
    unique = set(tips)
    if len(unique) != len(tips):
        raise ValueError("duplicate tip labels")
    return unique, has_lengths, text.count("(") > 1


def classify_quartet(split: str, quartet: tuple[str, str, str, str]) -> str:
    """Classify a split represented as A|B against fixed q1/q2/q3 order."""
    groups = [frozenset(x.split(",")) for x in split.split("|")]
    q = tuple(quartet)
    expected = [
        {frozenset((q[0], q[1])), frozenset((q[2], q[3]))},
        {frozenset((q[0], q[2])), frozenset((q[1], q[3]))},
        {frozenset((q[0], q[3])), frozenset((q[1], q[2]))},
    ]
    observed = set(groups)
    return {frozenset(expected[i]): f"q{i + 1}" for i in range(3)}.get(frozenset(observed), "unresolved")


def read_stratum_tokens(raw_dir: Path) -> dict[str, int]:
    """Count strata mentioned in deposited/source-side metadata if available."""
    counts = Counter()
    for path in [raw_dir / "README.md", raw_dir / "Orthogroups_filtered.txt"]:
        if not path.exists():
            continue
        text = path.read_text(errors="replace")
        for label in re.findall(r"\b(?:stratum|strata)[_ -]*([A-Za-z0-9]+)", text, flags=re.I):
            counts[label] += 1
    return dict(counts)


def audit_trees(raw_dir: Path) -> list[dict[str, str]]:
    path = raw_dir / "gene_trees"
    if not path.exists():
        return [{"resource": "gene_trees", "status": "not_downloaded", "n_records": "NA", "n_taxa": "NA", "branch_lengths": "NA", "resolved": "NA", "malformed": "NA", "stratum": "NA", "notes": "Dryad file endpoint returned 401/403"}]
    records = split_newicks(path.read_text(errors="replace"))
    rows = []
    for i, tree in enumerate(records, 1):
        try:
            taxa, lengths, resolved = parse_newick(tree)
            rows.append({"resource": "gene_trees", "status": "usable", "n_records": str(len(records)), "n_taxa": str(len(taxa)), "branch_lengths": str(lengths).lower(), "resolved": str(resolved).lower(), "malformed": "false", "stratum": "unassigned", "notes": f"record={i}; taxa={','.join(sorted(taxa))}"})
        except Exception as exc:
            rows.append({"resource": "gene_trees", "status": "malformed", "n_records": str(len(records)), "n_taxa": "NA", "branch_lengths": "NA", "resolved": "NA", "malformed": "true", "stratum": "unassigned", "notes": f"record={i}; {exc}"})
    return rows or [{"resource": "gene_trees", "status": "no_newick_records", "n_records": "0", "n_taxa": "NA", "branch_lengths": "NA", "resolved": "NA", "malformed": "NA", "stratum": "NA", "notes": "No semicolon-terminated Newick records found"}]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, default=RAW)
    args = parser.parse_args()
    raw_dir = args.raw_dir
    OUT.mkdir(parents=True, exist_ok=True)
    META.mkdir(parents=True, exist_ok=True)

    rows = audit_trees(raw_dir)
    with (OUT / "gene_tree_audit.tsv").open("w", newline="") as fh:
        fields = list(rows[0])
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)

    with (OUT / "stratum_tree_counts.tsv").open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["stratum", "n_trees", "status"])
        writer.writerow(["unassigned", "NA", "not_recoverable_without_gene-tree/stratum mapping"])

    with (OUT / "coordinate_mapping_audit.tsv").open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["resource", "identifier", "chromosome_or_scaffold", "start", "end", "gene_order", "stratum", "status", "evidence"])
        writer.writerow(["gene_trees", "orthogroup identifiers", "NA", "NA", "NA", "unknown", "unknown", "not_verified", "Dryad gene_trees unavailable; source pipeline maps OG to stratum but not deposited physical coordinates"])

    with (OUT / "preliminary_quartet_support.tsv").open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["stratum", "quartet", "q1_species", "q2_mating_type", "q3_third", "n_usable_trees", "status"])
        writer.writerow(["NA", "S1_a1|S1_a2|S2_a1|S2_a2", "NA", "NA", "NA", "0", "not_run_gene_trees_unavailable"])

    files = []
    for p in sorted(raw_dir.iterdir()) if raw_dir.exists() else []:
        if p.is_file():
            files.append((p.name, p.stat().st_size, sha256(p)))
    with (META / "source_manifest.tsv").open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["source", "doi_or_url", "file", "size_bytes", "sha256", "status", "notes"])
        writer.writerow(["Lucotte et al. 2025", "https://doi.org/10.1038/s41467-025-60222-5", "published study", "NA", "NA", "verified_metadata", "gene trees/species tree deposited in Dryad"])
        writer.writerow(["Dryad", "https://doi.org/10.5061/dryad.0k6djhb89", "gene_trees", "59490", "7ca1c37b8dd6deeda786e57f42a3041438f3bd248437b4c3968317f60b10c968", "metadata_verified_download_blocked", "file endpoint returned HTTP 401/403 in this environment"])
        writer.writerow(["Dryad", "https://doi.org/10.5061/dryad.0k6djhb89", "speciestree.tre", "1619", "ee253958e92034d7fb50c3284a7ab607a8e7c52307c37dda187675d21a4db722", "metadata_verified_download_blocked", "file endpoint returned HTTP 401/403 in this environment"])
        for name, size, digest in files:
            writer.writerow(["local", "", name, size, digest, "downloaded", "source-side audit input"])

    print(f"gene-tree audit: {OUT / 'gene_tree_audit.tsv'}")
    print("Dryad gene_trees available:", (RAW / "gene_trees").exists())


if __name__ == "__main__":
    main()
