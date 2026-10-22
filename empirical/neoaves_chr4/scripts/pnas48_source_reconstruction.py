#!/usr/bin/env python3
"""Audit the published 63K collapsed archive without inventing locus IDs.

The collapsed Stiller archive contains one bare Newick tree per line.  The
published 63K header supplies the genomic locus for each Gene index, while
the named archive supplies the same locus names explicitly.  This script
checks that relationship and writes a mapping audit.  It deliberately does
not create a claimed exact PNAS input when any source records remain
unresolved.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "empirical/neoaves_chr4"
EXT = BASE / "external"
STILLER = EXT / "stiller2024"
PNAS = EXT / "chr4avian"
RESULTS = BASE / "results/stage4e"
OUTDIR = BASE / "data/stage4e_pnas48"

EXPECTED = {
    "named_sha256": "cbf5a4f631fa92ef6da65c48a5aaf39b424fce94d2812c526eaed904dc61392b",
    "collapsed_sha256": "5b058690d31531b3eb4f287556597bf5aa7db380ddf74c84bc53f71b512497df",
    "header_rows": 63430,
    "pnas_excluded": 1431,
}
LEAF = re.compile(r"(?:\(|,)([^(),:]+):")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def header_rows() -> list[dict[str, str]]:
    archive = PNAS / "chr4avian-v1.0.0.zip"
    with zipfile.ZipFile(archive) as zf:
        member = next(n for n in zf.namelist() if n.endswith("genetreesupport/63K_trees.names_header.txt.xz"))
        raw = subprocess.check_output(["xz", "-cd"], input=zf.read(member))
    rows = []
    for line in raw.decode().splitlines()[1:]:
        gene, chrom, start, end, size, offset = line.split()
        rows.append({
            "gene": gene,
            "locus_id": f"{chrom}_{start}_{end}_{size}_start{offset}_fasta.treefile",
            "chromosome": chrom,
            "start": start,
            "end": end,
            "offset": offset,
        })
    return rows


def leaf_digest(tree: str) -> tuple[int, str]:
    taxa = sorted(LEAF.findall(tree))
    return len(taxa), hashlib.sha256("\0".join(taxa).encode()).hexdigest()


def named_index() -> dict[str, tuple[int, str]]:
    result: dict[str, tuple[int, str]] = {}
    with gzip.open(STILLER / "63430.named.gene.trees.gz", "rt") as handle:
        for line in handle:
            locus, tree = line.rstrip("\n").split("\t", 1)
            if locus in result:
                raise ValueError(f"duplicate named locus: {locus}")
            result[locus] = leaf_digest(tree)
    return result


def source_rows() -> list[dict[str, object]]:
    header = header_rows()
    named = named_index()
    if len(header) != EXPECTED["header_rows"]:
        raise ValueError(f"header row count {len(header)} != {EXPECTED['header_rows']}")
    if len(named) != EXPECTED["header_rows"]:
        raise ValueError(f"named row count {len(named)} != {EXPECTED['header_rows']}")
    header_ids = {r["locus_id"] for r in header}
    if header_ids != set(named):
        raise ValueError("header and named locus IDs differ")

    rows: list[dict[str, object]] = []
    with gzip.open(STILLER / "63430.aLRT-0.95-collapsed.gene.trees.gz", "rt") as handle:
        for i, line in enumerate(handle):
            meta = header[i]
            collapsed_n, collapsed_digest = leaf_digest(line.rstrip("\n"))
            named_n, named_digest = named[meta["locus_id"]]
            rows.append({
                "collapsed_record": i + 1,
                "gene": meta["gene"],
                "locus_id": meta["locus_id"],
                "chromosome": meta["chromosome"],
                "start": meta["start"],
                "end": meta["end"],
                "collapsed_n_taxa": collapsed_n,
                "named_n_taxa": named_n,
                "taxon_set_status": "match" if (collapsed_n, collapsed_digest) == (named_n, named_digest) else "unresolved",
                "mapping_status": "strongly_supported" if (collapsed_n, collapsed_digest) == (named_n, named_digest) else "unverified",
            })
    if len(rows) != EXPECTED["header_rows"]:
        raise ValueError(f"collapsed row count {len(rows)} != {EXPECTED['header_rows']}")
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def existing_j48_rows() -> list[dict[str, str]]:
    path = BASE / "results/stage4e_input_manifest.tsv"
    with path.open() as handle:
        return [r for r in csv.DictReader(handle, delimiter="\t") if r["scope"] == "jarvis48"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write the audit and reports")
    args = parser.parse_args()

    named_path = STILLER / "63430.named.gene.trees.gz"
    collapsed_path = STILLER / "63430.aLRT-0.95-collapsed.gene.trees.gz"
    if sha256(named_path) != EXPECTED["named_sha256"]:
        raise ValueError("named archive checksum changed")
    if sha256(collapsed_path) != EXPECTED["collapsed_sha256"]:
        raise ValueError("collapsed archive checksum changed")

    rows = source_rows()
    counts = {
        "n_rows": len(rows),
        "n_unique_locus_ids": len({r["locus_id"] for r in rows}),
        "taxon_set_matches": sum(r["taxon_set_status"] == "match" for r in rows),
        "taxon_set_unresolved": sum(r["taxon_set_status"] != "match" for r in rows),
        "published_outlier_loci": sum(
            r["chromosome"] == "chr4" and (
                25030000 <= int(r["start"]) <= 32670000
                or 33510000 <= int(r["start"]) <= 34470000
                or 44130000 <= int(r["start"]) <= 56810000
            )
            for r in rows
        ),
    }
    if counts["n_unique_locus_ids"] != EXPECTED["header_rows"]:
        raise ValueError("mapped locus IDs are not unique")
    if counts["published_outlier_loci"] != EXPECTED["pnas_excluded"]:
        raise ValueError("published outlier count changed")

    old = existing_j48_rows()
    mapping_audit = [
        {"component": "collapsed_archive_record_count", "status": "verified", "value": counts["n_rows"], "evidence": "local archive line count"},
        {"component": "header_gene_to_locus_mapping", "status": "strongly_supported", "value": counts["n_unique_locus_ids"], "evidence": "63K_trees.names_header.txt.xz plus exact named-archive locus set"},
        {"component": "collapsed_to_locus_taxon_set_mapping", "status": "strongly_supported_with_unresolved_rows", "value": f"{counts['taxon_set_matches']} matched; {counts['taxon_set_unresolved']} unresolved", "evidence": "collapsed Gene order compared with named tree at header-derived locus"},
        {"component": "published_outlier_mask", "status": "verified_for_header_loci", "value": counts["published_outlier_loci"], "evidence": "delta_quartet.R inclusive GalGal4 window-start predicate"},
        {"component": "exact_published_48_input", "status": "unverified", "value": "not constructed", "evidence": "317 collapsed records lack a matching named taxon set; source does not provide IDs in the collapsed file"},
        {"component": "existing_j48_inputs", "status": "preserved", "value": len(old), "evidence": "stage4e_input_manifest.tsv"},
    ]
    comparisons = []
    for row in old:
        comparisons.append({
            "treatment": row["treatment"],
            "old_input": row["relative_path"],
            "old_loci_or_trees": row["n_loci"],
            "old_sha256": row["sha256"],
            "corrected_input": "not_constructed_unverified_mapping",
            "comparison_status": "preserved_old_input; exact_published_counterpart_unavailable",
            "difference_category": "published_collapsed_source_identity_unresolved",
            "taxon_sets": "old_48_taxon_pruned; published_unresolved",
            "locus_sets": "old_named_manifest; published_317_unresolved",
            "gene_tree_counts": "old_recorded; published_exact_count_unresolved",
            "resolutions": "not_comparable_without_exact_locus_pairing",
            "branch_lengths": "not_comparable_without_exact_locus_pairing",
            "support_annotations": "published_collapsed_support_not_used_in_old_input",
            "genomic_filtering": "old_mask_verified; published_unresolved_rows_not_selected",
        })

    if args.write:
        write_tsv(RESULTS / "pnas48_locus_mapping_audit.tsv", rows)
        write_tsv(RESULTS / "pnas48_dataset_comparison.tsv", comparisons)
        write_tsv(RESULTS / "pnas48_mapping_summary.tsv", mapping_audit)
        report = f"""# PNAS 48-taxon source reconstruction

Date: 2026-10-22. No ASTRAL-IV inference or SLURM job was run.

## Determination

The original collapsed archive can be assigned to genomic loci with a strong
but incomplete provenance chain. The 63,430-row published header supplies the
`Gene` index and genomic coordinates. Its 63,430 derived locus IDs are all
present exactly once in `63430.named.gene.trees.gz`. Comparing collapsed record
`i` with the named tree for header `Gene=i` gives identical taxon sets for
**{counts['taxon_set_matches']:,}** records. This explains why the previous
same-row comparison was misleading: the named archive changes order after the
first block of loci, producing 31,066 apparent disagreements when the two
files are compared by physical line number.

The remaining **{counts['taxon_set_unresolved']}** records have different taxon
sets from their header-derived named locus. Their identities are retained as
unresolved. No locus ID is assigned from topology similarity or from a guessed
row correspondence. Therefore the complete published collapsed collection is
not an exact verifiable input for locus-filtered 48-taxon treatments.

## Source evidence

| source | checksum / version |
|---|---|
| `63430.aLRT-0.95-collapsed.gene.trees.gz` | `{EXPECTED['collapsed_sha256']}` |
| `63430.named.gene.trees.gz` | `{EXPECTED['named_sha256']}` |
| `master_table_gene_trees.txt` | `e6d01f4ecd9811f525da07a8702f811804998f135342d94a46a2f670247332f8` |
| PNAS repository | `17b408a622c9589fe9f6ca8a434c150974e39cd4` |
| PNAS Zenodo archive | `chr4avian-v1.0.0.zip`, `0dc24ebde4dd904115bd1e56735550e792d5cf8f7af437dbadf809108c313763` |
| Stiller Figshare | `10.6084/m9.figshare.25825408.v1` |
| Stiller ERDA | `https://sid.erda.dk/sharelink/ENhZODU9YE` |

The Stiller README explicitly identifies the collapsed aLRT-0.95 archive as
the ASTRAL input and describes the named archive as locus-name plus Newick
records. The PNAS `delta_quartet.R` script independently establishes the
inclusive GalGal4 window-start predicate used for outlier selection.

## Counts and exclusions

- Collapsed records: **{counts['n_rows']:,}**.
- Unique header-derived loci: **{counts['n_unique_locus_ids']:,}**.
- Taxon-set matches: **{counts['taxon_set_matches']:,}**.
- Unresolved collapsed records: **{counts['taxon_set_unresolved']:,}**.
- Published outlier loci under the frozen predicate: **{counts['published_outlier_loci']:,}**.
- Existing J48 inputs preserved: **{len(old)}**.

The 1,431-locus outlier mask is verified against the header-derived locus
coordinates and remains unchanged. This does not prove that the same 1,431
IDs can be selected from the unresolved collapsed records.

## Dataset decision

No `PNAS48_T0` or outlier-removed exact input is created. Doing so would
silently turn the 317 unresolved records into genomic loci. Existing `J48_*`
files remain the reconstructed named-tree re-collapse datasets. They are
appropriate for the MSRC extension after being labeled as reconstructed
design inputs, but they should not be described as exact PNAS input
reproductions.

The source archive is usable for a future **partial verified** dataset made of
the {counts['taxon_set_matches']:,} matched records, but that dataset would omit
317 loci and must be labeled as such. It is not produced by this audit.

## Published species-tree outputs

The available `score-columbea/` files are Figure S9 scoring outputs and the
available `score.tre` is a full published tree; the archive does not identify a
48-taxon full-data or outlier-removed species tree. No new species tree was
inferred here. The 363-taxon published reference remains a reference output,
not validation of a 48-taxon reproduction.

## Recommended next steps

Use the existing reconstructed-design `J48_T0`, `J48_T_PNAS`,
`J48_T_STRUCT_PRIMARY`, `J48_T_BLOCK_500K`, and `J48_NONCHR4` inputs for the
MSRC extension only, with their current checksums. For a strict PNAS replication,
first obtain a source-side file or author confirmation that links all 63,430
collapsed records to locus IDs, or explicitly accept the 317-locus omission as
a partial verified sensitivity. Then prepare Bridges-2 single-treatment jobs
without changing the 363-species inputs; do not submit them from this audit.
"""
        (RESULTS / "pnas48_source_reconstruction.md").write_text(report)
    print(json.dumps(counts, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
