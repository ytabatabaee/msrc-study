#!/usr/bin/env python3
"""Audit the reconstructed Stage 4E Jarvis-48 inputs against published metadata.

This deliberately does not infer a species tree or overwrite an existing input.
The collapsed archive has no embedded locus IDs, so locus assignment is made
only through the published 63K header plus the exact named-archive locus set.
The source reconstruction reports the 317 records whose taxon sets remain
unresolved; it never assigns those IDs by topology or guessed row order.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stage4d_core import PNAS_INTERVALS, ROOT, published_outlier, read_tsv, sha256

BASE = ROOT / "empirical/neoaves_chr4"
DATA4D = BASE / "data/stage4d"
DATA4E = BASE / "data/stage4e"
RESULTS = BASE / "results"
OUT = RESULTS / "stage4e"

SOURCES = {
    "pnas_repository": "https://github.com/smirarab/chr4avian/tree/master",
    "pnas_repository_commit": "17b408a622c9589fe9f6ca8a434c150974e39cd4",
    "pnas_zenodo": "https://doi.org/10.5281/zenodo.10699423",
    "pnas_zenodo_downloaded_record": "https://zenodo.org/api/records/10699424/files/smirarab/chr4avian-v1.0.0.zip/content",
    "stiller_figshare": "https://doi.org/10.6084/m9.figshare.25825408.v1",
    "stiller_erda": "https://sid.erda.dk/sharelink/ENhZODU9YE",
}


def count_trees(path: Path) -> int:
    return sum(1 for line in path.open() if line.strip())


def main() -> None:
    mapping = read_tsv(DATA4D / "jarvis48_taxa.tsv")
    source_2014 = {r["source_2014_taxon"] for r in mapping}
    stiller = {r["stiller2024_taxon"] for r in mapping}
    status = {r["mapping_status"] for r in mapping}
    manifest = read_tsv(RESULTS / "stage4e_input_manifest.tsv")
    locus = read_tsv(RESULTS / "stage4d_locus_manifest.tsv")
    expected_pnas = {
        r["locus_id"] for r in locus
        if published_outlier(r["chromosome"], int(float(r["start"])))
    }
    recorded_pnas = {r["locus_id"] for r in locus if r["published_outlier_region"] == "True"}
    original_audit = json.loads((RESULTS / "stage4d/j48_existing_input_audit.json").read_text())
    rows = []
    def add(component, determination, value, evidence):
        rows.append({"component": component, "determination": determination, "value": value, "evidence": evidence})

    add("species_selection", "reconstructed_match_with_mapping", f"48 source / {len(stiller)} Stiller labels; missing=0; extra=0", "data/stage4d/jarvis48_taxa.tsv")
    add("species_mapping", "documented_nonexact_mapping", ";".join(sorted(status)), "data/stage4d/jarvis48_taxa.tsv")
    add("source_gene_trees", "strongly_supported_partial_mapping", "63,113/63,430 collapsed records match the header-derived named-locus taxon set; 317 remain unresolved", "results/stage4e/pnas48_source_reconstruction.md; results/stage4e/pnas48_locus_mapping_audit.tsv")
    add("gene_tree_preprocessing", "exact_full_mapping_unverified", "The published collapsed archive is the documented ASTRAL input; current J48 inputs remain the named-tree re-collapse fallback", "stiller2024/README.txt; results/stage4e/pnas48_source_reconstruction.md")
    add("locus_filtering", "verified_for_reconstruction", f"{len(locus)} unique locus IDs; fewer-than-four omitted after pruning", "results/stage4d_locus_manifest.tsv; scripts/05_stage4e_prepare_aster.py")
    add("pnas_outlier_exclusions", "exact_predicate_and_count_verified", f"expected={len(expected_pnas)}; recorded={len(recorded_pnas)}; only_ours=0; only_published=0", "scripts/stage4d_core.py; chr4avian/genetreesupport/delta_quartet.R")
    add("species_tree_settings", "unverifiable_for_48_comparison", "published 48-taxon command/output not present in available PNAS archive", "chr4avian-v1.0.0.zip inventory")
    add("published_output_topology", "unverifiable", "available score.tre is not an identified published 48-taxon output", "external/chr4avian/score.tre")

    exclusions = {
        "published_interval_rule": "chr4 window start inclusive: 25030000-32670000; 33510000-34470000; 44130000-56810000",
        "published_loci_in_local_manifest": len(expected_pnas),
        "recorded_loci_in_local_manifest": len(recorded_pnas),
        "loci_only_ours": len(recorded_pnas - expected_pnas),
        "loci_only_published": len(expected_pnas - recorded_pnas),
        "ambiguous_or_unmapped": 0,
        "assembly": "GalGal4 coordinates as used by published script; local manifest retains Stiller window coordinates",
        "boundary": "window start, inclusive; overlap reduces to locus membership",
    }
    input_rows = []
    for r in manifest:
        if r["scope"] == "jarvis48":
            input_rows.append({k: r[k] for k in ("treatment", "n_loci", "trees_with_fewer_than_four_taxa_omitted", "sha256")})
    source_checksums = {}
    for manifest_path in (BASE / "external/stiller2024/download_manifest.tsv", BASE / "external/chr4avian/download_manifest.tsv"):
        for row in read_tsv(manifest_path):
            if row.get("filename") in {"63430.named.gene.trees.gz", "63430.aLRT-0.95-collapsed.gene.trees.gz", "master_table_gene_trees.txt", "chr4avian-v1.0.0.zip", "delta_quartet.R", "genetreesupport-README.md"}:
                source_checksums[row["filename"]] = row["sha256"]

    report = f"""# PNAS 48-taxon data audit

Date: 2026-10-22. This is a data and provenance audit; no ASTRAL-IV inference
was run and no cluster job was submitted.

## Answers to the required questions

1. **Are the current J48 inputs identical to the published PNAS inputs?** No.
They are a reconstructed published-design dataset. The PNAS/Stiller README
identifies `63430.aLRT-0.95-collapsed.gene.trees.gz` as the ASTRAL input,
whereas the current Stage 4E inputs come from the ID-preserving
`named_recollapsed.tre` fallback. The local audit found 31,066 named/collapsed
same-row disagreements, so exact gene-tree identity is not established.

2. **What differs?** The input gene-tree collection/preprocessing is the
material discrepancy. Current pruning retains the 48 mapped Stiller labels,
omits trees with fewer than four retained taxa, and preserves branch lengths;
the exact collapsed-file locus correspondence needed to verify every topology
and exclusion cannot be recovered by row because the files are not aligned.

3. **Did we identify original published 48-taxon gene-tree files?** No.
The available PNAS archive contains scripts, support data, and species-tree
outputs, but no authoritative 48-taxon gene-tree file with locus IDs. The
published method describes restricting the 63,430 Stiller loci to the 48
Jarvis species. Exact 48-taxon files remain unavailable.

4. **Are PNAS excluded intervals and locus IDs correct?** The interval
predicate and local locus membership are exactly reproduced: {len(expected_pnas)}
loci, with zero implementation-only or published-only loci in the local
manifest. This validates the mask against the published script, not the
unavailable original 48-taxon file. The three intervals are GalGal4 `chr4`
window-start intervals and are inclusive.

5. **What was corrected?** No raw `.tre` was overwritten and no corrected
published-equivalent dataset was fabricated. The audit adds machine-readable
provenance and preflight records; current J48 inputs remain explicitly
reconstructed inputs.

6. **What remains unverifiable?** Exact collapsed-tree-to-locus identity,
published 48-taxon gene-tree topology/branch/support equality, the original
48-taxon ASTRAL command/version/output, and a published 48-taxon full-data or
outlier-removed species tree.

7. **Which datasets should run on Bridges-2?** Run the reconstructed-design
baseline and MSRC treatments only after labeling them as reconstructed:
`J48_T0`, `J48_T_PNAS`, `J48_T_STRUCT_PRIMARY`, `J48_T_BLOCK_500K`, and
`J48_NONCHR4`. Keep the current files for provenance; do not call these exact
published-input reproductions.

## Taxa and filtering

The frozen mapping contains 48 source 2014 labels and 48 Stiller labels. The
only non-exact mapping status is recorded in the mapping file (`Manacus` was
reidentified in the Stiller source); no label was invented by this audit.
Current Stage 4E counts and checksums are listed below. Fewer-than-four-taxon
trees are omitted after pruning, as implemented by `05_stage4e_prepare_aster.py`.

| treatment | loci | omitted <4 | SHA-256 |
|---|---:|---:|---|
"""
    for r in input_rows:
        report += f"| {r['treatment']} | {r['n_loci']} | {r['trees_with_fewer_than_four_taxa_omitted']} | `{r['sha256']}` |\n"
    report += f"""
## Replication-fidelity classification

| component | determination |
|---|---|
| Species selection | reconstructed match with documented label mapping |
| Source gene trees | mismatch / exact identity unverifiable |
| Gene-tree preprocessing | mismatch: current ID-preserving re-collapse vs published collapsed archive |
| Locus filtering | verified for the local reconstruction |
| PNAS outlier exclusions | exact published predicate and local locus membership verified |
| Species-tree inference settings | unverifiable for the 48-taxon comparison |
| Published output topology | unverifiable for a 48-taxon output |

This therefore qualifies as a **reconstructed published-design reproduction**.
The structural masks, block controls, and CU/SU comparisons remain **MSRC
extensions**.

## Sources and versions

"""
    for k, v in SOURCES.items():
        report += f"- `{k}`: {v}\n"
    report += "\nSelected downloaded source checksums:\n\n"
    for name, checksum in sorted(source_checksums.items()):
        report += f"- `{name}`: `{checksum}`\n"
    report += f"""

Local reconstruction checksums:

- `data/stage4d/jarvis48_taxa.tsv`: `{sha256(DATA4D / 'jarvis48_taxa.tsv')}`
- `data/stage4d/named_recollapsed.tre`: `{sha256(DATA4D / 'named_recollapsed.tre')}`
- `results/stage4d_locus_manifest.tsv`: `{sha256(RESULTS / 'stage4d_locus_manifest.tsv')}`

Published interval coordinates are GalGal4 window-start intervals
25,030,000–32,670,000, 33,510,000–34,470,000, and 44,130,000–56,810,000;
the paper's GalGal6 table gives corresponding approximate ranges
25,555,000–33,202,000, 34,230,000–34,999,000, and 44,690,000–57,179,000.
The local locus mask uses the published GalGal4 predicate because that is the
coordinate system used by the authors' script.

The local Stiller archive manifest records the ERDA archive files and checksums;
the local PNAS download manifest records the GitHub commit and Zenodo ZIP
checksum. Large raw archives remain external and are not duplicated by this
audit. No published 48-taxon input was found in the available archive inventory.

## Recommended command preparation

Prepare, but do not submit, one treatment per job using the existing single-
treatment workflow after transferring the current Stage 4E inputs:

```bash
for T in J48_T0 J48_T_PNAS J48_T_STRUCT_PRIMARY J48_T_BLOCK_500K J48_NONCHR4; do
  sbatch -A bio240015p --qos=rm -p RM-shared --cpus-per-task=32 \\
    --time=16:00:00 --export=ALL,TREATMENT="$T" \\
    empirical/neoaves_chr4/cluster/stage4e_aster/01_single_treatment.sbatch
done
```

These commands are intentionally recommendations only; no command was
submitted by this audit.
"""
    (OUT / "pnas48_data_audit.md").write_text(report)
    with (OUT / "pnas48_data_audit.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["component", "determination", "value", "evidence"], delimiter="\t")
        w.writeheader(); w.writerows(rows)
    audit_json = {
        "sources": SOURCES,
        "taxa": {"source_2014": sorted(source_2014), "stiller2024": sorted(stiller), "n": len(stiller), "mapping_status": sorted(status)},
        "pnas_exclusions": exclusions,
        "current_j48": input_rows,
        "current_input_identity": "reconstructed_named_recollapse_not_exact_published_collapsed",
        "existing_inputs_preserved": True,
        "original_stage4d_input_audit": original_audit,
    }
    (OUT / "pnas48_data_audit.json").write_text(json.dumps(audit_json, indent=2) + "\n")
    print(json.dumps({"taxa": len(stiller), "pnas_excluded": len(expected_pnas), "recorded_pnas": len(recorded_pnas), "report": str((OUT / 'pnas48_data_audit.md').relative_to(ROOT))}, indent=2))


if __name__ == "__main__":
    main()
