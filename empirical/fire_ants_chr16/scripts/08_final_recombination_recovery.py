#!/usr/bin/env python3
"""Final targeted fire-ant recombination recovery audit.

This stage records the published Si_gnG/optical-chromosome correspondence
without treating it as a validated transfer to the frozen Si_gnGA axis.  It
therefore deliberately produces no cM/Mb estimates.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import zipfile
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "fire_ants_chr16"
EMP = ROOT / "empirical" / "fire_ants_chr16"
AGP_ZIP = DATA / "raw/recombination_sources/stolle_2019/assembly_corrections_AGP.zip"
AGP_NAME = "S.inv.AGP.txt"
MARKERS = DATA / "processed/recombination_marker_level_recovery.tsv"
WINDOWS = DATA / "processed/recombination_window_integration.tsv"
STAGE6C = EMP / "results/stage6c_manifest.json"
OUT_MAP = DATA / "processed/final_historical_chr16_scaffold_mapping.tsv"
OUT_WINDOWS = DATA / "processed/final_recombination_window_integration.tsv"
OUT_SOURCES = DATA / "metadata/final_recombination_source_manifest.tsv"
OUT_REPORT = EMP / "results/final_recombination_recovery_report.md"
OUT_MANIFEST = EMP / "results/final_recombination_recovery_manifest.json"

CLASSIFICATION = "QUANTITATIVE_RECOMBINATION_TRACK_NOT_RECOVERED"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_agp() -> list[dict[str, str]]:
    rows = []
    with zipfile.ZipFile(AGP_ZIP) as archive:
        with archive.open(AGP_NAME) as raw:
            for line in raw.read().decode().splitlines():
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                if fields[0] != "chr16" or fields[4] != "W":
                    continue
                rows.append(
                    {
                        "historical_assembly": "Si_gnG / GCF_000188075.1",
                        "historical_scaffold": fields[5],
                        "historical_linkage_group": "lg16 (optical chr16)",
                        "source_start": fields[6],
                        "source_end": fields[7],
                        "optical_chr16_start": fields[1],
                        "optical_chr16_end": fields[2],
                        "orientation": fields[8],
                        "target_assembly": "Si_gnGA / gng20170922wFex.fa",
                        "target_chromosome": "chr16",
                        "target_start": "",
                        "target_end": "",
                        "identity": "",
                        "coverage": "",
                        "mapping_quality": "historical_to_optical_only",
                        "mapping_status": "not_transferred_to_Si_gnGA",
                        "exclusion_reason": (
                            "AGP uses GCF_000188075.1 scaffold accessions; no validated "
                            "Si_gnG/Si_gnF-to-Si_gnGA correspondence or sequence alignment "
                            "was recovered."
                        ),
                        "source_url": (
                            "https://wurmlab.com/data/optical_mapping/"
                            "Stolle.et.al.degenerative.expansion.of.a.young.supergene."
                            "Solenopsis.assembly.corrections.AGP.zip"
                        ),
                    }
                )
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def make_window_file() -> tuple[int, int]:
    rows = read_tsv(WINDOWS)
    fields = [
        "window_id", "chromosome", "start", "end", "supergene_classification",
        "available_recombination_estimate", "recombination_estimate_type",
        "number_of_supporting_markers", "supporting_physical_interval",
        "distance_to_nearest_marker", "coverage_or_reliability_flag",
        "quartet_discordance_statistic", "topology_support_statistics",
        "branch_length_statistics",
    ]
    out = []
    for row in rows:
        out.append(
            {
                "window_id": row["window_id"],
                "chromosome": row["chromosome"],
                "start": row["start"],
                "end": row["end"],
                "supergene_classification": row["supergene_classification"],
                "available_recombination_estimate": "",
                "recombination_estimate_type": "missing_no_validated_target_interval",
                "number_of_supporting_markers": "0",
                "supporting_physical_interval": "",
                "distance_to_nearest_marker": "",
                "coverage_or_reliability_flag": "uncovered_quantitative_recombination",
                "quartet_discordance_statistic": row["quartet_discordance_D"],
                "topology_support_statistics": "q_species,q_haplotype,q_third from frozen Stage 4A",
                "branch_length_statistics": row["branch_length_statistics"],
            }
        )
    write_tsv(OUT_WINDOWS, out, fields)
    return len(rows), sum(row["chromosome"] == "chr16" for row in rows)


def write_sources() -> None:
    fields = ["source_id", "study", "data_type", "assembly_or_axis", "accession_or_identifier", "url", "retrieved_or_status", "usable_for_target_rates", "reason"]
    rows = [
        {"source_id": "wang2013", "study": "Wang et al. 2013", "data_type": "direct linkage map", "assembly_or_axis": "Si_gnF scaffolds", "accession_or_identifier": "nature11832; 27,934 rows in frozen audit", "url": "https://www.nature.com/articles/nature11832", "retrieved_or_status": "retrieved", "usable_for_target_rates": "no", "reason": "No marker sequences and no validated Si_gnF-to-Si_gnGA transfer."},
        {"source_id": "pracana2017", "study": "Pracana et al. 2017", "data_type": "ALLMAPS linkage-map placement method", "assembly_or_axis": "Si_gnG / GCA_000188075.1", "accession_or_identifier": "Dryad 10.5061/dryad.js509", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5485014/", "retrieved_or_status": "paper and supplement audited; ALLMAPS inputs not public in cited supplement", "usable_for_target_rates": "no", "reason": "Paper documents seven maps and linkage-group placement but does not provide a marker-transfer table to Si_gnGA."},
        {"source_id": "stolle2019_agp", "study": "Stolle et al. 2019", "data_type": "assembly-correction AGP", "assembly_or_axis": "Si_gnG/GCF scaffolds to optical chr16", "accession_or_identifier": "31 chr16 W-components; zip SHA256 56d2686ac7b072d0bf5e35f99a6c58d218230efdd9d8405ddd5ccfdc9e2d87dd", "url": "https://wurmlab.com/data/optical_mapping/", "retrieved_or_status": "retrieved", "usable_for_target_rates": "no", "reason": "Validates historical scaffold order/orientation only; target Si_gnGA equivalence remains unresolved."},
        {"source_id": "ross2018", "study": "Ross and Shoemaker 2018", "data_type": "pedigree pairwise fractions", "assembly_or_axis": "Si_gnH_C3 marker panel", "accession_or_identifier": "Dryad 10.5061/dryad.hh6342f; 15-locus panel", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6223060/", "retrieved_or_status": "article and supplementary PDFs retrieved; Dryad workbook stream HTTP 403", "usable_for_target_rates": "no", "reason": "Three supergene-linked loci support pairwise fractions but no validated physical interval or machine-readable numeric table was recovered."},
        {"source_id": "yan2020", "study": "Yan et al. 2020", "data_type": "LD and inversion evidence", "assembly_or_axis": "Si_gnH/chr16", "accession_or_identifier": "PRJNA421367", "url": "https://www.nature.com/articles/s41559-019-1081-1", "retrieved_or_status": "supplementary files audited", "usable_for_target_rates": "no", "reason": "LD and historical gene flux are not direct cM/Mb measurements."},
    ]
    write_tsv(OUT_SOURCES, rows, fields)


def write_report(marker_rows: int, wang_rows: int, mapped: int, windows: int, chr16_windows: int, agp_rows: int) -> None:
    OUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    text = f"""# Final targeted fire-ant recombination recovery

## Classification

`{CLASSIFICATION}`

The final targeted attempt recovered a validated historical assembly placement
but no validated marker-level physical coordinates on the frozen `Si_gnGA`
axis. No cM/Mb estimate, interval-level rate, window-level rate, association
test, or recombination figure is therefore produced.

## Frozen analysis audit

The prior recovery state was audited at commit `ed101dcc5040ca3ab3fed5281b49a74bda9fe5cb`.

- Reference: Stolle `Si_gnGA`, `gng20170922wFex.fa`; target `chr16`.
- Frozen supergene: `chr16:11,680,438–27,917,498`.
- Existing windows: {windows}; chromosome-16 windows: {chr16_windows}.
- Wang marker rows retained from Stage 7: {wang_rows}; total marker/pairwise rows: {marker_rows}; target chr16 mappings: {mapped}.
- Existing Stage 4A–6C topology, discordance, and branch-length outputs were not modified.

## What the targeted correspondence recovered

Pracana et al. (2017) used seven Wang linkage maps to order and orient
scaffolds from the older `Si_gnG`/`GCA_000188075.1` reference with ALLMAPS.
The cited public supplement contains methods and population-genetic tables,
not the ALLMAPS input files, scaffold-order table, or an AGP. The cited Dryad
record is `10.5061/dryad.js509`.

Stolle et al. (2019) provide the closest recoverable bridge: the archived
`S.inv.AGP.txt` orders and orients {agp_rows} `GCF_000188075.1` scaffold
components on an optical `chr16` pseudomolecule. This is preserved in
`data/fire_ants_chr16/processed/final_historical_chr16_scaffold_mapping.tsv`.
The table deliberately leaves Si_gnGA coordinates blank. Its scaffold names
are accessioned `NW_...` identifiers, whereas every Wang row uses a
`Si_gnF.scaffold..._nt...` identifier. The available records do not establish
that identifier correspondence, and no sequence-based transfer to
`gng20170922wFex.fa` was available. Optical chromosome `chr16` was therefore
not silently treated as the current Si_gnGA `chr16`.

## Wang et al. (2013)

The Wang tables contain linkage-group labels, marker identifiers, original
Si_gnF scaffold/offset strings, and family-specific cM positions. They contain
no marker sequences. The historical AGP can establish an optical linkage-group
context, but it cannot lift the Wang Si_gnF marker offsets without a validated
Si_gnF-to-GCF/Si_gnGA correspondence. Consequently, zero markers were newly
mapped to target chr16, including zero inside and zero outside the supergene.
The earlier zero-marker result was therefore an assembly/identifier
correspondence failure, rather than evidence that Wang measured no chr16
markers.

## Ross and Shoemaker (2018)

Ross and Shoemaker did not publish a 2,796-marker RAD-seq map. The 2,796
RAD-marker figure belongs to Wang's P034 linkage map. Ross and Shoemaker used
a 15-locus pedigree panel and reported pairwise segregation/recombination
measurements for the three supergene-linked loci `C294`, `Gp-9`, and `i_126`
among polygyne queen progenies. These measurements remain pairwise-only:
their physical positions were not validated on Si_gnGA and the public raw
workbook could not be retrieved from Dryad in this environment. They cannot
be converted to cM/Mb or assigned to windows without unsupported assumptions.

## Quantitative estimates and integration

No adjacent marker pair has both validated Si_gnGA physical coordinates and
genetic positions. Therefore the rate formula was not applied. No marker
spacing, interval length, confidence bound, genotype-specific map, or
SB/SB-versus-SB/Sb comparison is reported. All {windows} window-level
recombination values remain missing; {chr16_windows} chr16 windows have
zero supporting markers. Quartet discordance, topology support, and branch
lengths are retained as existing measurements but are not tested against
recombination because the required spatially resolved quantitative covariate
is absent. No continuous track or association figure is scientifically
supported.

## Remaining blockers

1. Original Wang marker sequences or a validated Si_gnF-to-GCF/Si_gnGA
   scaffold correspondence is unavailable.
2. Pracana's described ALLMAPS inputs and intermediate scaffold-order files
   are not present in the cited supplement/data record recovered here.
3. Ross's pedigree panel supplies useful pairwise evidence but not a
   validated physical map on the frozen axis; its raw workbook was not
   accessible.
4. The historical optical chr16 AGP is an assembly correspondence artifact,
   not a direct meiotic recombination map.

The final-recovery source manifest, historical mapping table, explicit
missing-value window table, and regression checks are self-contained in this
repository. Execution: `python3 empirical/fire_ants_chr16/scripts/08_final_recombination_recovery.py`.
"""
    OUT_REPORT.write_text(text)


def run() -> None:
    if not AGP_ZIP.exists():
        raise FileNotFoundError(AGP_ZIP)
    agp_rows = parse_agp()
    fields = list(agp_rows[0])
    write_tsv(OUT_MAP, agp_rows, fields)
    marker_rows = read_tsv(MARKERS)
    wang_rows = sum(row.get("source_study") == "Wang et al. 2013" for row in marker_rows)
    mapped = sum(row.get("mapping_status") == "mapped_to_target_chr16" for row in marker_rows)
    windows, chr16_windows = make_window_file()
    write_sources()
    write_report(len(marker_rows), wang_rows, mapped, windows, chr16_windows, len(agp_rows))
    manifest = {
        "classification": CLASSIFICATION,
        "historical_agp_chr16_components": len(agp_rows),
        "wang_marker_rows": len(marker_rows),
        "wang_markers_mapped_to_target_chr16": mapped,
        "windows": windows,
        "chr16_windows": chr16_windows,
        "target_rate_intervals": 0,
        "stage6c_manifest_sha256": sha256(STAGE6C),
        "agp_zip_sha256": sha256(AGP_ZIP),
        "outputs": [str(path.relative_to(ROOT)) for path in (OUT_MAP, OUT_WINDOWS, OUT_SOURCES, OUT_REPORT)],
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


class FinalRecoveryTests(unittest.TestCase):
    def test_historical_agp_is_chr16_and_not_target_transfer(self):
        rows = parse_agp()
        self.assertEqual(len(rows), 31)
        self.assertTrue(all(row["historical_scaffold"].startswith("NW_") for row in rows))
        self.assertTrue(all(row["target_start"] == "" for row in rows))
        self.assertTrue(all(row["mapping_status"] == "not_transferred_to_Si_gnGA" for row in rows))

    def test_frozen_marker_and_window_counts(self):
        markers = read_tsv(MARKERS)
        self.assertEqual(len(markers), 27937)
        self.assertEqual(sum(row.get("mapping_status") == "mapped_to_target_chr16" for row in markers), 0)
        windows = read_tsv(WINDOWS)
        self.assertEqual(len(windows), 213)
        self.assertEqual(sum(row["chromosome"] == "chr16" for row in windows), 96)

    def test_classification_and_frozen_stage_manifest(self):
        self.assertEqual(CLASSIFICATION, "QUANTITATIVE_RECOMBINATION_TRACK_NOT_RECOVERED")
        manifest = json.loads(STAGE6C.read_text())
        for path, expected in manifest["input_checksums"].items():
            candidate = ROOT / path
            if candidate.exists():
                self.assertEqual(sha256(candidate), expected, path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-tests", action="store_true")
    args = parser.parse_args()
    if args.run_tests:
        unittest.main(argv=[__file__], exit=False)
    else:
        run()


if __name__ == "__main__":
    main()
