#!/usr/bin/env python3
"""Stage 1A audit for the published shared-PIP analysis of CM027530.1.

This script reconstructs source objects, representative samples, shared-region
coordinates, and workflow provenance from public small metadata/code artifacts.
It does not download reads, process BAMs, infer trees, or inspect topology.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "songbirds_pips"
EMP = ROOT / "empirical" / "songbirds_pips"
RAW = DATA / "raw" / "oup_supplement"
PUB = DATA / "raw" / "public_metadata"
FIG = DATA / "raw" / "figshare_27284553"
XLSX = RAW / "DatasetS1_S4_S5.xlsx"
TODAY = "2026-10-02"

NS = {
    "a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}

FOCAL_SHARED = "CM027530.1"
FOCAL_ID = "shared_CM027530.1"
REFERENCE = "GCA_001746935.2"
REFERENCE_SPECIES = "Setophaga coronata"


def ensure_dirs() -> None:
    for p in [DATA / "metadata", DATA / "processed", EMP / "results"]:
        p.mkdir(parents=True, exist_ok=True)


def sha256(path: Path) -> str:
    if not path.exists() or path.is_dir():
        return ""
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fields})


def xlsx_sheets(path: Path) -> dict[str, list[dict[str, str]]]:
    with zipfile.ZipFile(path) as z:
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rid_to_target = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}
        strings = []
        if "xl/sharedStrings.xml" in z.namelist():
            ss = ET.fromstring(z.read("xl/sharedStrings.xml"))
            for si in ss.findall("a:si", NS):
                strings.append("".join(t.text or "" for t in si.iter(f"{{{NS['a']}}}t")))

        def cell_value(c: ET.Element) -> str:
            t = c.attrib.get("t")
            v = c.find("a:v", NS)
            if t == "inlineStr":
                is_ = c.find("a:is", NS)
                return "" if is_ is None else "".join(x.text or "" for x in is_.iter(f"{{{NS['a']}}}t"))
            if v is None:
                return ""
            val = v.text or ""
            return strings[int(val)] if t == "s" else val

        out = {}
        for s in wb.find("a:sheets", NS):
            name = s.attrib["name"]
            rid = s.attrib[f"{{{NS['r']}}}id"]
            target = rid_to_target[rid]
            if not target.startswith("worksheets/"):
                continue
            ws = ET.fromstring(z.read("xl/" + target))
            rows = []
            for row in ws.findall(".//a:sheetData/a:row", NS):
                vals = []
                last_col = 0
                for c in row.findall("a:c", NS):
                    ref = c.attrib.get("r", "")
                    letters = re.sub(r"[^A-Z]", "", ref)
                    col = 0
                    for ch in letters:
                        col = col * 26 + (ord(ch) - 64)
                    while last_col + 1 < col:
                        vals.append("")
                        last_col += 1
                    vals.append(cell_value(c))
                    last_col = col
                rows.append(vals)
            if not rows:
                out[name] = []
                continue
            header = [str(x).strip() for x in rows[0]]
            data = []
            for vals in rows[1:]:
                if not any(vals):
                    continue
                data.append({header[i]: vals[i] if i < len(vals) else "" for i in range(len(header))})
            out[name] = data
        return out


def intish(x: str) -> int:
    return int(round(float(x)))


def main() -> None:
    ensure_dirs()
    sheets = xlsx_sheets(XLSX)
    d1 = sheets["DatasetS1"]
    d4 = sheets["DatasetS4"]
    d5 = sheets["DatasetS5"]

    d4_index = {(r["Specimen.number"], r["Species"], r["PIP.name"], r["PIP.genotype"]): r for r in d4}
    focal_d1 = [r for r in d1 if r["Chromosome.accession"] == FOCAL_SHARED and r["PIP.name"] == FOCAL_SHARED]
    focal_d5 = [r for r in d5 if r["Shared.PIP.tree"] == FOCAL_SHARED]

    starts = [intish(r["PIP.start"]) for r in focal_d1]
    ends = [intish(r["PIP.end"]) for r in focal_d1]
    shared_start = max(starts)
    shared_end = min(ends)
    shared_len = shared_end - shared_start
    left_comp_end = min(intish(r["comp_break_right"]) for r in focal_d1)
    right_comp_start = max(intish(r["comp_break_left"]) for r in focal_d1)

    source_rows = []

    def add_source(source: str, description: str, url: str, version: str, path: Path | None, role: str, notes: str) -> None:
        source_rows.append({
            "source": source,
            "description": description,
            "url": url,
            "version_or_commit": version,
            "local_path": "" if path is None else str(path.relative_to(ROOT)),
            "sha256": "" if path is None else sha256(path),
            "analysis_role": role,
            "notes": notes,
        })

    add_source("primary_article_methods", "Pegan and Winger 2025 article Methods and Results text for shared-PIP analysis.", "https://doi.org/10.1093/gbe/evaf205", "published 2025-11-03; corrected/typeset 2025-11-19", None, "describes shared-PIP overlap, representative selection, genetic-distance calculation, and neighbor-joining plotting", f"Accessed {TODAY}.")
    add_source("supplement_workbook", "OUP supplementary workbook containing Dataset S1, Dataset S4, and Dataset S5.", "https://academic.oup.com/gbe/article/17/11/evaf205/8313116#supplementary-data", "DatasetS1_S4_S5.xlsx", XLSX, "source for species-specific PIP coordinates, genotype labels, sample metadata, and 22 shared-PIP representatives", "Dataset S5 is the fixed published representative set used here.")
    add_source("supplement_pdf", "OUP supplementary information PDF.", "https://academic.oup.com/gbe/article/17/11/evaf205/8313116#supplementary-data", "SI_GBE_6Oct25.pdf", RAW / "SI_GBE_6Oct25.pdf", "supplementary methods/figures; no reusable CM027530.1 intermediate matrix found locally", "Full text extraction was not required for coordinate/sample reconstruction.")
    for name, role in [
        ("README.txt", "states that large sequence datasets/intermediates are not included and raw reads require preprocessing"),
        ("1_Bioinformatic_code/B3a_SAFs_for_genetrees.txt", "documents representative selection and ANGSD SAF generation for shared-PIP and COMP regions"),
        ("1_Bioinformatic_code/B3b_genetic_distance_for_genetrees.txt", "documents pairwise winsfs 2D-SFS and dxy calculation"),
        ("1_Bioinformatic_code/B3c_plot_genetrees.txt", "documents neighbor-joining plotting from concatenated dxy matrices"),
        ("1_Bioinformatic_code/draw_pip_inds.R", "documents random representative selection by species and non-heterozygous genotype"),
        ("1_Bioinformatic_code/dxy_sfs.py", "converts 2D-SFS values to genetic distance/dxy"),
    ]:
        add_source("oup_author_code_" + Path(name).stem, "Packaged OUP author code: " + name, "https://academic.oup.com/gbe/article/17/11/evaf205/8313116#supplementary-data", "packaged supplement code, not a moving branch", RAW / "code" / "code" / name, role, "Small text/code file frozen locally.")
    add_source("figshare_metadata", "Figshare API metadata for preprocessing-code repository cited by the OUP README.", "https://api.figshare.com/v2/articles/27284553", "figshare article 27284553, version 1", PUB / "figshare_27284553.json", "inventory of public preprocessing scripts and absence of shared-PIP intermediates", "Metadata lists code and small source data only; no CM027530.1 VCF/BAM/SAF/dxy/Newick files.")
    for name in ["1.Sample_Processing", "2.Realign_around_indels", "3.Print_SNP_GLs_with_ANGSD", "4.Filter_SNPs_with_ngsParalog", "README.md"]:
        add_source("figshare_" + name, "Figshare preprocessing source object cited as prerequisite by OUP supplement README.", "https://doi.org/10.6084/m9.figshare.27284553", "version 1", FIG / name, "raw-read preprocessing prerequisite; not sufficient to reproduce shared-PIP analysis without reads/BAMs", "Downloaded small code file only.")
    add_source("ncbi_prjna1130443_metadata", "NCBI BioProject metadata for public raw reads.", "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1130443", "NCBI ESummary downloaded 2026-10-02", PUB / "ncbi_prjna1130443.xml", "raw-read archive metadata", "BioProject data type is raw sequence reads.")
    add_source("ncbi_prjna1043688_metadata", "NCBI BioProject metadata for public raw reads.", "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1043688", "NCBI ESummary downloaded 2026-10-02", PUB / "ncbi_prjna1043688.xml", "raw-read archive metadata", "BioProject data type is raw sequence reads.")
    add_source("dataset_s3_plots", "Article-referenced Dataset S3 shared-PIP/COMP tree plots.", "https://academic.oup.com/gbe/article/17/11/evaf205/8313116#supplementary-data", "referenced by article; reusable tree/matrix object not found in downloaded OUP ZIP", None, "published visual output only; not used for reconstruction or candidate selection", "No Newick, distance matrix, or machine-readable topology output was found locally.")
    write_tsv(DATA / "metadata" / "stage1a_shared_pip_manifest.tsv", source_rows, ["source", "description", "url", "version_or_commit", "local_path", "sha256", "analysis_role", "notes"])

    sample_rows = []
    missing = []
    for r in sorted(focal_d5, key=lambda x: (x["Species"], x["PIP.geno"], x["Specimen.number"])):
        key = (r["Specimen.number"], r["Species"], r["PIP.name"], r["PIP.geno"])
        meta = d4_index.get(key)
        if meta is None:
            missing.append(key)
            meta = {}
        sample_rows.append({
            "sample_id": r["Specimen.number"],
            "species": r["Species"],
            "pip_genotype": r["PIP.geno"],
            "population": "not_provided; latitude_longitude_available" if meta else "not_provided",
            "sra_accession": r["SRA.Accession"],
            "reference_assembly": f"{REFERENCE} ({REFERENCE_SPECIES})",
            "selection_source": "DatasetS5 Shared.PIP.tree=CM027530.1",
            "notes": f"published representative; PIP.name={r['PIP.name']}; voucher_museum={meta.get('Voucher.museum','')}; lat={meta.get('Latitude','')}; lon={meta.get('Longitude','')}",
        })
    write_tsv(DATA / "processed" / "stage1a_shared_pip_samples.tsv", sample_rows, ["sample_id", "species", "pip_genotype", "population", "sra_accession", "reference_assembly", "selection_source", "notes"])

    interval_rows = []
    for r in sorted(focal_d1, key=lambda x: x["Species"]):
        interval_rows.append({
            "interval_id": f"{r['Species'].replace(' ', '_')}|{FOCAL_SHARED}",
            "interval_type": "species_specific_PIP_from_DatasetS1",
            "chromosome": FOCAL_SHARED,
            "start": intish(r["PIP.start"]),
            "end": intish(r["PIP.end"]),
            "length_bp": intish(r["PIP.end"]) - intish(r["PIP.start"]),
            "reference_assembly": f"{r['Reference.genome.assembly']} ({r['Reference.genome.species']})",
            "species": r["Species"],
            "derivation_rule": "author species-specific PIP coordinates from Dataset S1",
            "provenance": "DatasetS1",
            "notes": f"LD.supported={r['LD.supported']}; comp_break_right={r['comp_break_right']}; comp_break_left={r['comp_break_left']}",
        })
    interval_rows.append({
        "interval_id": FOCAL_ID,
        "interval_type": "derived_common_overlap",
        "chromosome": FOCAL_SHARED,
        "start": shared_start,
        "end": shared_end,
        "length_bp": shared_len,
        "reference_assembly": f"{REFERENCE} ({REFERENCE_SPECIES})",
        "species": ";".join(sorted({r["Species"] for r in focal_d1})),
        "derivation_rule": "max(PIP.start) and min(PIP.end) across the 11 DatasetS1 species-specific CM027530.1 PIPs; length uses source end-start convention",
        "provenance": "DatasetS1 plus article Methods statement that shared PIP is largest region overlapping all identified PIPs",
        "notes": f"all_11_species_mapped_to_same_reference=yes; common_left_COMP_candidate=0-{left_comp_end}; common_right_COMP_candidate={right_comp_start}+; exact author shared COMP start/end not present in DatasetS1/S5",
    })
    write_tsv(DATA / "processed" / "stage1a_shared_pip_interval.tsv", interval_rows, ["interval_id", "interval_type", "chromosome", "start", "end", "length_bp", "reference_assembly", "species", "derivation_rule", "provenance", "notes"])

    workflow_rows = [
        {"step": "1", "software": "ANGSD / bwa / samtools / GATK / ngsParalog upstream preprocessing", "version_if_known": "versions in Figshare/OUP methods; exact tool versions vary by script", "input": "NCBI SRA raw reads and reference genomes", "output": "processed BAMs and genotype-likelihood/BEAGLE-style intermediates", "parameters": "see Figshare scripts 1-4 and OUP README", "source_file_or_method": "OUP README; Figshare 27284553"},
        {"step": "2", "software": "manual metadata construction plus R draw_pip_inds.R", "version_if_known": "R 4.2.0; IRanges 2.36.0; dplyr 1.1.4", "input": "shared_pips_meta, genotype table from B1, all-samples table, BAM list", "output": "selected individuals CSV for each shared PIP", "parameters": "select one individual per non-AB genotype per species; if no PIP in species, select one random individual; selection is random in the script", "source_file_or_method": "B3a_SAFs_for_genetrees.txt; draw_pip_inds.R; DatasetS5 records fixed published choices"},
        {"step": "3", "software": "ANGSD", "version_if_known": "0.941-6-g67b6b3b; htslib 1.15-1-gf307dfb", "input": "selected individual BAMs, reference genome, shared PIP BED, optional COMP BED", "output": "per-sample .saf.idx for PIP and COMP", "parameters": "-doSaf 1 -GL 2 -minMapQ 30 -minQ 30 -anc $REF -r $REG -sites $PSITES/$CSITES", "source_file_or_method": "B3a_SAFs_for_genetrees.txt"},
        {"step": "4", "software": "winsfs", "version_if_known": "0.7.0", "input": "pairwise .saf.idx files", "output": "pairwise 2D-SFS files and folded/trimmed SFS files", "parameters": "winsfs SAF1 SAF2; winsfs view -f", "source_file_or_method": "B3b_genetic_distance_for_genetrees.txt"},
        {"step": "5", "software": "dxy_sfs.py", "version_if_known": "Python 3.10.4 in author script", "input": "trimmed folded 2D-SFS", "output": "pairwise dxy/genetic distance values", "parameters": "-m 1 -n 1 -s <sfs> -o <out>", "source_file_or_method": "B3b_genetic_distance_for_genetrees.txt; dxy_sfs.py"},
        {"step": "6", "software": "R ape::nj plus plotting libraries", "version_if_known": "R 4.2.0; ape 5.7-1; reshape2 1.4.4; tidyverse 2.0.0; ggtree 3.10.0; phylobase 0.8.10", "input": "concatenated pairwise dxy CSV with samp1,samp2,sp1,sp2,geno1,geno2,gendist", "output": "neighbor-joining tree plot object for PIP and COMP", "parameters": "acast distance matrix; missing values set to 1; negative branch lengths reset to 0 for plotting", "source_file_or_method": "B3c_plot_genetrees.txt"},
    ]
    write_tsv(EMP / "results" / "stage1a_published_workflow.tsv", workflow_rows, ["step", "software", "version_if_known", "input", "output", "parameters", "source_file_or_method"])

    with (EMP / "results" / "stage1a_published_workflow.md").open("w") as f:
        f.write("# Stage 1A published shared-PIP workflow reconstruction\n\n")
        f.write("Focal published object: `shared_CM027530.1`, an author-reported shared PIP on `CM027530.1` among 11 warbler species mapped to `GCA_001746935.2` (*Setophaga coronata*). Stage 1A reconstructed the workflow from public small metadata/code only and did not run ANGSD, winsfs, NJ, or any tree inference.\n\n")
        f.write("| step | software | version | input | output | parameters | source |\n")
        f.write("|---|---|---|---|---|---|---|\n")
        for r in workflow_rows:
            f.write("| {step} | {software} | {version_if_known} | {input} | {output} | {parameters} | {source_file_or_method} |\n".format(**{k: str(v).replace("|", "/") for k, v in r.items()}))
        f.write("\nThe published object is best described as a neighbor-joining genetic-distance dendrogram/tree built from pairwise `dxy` estimates derived from ANGSD SAF and winsfs 2D-SFS calculations. The public supplement does not provide the SAF files, pairwise SFS files, dxy matrices, or Newick trees needed to reproduce the figure directly without reprocessing BAMs or reads.\n")

    intermediate_rows = [
        {"data_object": "VCF/BCF for CM027530.1 22 samples", "available": "no", "public_source_checked": "OUP supplement; Figshare 27284553; NCBI BioProject metadata; article Methods", "local_path": "", "usable_without_raw_reprocessing": "no", "notes": "No VCF/BCF found."},
        {"data_object": "BEAGLE genotype likelihoods / ANGSD .beagle.gz", "available": "no", "public_source_checked": "OUP supplement scripts and Figshare preprocessing metadata", "local_path": "", "usable_without_raw_reprocessing": "no", "notes": "Scripts refer to BEAGLE inputs, but files are not packaged."},
        {"data_object": "ANGSD SAF / .saf.idx", "available": "no", "public_source_checked": "OUP B3a/B3b scripts; supplement ZIP", "local_path": "", "usable_without_raw_reprocessing": "no", "notes": "B3a documents how to create per-sample SAF from BAMs."},
        {"data_object": "winsfs pairwise SFS files", "available": "no", "public_source_checked": "OUP B3b script; supplement ZIP", "local_path": "", "usable_without_raw_reprocessing": "no", "notes": "B3b documents generation; files absent."},
        {"data_object": "pairwise dxy/genetic-distance matrices", "available": "no", "public_source_checked": "OUP B3c script; supplement ZIP", "local_path": "", "usable_without_raw_reprocessing": "no", "notes": "B3c expects precomputed .dxy CSV files; they are not public in the downloaded supplement."},
        {"data_object": "Newick/NJ tree files", "available": "no", "public_source_checked": "OUP supplement; article Dataset S3 reference", "local_path": "", "usable_without_raw_reprocessing": "no", "notes": "Only plotted shared-PIP/COMP results are referenced; no reusable Newick files found."},
        {"data_object": "BAM/CRAM", "available": "not_publicly_packaged", "public_source_checked": "OUP README and NCBI BioProjects", "local_path": "", "usable_without_raw_reprocessing": "no_unless_authors_provide_BAMs", "notes": "Author B3a requires BAM list. Public NCBI objects are raw sequence reads."},
        {"data_object": "raw sequence reads", "available": "yes", "public_source_checked": "NCBI BioProjects PRJNA1043688 and PRJNA1130443; Dataset S4/S5 sample accessions", "local_path": "", "usable_without_raw_reprocessing": "no", "notes": "Do not download/process during Stage 1A."},
        {"data_object": "reference genome", "available": "yes", "public_source_checked": "Dataset S1", "local_path": "", "usable_without_raw_reprocessing": "partly", "notes": "Reference assembly is GCA_001746935.2; downloading/indexing would be needed for ANGSD."},
        {"data_object": "representative sample/genotype metadata", "available": "yes", "public_source_checked": "Dataset S5", "local_path": str((DATA / "processed" / "stage1a_shared_pip_samples.tsv").relative_to(ROOT)), "usable_without_raw_reprocessing": "metadata_only", "notes": "Exactly 22 representatives reconstructed."},
        {"data_object": "shared-PIP interval", "available": "yes_derived", "public_source_checked": "Dataset S1; article Methods", "local_path": str((DATA / "processed" / "stage1a_shared_pip_interval.tsv").relative_to(ROOT)), "usable_without_raw_reprocessing": "metadata_only", "notes": "Derived common overlap; exact shared COMP interval absent."},
    ]
    write_tsv(EMP / "results" / "stage1a_reusable_intermediates.tsv", intermediate_rows, ["data_object", "available", "public_source_checked", "local_path", "usable_without_raw_reprocessing", "notes"])

    species_counts = Counter(r["species"] for r in sample_rows)
    geno_counts = Counter(r["pip_genotype"] for r in sample_rows)
    validation = {
        "sample_count": len(sample_rows),
        "species_count": len(species_counts),
        "genotype_counts": dict(geno_counts),
        "one_AA_and_one_BB_per_species": all(species_counts[sp] == 2 for sp in species_counts) and all(
            Counter(r["pip_genotype"] for r in sample_rows if r["species"] == sp) == Counter({"AA": 1, "BB": 1}) for sp in species_counts
        ),
        "missing_datasetS4_joins": [list(x) for x in missing],
        "shared_interval_start": shared_start,
        "shared_interval_end": shared_end,
        "shared_interval_length_bp_end_minus_start": shared_len,
        "can_reproduce_without_raw_read_or_BAM_processing": False,
    }
    with (EMP / "results" / "stage1a_validation.json").open("w") as f:
        json.dump(validation, f, indent=2, sort_keys=True)

    with (EMP / "results" / "stage1a_report.md").open("w") as f:
        f.write("# Songbirds PIPs Stage 1A: shared_CM027530.1 reproducibility audit\n\n")
        f.write("Stage 1A asked whether the published shared-PIP analysis for `CM027530.1` can be reproduced from public intermediate data or author code without raw-read reprocessing. The answer is no: public sources provide the representative sample set, coordinates, and workflow code, but not the BAMs, SAFs, pairwise SFS files, dxy matrices, or Newick/tree objects required to regenerate the published neighbor-joining analysis directly.\n\n")
        f.write("## Frozen representative set\n\n")
        f.write(f"- Samples: {len(sample_rows)}\n")
        f.write(f"- Species: {len(species_counts)}\n")
        f.write(f"- Genotypes: AA={geno_counts.get('AA', 0)}, BB={geno_counts.get('BB', 0)}\n")
        f.write(f"- Design: one AA and one BB representative per species: `{validation['one_AA_and_one_BB_per_species']}`\n")
        f.write("- Selection source: Dataset S5, not a new Stage-1A random draw. The author code contains a random picker, but the published representatives are fixed in Dataset S5.\n\n")
        f.write("## Shared interval\n\n")
        f.write(f"- Reference: `{REFERENCE}` (*{REFERENCE_SPECIES}*)\n")
        f.write(f"- Chromosome/scaffold: `{FOCAL_SHARED}`\n")
        f.write(f"- Derived shared overlap: `{FOCAL_SHARED}:{shared_start}-{shared_end}`\n")
        f.write(f"- Length: `{shared_len}` bp under the source end-minus-start convention used by Dataset S1.\n")
        f.write("- Rule: `max(PIP.start)` and `min(PIP.end)` across the 11 species-specific `CM027530.1` PIPs in Dataset S1, matching the article Methods statement that shared PIP regions are the largest possible overlap among species-level PIPs.\n")
        f.write("- All 11 species-specific PIPs are mapped to `GCA_001746935.2` (*Setophaga coronata*).\n")
        f.write("- Exact author shared-COMP interval was not present in Dataset S1, Dataset S5, or packaged B3 code; B3a expects it in a user-provided `shared_pips_meta` file.\n\n")
        f.write("## Reusable intermediates\n\n")
        f.write("No reusable public CM027530.1 VCF/BCF, BEAGLE genotype-likelihood file, SAF, pairwise SFS, dxy matrix, BAM/CRAM, or Newick tree was found. Figshare 27284553 is a preprocessing/code repository and contains no CM027530.1 shared-PIP intermediate data. NCBI BioProjects list raw sequence reads.\n\n")
        f.write("## Cheapest Stage 1B route\n\n")
        f.write("First try to obtain author intermediate BAM/SAF/dxy files or the `shared_pips_meta` table for `CM027530.1`. If those are unavailable, a limited selected-sample regional reprocessing workflow is required: retrieve the 22 samples' raw reads or BAMs if authors can provide them, align/index to `GCA_001746935.2`, run ANGSD SAF over the shared PIP and matched COMP region, compute pairwise winsfs/dxy, and only then reconstruct NJ/local-tree outputs. Do not run this automatically from Stage 1A.\n")

    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
