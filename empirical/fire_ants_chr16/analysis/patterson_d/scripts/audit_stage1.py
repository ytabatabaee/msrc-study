#!/usr/bin/env python3
"""Audit the published Solenopsis VCF, sample metadata, and BUSCO-window mapping.

Run from repository root. The optional --scan-vcf streams the 1.5 GB file and
maps variant positions into BUSCO gene bodies from Supplementary Data 2.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import sys
import zipfile
from bisect import bisect_right
from collections import Counter, defaultdict
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / "empirical/fire_ants_chr16/analysis/patterson_d"
RESULTS = BASE / "results"
VCF = ROOT / "data/fire_ants_chr16/raw/genotypes/gt.vcf.gz"
SUPP1 = ROOT / "data/fire_ants_chr16/raw/upstream/nature_communications_2022/41467_2022_28806_MOESM4_ESM_supplementary_data_1.xlsx"
SUPP2 = ROOT / "data/fire_ants_chr16/raw/upstream/nature_communications_2022/41467_2022_28806_MOESM5_ESM_supplementary_data_2.xlsx"
WINDOWS = ROOT / "data/fire_ants_chr16/processed/stage1_window_index.tsv"
TIP_GROUPS = ROOT / "data/fire_ants_chr16/processed/stage6c_individual_to_group.tsv"
TREE_TIPS = ROOT / "data/fire_ants_chr16/processed/stage6c_tree_tip_inventory.tsv"
TREE = ROOT / "data/fire_ants_chr16/processed/stage3_background_tree.nwk"
RANDOM_ALLELE_SOURCE_IDS = [
    "USP3-1-bigB-m", "AR223-1-bigB-m", "Copa2-1-bigB-m", "Par1-1-bigB-m",
    "AR112-1-bigB-p", "AR142-1-bigB-m", "SS1-2-bigB-m",
]

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

def spreadsheet_rows(path):
    """Read the first worksheet of an xlsx using only the Python standard library."""
    with zipfile.ZipFile(str(path)) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared = ["".join(t.text or "" for t in si.findall(".//m:t", NS))
                      for si in root.findall("m:si", NS)]
        sheet = ET.fromstring(z.read("xl/worksheets/sheet1.xml"))
        rows = []
        for row in sheet.findall(".//m:sheetData/m:row", NS):
            out = {}
            for cell in row.findall("m:c", NS):
                value = cell.find("m:v", NS)
                value = value.text if value is not None else ""
                if cell.attrib.get("t") == "s" and value:
                    value = shared[int(value)]
                col = re.match(r"[A-Z]+", cell.attrib["r"]).group(0)
                out[col] = value
            rows.append(out)
        return rows

def read_tsv(path):
    with open(str(path), newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path, rows, fieldnames=None):
    if not rows and not fieldnames:
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(str(path), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({k: clean(row.get(k, "")) for k in fieldnames})

def clean(value):
    if value is None:
        return ""
    return str(value).replace("\t", " ").replace("\r", " ").replace("\n", " ")

def sha256(path):
    h = hashlib.sha256()
    with open(str(path), "rb") as f:
        while True:
            b = f.read(4 * 1024 * 1024)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def parse_vcf_header(path):
    ref = ""
    contigs = []
    with gzip.open(str(path), "rb") as f:
        for line in f:
            if line.startswith(b"##reference="):
                ref = line.decode("utf-8", "replace").strip().split("=", 1)[1]
            elif line.startswith(b"##contig=<ID="):
                m = re.match(br"##contig=<ID=([^,>]+)", line)
                if m:
                    contigs.append(m.group(1).decode("utf-8", "replace"))
            elif line.startswith(b"#CHROM\t"):
                columns = line.decode("utf-8", "replace").rstrip("\r\n").split("\t")
                return ref, contigs, columns[9:]
    raise ValueError("VCF lacks #CHROM header")

def normalize_species(x):
    return (x or "").replace("\u00a0", " ").strip()

def group_from_metadata(species, variant):
    species = normalize_species(species)
    v = (variant or "").strip()
    if species in ("S.invicta", "S.macdonaghi") and v in ("SB", "Sb"):
        return "invicta/macdonaghi_" + v
    if species == "S.richteri" and v in ("SB", "Sb"):
        return "richteri_" + v
    return species.replace("S.", "").replace(" ", "")

def parse_data1(path):
    rows = spreadsheet_rows(path)
    header = {v: k for k, v in rows[0].items()}
    result = {}
    for r in rows[1:]:
        sample = r.get(header["Sample name"], "").strip()
        if not sample:
            continue
        result[sample] = {
            "species": normalize_species(r.get(header["Species"], "")),
            "country": r.get(header["Country"], "").strip(),
            "latitude": r.get(header["Latitude"], ""),
            "longitude": r.get(header["Longitude"], ""),
            "colony": r.get(header["Colony Identity"], ""),
            "caste": r.get(header["Caste"], "").strip(),
            "supergene_variant": r.get(header["Supergene Variant"], "").strip(),
            "source": r.get(header["Source"], ""),
            "bioproject": r.get(header["BioProject"], ""),
            "sra_identifier": r.get(header["SRA identifier"], ""),
            "read_count_millions": r.get(header["Number of sequence reads (M)"], ""),
            "genome_coverage": r.get(header["Genome coverage"], ""),
            "mapping_rate": r.get(header["Map Reads %"], ""),
            "rflp_result_colony": r.get(header["RFLP result (colony)"], ""),
            "gp9_variant_1": r.get(header["Gp-9 variant 1"], ""),
            "gp9_variant_2": r.get(header["Gp-9 variant 2"], ""),
            "gp9_variant_3": r.get(header["Gp-9 variant 3"], ""),
            "fraction_missing": r.get(header["Fraction of sites missing"], ""),
            "used_for_twisst": r.get(header["Used for Twisst"], ""),
            "used_for_bbbaa_bbaba": r.get(header["Used for BBBAA-BBABA (invicta and richteri)"], ""),
        }
    return result

def parse_data2(path):
    rows = spreadsheet_rows(path)
    header = {v: k for k, v in rows[0].items()}
    genes = {}
    for r in rows[1:]:
        gene = r.get(header["BUSCO_gene_ID_Hymenoptera_ODB10"], "").strip()
        if not gene:
            continue
        genes[gene] = {
            "genomic_region": r.get(header["Genomic_region"], ""),
            "scaffold": r.get(header["Scaffold"], ""),
            "start": int(float(r.get(header["Position_start"], "0"))),
            "end": int(float(r.get(header["Position_end"], "0"))),
            "profile_status": r.get(header["BUSCO_profile_status"], ""),
        }
    return genes

def random_treatment_matches(samples):
    matches = defaultdict(list)
    for legacy in RANDOM_ALLELE_SOURCE_IDS:
        c = [s for s in samples if s == legacy or s.startswith(legacy + "-")]
        matches[legacy] = c
    return matches

def build_and_write_audits(scan_vcf=False):
    RESULTS.mkdir(parents=True, exist_ok=True)
    data1 = parse_data1(SUPP1)
    data2 = parse_data2(SUPP2)
    vcf_ref, vcf_contigs, vcf_samples = parse_vcf_header(VCF)
    windows = read_tsv(WINDOWS)
    tree_map = {r["sample"]: r for r in read_tsv(TIP_GROUPS)}
    tips = {r["tree_tip"] for r in read_tsv(TREE_TIPS)}
    random_matches = random_treatment_matches(vcf_samples)
    random_samples = {s for vals in random_matches.values() for s in vals}

    # Source and integrity inventory.
    inventory = []
    for path, source, url, fmt in [
        (VCF, "Stolle et al. 2022 public genotype VCF", "https://github.com/wurmlab/wurmlab.github.io/raw/master/data/supergene_introgression/gt.vcf.gz/gt.vcf.gz", "BGZF/gzip VCFv4.2"),
        (SUPP1, "Stolle et al. 2022 Supplementary Data 1", "https://pmc.ncbi.nlm.nih.gov/articles/instance/8917144/bin/41467_2022_28806_MOESM4_ESM.xlsx", "XLSX sample metadata"),
        (SUPP2, "Stolle et al. 2022 Supplementary Data 2", "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-28806-7/MediaObjects/41467_2022_28806_MOESM5_ESM.xlsx", "XLSX BUSCO/scaffold coordinate metadata"),
    ]:
        inventory.append({"local_path": str(path.relative_to(ROOT)), "source": source, "url": url,
                          "format": fmt, "bytes": path.stat().st_size, "sha256": sha256(path)})
    write_tsv(RESULTS / "genotype_source_inventory.tsv", inventory)

    expected_sha = "e782b2e523357ee09ffe81ce0cc2c6b25f733159ba965c638cc03c6785306775"
    expected_bytes = 1503782377
    vcf_bytes = VCF.stat().st_size
    vcf_sha = inventory[0]["sha256"]
    # gzip integrity is performed separately with gzip -t and recorded in this source audit.
    validation = [{
        "source_url": inventory[0]["url"], "lfs_pointer_url": "https://raw.githubusercontent.com/wurmlab/wurmlab.github.io/master/data/supergene_introgression/gt.vcf.gz/gt.vcf.gz",
        "lfs_pointer_bytes": 135, "lfs_expected_bytes": expected_bytes, "actual_bytes": vcf_bytes,
        "bytes_match": vcf_bytes == expected_bytes, "expected_sha256": expected_sha,
        "actual_sha256": vcf_sha, "sha256_match": vcf_sha == expected_sha,
        "gzip_integrity": "PASS (gzip -t)", "vcf_version": "VCFv4.2", "reference_header": vcf_ref,
        "sample_count": len(vcf_samples), "contig_header_count": len(vcf_contigs),
        "remote_tabix_index": "not found (HTTP 404 for published .tbi endpoint)",
        "local_tabix_index": "absent", "local_tabix_tools": "bcftools/tabix/bgzip not installed",
    }]
    write_tsv(RESULTS / "download_validation.tsv", validation)

    # Mapping to the sample supplement and to the exact 267-tip local-tree sample set.
    missing_meta = sorted(set(vcf_samples) - set(data1))
    unexpected_supp = sorted(set(data1) - set(vcf_samples))
    if missing_meta:
        raise ValueError("VCF sample IDs lack exact Supplementary Data 1 matches: " + ",".join(missing_meta[:10]))
    if len(tips) != 267 or not tips.issubset(set(vcf_samples)):
        raise ValueError("Expected exact 267 local-tree tips all present in VCF")
    if len(set(vcf_samples)) != len(vcf_samples):
        raise ValueError("Duplicate sample ID in VCF header")
    random_unresolved = [x for x, y in random_matches.items() if len(y) != 1]

    sample_rows = []
    map_rows = []
    for sid in vcf_samples:
        m = data1[sid]
        tree = tree_map.get(sid, {})
        vcf_group = tree.get("astral_group", "") or group_from_metadata(m["species"], m["supergene_variant"])
        miss = float(m["fraction_missing"]) if m["fraction_missing"] else float("nan")
        random_flag = sid in random_samples
        pooled = m["caste"].lower() == "workers"
        twisst = m["used_for_twisst"] == "1"
        map_rows.append({
            "vcf_sample_id": sid, "supplement_sample_id": sid, "supplement_match": "exact",
            "local_tree_tip": int(sid in tips), "twisst_selected": int(twisst),
            "twisst_group": vcf_group, "species": m["species"], "country": m["country"],
            "latitude": m["latitude"], "longitude": m["longitude"], "colony": m["colony"],
            "caste": m["caste"], "supergene_variant_source": m["supergene_variant"],
            "source_sample_missing_fraction": m["fraction_missing"], "source_bioproject": m["bioproject"],
            "source_sra_identifier": m["sra_identifier"], "random_allele_treatment": int(random_flag),
            "pooled_worker_sample": int(pooled),
            "primary_haploid_candidate": int(sid in tips and m["caste"].lower() == "male" and math.isfinite(miss) and miss < 0.05 and not random_flag),
            "group_assignment_rule": "Supplementary Data 1 Supergene Variant + Species; cross-check stage6c mapping",
        })
        sample_rows.append({
            "vcf_sample_id": sid, "species": m["species"], "country": m["country"],
            "latitude": m["latitude"], "longitude": m["longitude"], "colony": m["colony"],
            "caste": m["caste"], "supergene_variant": m["supergene_variant"],
            "source": m["source"], "bioproject": m["bioproject"], "sra_identifier": m["sra_identifier"],
            "read_count_millions": m["read_count_millions"], "genome_coverage": m["genome_coverage"],
            "mapping_rate": m["mapping_rate"], "fraction_missing": m["fraction_missing"],
            "rflp_result_colony": m["rflp_result_colony"], "gp9_variant_1": m["gp9_variant_1"],
            "gp9_variant_2": m["gp9_variant_2"], "gp9_variant_3": m["gp9_variant_3"],
            "local_tree_tip": int(sid in tips), "twisst_selected": int(twisst), "twisst_group": vcf_group,
            "pooled_worker_sample": int(pooled), "random_allele_treatment": int(random_flag),
            "sample_quality_flag": "exclude_pool" if pooled else ("exclude_random_allele" if random_flag else ("eligible_if_missing_lt_0.05" if sid in tips else "not_in_267_tip_tree")),
        })
    write_tsv(RESULTS / "vcf_sample_inventory.tsv", sample_rows)
    write_tsv(RESULTS / "sample_mapping_audit.tsv", map_rows)

    # Exact BUSCO gene identity bridge: VCF scaffold-position -> Supplementary Data 2 BUSCO ID -> frozen tree window.
    gene_window = {}
    for w in windows:
        for i in range(1, 5):
            gene = w["gene%d" % i]
            if gene in gene_window:
                raise ValueError("A BUSCO gene appears in more than one frozen window: " + gene)
            gene_window[gene] = w
    absent_genes = sorted(set(gene_window) - set(data2))
    if absent_genes:
        raise ValueError("Window BUSCO IDs absent from Supplementary Data 2: " + ",".join(absent_genes[:10]))
    if len(gene_window) != 852 or len(set(gene_window) & set(data2)) != 852:
        raise ValueError("Expected 852 unique BUSCO identifiers across 213 windows")
    contig_set = set(vcf_contigs)
    gene_rows = []
    intervals = defaultdict(list)
    for gene, w in gene_window.items():
        g = data2[gene]
        in_vcf = g["scaffold"] in contig_set
        if g["start"] < 1 or g["end"] < g["start"]:
            raise ValueError("Invalid Supplementary Data 2 coordinates for " + gene)
        if not in_vcf:
            raise ValueError("BUSCO scaffold missing from VCF contig header: " + g["scaffold"])
        intervals[g["scaffold"]].append((g["start"], g["end"], gene, w["window_index"]))
        gene_rows.append({
            "window_index": w["window_index"], "tree_chromosome": w["chrom"], "tree_region": w["region"],
            "tree_window_start": w["start"], "tree_window_end": w["end"], "tree_window_mid": w["mid"],
            "gene_slot": "gene%d" % next(i for i in range(1, 5) if w["gene%d" % i] == gene),
            "busco_gene_id": gene, "supplement_genomic_region": g["genomic_region"],
            "vcf_contig": g["scaffold"], "vcf_gene_start_1based": g["start"],
            "vcf_gene_end_1based_inclusive": g["end"], "vcf_contig_present": int(in_vcf),
            "mapping_rule": "exact BUSCO gene ID in frozen four-gene window; VCF contig/position from Supplementary Data 2; no chromosome-scale offset/strand transform",
        })
    write_tsv(RESULTS / "window_gene_coordinate_mapping.tsv", gene_rows)

    scan_summary = None
    if scan_vcf:
        scan_summary = scan_records(intervals)
        write_tsv(RESULTS / "vcf_record_scan_by_contig.tsv", scan_summary["contig_rows"])
        write_tsv(RESULTS / "window_vcf_coverage.tsv", scan_summary["window_rows"])
        (RESULTS / "vcf_scan_summary.json").write_text(json.dumps(scan_summary["summary"], indent=2) + "\n")

    # Tabular info for reuse by the second-stage analysis.
    target_scaffolds = {data2[g]["scaffold"] for g in gene_window}
    prior_scan_path = RESULTS / "vcf_scan_summary.json"
    prior_scan_complete = False
    if prior_scan_path.exists():
        try:
            prior_scan_complete = int(json.loads(prior_scan_path.read_text()).get("vcf_variant_records", 0)) > 0
        except (ValueError, OSError):
            prior_scan_complete = False
    manifest = {
        "vcf_path": str(VCF.relative_to(ROOT)), "vcf_bytes": vcf_bytes, "vcf_sha256": vcf_sha,
        "vcf_sha256_expected": expected_sha, "vcf_sha256_match": vcf_sha == expected_sha,
        "vcf_reference_header": vcf_ref, "vcf_format": "VCFv4.2 BGZF/gzip",
        "vcf_contig_count": len(vcf_contigs), "vcf_sample_count": len(vcf_samples),
        "supplementary_data_1_rows": len(data1), "supplementary_rows_not_in_vcf": len(unexpected_supp),
        "vcf_samples_without_supplement_match": len(missing_meta),
        "tree_tip_count": len(tips), "tree_tips_exactly_in_vcf": len(tips & set(vcf_samples)),
        "random_allele_source_ids": RANDOM_ALLELE_SOURCE_IDS,
        "random_allele_matches": dict(random_matches), "random_allele_source_ids_unresolved": random_unresolved,
        "supplementary_data_2_busco_rows": len(data2), "frozen_windows": len(windows),
        "frozen_window_gene_ids": len(gene_window), "window_genes_in_supplementary_data_2": len(set(gene_window) & set(data2)),
        "window_gene_scaffolds_expected": len(target_scaffolds),
        "window_gene_scaffolds_in_vcf_header": len(target_scaffolds & contig_set),
        "exact_gene_identity_mapping": True, "vcf_index_present": False,
        "record_scan_complete": bool(scan_vcf or prior_scan_complete),
    }
    (RESULTS / "stage1_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return {"data1": data1, "data2": data2, "vcf_samples": vcf_samples,
            "vcf_contigs": vcf_contigs, "vcf_ref": vcf_ref, "windows": windows,
            "tree_map": tree_map, "tips": tips, "random_samples": random_samples,
            "random_matches": random_matches, "manifest": manifest, "scan": scan_summary}

def scan_records(gene_intervals):
    """Stream full VCF; assign positions only within exact BUSCO gene bodies."""
    by_contig = {}
    for contig, ivals in gene_intervals.items():
        ivals.sort(key=lambda z: (z[0], z[1]))
        starts = [x[0] for x in ivals]
        prefix_max_end = []
        best = 0
        for x in ivals:
            best = max(best, x[1])
            prefix_max_end.append(best)
        by_contig[contig] = (starts, prefix_max_end, ivals)
    contig_counts = Counter()
    window_counts = Counter()
    window_gene_hits = defaultdict(set)
    n_records = n_gene_body = n_target_mapped = n_ambiguous = n_unmapped_contig = 0
    with gzip.open(str(VCF), "rb") as f:
        for line in f:
            if not line or line[0] == 35:  # #
                continue
            n_records += 1
            a = line.find(b"\t")
            b = line.find(b"\t", a + 1)
            contig = line[:a].decode("utf-8", "replace")
            try:
                pos = int(line[a + 1:b])
            except ValueError:
                continue
            contig_counts[contig] += 1
            info = by_contig.get(contig)
            if info is None:
                n_unmapped_contig += 1
                continue
            starts, prefix_max_end, ivals = info
            j = bisect_right(starts, pos) - 1
            hits = []
            while j >= 0 and prefix_max_end[j] >= pos:
                if ivals[j][0] <= pos <= ivals[j][1]:
                    hits.append(ivals[j])
                j -= 1
            if hits:
                n_gene_body += 1
                if len(hits) > 1:
                    n_ambiguous += 1
                mapped = [x for x in hits if x[3]]
                if mapped:
                    n_target_mapped += 1
                    for hit in mapped:
                        window_counts[hit[3]] += 1
                        window_gene_hits[hit[3]].add(hit[2])
    contig_rows = [{"contig": c, "variant_records": n} for c, n in sorted(contig_counts.items())]
    window_rows = []
    for w in read_tsv(WINDOWS):
        wi = w["window_index"]
        window_rows.append({
            "window_index": wi, "chromosome": w["chrom"], "region": w["region"],
            "window_start": w["start"], "window_end": w["end"], "window_mid": w["mid"],
            "busco_gene_count": 4, "busco_genes_with_vcf_records": len(window_gene_hits[wi]),
            "vcf_variant_records_in_busco_gene_bodies": window_counts[wi],
        })
    summary = {
        "vcf_variant_records": n_records, "vcf_contigs_with_records": len(contig_counts),
        "variant_records_overlapping_any_of_852_tree_window_busco_genes": n_gene_body,
        "variant_records_in_213_tree_window_busco_genes": n_target_mapped,
        "variant_records_overlapping_multiple_busco_genes": n_ambiguous,
        "records_on_contigs_without_window_gene_intervals": n_unmapped_contig,
        "assignment_scope": "exact gene-body intervals; no 1-kb flanks or unsupported scaffold-to-chromosome coordinate lift",
    }
    return {"summary": summary, "contig_rows": contig_rows, "window_rows": window_rows}

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-vcf", action="store_true", help="stream the complete compressed VCF to count records and map BUSCO gene positions")
    args = parser.parse_args()
    data = build_and_write_audits(scan_vcf=args.scan_vcf)
    print(json.dumps(data["manifest"], indent=2))
    if args.scan_vcf:
        print(json.dumps(data["scan"]["summary"], indent=2))
