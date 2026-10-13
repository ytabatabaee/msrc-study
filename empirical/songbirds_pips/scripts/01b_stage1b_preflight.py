#!/usr/bin/env python3
"""Stage 1B preflight and failure-safe scaffolding for shared_CM027530.1.

This script validates the frozen Stage-1A design and records whether the local
environment can execute the requested raw-read pilot. It deliberately does not
substitute samples, alter the frozen shared interval, or fabricate downstream
ANGSD/winsfs/NJ outputs when required executables are missing.
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "songbirds_pips"
EMP = ROOT / "empirical" / "songbirds_pips"

SAMPLES = DATA / "processed" / "stage1a_shared_pip_samples.tsv"
INTERVALS = DATA / "processed" / "stage1a_shared_pip_interval.tsv"

RAW_DIR = DATA / "raw" / "stage1b_shared_pip"
INTERMEDIATE = DATA / "intermediate" / "stage1b_shared_pip"
PROCESSED = DATA / "processed"
RESULTS = EMP / "results"
FIGURES = EMP / "figures"
LOGS = INTERMEDIATE / "logs"

FOCAL_ID = "shared_CM027530.1"
FOCAL_CHROM = "CM027530.1"
FOCAL_START = 574805
FOCAL_END = 874807
FOCAL_LENGTH = 300002
REFERENCE = "GCA_001746935.2"
EXPECTED_N = 22
EXPECTED_SPECIES = 11
EXPECTED_GENOTYPES = {"AA": 11, "BB": 11}
EXPECTED_PAIRS = 231

REQUIRED_TOOLS = [
    "prefetch",
    "fasterq-dump",
    "bwa",
    "samtools",
    "angsd",
    "winsfs",
]
OPTIONAL_TOOLS = [
    "AdapterRemoval",
    "fastp",
    "bam",
    "picard",
    "gatk",
    "java",
    "Rscript",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def tool_version(tool: str, path: str | None) -> str:
    if path is None:
        return "not_found"
    commands = {
        "prefetch": [path, "--version"],
        "fasterq-dump": [path, "--version"],
        "bwa": [path],
        "samtools": [path, "--version"],
        "angsd": [path, "-version"],
        "winsfs": [path, "--version"],
        "AdapterRemoval": [path, "--version"],
        "fastp": [path, "--version"],
        "bam": [path, "--version"],
        "picard": [path, "--version"],
        "gatk": [path, "--version"],
        "java": [path, "-version"],
        "Rscript": [path, "--version"],
    }
    try:
        proc = subprocess.run(
            commands.get(tool, [path, "--version"]),
            check=False,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except Exception as exc:  # pragma: no cover - defensive provenance capture
        return f"version_check_failed: {exc}"
    out = (proc.stdout + proc.stderr).strip().splitlines()
    return " | ".join(out[:3]) if out else f"exit_{proc.returncode}_no_output"


def validate_design(samples: list[dict[str, str]], intervals: list[dict[str, str]]) -> dict[str, object]:
    genotype_counts = Counter(row["pip_genotype"] for row in samples)
    species_counts = Counter(row["species"] for row in samples)
    assembly_ok = all(REFERENCE in row["reference_assembly"] for row in samples)
    accessions_present = all(row.get("sra_accession", "").strip() for row in samples)
    one_aa_bb = all(
        Counter(row["pip_genotype"] for row in samples if row["species"] == species) == {"AA": 1, "BB": 1}
        for species in species_counts
    )
    frozen_interval = [row for row in intervals if row["interval_id"] == FOCAL_ID]
    interval_ok = False
    interval_record = {}
    if len(frozen_interval) == 1:
        interval_record = frozen_interval[0]
        interval_ok = (
            interval_record.get("chromosome") == FOCAL_CHROM
            and int(interval_record.get("start", -1)) == FOCAL_START
            and int(interval_record.get("end", -1)) == FOCAL_END
            and int(interval_record.get("length_bp", -1)) == FOCAL_LENGTH
            and REFERENCE in interval_record.get("reference_assembly", "")
        )
    checks = {
        "exactly_22_samples": len(samples) == EXPECTED_N,
        "eleven_species": len(species_counts) == EXPECTED_SPECIES,
        "aa_11_bb_11": dict(genotype_counts) == EXPECTED_GENOTYPES,
        "one_AA_and_one_BB_per_species": one_aa_bb,
        "all_accessions_present": accessions_present,
        "accession_namespace": "SAMN/BioSample values in frozen sra_accession column; SRR run accession resolution not yet performed",
        "reference_assembly_GCA_001746935_2": assembly_ok,
        "focal_interval_matches_stage1a": interval_ok,
        "sample_count": len(samples),
        "species_count": len(species_counts),
        "genotype_counts": dict(genotype_counts),
        "interval": interval_record,
    }
    checks["design_valid"] = all(
        bool(checks[key])
        for key in [
            "exactly_22_samples",
            "eleven_species",
            "aa_11_bb_11",
            "one_AA_and_one_BB_per_species",
            "all_accessions_present",
            "reference_assembly_GCA_001746935_2",
            "focal_interval_matches_stage1a",
        ]
    )
    return checks


def write_logs(samples: list[dict[str, str]], tool_rows: list[dict[str, str]]) -> None:
    LOGS.mkdir(parents=True, exist_ok=True)
    command_log = LOGS / "stage1b_planned_commands.sh"
    accessions = " ".join(row["sra_accession"] for row in samples)
    sample_ids = " ".join(row["sample_id"] for row in samples)
    command_log.write_text(
        "\n".join(
            [
                "#!/usr/bin/env bash",
                "set -euo pipefail",
                "",
                "# Stage 1B planned commands. Not executed in the current run because required tools were missing.",
                "# Frozen BioSample accessions from stage1a_shared_pip_samples.tsv:",
                f"# {accessions}",
                "# Frozen sample IDs:",
                f"# {sample_ids}",
                "",
                "# 1. Resolve SAMN BioSample accessions to run-level SRR accessions before prefetch.",
                "# 2. Download only those resolved SRR accessions with prefetch/fasterq-dump.",
                "# prefetch <SRR> --output-directory data/songbirds_pips/raw/stage1b_shared_pip/sra/",
                "# fasterq-dump <SRR> --split-files --outdir data/songbirds_pips/raw/stage1b_shared_pip/fastq/",
                "",
                "# 3. Obtain GCA_001746935.2 Setophaga coronata reference FASTA and build indices.",
                "# bwa index <GCA_001746935.2.fa>",
                "# samtools faidx <GCA_001746935.2.fa>",
                "",
                "# 4. Align with the author-matching mapper/filtering where available.",
                "# bwa mem -t <threads> <ref.fa> <R1.fastq.gz> <R2.fastq.gz> | samtools sort -@ <threads> -o <sample>.bam",
                "# samtools index <sample>.bam",
                "# samtools view -b <sample>.bam CM027530.1:574805-874807 > <sample>.CM027530.1_574805_874807.bam",
                "",
                "# 5. Generate per-sample SAFs for the frozen interval.",
                "# angsd -i <regional.bam> -out <sample> -doSaf 1 -GL 2 -minMapQ 30 -minQ 30 -anc <ref.fa> -r CM027530.1:574805-874807 -sites <sites>",
                "",
                "# 6. Run all 231 pairwise winsfs analyses and author dxy_sfs.py.",
                "# winsfs <sample1.saf.idx> <sample2.saf.idx> > <pair>.sfs",
                "# python dxy_sfs.py -m 1 -n 1 -s <pair>.sfs -o <pair>.dxy",
                "",
                "# 7. Build distance matrix and ape::nj tree only after dxy is complete.",
                "",
            ]
        )
        + "\n"
    )
    write_tsv(
        LOGS / "stage1b_preflight_tool_versions.tsv",
        tool_rows,
        ["tool", "path", "version", "required_for_stage1b"],
    )


def append_once(path: Path, marker: str, text: str) -> None:
    existing = path.read_text() if path.exists() else ""
    if marker not in existing:
        with path.open("a") as handle:
            if existing and not existing.endswith("\n"):
                handle.write("\n")
            handle.write(text)


def main() -> None:
    for directory in [
        RAW_DIR / "sra",
        RAW_DIR / "fastq",
        INTERMEDIATE / "reference",
        INTERMEDIATE / "bam",
        INTERMEDIATE / "regional_bam",
        INTERMEDIATE / "saf",
        INTERMEDIATE / "sfs",
        LOGS,
        PROCESSED,
        RESULTS,
        FIGURES,
    ]:
        directory.mkdir(parents=True, exist_ok=True)

    samples = read_tsv(SAMPLES)
    intervals = read_tsv(INTERVALS)
    design = validate_design(samples, intervals)

    tool_rows: list[dict[str, str]] = []
    missing_required: list[str] = []
    for tool in REQUIRED_TOOLS + OPTIONAL_TOOLS:
        path = shutil.which(tool)
        if tool in REQUIRED_TOOLS and path is None:
            missing_required.append(tool)
        tool_rows.append(
            {
                "tool": tool,
                "path": path or "not_found",
                "version": tool_version(tool, path),
                "required_for_stage1b": "yes" if tool in REQUIRED_TOOLS else "optional_or_author_pipeline",
            }
        )

    blocked = bool(missing_required)
    timestamp = datetime.now(timezone.utc).isoformat()

    pairs = list(combinations(samples, 2))
    pair_rows = []
    dxy_rows = []
    for left, right in pairs:
        pair_id = f"{left['sample_id']}__{right['sample_id']}"
        note = "not_run_preflight_blocked_missing_tools" if blocked else "pending_execution"
        pair_rows.append(
            {
                "pair_id": pair_id,
                "sample1": left["sample_id"],
                "sample2": right["sample_id"],
                "species1": left["species"],
                "species2": right["species"],
                "genotype1": left["pip_genotype"],
                "genotype2": right["pip_genotype"],
                "sfs_path": "",
                "status": note,
                "runtime_seconds": "",
                "notes": ",".join(missing_required) if blocked else "",
            }
        )
        dxy_rows.append(
            {
                "sample1": left["sample_id"],
                "sample2": right["sample_id"],
                "species1": left["species"],
                "species2": right["species"],
                "genotype1": left["pip_genotype"],
                "genotype2": right["pip_genotype"],
                "dxy": "",
                "status": note,
            }
        )

    qc_rows = [
        {
            "sample": row["sample_id"],
            "species": row["species"],
            "pip_genotype": row["pip_genotype"],
            "mapped_reads_region": "",
            "mean_depth": "",
            "median_depth": "",
            "breadth_1x": "",
            "breadth_3x": "",
            "notes": "not_run_preflight_blocked_missing_tools" if blocked else "pending_execution",
        }
        for row in samples
    ]

    matrix_fieldnames = ["sample"] + [row["sample_id"] for row in samples]
    matrix_rows = []
    for row in samples:
        matrix_row = {"sample": row["sample_id"]}
        for col in samples:
            matrix_row[col["sample_id"]] = "0" if row["sample_id"] == col["sample_id"] else "NA"
        matrix_rows.append(matrix_row)

    signal_rows = [
        {
            "metric": "mean_distance_same_genotype_different_species",
            "value": "",
            "status": "not_run_preflight_blocked_missing_tools" if blocked else "pending_execution",
            "notes": "requires complete dxy matrix",
        },
        {
            "metric": "mean_distance_different_genotype_same_species",
            "value": "",
            "status": "not_run_preflight_blocked_missing_tools" if blocked else "pending_execution",
            "notes": "requires complete dxy matrix",
        },
        {
            "metric": "mean_delta_i",
            "value": "",
            "status": "not_run_preflight_blocked_missing_tools" if blocked else "pending_execution",
            "notes": "requires complete dxy matrix",
        },
    ]

    write_tsv(
        RESULTS / "stage1b_focal_region_qc.tsv",
        qc_rows,
        [
            "sample",
            "species",
            "pip_genotype",
            "mapped_reads_region",
            "mean_depth",
            "median_depth",
            "breadth_1x",
            "breadth_3x",
            "notes",
        ],
    )
    write_tsv(
        RESULTS / "stage1b_pairwise_sfs_status.tsv",
        pair_rows,
        [
            "pair_id",
            "sample1",
            "sample2",
            "species1",
            "species2",
            "genotype1",
            "genotype2",
            "sfs_path",
            "status",
            "runtime_seconds",
            "notes",
        ],
    )
    write_tsv(
        PROCESSED / "stage1b_shared_pip_dxy.tsv",
        dxy_rows,
        ["sample1", "sample2", "species1", "species2", "genotype1", "genotype2", "dxy", "status"],
    )
    write_tsv(PROCESSED / "stage1b_shared_pip_dxy_matrix.tsv", matrix_rows, matrix_fieldnames)
    write_tsv(RESULTS / "stage1b_shared_pip_signal.tsv", signal_rows, ["metric", "value", "status", "notes"])

    design_path = RESULTS / "stage1b_design_validation.json"
    checks = {
        **design,
        "expected_pair_count": EXPECTED_PAIRS,
        "enumerated_pair_count": len(pairs),
        "pair_count_231": len(pairs) == EXPECTED_PAIRS,
        "no_background_control_analysis_run": True,
        "no_candidate_set_changes": True,
        "no_tree_or_distance_results_fabricated": True,
        "required_tools_present": not missing_required,
        "missing_required_tools": missing_required,
        "stage1b_executable": design["design_valid"] and not missing_required,
    }
    design_path.write_text(json.dumps(checks, indent=2, sort_keys=True) + "\n")

    write_logs(samples, tool_rows)

    summary = RESULTS / "stage1b_shared_pip_signal_summary.txt"
    summary.write_text(
        "\n".join(
            [
                "Stage 1B shared_CM027530.1 signal summary",
                f"Timestamp UTC: {timestamp}",
                "",
                "No genotype-vs-species distance statistic was calculated in this run.",
                "The frozen 22-sample design and 231 expected pairs validate, but the raw-read pilot did not proceed because required executables are missing from PATH: "
                + ", ".join(missing_required)
                + ".",
                "Therefore there is no complete dXY matrix, no NJ tree, and no support for claiming reproduction of the published AA/BB-across-species signal.",
                "",
            ]
        )
    )

    reproduction = RESULTS / "stage1b_reproduction_report.md"
    reproduction.write_text(
        "\n".join(
            [
                "# Stage 1B reproduction report: shared_CM027530.1",
                "",
                "Stage 1B did not reach raw-read processing, ANGSD SAF generation, winsfs pairwise SFS estimation, dXY calculation, or NJ plotting in this environment.",
                "",
                "The frozen Stage-1A design validated: 22 representatives, 11 species, AA = 11, BB = 11, one AA and one BB representative per species, and the frozen interval `CM027530.1:574805-874807` on `GCA_001746935.2`.",
                "",
                "Execution is blocked because the required command-line tools were absent from PATH: "
                + ", ".join(missing_required)
                + ". `Rscript` is present, but the upstream raw-read/SFS stack is not.",
                "",
                "No samples were substituted, no focal interval was changed, and no background/control regions were analyzed.",
                "",
                "Because no complete dXY matrix or NJ tree exists, the published AA/BB-across-species signal has not yet been reproduced. The exact topology cannot be compared, and no anomalous species/sample can be identified from Stage 1B output.",
                "",
                "The frozen sample accession field contains SAMN BioSample identifiers rather than run-level SRR accessions. Before downloading reads, the next executable run must resolve each SAMN accession to its public run accession(s) and record that mapping.",
                "",
            ]
        )
    )

    preanalysis = RESULTS / "stage1b_stage1c_preanalysis_plan.md"
    preanalysis.write_text(
        "\n".join(
            [
                "# Stage 1C preanalysis plan placeholder",
                "",
                "Stage 1C is not justified by the present Stage 1B run because the focal raw-read pilot did not execute and the published shared-PIP signal was not reproduced.",
                "",
                "The intended Stage 1C question remains: is genotype-across-species structure at `shared_CM027530.1` stronger than in matched non-PIP regions?",
                "",
                "A final control-window rule must be frozen only after Stage 1B produces usable focal-region data, stable 231-pair SFS estimation, and a complete or nearly complete dXY matrix. Until then, no background windows should be selected, no control NJ trees should be inspected, and no Stage 1C result should be generated.",
                "",
                "Provisional independent eligibility constraints, to be revisited after successful Stage 1B:",
                "",
                "- same scaffold `CM027530.1`;",
                "- 300 kb windows outside all known PIP intervals on the scaffold;",
                "- exclude low-mappability or low-callability windows using thresholds defined from coverage/QC, before inspecting any control distance tree;",
                "- target 10-50 controls only if computationally manageable;",
                "- primary statistic: focal genotype-vs-species distance contrast compared with the distribution from eligible controls.",
                "",
            ]
        )
    )

    manifest_files = [
        EMP / "scripts" / "01b_stage1b_preflight.py",
        EMP / "scripts" / "01b_validate_stage1b.py",
        design_path,
        RESULTS / "stage1b_focal_region_qc.tsv",
        RESULTS / "stage1b_pairwise_sfs_status.tsv",
        RESULTS / "stage1b_shared_pip_signal.tsv",
        summary,
        reproduction,
        preanalysis,
        PROCESSED / "stage1b_shared_pip_dxy.tsv",
        PROCESSED / "stage1b_shared_pip_dxy_matrix.tsv",
        LOGS / "stage1b_planned_commands.sh",
        LOGS / "stage1b_preflight_tool_versions.tsv",
    ]
    manifest = {
        "stage": "1B",
        "pip_id": FOCAL_ID,
        "created_utc": timestamp,
        "status": "blocked_preflight_missing_required_tools" if blocked else "ready_for_execution",
        "frozen_inputs": {
            "samples_tsv": str(SAMPLES.relative_to(ROOT)),
            "interval_tsv": str(INTERVALS.relative_to(ROOT)),
            "reference_accession": REFERENCE,
            "interval": f"{FOCAL_CHROM}:{FOCAL_START}-{FOCAL_END}",
            "source_convention_length_bp": FOCAL_LENGTH,
        },
        "design_validation": checks,
        "tool_versions": tool_rows,
        "compute_storage": {
            "raw_fastq_downloaded": False,
            "all_22_samples_downloaded_successfully": False,
            "angsd_saf_completed_samples": 0,
            "winsfs_succeeded_pairs": 0,
            "dxy_matrix_complete": False,
            "nj_tree_created": False,
            "stage1c_justified": False,
        },
        "outputs": {},
        "not_created_because_blocked": [
            "empirical/songbirds_pips/results/stage1b_shared_pip_nj.nwk",
            "empirical/songbirds_pips/figures/songbirds_stage1b_shared_pip_nj.pdf",
            "empirical/songbirds_pips/figures/songbirds_stage1b_shared_pip_nj.png",
        ],
        "deviations_from_published_pipeline": [
            "No raw-read processing was launched because required executables were absent.",
            "Run-level SRR accessions were not resolved from the frozen SAMN BioSample accessions.",
        ],
    }
    for path in manifest_files:
        manifest["outputs"][str(path.relative_to(ROOT))] = {"sha256": sha256(path), "bytes": path.stat().st_size}

    manifest_path = RESULTS / "stage1b_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    checksum_rows = []
    for rel, meta in manifest["outputs"].items():
        checksum_rows.append({"path": rel, "sha256": meta["sha256"], "bytes": meta["bytes"]})
    checksum_rows.append(
        {
            "path": str(manifest_path.relative_to(ROOT)),
            "sha256": sha256(manifest_path),
            "bytes": manifest_path.stat().st_size,
        }
    )
    write_tsv(RESULTS / "stage1b_checksums.tsv", checksum_rows, ["path", "sha256", "bytes"])

    print(json.dumps({"status": manifest["status"], "missing_required_tools": missing_required}, indent=2))


if __name__ == "__main__":
    main()
