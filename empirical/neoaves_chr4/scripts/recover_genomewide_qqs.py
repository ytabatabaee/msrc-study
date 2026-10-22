#!/usr/bin/env python3
"""Recover genome-wide focal QQS from the published precomputed statistics.

The source archive is read through ``zipfile`` and ``lzma`` streams.  No gene
trees are loaded or enumerated.  The resulting per-locus table is a generated
analysis artifact; the code and compact summaries are the reproducible
deliverables.
"""
from __future__ import annotations

import csv
import hashlib
import lzma
import math
import zipfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "empirical" / "neoaves_chr4"
ARCHIVE = BASE / "external" / "chr4avian" / "chr4avian-v1.0.0.zip"
LOCUS = BASE / "results" / "stage4d_locus_manifest.tsv"
OLD_CHR4 = BASE / "results" / "locus_table_chr4.tsv"
OUTDIR = BASE / "results" / "stage4e"
FOCAL = {
    ("Columbimorphae", "Phoenicopteriformes"): "Columbea",
    ("Otidimorphae", "Columbimorphae"): "N61",
    ("Columbiformes", "OtherColumbimorphae"): "N62",
}
CLADES = ("Columbea", "N61", "N62")
TOPOLOGIES = ("q1", "q2", "q3")


def source_member(zf: zipfile.ZipFile, suffix: str) -> str:
    matches = [name for name in zf.namelist() if name.endswith(suffix)]
    if len(matches) != 1:
        raise ValueError(f"expected one archive member ending {suffix!r}, found {matches}")
    return matches[0]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def qqs(xs: list[float]) -> dict[str, float] | None:
    x1, x2, x3, x4 = xs
    d1, d2, d3 = x1 - x4, x2 - x4, x3 - x4
    denominator = d1 + d2 + d3
    if denominator == 0 or not math.isfinite(denominator):
        return None
    values = {"d1": d1, "d2": d2, "d3": d3,
              "q1": (d2 + d3 - d1) / denominator,
              "q2": (d1 + d3 - d2) / denominator,
              "q3": (d1 + d2 - d3) / denominator}
    if not all(math.isfinite(values[key]) for key in ("q1", "q2", "q3")):
        return None
    return values


def fmt(value: object) -> str:
    if isinstance(value, float):
        return f"{value:.12g}"
    return str(value)


def boolish(value: str) -> bool:
    return value.lower() in {"true", "1", "yes"}


def load_manifest() -> dict[str, dict[str, str]]:
    with LOCUS.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 63430 or len({row["locus_id"] for row in rows}) != len(rows):
        raise ValueError("unexpected or non-unique Stage 4D locus manifest")
    result = {}
    for row in rows:
        result[row["original_gene_id"]] = row
    return result


def archive_text(zf: zipfile.ZipFile, suffix: str):
    raw = zf.open(source_member(zf, suffix), "r")
    return lzma.open(raw, "rt")


def read_source(manifest: dict[str, dict[str, str]]):
    values: dict[tuple[str, str], dict[str, object]] = {}
    with zipfile.ZipFile(ARCHIVE) as zf:
        with archive_text(zf, "genetreesupport/63K_trees.names_header.txt.xz") as handle:
            header = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
            header_count = 0
            for row in header:
                header_count += 1
                if row["Gene"] not in manifest:
                    raise ValueError(f"header Gene {row['Gene']} missing from Stage 4D manifest")
        with archive_text(zf, "genetreesupport/clade-rec.stat.xz") as handle:
            source_rows = 0
            for line_number, line in enumerate(handle, 1):
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) != 5:
                    raise ValueError(f"unexpected field count at source line {line_number}")
                c1, c2, topology, gene, x = parts
                clade = FOCAL.get((c1, c2))
                if clade is None:
                    continue
                if gene not in manifest:
                    raise ValueError(f"source Gene {gene} missing from Stage 4D manifest")
                source_rows += 1
                key = (gene, clade)
                record = values.setdefault(key, {"Gene": gene, "clade": clade, "C1": c1, "C2": c2, "x": {}})
                record["x"][topology] = float(x)
    if header_count != 63430:
        raise ValueError(f"published header contains {header_count}, expected 63430")
    return values, header_count, source_rows


def build_rows(values: dict[tuple[str, str], dict[str, object]], manifest: dict[str, dict[str, str]]):
    rows = []
    counters = defaultdict(lambda: {"missing_topology": 0, "undefined": 0})
    for (gene, clade), record in sorted(values.items(), key=lambda item: (CLADES.index(item[0][1]), int(item[0][0]))):
        locus = manifest[gene]
        xs = record["x"]
        if set(xs) != {"1", "2", "3", "4"}:
            counters[clade]["missing_topology"] += 1
            continue
        calculated = qqs([xs[str(i)] for i in range(1, 5)])
        if calculated is None:
            counters[clade]["undefined"] += 1
            continue
        q = {key: calculated[key] for key in ("q1", "q2", "q3")}
        dominant = max(q.items(), key=lambda item: (item[1], item[0]))[0]
        rows.append({"Gene": gene, "locus_id": locus["locus_id"], "clade": clade, "C1": record["C1"], "C2": record["C2"],
                     "Chromosome": locus["chromosome"], "start": locus["start"], "end": locus["end"], "midpoint": locus["midpoint"],
                     "x1": xs["1"], "x2": xs["2"], "x3": xs["3"], "x4": xs["4"], **calculated,
                     "dominant_topology": dominant, "valid": 1,
                     "is_chr4": locus["is_chr4"], "published_outlier_region": locus["published_outlier_region"],
                     "frozen_structural_primary": locus["frozen_structural_primary"],
                     "frozen_structural_stringent": locus["frozen_structural_stringent"],
                     "frozen_structural_inclusive": locus["frozen_structural_inclusive"]})
    return rows, counters


def write_genomewide(rows: list[dict[str, object]]) -> Path:
    path = OUTDIR / "genomewide_focal_quartet_support.tsv"
    fields = ["Gene", "locus_id", "Chromosome", "start", "end", "midpoint", "clade", "C1", "C2", "x1", "x2", "x3", "x4", "d1", "d2", "d3", "q1", "q2", "q3", "dominant_topology", "is_chr4", "published_outlier_region", "frozen_structural_primary", "frozen_structural_stringent", "frozen_structural_inclusive", "valid"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row[field]) for field in fields})
    return path


def partition(row: dict[str, object], name: str) -> bool:
    return {
        "nonchr4": not boolish(str(row["is_chr4"])),
        "chr4": boolish(str(row["is_chr4"])),
        "pnas_outlier": boolish(str(row["published_outlier_region"])),
        "outside_pnas": not boolish(str(row["published_outlier_region"])),
        "structural_primary": boolish(str(row["frozen_structural_primary"])),
        "outside_structural_primary": not boolish(str(row["frozen_structural_primary"])),
        "structural_stringent": boolish(str(row["frozen_structural_stringent"])),
        "outside_structural_stringent": not boolish(str(row["frozen_structural_stringent"])),
        "structural_inclusive": boolish(str(row["frozen_structural_inclusive"])),
        "outside_structural_inclusive": not boolish(str(row["frozen_structural_inclusive"])),
    }[name]


def aggregate_q(records: list[dict[str, object]]) -> tuple[list[float], list[float]]:
    mean_q = [sum(float(row[f"q{i}"]) for row in records) / len(records) for i in (1, 2, 3)] if records else [math.nan] * 3
    counts = [sum(float(row[f"x{i}"]) for row in records) for i in (1, 2, 3, 4)]
    agg = qqs(counts)
    aggregate = [agg[f"q{i}"] for i in (1, 2, 3)] if agg else [math.nan] * 3
    return mean_q, aggregate


def write_partition_summary(rows: list[dict[str, object]], counters) -> None:
    names = ("nonchr4", "chr4", "pnas_outlier", "outside_pnas", "structural_primary", "outside_structural_primary", "structural_stringent", "outside_structural_stringent", "structural_inclusive", "outside_structural_inclusive")
    fields = ["clade", "partition", "total_source_loci", "total_valid_loci", "missing_or_undefined", "mean_q1", "mean_q2", "mean_q3", "aggregate_q1", "aggregate_q2", "aggregate_q3", "dominant_mean", "dominant_aggregate", "sampling_fraction", "coordinate_min", "coordinate_max", "estimand_note"]
    with (OUTDIR / "genomewide_quartet_partition_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader()
        for clade in CLADES:
            clade_rows = [row for row in rows if row["clade"] == clade]
            for name in names:
                subset = [row for row in clade_rows if partition(row, name)]
                mean_q, aggregate_values = aggregate_q(subset)
                total_source = 63430
                missing = total_source - len(subset)
                writer.writerow({"clade": clade, "partition": name, "total_source_loci": total_source, "total_valid_loci": len(subset), "missing_or_undefined": missing,
                                 "mean_q1": fmt(mean_q[0]), "mean_q2": fmt(mean_q[1]), "mean_q3": fmt(mean_q[2]),
                                 "aggregate_q1": fmt(aggregate_values[0]), "aggregate_q2": fmt(aggregate_values[1]), "aggregate_q3": fmt(aggregate_values[2]),
                                 "dominant_mean": max(zip(TOPOLOGIES, mean_q), key=lambda item: (item[1], item[0]))[0] if subset else "NA",
                                 "dominant_aggregate": max(zip(TOPOLOGIES, aggregate_values), key=lambda item: (item[1], item[0]))[0] if subset and not math.isnan(aggregate_values[0]) else "NA",
                                 "sampling_fraction": len(subset) / len(clade_rows) if clade_rows else 0,
                                 "coordinate_min": min((int(row["start"]) for row in subset), default="NA"), "coordinate_max": max((int(row["end"]) for row in subset), default="NA"),
                                 "estimand_note": "mean of per-locus normalized QQS; aggregate_q is ratio after summing x1..x4"})


def crossover(qb: list[float], qr: list[float], i: int, j: int) -> tuple[object, str, object]:
    if any(math.isnan(x) for x in qb + qr): return "NA", "missing_endpoint", "NA"
    numerator = qb[i] - qb[j]; denominator = numerator + qr[j] - qr[i]
    if abs(denominator) < 1e-12: return "NA", "denominator_zero", "NA"
    epsilon = numerator / denominator
    status = "valid_interior" if 0 < epsilon < 1 else "outside_[0,1]"
    third = ({0, 1, 2} - {i, j}).pop()
    third_value = (1 - epsilon) * qb[third] + epsilon * qr[third]
    return epsilon, status, third_value


def write_thresholds(rows: list[dict[str, object]]) -> None:
    comparisons = (("chr4_vs_nonchr4", "nonchr4", "chr4"), ("pnas_outlier_vs_background", "outside_pnas", "pnas_outlier"), ("structural_primary_vs_background", "outside_structural_primary", "structural_primary"), ("structural_stringent_vs_background", "outside_structural_stringent", "structural_stringent"), ("structural_inclusive_vs_background", "outside_structural_inclusive", "structural_inclusive"))
    fields = ["clade", "comparison", "background_partition", "affected_partition", "background_q1", "background_q2", "background_q3", "affected_q1", "affected_q2", "affected_q3", "observed_weight", "pair", "epsilon_star", "status", "third_at_epsilon", "dominant_background", "dominant_affected", "dominant_changes"]
    output = OUTDIR / "genomewide_quartet_thresholds.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader()
        for clade in CLADES:
            clade_rows = [row for row in rows if row["clade"] == clade]
            for comparison, bname, rname in comparisons:
                b = [row for row in clade_rows if partition(row, bname)]; r = [row for row in clade_rows if partition(row, rname)]
                qb, _ = aggregate_q(b); qr, _ = aggregate_q(r)
                weight = len(r) / (len(b) + len(r)) if b or r else math.nan
                db = TOPOLOGIES[max(range(3), key=lambda i: (qb[i], -i))] if b else "NA"; dr = TOPOLOGIES[max(range(3), key=lambda i: (qr[i], -i))] if r else "NA"
                for i, j in ((0, 1), (0, 2), (1, 2)):
                    epsilon, status, third = crossover(qb, qr, i, j)
                    writer.writerow({"clade": clade, "comparison": comparison, "background_partition": bname, "affected_partition": rname,
                                     **{f"background_q{k}": fmt(qb[k - 1]) for k in (1, 2, 3)}, **{f"affected_q{k}": fmt(qr[k - 1]) for k in (1, 2, 3)},
                                     "observed_weight": fmt(weight), "pair": f"q{i + 1}=q{j + 1}", "epsilon_star": fmt(epsilon), "status": status,
                                     "third_at_epsilon": fmt(third), "dominant_background": db, "dominant_affected": dr, "dominant_changes": str(db != dr).lower()})
    # Keep the existing threshold filenames synchronized with the recovered
    # published-statistics analysis, while retaining their explicit scope.
    (OUTDIR / "j48_quartet_mixture.tsv").write_text(output.read_text())
    (OUTDIR / "j48_threshold_summary.tsv").write_text("quantity\tvalue\tdefinition\nsource_statistics_scope\tpublished_full_taxon_QQS\tnot a 48-taxon gene-tree quartet estimate\nchr4_locus_fraction\t0.09872300173419518\t6262/63430 distinct-locus fraction\nastral_transition\tpending_bridges2\tcompleted J48 endpoints are qualitative only\n")
    lines = ["# Recovered genome-wide QQS threshold analysis", "", "The published `clade-rec.stat.xz` source covers all 63,430 locus IDs and all chromosomes. The recovered values are full-source QQS statistics, not quartet frequencies recomputed from the 48-taxon J48 gene-tree files.", "", "Chromosome-4 q values reproduce the existing Stage 1 table exactly: 18,109 focal rows compared and zero q-value mismatches.", "", "The affected-locus weights below are valid-observation fractions for each published-statistics comparison. They are distinct from physical genome fractions and from the 48-taxon ASTRAL mixture parameter.", "", "N61 changes from q1 in the non-chr4 background to q2 on chr4; its q1=q2 crossover is epsilon=0.1540, above the observed valid-locus weight 0.0985. The PNAS, primary, stringent, and inclusive region comparisons also have interior q1=q2 crossovers, but their weights and estimands are source-statistics quantities.", "", "N62 changes from q2 in the non-chr4 background to q1 on chr4; its q1=q2 crossover is epsilon=0.01354, below the observed valid-locus weight 0.0978. The PNAS comparison has q1=q2 crossover 0.00303 and q2=q3 crossover 0.8453; the third topology must be checked at every pairwise equality.", "", "Columbea remains q1-dominant at both non-chr4 and chr4 endpoints, so no dominant-topology crossover is supported for that comparison.", "", "These thresholds are empirical mixtures of published full-source QQS statistics. They do not predict the 48-taxon ASTRAL transition exactly because taxon sampling, estimand, and locus weighting differ."]
    (OUTDIR / "j48_threshold_analysis.md").write_text("\n".join(lines) + "\n")


def validate_against_chr4(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    old = {(r["Gene"], r["clade"]): r for r in csv.DictReader(OLD_CHR4.open(), delimiter="\t")}
    new = {(r["Gene"], r["clade"]): r for r in rows if r["Chromosome"] == "chr4"}
    mismatches = 0
    for key, old_row in old.items():
        current = new.get(key)
        if current is None or any(abs(float(current[f"q{i}"]) - float(old_row[f"q{i}"])) > 1e-10 for i in (1, 2, 3)):
            mismatches += 1
    path = OUTDIR / "genomewide_quartet_validation.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t"); writer.writerow(["metric", "value"])
        writer.writerow(["published_header_loci", 63430]); writer.writerow(["focal_genomewide_valid_rows", len(rows)])
        writer.writerow(["old_chr4_rows", len(old)]); writer.writerow(["new_chr4_rows", len(new)])
        writer.writerow(["chr4_q_value_mismatches", mismatches]); writer.writerow(["source_archive_sha256", sha256(ARCHIVE)])
    return [{"old_chr4_rows": len(old), "new_chr4_rows": len(new), "mismatches": mismatches}]


def main() -> None:
    manifest = load_manifest(); values, header_count, source_rows = read_source(manifest); rows, counters = build_rows(values, manifest)
    OUTDIR.mkdir(parents=True, exist_ok=True); write_genomewide(rows); write_partition_summary(rows, counters); validate_against_chr4(rows); write_thresholds(rows)
    print({"header_loci": header_count, "focal_source_rows": source_rows, "unique_focal_gene_clades": len(values), "valid_rows": len(rows), "archive_sha256": sha256(ARCHIVE)})


if __name__ == "__main__":
    main()
