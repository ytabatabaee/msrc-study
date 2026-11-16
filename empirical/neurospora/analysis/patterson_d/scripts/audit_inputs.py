#!/usr/bin/env python3
"""Audit deposited Neurospora genotype/window inputs; no D is calculated."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
OUT = ROOT / "empirical/neurospora/analysis/patterson_d/results"
GENO = ROOT / "data/neurospora/raw/neurospora/geno/cras_sit_hist_tet_L10.min20MAC1.geno.gz"
WINDOWS = ROOT / "data/neurospora/raw/neurospora/trees/cras_sit_hist_tet_L10.min20MAC1.phyml_GTR_bionj_w50sites.data.tsv.gz"
ARCHIVE = ROOT / "data/neurospora/raw/neurospora.tar.gz"

def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    windows: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with gzip.open(WINDOWS, "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            windows[row["scaffold"]].append((int(row["start"]), int(row["end"])))
    ix = Counter()
    assigned = Counter()
    call_hist = Counter()
    chrom_hist = Counter()
    nrows = 0
    extras = []
    with gzip.open(GENO, "rt") as f:
        header = f.readline().strip().split()
        samples = header[2:]
        for line in f:
            fields = line.split()
            chrom, pos = fields[0], int(fields[1])
            nrows += 1
            chrom_hist[chrom] += 1
            call_hist[sum(x != "N" for x in fields[2:])] += 1
            k = ix[chrom]
            ws = windows[chrom]
            while k < len(ws) and ws[k][1] < pos:
                k += 1
            ix[chrom] = k
            if k < len(ws) and ws[k][0] <= pos <= ws[k][1]:
                assigned[(chrom, k)] += 1
            else:
                extras.append((chrom, pos))
    window_hist = Counter(assigned.values())
    with (OUT / "genotype_source_inventory.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["file", "bytes", "md5", "format", "records", "taxa_or_note"])
        w.writerow([str(ARCHIVE.relative_to(ROOT)), ARCHIVE.stat().st_size, md5(ARCHIVE), "tar.gz", "archive", "Zenodo record 4946041"])
        w.writerow([str(GENO.relative_to(ROOT)), GENO.stat().st_size, md5(GENO), "gzip TSV; haploid nucleotide characters", nrows, ",".join(samples)])
        w.writerow([str(WINDOWS.relative_to(ROOT)), WINDOWS.stat().st_size, md5(WINDOWS), "gzip TSV; 1-based inclusive source intervals", sum(map(len, windows.values())), "50 variable sites/window"])
    summary = {
        "archive_md5_expected_and_observed": md5(ARCHIVE),
        "archive_bytes": ARCHIVE.stat().st_size,
        "genotype_rows": nrows,
        "sample_count": len(samples),
        "sample_ids": samples,
        "windows": sum(map(len, windows.values())),
        "windows_by_linkage_group": {k: len(v) for k, v in sorted(windows.items())},
        "genotype_rows_by_linkage_group": dict(sorted(chrom_hist.items())),
        "genotype_rows_per_window_histogram": dict(sorted(Counter(assigned.values()).items())),
        "windows_with_genotype_rows": len(assigned),
        "rows_outside_window_intervals": len(extras),
        "first_outside_rows": extras[:10],
        "called_sample_count_histogram": dict(sorted(call_hist.items())),
    }
    (OUT / "coordinate_audit.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
