#!/usr/bin/env python3
"""Audit the Stiller alignment archive and its exact Stage-4D locus mapping."""
from __future__ import annotations
import argparse, csv, hashlib, json, os, re, sys, tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
EXT = ROOT / "empirical/neoaves_chr4/external/stiller2024"
MANIFEST = ROOT / "empirical/neoaves_chr4/results/stage4d_locus_manifest.tsv"
OUT = Path(__file__).resolve().parents[1] / "results"
ARCHIVE = EXT / "63430.alns.tar.gz"
URL = "https://sid.erda.dk/share_redirect/ENhZODU9YE/01_alignments_and_gene_trees/intergenic_regions/63430.alns.tar.gz"

def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""): h.update(block)
    return h.hexdigest()

def stem(name):
    base = Path(name).name
    if base.startswith("._"): return ""
    base = re.sub(r"\.(?:fasta|fas|fa|aln)$", "", base, flags=re.I)
    base = re.sub(r"_fasta\.treefile$", "", base)
    base = base.replace(".1k.start", "_1k_start")
    return base

def main(argv=None):
    global OUT
    ap = argparse.ArgumentParser(description="Audit exact Stiller alignment-to-manifest mapping")
    ap.add_argument("--output-dir", type=Path, default=OUT)
    OUT = ap.parse_args(argv).output_dir
    OUT.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(MANIFEST.open(), delimiter="\t"))
    chr4 = [r for r in rows if r["chromosome"] == "chr4"]
    manifest_stems = {stem(r["locus_id"]): r for r in chr4}
    members = []
    matched = {}
    if ARCHIVE.exists():
        with tarfile.open(ARCHIVE, "r:gz") as tf:
            for m in tf:
                if m.isfile():
                    s = stem(m.name)
                    members.append(m.name)
                    if s in manifest_stems: matched[s] = m.name
    audit = {
        "archive": str(ARCHIVE), "source_url": URL,
        "archive_exists": ARCHIVE.exists(),
        "archive_bytes": ARCHIVE.stat().st_size if ARCHIVE.exists() else 0,
        "archive_sha256": sha256(ARCHIVE) if ARCHIVE.exists() else "",
        "manifest_rows": len(rows), "chr4_manifest_rows": len(chr4),
        "archive_file_members": len(members), "exact_chr4_matches": len(matched),
        "missing_chr4_alignments": len(chr4) - len(matched),
        "example_members": members[:5],
    }
    (OUT / "input_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    with (OUT / "locus_alignment_mapping.tsv").open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["locus_id", "chromosome", "start", "end", "midpoint", "alignment_member", "mapping_status"])
        for r in chr4:
            s = stem(r["locus_id"]); w.writerow([r["locus_id"], r["chromosome"], r["start"], r["end"], r["midpoint"], matched.get(s, ""), "exact" if s in matched else "missing"])
    print(json.dumps(audit, indent=2))

if __name__ == "__main__": main()
