#!/usr/bin/env python3
"""Make the current partial MSRC-threshold figure without new inference."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "empirical" / "neoaves_chr4"
RESULTS = BASE / "results" / "stage4e"
FIGURES = BASE / "figures"


def rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    thresholds = rows(RESULTS / "genomewide_quartet_thresholds.tsv")
    endpoints = {}
    crossovers = {}
    for row in thresholds:
        if row["comparison"] == "chr4_vs_nonchr4":
            endpoints[row["clade"]] = ([float(row[f"background_q{i}"]) for i in (1, 2, 3)], [float(row[f"affected_q{i}"]) for i in (1, 2, 3)])
            if row["pair"] == "q1=q2":
                crossovers[row["clade"]] = float(row["epsilon_star"]) if row["status"] == "valid_interior" else None
    def esc(value):
        return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="720" viewBox="0 0 1100 720">', '<rect width="1100" height="720" fill="white"/><style>text{font-family:sans-serif;fill:#222}.title{font-size:20px;font-weight:bold}.body{font-size:14px}.small{font-size:11px}</style>']
    parts += ['<rect x="30" y="30" width="500" height="300" fill="#fafafa" stroke="#555"/><text x="50" y="65" class="title">A  Genomic architecture</text>', '<line x1="70" y1="180" x2="490" y2="180" stroke="#222" stroke-width="8"/>']
    for start, end in ((25.03, 32.67), (33.51, 34.47), (44.13, 56.81)):
        x1, x2 = 70 + start / 60 * 420, 70 + end / 60 * 420
        parts.append(f'<line x1="{x1:.1f}" y1="180" x2="{x2:.1f}" y2="180" stroke="#d95f02" stroke-width="20"/>')
    parts += ['<text x="75" y="215" class="small">0</text><text x="455" y="215" class="small">60 Mb</text>', '<text x="55" y="275" class="small">Orange: published GalGal4 PNAS intervals</text>', '<text x="55" y="298" class="small">Structural breakpoint/region coordinates remain a separate track.</text>']
    parts += ['<rect x="570" y="30" width="500" height="300" fill="#fafafa" stroke="#555"/><text x="590" y="65" class="title">B  Published QQS mixture</text>']
    colors = {"q1": "#222222", "q2": "#1f78b4", "q3": "#e31a1c"}
    for row_index, clade in enumerate(("Columbea", "N61", "N62")):
        y0 = 105 + row_index * 65
        b, affected = endpoints[clade]
        for qi, qname in enumerate(("q1", "q2", "q3")):
            y1 = y0 + 42 - b[qi] * 38; y2 = y0 + 42 - affected[qi] * 38
            parts.append(f'<line x1="620" y1="{y1:.1f}" x2="1010" y2="{y2:.1f}" stroke="{colors[qname]}" stroke-width="2"/>')
        observed = {"Columbea": 0.0984, "N61": 0.0985, "N62": 0.0978}[clade]
        xo = 620 + observed * 390
        parts.append(f'<line x1="{xo:.1f}" y1="{y0}" x2="{xo:.1f}" y2="{y0+45}" stroke="#888" stroke-dasharray="3,2"/><text x="590" y="{y0+12}" class="small">{clade}</text>')
        if crossovers[clade] is not None:
            xc = 620 + crossovers[clade] * 390
            parts.append(f'<line x1="{xc:.1f}" y1="{y0}" x2="{xc:.1f}" y2="{y0+45}" stroke="#6a3d9a"/><text x="{min(xc+3,1000):.1f}" y="{y0+12}" class="small">ε*={crossovers[clade]:.3f}</text>')
    parts += ['<text x="590" y="300" class="small">Black/blue/red: q1/q2/q3; dashed: observed valid-locus weight</text>', '<text x="590" y="316" class="small">Full-source QQS; not a 48-taxon gene-tree estimate.</text>']
    parts += ['<rect x="30" y="370" width="500" height="300" fill="#fafafa" stroke="#555"/><text x="50" y="405" class="title">C  Completed 48-taxon ASTRAL</text>', '<text x="55" y="465" class="body">Baseline: T0, structural primary</text><text x="55" y="500" class="body">Alternative: PNAS, 500-kb blocks, non-chr4</text><text x="55" y="535" class="body">Between-group unrooted RF = 6</text><text x="55" y="590" class="small">These are discrete deletion treatments, not a continuous epsilon axis.</text>']
    parts += ['<rect x="570" y="370" width="500" height="300" fill="#fafafa" stroke="#555"/><text x="590" y="405" class="title">D  Future empirical ASTRAL transition</text><text x="595" y="470" class="body">Pending Bridges-2 mixture grid</text><text x="595" y="515" class="small">No intermediate species-tree results or significance values are plotted.</text><text x="595" y="545" class="small">See j48_threshold_experiment_manifest.tsv.</text>', '</svg>']
    svg = "".join(parts)
    svg_path = FIGURES / "neoaves_msrc_threshold.svg"
    svg_path.write_text(svg)


if __name__ == "__main__":
    main()
