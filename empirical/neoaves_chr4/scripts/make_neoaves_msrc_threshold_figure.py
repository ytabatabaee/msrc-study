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
    svg = '''<svg xmlns="http://www.w3.org/2000/svg" width="1100" height="720" viewBox="0 0 1100 720">
<rect width="1100" height="720" fill="white"/><style>text{font-family:sans-serif;fill:#222} .title{font-size:20px;font-weight:bold}.body{font-size:16px}.small{font-size:13px}</style>
<rect x="30" y="30" width="500" height="300" fill="#fafafa" stroke="#555"/><text x="50" y="65" class="title">A  Genomic architecture</text>
<line x1="70" y1="180" x2="490" y2="180" stroke="#222" stroke-width="8"/><line x1="245" y1="180" x2="298" y2="180" stroke="#d95f02" stroke-width="20"/><line x1="305" y1="180" x2="336" y2="180" stroke="#d95f02" stroke-width="20"/><line x1="379" y1="180" x2="468" y2="180" stroke="#d95f02" stroke-width="20"/>
<text x="75" y="215" class="small">0</text><text x="455" y="215" class="small">60 Mb</text><text x="55" y="275" class="small">Orange: published GalGal4 PNAS intervals</text><text x="55" y="298" class="small">Structural tracks withheld pending coordinate-specific source.</text>
<rect x="570" y="30" width="500" height="300" fill="#fafafa" stroke="#555"/><text x="590" y="65" class="title">B  Quartet mixture</text><text x="590" y="125" class="body">Background QQS endpoint: unavailable</text><text x="590" y="160" class="body">Affected chr4 endpoint: available</text><text x="590" y="215" class="small">Linear mixture curves withheld because Q_B</text><text x="590" y="238" class="small">cannot be estimated from the cached summary.</text><text x="590" y="285" class="small">No crossover or uncertainty band is claimed.</text>
<rect x="30" y="370" width="500" height="300" fill="#fafafa" stroke="#555"/><text x="50" y="405" class="title">C  ASTRAL topology transition</text><text x="55" y="470" class="body">Pending Bridges-2 coarse mixture grid</text><text x="55" y="515" class="small">No intermediate species-tree results were fabricated.</text><text x="55" y="545" class="small">Completed endpoint RF comparisons are reported separately.</text>
<rect x="570" y="370" width="500" height="300" fill="#fafafa" stroke="#555"/><text x="590" y="405" class="title">D  Specificity controls</text><text x="595" y="470" class="body">Pending deterministic matched controls</text><text x="595" y="515" class="small">No significance values or uncertainty bands are claimed.</text><text x="595" y="545" class="small">See the Bridges-2 experiment manifest.</text>
</svg>'''
    svg_path = FIGURES / "neoaves_msrc_threshold.svg"
    svg_path.write_text(svg)


if __name__ == "__main__":
    main()
