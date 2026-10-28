"""Bounded exact-Moran geometry audit for the conditional quartet image."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .full_species_tree import full_species_tree_quartet_probabilities
from .unbalanced_species_tree import unbalanced_species_tree_conditional_quartet_probabilities

PATTERNS = ("0011", "0101", "0110")
COLORS = {"0011": "#0072B2", "0101": "#D55E00", "0110": "#009E73"}


def xy(q):
    return 0.5 * (2 * q[0] + q[2]), np.sqrt(3.0) * q[0] / 2.0


def _balanced(z, root, t1, rho=.1, sister=(.2, .3)):
    out = full_species_tree_quartet_probabilities(
        M_root=4, M_left=4, M_right=4, M_terminal=(4,) * 4, rho=rho,
        root_length=root, sister_lengths=sister,
        terminal_lengths=(t1, .1, .2, .3), pi_old=np.eye(5)[2],
    )
    return out.conditional_quartet_probabilities[int(z, 2)]


def _pectinate(z, root, t1, rho=.1, intermediate=.2, sister=.3):
    out = unbalanced_species_tree_conditional_quartet_probabilities(
        M_root=4, M_intermediate=4, M_sister=4, M_terminal=(4,) * 4, rho=rho,
        root_length=root, intermediate_length=intermediate, sister_length=sister,
        terminal_lengths=(t1, .1, .2, .3), pi_old=np.eye(5)[2],
    )
    return out.conditional_quartet_probabilities[int(z, 2)]


def _svg(output: Path, rows):
    width, height = 920, 470
    panels = {"balanced": (35, 45), "pectinate": (485, 45)}
    scale = 330

    def point(q, origin):
        x, y = xy(q)
        return origin[0] + scale * x, origin[1] + scale * (0.8660254038 - y)

    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif} .small{font-size:12px}</style>']
    for shape, origin in panels.items():
        tri = [point((1, 0, 0), origin), point((0, 1, 0), origin), point((0, 0, 1), origin)]
        parts.append(f'<text x="{origin[0]+150}" y="25" text-anchor="middle" font-size="20">{shape}</text>')
        parts.append('<polygon points="{}" fill="none" stroke="#555" stroke-width="1.2"/>'.format(" ".join(f"{x:.2f},{y:.2f}" for x,y in tri)))
        a, b = point((1, 0, 0), origin), point((.25, .75, 0), origin)
        parts.append(f'<line x1="{a[0]:.2f}" y1="{a[1]:.2f}" x2="{b[0]:.2f}" y2="{b[1]:.2f}" stroke="#999" stroke-dasharray="5,4"/>')
        for row in rows:
            if row[0] != shape:
                continue
            x, y = point(row[5:8], origin)
            parts.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.2" fill="{COLORS[row[1]]}" fill-opacity="0.55"/>')
        h1, h2 = point((.8,.1,.1), origin), point((.1,.8,.1), origin)
        parts.append(f'<line x1="{h1[0]:.2f}" y1="{h1[1]:.2f}" x2="{h2[0]:.2f}" y2="{h2[1]:.2f}" stroke="#222" stroke-dasharray="3,3"/>')
    parts += ['<text x="35" y="445" class="small">Blue/orange/green: exact sampled images for 0011/0101/0110; black: illustrative hybridization line</text>', '</svg>']
    (output / "conditional_moran_global_geometry.svg").write_text("\n".join(parts))


def run(output: Path):
    output.mkdir(parents=True, exist_ok=True)
    settings = (("base", 1.0, .2, .1, .2, .3), ("root_long", 8.0, .2, .1, .2, .3),
                ("terminal_long", 1.0, 1.0, .1, .2, .3), ("switching", 1.0, .2, .5, .2, .3),
                ("ancestral_long", 1.0, .2, .1, .8, .8))
    rows = []
    for shape in ("balanced", "pectinate"):
        for pattern in PATTERNS:
            for label, root, t1, rho, first, second in settings:
                q = (_balanced(pattern, root, t1, rho, (first, second)) if shape == "balanced"
                     else _pectinate(pattern, root, t1, rho, first, second))
                rows.append((shape, pattern, label, root, t1, *q, *xy(q)))
    with (output / "global_geometry_samples.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("shape", "pattern", "family", "T_root", "T_terminal1", "q1", "q2", "q3", "x", "y"))
        writer.writerows(rows)
    with (output / "hybridization_line.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("gamma", "q1", "q2", "q3"))
        for gamma in np.linspace(0, 1, 21):
            writer.writerow((gamma, *((1-gamma)*np.array([.8,.1,.1]) + gamma*np.array([.1,.8,.1]))))
    _svg(output, rows)


if __name__ == "__main__":
    run(Path(__file__).parent / "results" / "conditional_moran_rigorous_coverage")
