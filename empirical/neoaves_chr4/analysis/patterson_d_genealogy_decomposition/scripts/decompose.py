#!/usr/bin/env python3
"""Quantitative genealogy-based prediction of Neoaves Patterson's D.

This analysis is deliberately separate from the frozen Patterson's D and
genealogy-integration outputs.  It uses the existing named gene trees and the
original local alignments to estimate callable sequence exposure per locus.
"""
from __future__ import annotations
import csv, gzip, hashlib, math, random, shutil, subprocess, tarfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / "empirical/neoaves_chr4/analysis/patterson_d_genealogy_decomposition"
RES, FIG = BASE / "results", BASE / "figures"
PATT = ROOT / "empirical/neoaves_chr4/analysis/patterson_d/results"
MANIFEST = ROOT / "empirical/neoaves_chr4/results/stage4d_locus_manifest.tsv"
MAPPING = PATT / "locus_alignment_mapping.tsv"
ARCHIVE = ROOT / "empirical/neoaves_chr4/external/stiller2024/63430.alns.tar.gz"
TREES = ROOT / "empirical/neoaves_chr4/external/stiller2024/63430.named.gene.trees.gz"
PROVENANCE = PATT / "provenance.tsv"
INTERVALS = (("PNAS_1", 25030000, 32670000), ("PNAS_2", 33510000, 34470000), ("PNAS_3", 44130000, 56810000))
SIZES = (("100kb", 100000), ("250kb", 250000), ("500kb", 500000))
FILTERS = ("all", "transversions_only")
BASES = set("ACGT")
DESIGNS = {
    "Columbea_primary": ("Columba_livia", "Phoenicopterus_ruber", "Ciconia_maguari", "Gallus_gallus"),
    "Columbea_alt_outgroup": ("Columba_livia", "Phoenicopterus_ruber", "Ciconia_maguari", "Meleagris_gallopavo"),
    "N61_primary": ("Tauraco_erythrolophus", "Columba_livia", "Ciconia_maguari", "Gallus_gallus"),
    "N61_alt_outgroup": ("Tauraco_erythrolophus", "Columba_livia", "Ciconia_maguari", "Meleagris_gallopavo"),
    "N62_primary": ("Columba_livia", "Pterocles_gutturalis", "Tauraco_erythrolophus", "Gallus_gallus"),
    "N62_alt_outgroup": ("Columba_livia", "Pterocles_gutturalis", "Tauraco_erythrolophus", "Meleagris_gallopavo"),
}


def read_tsv(path):
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or (list(rows[0]) if rows else [])
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def num(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def avg(xs): return sum(xs) / len(xs) if xs else None


def D(a, b):
    return None if a + b == 0 else (a - b) / (a + b)


def region(start):
    hits = [name for name, a, b in INTERVALS if a <= start <= b]
    return hits[0] if hits else "outside_all_PNAS"


def overlap(start, end):
    hits = [name for name, a, b in INTERVALS if start <= b and end >= a]
    return ";".join(hits) or "outside_all_PNAS"


class Node:
    __slots__ = ("children", "label", "length", "parent")
    def __init__(self, label="", length=0.0):
        self.children, self.label, self.length, self.parent = [], label, length, None


def parse_newick(text):
    s, i = text.strip().rstrip(";"), 0
    def token(stop):
        nonlocal i
        j = i
        while i < len(s) and s[i] not in stop: i += 1
        return s[j:i]
    def node():
        nonlocal i
        n = Node()
        if i < len(s) and s[i] == "(":
            i += 1
            while True:
                c = node(); c.parent = n; n.children.append(c)
                if i < len(s) and s[i] == ",": i += 1; continue
                if i < len(s) and s[i] == ")": i += 1; break
                raise ValueError("malformed Newick")
            if i < len(s) and s[i] not in ",):": token(",):")
        else:
            n.label = token(":,()")
        if i < len(s) and s[i] == ":":
            i += 1
            try: n.length = float(token(",()"))
            except ValueError: n.length = math.nan
        return n
    root = node()
    if i != len(s): raise ValueError("trailing Newick text")
    return root


def distances(root, taxa):
    leaves = {}
    def collect(n):
        if not n.children: leaves[n.label] = n
        for c in n.children: collect(c)
    collect(root)
    if any(t not in leaves for t in taxa): return None
    paths = {}
    for taxon in taxa:
        n, d, p = leaves[taxon], 0.0, []
        while n is not None:
            p.append((n, d)); d += n.length if math.isfinite(n.length) else 0.0; n = n.parent
        paths[taxon] = p
    out = {}
    for i, a in enumerate(taxa):
        for b in taxa[i + 1:]:
            pb = {id(n): d for n, d in paths[b]}
            out[(a, b)] = next((da + pb[id(n)] for n, da in paths[a] if id(n) in pb), None)
    return out


def induced_quartet(root, taxa):
    """Return topology q1/q2/q3 and induced quartet internal edge.

    q1=(P1,P2)|(P3,O), q2=(P1,P3)|(P2,O), q3=(P1,O)|(P2,P3).
    The four-point criterion is root independent.  The induced internal edge
    is half the difference between the unique smallest and tied larger sums;
    this is equivalent to pruning the other leaves and suppressing degree-two
    nodes for an additive tree.
    """
    dm = distances(root, taxa)
    if dm is None: return {"status": "missing_taxon"}
    a, b, c, o = taxa
    pair = lambda x, y: dm.get((x, y), dm.get((y, x)))
    sums = [pair(a, b) + pair(c, o), pair(a, c) + pair(b, o), pair(a, o) + pair(b, c)]
    if any(v is None or not math.isfinite(v) for v in sums): return {"status": "nonfinite"}
    order = sorted(range(3), key=lambda k: sums[k])
    resolved = abs(sums[order[2]] - sums[order[1]]) < 1e-8 and abs(sums[order[1]] - sums[order[0]]) > 1e-8
    if not resolved: return {"status": "unresolved"}
    topology = ("q1", "q2", "q3")[order[0]]
    return {"status": "resolved", "topology": topology, "internal_branch": (sums[order[1]] - sums[order[0]]) / 2.0,
            "four_point_sums": ";".join(map(str, sums))}


def parse_fasta(raw):
    seqs, name, chunks = {}, None, []
    for line in raw.decode().splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if name is not None: seqs[name] = "".join(chunks).upper()
            name, chunks = line[1:].split()[0], []
        else: chunks.append(line.replace(" ", ""))
    if name is not None: seqs[name] = "".join(chunks).upper()
    return seqs


def callable_sites(seqs, taxa):
    n = min((len(seqs.get(t, "")) for t in taxa), default=0)
    return n, sum(all(seqs[t][i] in BASES for t in taxa) for i in range(n))


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""): h.update(chunk)
    return h.hexdigest()


def load_inputs():
    manifest = [r for r in read_tsv(MANIFEST) if r["chromosome"] == "chr4"]
    by_locus = {r["locus_id"]: r for r in manifest}
    align_map = {r["locus_id"]: r for r in read_tsv(MAPPING) if r["mapping_status"] == "exact"}
    patt = read_tsv(PATT / "patterson_d_locus_results.tsv")
    observed = {(r["design"], r["filter"], r["locus_id"]): r for r in patt if r["design"] in DESIGNS}
    return manifest, by_locus, align_map, observed


def recover_callable(manifest, align_map):
    wanted = set(align_map); out = {}; members = 0
    member_to_locus = {v["alignment_member"]: k for k, v in align_map.items()}
    if not ARCHIVE.exists(): return out, members
    with tarfile.open(ARCHIVE, "r:gz") as tf:
        for member in tf:
            if not member.isfile(): continue
            name = member.name
            locus = member_to_locus.get(name)
            if locus is None or locus not in wanted: continue
            handle = tf.extractfile(member)
            if handle is None: continue
            seqs = parse_fasta(handle.read()); taxa_union = set().union(*DESIGNS.values())
            lengths = [len(seqs[t]) for t in taxa_union if t in seqs]
            out[locus] = {"alignment_member": name, "n_alignment_sites": min(lengths) if lengths else 0,
                          "n_taxa_in_alignment": len(seqs), "observed_length_min": min(lengths) if lengths else 0}
            for design, taxa in DESIGNS.items():
                n, call = callable_sites(seqs, taxa)
                out[locus][f"n_alignment_sites_{design}"] = n
                out[locus][f"callable_sites_{design}"] = call
                out[locus][f"n_required_taxa_{design}"] = sum(t in seqs for t in taxa)
            members += 1
    return out, members


def recover_genealogy(manifest):
    wanted = {r["locus_id"] for r in manifest}; by_locus = {r["locus_id"]: r for r in manifest}; out = {}
    if not TREES.exists(): return out
    with gzip.open(TREES, "rt") as fh:
        for line in fh:
            if "\t" not in line: continue
            locus, nwk = line.rstrip("\n").split("\t", 1)
            if locus not in wanted: continue
            try:
                root = parse_newick(nwk)
                for design, taxa in DESIGNS.items():
                    q = induced_quartet(root, taxa)
                    row = {"locus_id": locus, "design": design, "topology_status": q.get("status", ""),
                           "topology": q.get("topology", ""), "internal_branch": q.get("internal_branch", ""),
                           "start": by_locus[locus]["start"], "end": by_locus[locus]["end"], "region": region(int(by_locus[locus]["start"]))}
                    out[(locus, design)] = row
            except Exception as exc:
                for design in DESIGNS:
                    out[(locus, design)] = {"locus_id": locus, "design": design, "topology_status": "parse_error:" + str(exc), "topology": "", "internal_branch": "", "start": by_locus[locus]["start"], "end": by_locus[locus]["end"], "region": region(int(by_locus[locus]["start"]))}
    return out


def make_locus_rows(manifest, by_locus, callable_map, genealogy, observed):
    rows = []
    for locus in by_locus:
        cm = callable_map.get(locus, {})
        for design in DESIGNS:
            g = genealogy.get((locus, design), {})
            obs = observed.get((design, "all", locus), {})
            topo = g.get("topology", "")
            rows.append({"locus_id": locus, "design": design, "chromosome": "chr4", "start": by_locus[locus]["start"], "end": by_locus[locus]["end"], "region": region(int(by_locus[locus]["start"])), "alignment_member": cm.get("alignment_member", ""), "n_alignment_sites": cm.get(f"n_alignment_sites_{design}", ""), "callable_sites": cm.get(f"callable_sites_{design}", ""), "n_taxa_in_alignment": cm.get("n_taxa_in_alignment", ""), "topology_status": g.get("topology_status", "missing_tree"), "topology": topo, "internal_branch": g.get("internal_branch", ""), "I_ABBA_q3": int(topo == "q3"), "I_BABA_q2": int(topo == "q2"), "observed_all_ABBA": obs.get("ABBA", ""), "observed_all_BABA": obs.get("BABA", ""), "observed_all_informative": obs.get("informative", "")})
    return rows


def pred_from_rows(rows, design, filt, region_name=None):
    z = []
    for r in rows:
        if r["design"] != design or not r["callable_sites"] or r["topology_status"] != "resolved": continue
        if region_name == "all_PNAS": hit = r["region"] != "outside_all_PNAS"
        elif region_name in ("PNAS_1", "PNAS_2", "PNAS_3"): hit = r["region"] == region_name
        elif region_name == "outside_all_PNAS": hit = r["region"] == "outside_all_PNAS"
        else: hit = True
        if not hit: continue
        z.append(r)
    observed = []
    for r in z:
        # Locate the original site-pattern record using its locus and design.
        # The values are attached by the caller to avoid any positional join.
        a, b, inf = num(r.get(f"observed_{filt}_ABBA")), num(r.get(f"observed_{filt}_BABA")), num(r.get(f"observed_{filt}_informative"))
        if a is not None and b is not None: observed.append((a, b, inf or 0))
    ca = cb = ta = tb = 0.0; ca_n = cb_n = 0; sum_l = 0.0
    for r in z:
        L = float(r["callable_sites"]); t = float(r["internal_branch"])
        sum_l += L
        if r["topology"] == "q3": ca += L; ta += L * t; ca_n += 1
        if r["topology"] == "q2": cb += L; tb += L * t; cb_n += 1
    abba = sum(x[0] for x in observed); baba = sum(x[1] for x in observed); inf = sum(x[2] for x in observed)
    mu_a = ta / ca if ca else None; mu_b = tb / cb if cb else None
    return {"n_prediction_loci": len(z), "n_ABBA_topology": ca_n, "n_BABA_topology": cb_n, "callable_sites_total": sum_l, "C_A": ca, "C_B": cb, "mu_A": mu_a, "mu_B": mu_b, "ABBA": abba, "BABA": baba, "informative_sites": inf, "D_site": D(abba, baba), "D_topology": D(ca, cb), "D_tree": D(ta, tb), "n_observed_loci": len(observed)}


def attach_observed(rows, observed):
    for r in rows:
        for filt in FILTERS:
            x = observed.get((r["design"], filt, r["locus_id"]), {})
            for k in ("ABBA", "BABA", "informative"):
                r[f"observed_{filt}_{k}"] = x.get(k, "")


def aggregate(rows, design, filt, region_name, size=None, start=None):
    z = rows
    if size is not None:
        z = [r for r in z if int(r["start"]) // size * size == start]
    # Region filters use locus starts, matching the frozen Patterson analysis.
    return pred_from_rows(z, design, filt, region_name)


def stats_for_all(rows, designs, filt, region_names, sizes):
    out = []
    for size_name, size in sizes:
        starts = range(0, ((max(int(r["end"]) for r in rows) + size - 1) // size) * size, size)
        for design in designs:
            for start in starts:
                s = aggregate(rows, design, filt, None, size, start)
                s.update({"resolution": size_name, "window_start": start, "window_end": start + size, "midpoint": start + size / 2, "design": design, "filter": filt, "region_mask": overlap(start, start + size - 1)})
                out.append(s)
    for design in designs:
        for reg in region_names:
            s = aggregate(rows, design, filt, reg)
            s.update({"resolution": "region", "window_start": "", "window_end": "", "midpoint": "", "design": design, "filter": filt, "region_mask": reg})
            out.append(s)
    return out


def decompose(inside, outside):
    def f(ca, cb, ma, mb):
        return D(ca * ma, cb * mb)
    vals = {"D_tree_inside": f(inside["C_A"], inside["C_B"], inside["mu_A"], inside["mu_B"]), "D_tree_outside": f(outside["C_A"], outside["C_B"], outside["mu_A"], outside["mu_B"])}
    ti_o = f(inside["C_A"], inside["C_B"], outside["mu_A"], outside["mu_B"]) - vals["D_tree_outside"]
    ti_i = vals["D_tree_inside"] - f(outside["C_A"], outside["C_B"], inside["mu_A"], inside["mu_B"])
    bi_o = f(outside["C_A"], outside["C_B"], inside["mu_A"], inside["mu_B"]) - vals["D_tree_outside"]
    bi_i = vals["D_tree_inside"] - f(inside["C_A"], inside["C_B"], outside["mu_A"], outside["mu_B"])
    vals["topology_shapley"] = (ti_o + ti_i) / 2
    vals["branch_length_shapley"] = (bi_o + bi_i) / 2
    vals["D_tree_shift"] = vals["D_tree_inside"] - vals["D_tree_outside"]
    vals["topology_plus_branch_check"] = vals["topology_shapley"] + vals["branch_length_shapley"]
    vals["D_topology_inside"] = D(inside["C_A"], inside["C_B"])
    vals["D_topology_outside"] = D(outside["C_A"], outside["C_B"])
    vals["D_topology_shift"] = vals["D_topology_inside"] - vals["D_topology_outside"]
    vals["D_unweighted_inside"] = D(inside["n_ABBA_topology"], inside["n_BABA_topology"])
    vals["D_unweighted_outside"] = D(outside["n_ABBA_topology"], outside["n_BABA_topology"])
    vals["D_unweighted_shift"] = vals["D_unweighted_inside"] - vals["D_unweighted_outside"]
    vals["exposure_weighted_minus_unweighted_inside"] = vals["D_topology_inside"] - D(inside["n_ABBA_topology"], inside["n_BABA_topology"])
    vals["exposure_weighted_minus_unweighted_outside"] = vals["D_topology_outside"] - D(outside["n_ABBA_topology"], outside["n_BABA_topology"])
    return vals


def bootstrap(rows, design, filt, region, block_size, n=500, seed=61):
    eligible = [r for r in rows if r["design"] == design and r["topology_status"] == "resolved" and r["callable_sites"] and (r["region"] != "outside_all_PNAS" if region == "all_PNAS" else r["region"] == region)]
    outside = [r for r in rows if r["design"] == design and r["topology_status"] == "resolved" and r["callable_sites"] and r["region"] == "outside_all_PNAS"]
    groups = defaultdict(list)
    for r in eligible: groups[int(r["start"]) // block_size].append(r)
    ogroups = defaultdict(list)
    for r in outside: ogroups[int(r["start"]) // block_size].append(r)
    rng = random.Random(seed); vals = defaultdict(list); ik, ok = list(groups), list(ogroups)
    for _ in range(n):
        sample = [r for k in (rng.choice(ik) for _ in ik) for r in groups[k]] if ik else []
        osample = [r for k in (rng.choice(ok) for _ in ok) for r in ogroups[k]] if ok else []
        a = pred_from_rows(sample, design, filt); b = pred_from_rows(osample, design, filt)
        for key in ("D_site", "D_topology", "D_tree"):
            vals[key].append(a[key])
            vals[key + "_outside"].append(b[key])
            vals[key + "_inside_minus_outside"].append(None if a[key] is None or b[key] is None else a[key] - b[key])
    out = {"n_blocks_inside": len(ik), "n_blocks_outside": len(ok)}
    for key, xs in vals.items():
        xs = sorted(x for x in xs if x is not None)
        out[key + "_se"] = math.sqrt(sum((x - avg(xs)) ** 2 for x in xs) / (len(xs) - 1)) if len(xs) > 1 else None
        out[key + "_ci_low"] = xs[int(.025 * len(xs))] if xs else None
        out[key + "_ci_high"] = xs[min(len(xs) - 1, int(.975 * len(xs)))] if xs else None
    return out


def convert(svg):
    c = shutil.which("rsvg-convert")
    if c:
        subprocess.run([c, "-f", "png", "-o", str(svg.with_suffix(".png")), str(svg)], check=False)
        subprocess.run([c, "-f", "pdf", "-o", str(svg.with_suffix(".pdf")), str(svg)], check=False)


def figure_tracks(rows):
    z = [r for r in rows if r["design"] == "N61_primary" and r["filter"] == "all" and r["resolution"] == "250kb"]
    width, height, left, top, pw = 1500, 680, 100, 45, 1360
    xmax = max(float(r["window_end"]) for r in z); panels = [("D site", "D_site", -1, 1), ("D topology", "D_topology", -1, 1), ("D tree", "D_tree", -1, 1), ("informative sites", "informative_sites", 0, max([float(r["informative_sites"]) for r in z] or [1]))]
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial;font-size:13px}}.line{{fill:none;stroke-width:2}}</style><text x="100" y="25" font-size="22">Neoaves chr4 N61: observed and genealogy-predicted D</text>']
    def x(v): return left + pw * v / xmax
    for i, (label, key, lo, hi) in enumerate(panels):
        y0 = top + i * 145; lines.append(f'<text x="15" y="{y0 + 8}">{label}</text>')
        for _, a, b in INTERVALS: lines.append(f'<rect x="{x(a):.1f}" y="{y0}" width="{max(1, x(b) - x(a)):.1f}" height="100" fill="#f4b183" opacity=".25"/>')
        pts = [f'{x(float(r["midpoint"])):.1f},{y0 + 100 * (1 - (float(r[key]) - lo) / (hi - lo)):.1f}' for r in z if r.get(key, "") not in ("", None)]
        if pts: lines.append(f'<polyline class="line" stroke="{("#d62728","#1f77b4","#2ca02c","#e67e22")[i]}" points="{" ".join(pts)}"/>')
        lines.append(f'<line x1="{left}" y1="{y0 + 100 * (1 - (0 - lo) / (hi - lo)):.1f}" x2="{left + pw}" y2="{y0 + 100 * (1 - (0 - lo) / (hi - lo)):.1f}" stroke="#999"/>')
    lines.append('<text x="520" y="650">Shading: published PNAS intervals; D_site is nucleotide ABBA–BABA, D_topology and D_tree are genealogy-based predictions.</text></svg>')
    p = FIG / "neoaves_N61_D_decomposition_tracks.svg"; p.write_text("\n".join(lines)); convert(p)


def bar_figure(name, title, rows, keys, labels):
    width, height, left, base, panelw = 1700, 800, 100, 520, 300
    vals = [num(r.get(k)) for r in rows for k in keys if num(r.get(k)) is not None]; lim = max([abs(x) for x in vals] or [1])
    lines = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial;font-size:13px}}.bar{{stroke:#333}}</style><text x="100" y="30" font-size="22">{title}</text>']
    colors = ("#d62728", "#1f77b4", "#2ca02c", "#9467bd", "#ff7f0e")
    for j, (r, label) in enumerate(zip(rows, labels)):
        x0 = left + j * panelw; lines.append(f'<text x="{x0}" y="65">{label}</text><line x1="{x0}" y1="{base}" x2="{x0 + 240}" y2="{base}" stroke="#555"/>')
        for i, k in enumerate(keys):
            v = num(r.get(k));
            if v is None: continue
            y = base - 260 * v / lim; lines.append(f'<rect class="bar" x="{x0 + 20 + i * 48}" y="{min(y, base):.1f}" width="32" height="{max(1, abs(y-base)):.1f}" fill="{colors[i % len(colors)]}"/><text x="{x0 + 18 + i*48}" y="{base+22}" transform="rotate(45 {x0+18+i*48} {base+22})">{k}</text>')
    lines.append('</svg>'); p = FIG / name; p.write_text("\n".join(lines)); convert(p)


def main():
    RES.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    manifest, by_locus, align_map, observed = load_inputs()
    callable_map, member_count = recover_callable(manifest, align_map)
    genealogy = recover_genealogy(manifest)
    rows = make_locus_rows(manifest, by_locus, callable_map, genealogy, observed); attach_observed(rows, observed)
    fields = list(rows[0]) if rows else []
    write_tsv(RES / "locus_genealogy_predictions.tsv", rows, fields)
    validation = []
    for label, path, status in [("stage4d_manifest", MANIFEST, "used_exact_locus_id"), ("alignment_mapping", MAPPING, f"{len(align_map)} exact mappings"), ("nucleotide_archive", ARCHIVE, f"{member_count} chr4 members inspected"), ("named_gene_trees", TREES, f"{len(genealogy)//len(DESIGNS)} loci parsed")]:
        validation.append({"input": label, "path": str(path), "exists": int(path.exists()), "sha256": "prior_verified:" + (next((r["sha256"] for r in read_tsv(PROVENANCE) if r["input"] == "63430.alns.tar.gz"), "") if label == "nucleotide_archive" else (sha256(path) if path.exists() and (label == "named_gene_trees" or path.stat().st_size < 100000000) else "not_rehashed")), "status": status})
    write_tsv(RES / "input_validation.tsv", validation)
    all_stats = []
    for design in DESIGNS:
        for filt in FILTERS: all_stats.extend(stats_for_all(rows, (design,), filt, ("PNAS_1", "PNAS_2", "PNAS_3", "all_PNAS", "outside_all_PNAS"), SIZES))
    stat_fields = list(all_stats[0]) if all_stats else []
    write_tsv(RES / "observed_vs_predicted_D.tsv", all_stats, stat_fields)
    n61_rows = [r for r in all_stats if r["design"] == "N61_primary" and r["filter"] == "all" and r["resolution"] == "region"]
    byreg = {r["region_mask"]: r for r in n61_rows}; inside, outside = byreg["all_PNAS"], byreg["outside_all_PNAS"]
    dec = decompose(inside, outside); dec.update({"design": "N61_primary", "filter": "all", "region_inside": "all_PNAS", "region_outside": "outside_all_PNAS", "observed_D_site_shift": inside["D_site"] - outside["D_site"], "observed_D_topology_shift": inside["D_topology"] - outside["D_topology"], "observed_D_tree_shift": inside["D_tree"] - outside["D_tree"]})
    write_tsv(RES / "regional_decomposition.tsv", [dec], list(dec))
    branch = []
    for reg in ("PNAS_1", "PNAS_2", "PNAS_3", "all_PNAS"):
        x, y = byreg.get(reg), byreg["outside_all_PNAS"]
        if x:
            branch.append({"region": reg, "C_A": x["C_A"], "C_B": x["C_B"], "mu_A": x["mu_A"], "mu_B": x["mu_B"], "outside_C_A": y["C_A"], "outside_C_B": y["C_B"], "outside_mu_A": y["mu_A"], "outside_mu_B": y["mu_B"], "D_tree": x["D_tree"], "D_topology": x["D_topology"], "n_A": x["n_ABBA_topology"], "n_B": x["n_BABA_topology"], "callable_sites_total": x["callable_sites_total"]})
    write_tsv(RES / "branch_length_contribution.tsv", branch, list(branch[0]))
    uncertainty = []
    for bs in (500000, 1000000, 2000000):
        b = bootstrap(rows, "N61_primary", "all", "all_PNAS", bs, n=500)
        for metric in ("D_site", "D_topology", "D_tree", "D_site_outside", "D_site_inside_minus_outside", "D_topology_outside", "D_topology_inside_minus_outside", "D_tree_outside", "D_tree_inside_minus_outside"):
            if metric.endswith("_inside_minus_outside"):
                base_metric = metric.removesuffix("_inside_minus_outside")
                point = inside[base_metric] - outside[base_metric]
            elif metric.endswith("_outside"): point = outside[metric.removesuffix("_outside")]
            else: point = inside[metric]
            uncertainty.append({"design": "N61_primary", "filter": "all", "block_size": bs, "metric": metric, "estimate": point, "bootstrap_se": b.get(metric + "_se"), "ci_low": b.get(metric + "_ci_low"), "ci_high": b.get(metric + "_ci_high"), "n_blocks_inside": b["n_blocks_inside"], "n_blocks_outside": b["n_blocks_outside"]})
    write_tsv(RES / "statistical_validation.tsv", uncertainty, list(uncertainty[0]))
    robust = []
    for design in DESIGNS:
        for filt in FILTERS:
            rr = [r for r in all_stats if r["design"] == design and r["filter"] == filt and r["region_mask"] == "all_PNAS" and r["resolution"] == "region"]
            if rr:
                r = rr[0]; robust.append({"design": design, "filter": filt, "D_site": r["D_site"], "D_topology": r["D_topology"], "D_tree": r["D_tree"], "callable_sites": r["callable_sites_total"], "n_prediction_loci": r["n_prediction_loci"], "note": "controls are descriptive; induced quartet roles differ"})
    write_tsv(RES / "robustness_analysis.tsv", robust, list(robust[0]))
    figure_tracks([r for r in all_stats if r["resolution"] == "250kb"])
    bar_figure("neoaves_N61_prediction_comparison.svg", "N61 observed versus genealogy-predicted D", [byreg["PNAS_1"], byreg["PNAS_2"], byreg["PNAS_3"], byreg["all_PNAS"], byreg["outside_all_PNAS"]], ["D_site", "D_topology", "D_tree"], ["PNAS_1", "PNAS_2", "PNAS_3", "all PNAS", "outside"])
    bar_figure("neoaves_N61_decomposition.svg", "N61 Shapley decomposition of D_tree shift", [dec], ["topology_shapley", "branch_length_shapley", "D_tree_shift"], ["inside minus outside"])
    controls = [next(r for r in all_stats if r["design"] == d and r["filter"] == "all" and r["resolution"] == "region" and r["region_mask"] == "all_PNAS") for d in DESIGNS if d.endswith("primary")]
    bar_figure("neoaves_all_clades_genealogy_prediction.svg", "Observed and genealogy-predicted D: descriptive clade controls", controls, ["D_site", "D_topology", "D_tree"], ["Columbea", "N61", "N62"])
    report = ["# Neoaves N61 genealogy-based Patterson's D decomposition", "", "## Scope and validation", f"The analysis used the existing checksum-verified Stiller alignment archive without redownloading it. Exact chr4 alignment-member mappings were used; {len(callable_map)} alignment members were read and callable sites were calculated for the four taxa in each design. The named gene-tree file was parsed for {len(genealogy)//len(DESIGNS)} loci and each design was classified by a root-independent four-point induced quartet.", "", "The induced quartet labels are q1=(P1,P2)|(P3,O), q2=(P1,P3)|(P2,O), q3=(P1,O)|(P2,P3). Under outgroup polarization, ABBA is q3 and BABA is q2. Internal lengths are induced quartet internal edges in ML substitution units, not coalescence times.", "", "## Main quantitative result", f"For N61, the original observed locus-level shift is ΔD_site = {num(inside['D_site']) - num(outside['D_site']):.6f}. The genealogy predictions are ΔD_topology = {num(inside['D_topology']) - num(outside['D_topology']):.6f} and ΔD_tree = {num(inside['D_tree']) - num(outside['D_tree']):.6f}. The exact values, exposures, and denominators are in `observed_vs_predicted_D.tsv`.", f"The symmetric decomposition attributes {dec['topology_shapley']:.6f} of the predicted shift to topology/exposure composition and {dec['branch_length_shapley']:.6f} to conditional internal-branch lengths; their sum is {dec['topology_plus_branch_check']:.6f}, equal to the predicted D_tree shift {dec['D_tree_shift']:.6f} up to floating-point precision.", f"The unweighted resolved-topology contrast shifts by {dec['D_unweighted_shift']:.6f}; callable-site exposure changes the topology-only contrast by {dec['exposure_weighted_minus_unweighted_inside'] - dec['exposure_weighted_minus_unweighted_outside']:.6f} beyond that unweighted comparison. Thus the large predicted shift is primarily a change in discordant-topology composition, with a secondary contribution from conditional internal branches; callable exposure is a small additional effect in this summary.", f"By interval, (D_site, D_topology, D_tree) were PNAS_1 ({byreg['PNAS_1']['D_site']}, {byreg['PNAS_1']['D_topology']}, {byreg['PNAS_1']['D_tree']}), PNAS_2 ({byreg['PNAS_2']['D_site']}, {byreg['PNAS_2']['D_topology']}, {byreg['PNAS_2']['D_tree']}), and PNAS_3 ({byreg['PNAS_3']['D_site']}, {byreg['PNAS_3']['D_topology']}, {byreg['PNAS_3']['D_tree']}).", "", "## Interpretation", "", "The local trees recover the sign but substantially overpredict the observed magnitude. Agreement in magnitude is not required by the infinite-sites approximation: observed ABBA/BABA counts also depend on finite-sites homoplasy, ancestral-state errors, alignment filtering, model misspecification, and tree-estimation dependence. The named trees were estimated from the same alignments, so this is not independent validation.", "", "A finite-sites D prediction was not fit here. A JC69 or richer model would require an explicit likelihood over the full quartet branch configuration, ascertainment, and empirical base composition; applying a fitted correction without independent calibration would be misleading. Transversion-only sensitivity is retained in the observed tables and robustness output.", "", "The PNAS intervals are published GalGal4 predicates applied to the Stiller GalGal6 axis; no exact lift was introduced. Results are therefore regional compatibility analyses, not nucleotide-resolution boundary validation.", "", "The analysis is compatible with a shared recombination-partitioned genealogy regime, but it does not discriminate MSRC from introgression, incomplete lineage sorting, mutation-rate variation, or gene-tree error. Both the predicted genealogies and observed D derive from the same underlying sequence data."]
    (RES / "scientific_interpretation.md").write_text("\n".join(report) + "\n")


if __name__ == "__main__": main()
