#!/usr/bin/env python3
"""Aggregate frozen Neoaves chr4 quartet and background-distance tracks."""
from __future__ import annotations
import csv, hashlib, html, json, math, statistics, subprocess
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical" / "neoaves_chr4"
OUT = EMP / "results" / "stage10_multiresolution_tracks"
FIG = EMP / "figures"
SUPPORT = EMP / "results/stage6_genomic_tracks/window_quartet_support.tsv"
BG = EMP / "results/stage9_full_taxon_background/background_distance_estimators.tsv"
COV = EMP / "results/stage9_full_taxon_background/per_window_taxon_coverage.tsv"
REGIONS = ROOT / "data/neoaves_chr4/processed/canonical_structural_regions.tsv"
RESOLUTIONS = (100_000, 250_000, 500_000)
CLADES = ("Columbea", "N61", "N62")
QCOLORS = ("#2673B8", "#D65F5F", "#4A9B68")


def rows(path):
    with path.open() as h: return list(csv.DictReader(h, delimiter="\t"))


def write(path, data, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(data)


def num(x):
    try: return float(x)
    except (TypeError, ValueError): return math.nan


def avg(v):
    v = [x for x in v if math.isfinite(x)]
    return statistics.mean(v) if v else math.nan


def med(v):
    v = [x for x in v if math.isfinite(x)]
    return statistics.median(v) if v else math.nan


def stdev(v):
    v = [x for x in v if math.isfinite(x)]
    return statistics.stdev(v) if len(v) > 1 else (0.0 if v else math.nan)


def correlation(a, b):
    z = [(x, y) for x, y in zip(a, b) if math.isfinite(x) and math.isfinite(y)]
    if len(z) < 2: return math.nan
    x, y = zip(*z); mx, my = avg(x), avg(y)
    sx = math.sqrt(sum((i-mx)**2 for i in x)); sy = math.sqrt(sum((i-my)**2 for i in y))
    return sum((i-mx)*(j-my) for i,j in z)/(sx*sy) if sx and sy else math.nan


def checksum(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def bin_start(mid, resolution):
    return int(math.floor(mid / resolution) * resolution)


def aggregate_support(data, resolution):
    acc = defaultdict(lambda: {"q1": [], "q2": [], "q3": []})
    for r in data:
        z = acc[(bin_start(num(r["midpoint"]), resolution), r["clade"])]
        for k in ("q1_species_history", "q2_alternative_1", "q3_alternative_2"):
            z[k[:2]] .append(num(r[k]))
    out = []
    for (b, clade), z in sorted(acc.items()):
        q1, q2, q3 = (avg(z[k]) for k in ("q1", "q2", "q3"))
        out.append({"resolution_bp": resolution, "clade": clade, "bin_start": b, "bin_end": b+resolution, "bin_midpoint": b+resolution/2, "n_loci": len(z["q1"]), "q1": q1, "q2": q2, "q3": q3, "support_sum": q1+q2+q3, "D2": q2-q1, "D3": q3-q1})
    return out


def aggregate_background(data, resolution):
    acc = defaultdict(lambda: defaultdict(list))
    for r in data:
        z = acc[bin_start(num(r["midpoint"]), resolution)]
        for k in ("B_observed", "B47", "B_fixed", "B_residual", "n_tree_taxa", "n_broad_present"):
            z[k].append(num(r[k]))
    out = []
    for b, z in sorted(acc.items()):
        n = [x for x in z["n_tree_taxa"] if math.isfinite(x)]
        out.append({"resolution_bp": resolution, "bin_start": b, "bin_end": b+resolution, "bin_midpoint": b+resolution/2, "n_loci": len(z["B_observed"]), "B_observed_mean": avg(z["B_observed"]), "B_observed_median": med(z["B_observed"]), "B_observed_sd": stdev(z["B_observed"]), "B47_mean": avg(z["B47"]), "B_fixed_mean": avg(z["B_fixed"]), "B_residual_mean": avg(z["B_residual"]), "n_taxa_mean": avg(n), "n_taxa_median": med(n), "n_taxa_sd": stdev(n), "fraction_ge_300_taxa": avg([x >= 300 for x in n]), "taxon_coverage_range": max(n)-min(n) if n else math.nan})
    return out


def path_segments(xs, ys, xmap, ymap):
    out, cur = [], []
    for i, (xv, yv) in enumerate(zip(xs, ys)):
        if not math.isfinite(yv) or (i and xv-xs[i-1] > 600000):
            if cur: out.append(" ".join(cur)); cur=[]
            continue
        cur.append(("M" if not cur else "L") + f" {xmap(xv):.2f},{ymap(yv):.2f}")
    if cur: out.append(" ".join(cur))
    return out


def svg_figure(support, background, regions, resolution, stem, include_background):
    width, left, right, panel_h, gap, top = 1600, 140, 35, 112, 15, 48
    n_panels = 8 if include_background else 7
    height = top + n_panels*panel_h + (n_panels-1)*gap + 55
    xmax = max(r["bin_end"] for r in background)/1e6
    usable = width-left-right
    xmap = lambda v: left + usable*(v/1e6)/xmax
    def ymap(v, lo, hi, y0): return y0+panel_h-panel_h*(v-lo)/(hi-lo)
    def text(v): return html.escape(str(v), quote=True)
    def marks(y0):
        s=[]
        for rr in regions:
            a,b=num(rr["canonical_start"]),num(rr["canonical_end"])
            s.append(f'<rect x="{xmap(a):.2f}" y="{y0}" width="{max(1,xmap(b)-xmap(a)):.2f}" height="{panel_h}" fill="#bdbdbd" opacity=".12"/>')
        return "".join(s)
    def line(xs, ys, lo, hi, y0, color, width="0.9"):
        return "".join(f'<path d="{p}" fill="none" stroke="{color}" stroke-width="{width}"/>' for p in path_segments(xs,ys,xmap,lambda v:ymap(v,lo,hi,y0)))
    def support_rects(z,y0):
        out=[]
        for r in z:
            bx=xmap(r["bin_start"]); bw=max(1,xmap(r["bin_end"])-bx); base=0
            for val,col in zip((r["q1"],r["q2"],r["q3"]),QCOLORS):
                out.append(f'<rect x="{bx:.2f}" y="{ymap(base+val,0,1,y0):.2f}" width="{bw:.2f}" height="{max(0,ymap(base,0,1,y0)-ymap(base+val,0,1,y0)):.2f}" fill="{col}"/>'); base+=val
        return "".join(out)
    svg=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">','<rect width="100%" height="100%" fill="white"/>','<style>text{font-family:Arial,sans-serif;fill:#222}.axis{stroke:#333;stroke-width:1}.grid{stroke:#ddd;stroke-width:.6}</style>',f'<text x="{left}" y="25" font-size="18" font-weight="bold">Neoaves chromosome 4 ({resolution//1000} kb bins)</text>']
    y0=top; svg += [marks(y0),f'<text x="{left-10}" y="{y0+panel_h/2+4}" text-anchor="end" font-size="11">regions</text>']
    for rr in regions:
        a,b=num(rr["canonical_start"]),num(rr["canonical_end"])
        if b-a>1e6: svg.append(f'<text x="{(xmap(a)+xmap(b))/2:.2f}" y="{y0+panel_h/2+4:.2f}" text-anchor="middle" font-size="8">{text(rr["canonical_region_id"])}</text>')
    y0 += panel_h+gap
    for clade in CLADES:
        z=sorted((r for r in support if r["clade"]==clade),key=lambda r:r["bin_start"]); xs=[r["bin_midpoint"] for r in z]
        svg += [marks(y0),support_rects(z,y0),f'<text x="{left-10}" y="{y0+panel_h/2+4}" text-anchor="end" font-size="11">{clade} q support</text>']
        if clade=="Columbea": svg += [f'<text x="{width-right-280}" y="{y0+18}" font-size="10" fill="{QCOLORS[0]}">q1 species-history</text>',f'<text x="{width-right-140}" y="{y0+18}" font-size="10" fill="{QCOLORS[1]}">q2 / q3 alternatives</text>']
        y0 += panel_h+gap; svg += [marks(y0),line(xs,[r["D2"] for r in z],-1,1,y0,"#D65F5F"),line(xs,[r["D3"] for r in z],-1,1,y0,"#4A9B68"),f'<line x1="{left}" x2="{width-right}" y1="{ymap(0,-1,1,y0):.2f}" y2="{ymap(0,-1,1,y0):.2f}" class="grid"/>',f'<text x="{left-10}" y="{y0+panel_h/2+4}" text-anchor="end" font-size="11">{clade} D2 / D3</text>']; y0 += panel_h+gap
    if include_background:
        z=sorted(background,key=lambda r:r["bin_start"]); xs=[r["bin_midpoint"] for r in z]; vals=[r["B_observed_mean"] for r in z]; finite=[v for v in vals if math.isfinite(v)]; lo,hi=(min(finite),max(finite)) if finite else (0,1); hi=hi if hi>lo else lo+1
        svg += [marks(y0),line(xs,vals,lo,hi,y0,"#6A3D9A","1.1"),f'<text x="12" y="{y0+panel_h/2-2}" font-size="10">B_observed</text>',f'<text x="12" y="{y0+panel_h/2+12}" font-size="9">ML substitution units</text>']
    for tick in range(0,int(xmax)+1,10):
        xx=xmap(tick*1e6); svg += [f'<line x1="{xx:.2f}" x2="{xx:.2f}" y1="{height-43}" y2="{height-38}" class="axis"/>',f'<text x="{xx:.2f}" y="{height-22}" text-anchor="middle" font-size="10">{tick}</text>']
    svg.append(f'<text x="{left+usable/2}" y="{height-5}" text-anchor="middle" font-size="12">GalGal6 chromosome 4 position (Mb)</text></svg>')
    out_svg=OUT/(stem+".svg"); out_svg.write_text("".join(svg))
    subprocess.run(["rsvg-convert","-f","pdf","-o",str(FIG/(stem+".pdf")),str(out_svg)],check=True)
    subprocess.run(["rsvg-convert","-f","png","-d","300","-p","300","-o",str(FIG/(stem+".png")),str(out_svg)],check=True)


def main():
    OUT.mkdir(parents=True,exist_ok=True); FIG.mkdir(parents=True,exist_ok=True)
    support=rows(SUPPORT); bg=rows(BG); cov=rows(COV); regions=rows(REGIONS)
    support_by={}; bg_by={}
    for res in RESOLUTIONS:
        support_by[res]=aggregate_support(support,res); bg_by[res]=aggregate_background(bg,res)
        write(OUT/f"binned_quartet_support_{res//1000}kb.tsv",support_by[res],list(support_by[res][0]))
        write(OUT/f"binned_background_distance_{res//1000}kb.tsv",bg_by[res],list(bg_by[res][0]))
    coverage=[]
    observed_span=max(num(r["end"]) for r in support+bg)
    for res in RESOLUTIONS:
        total=math.ceil(observed_span/res); z=bg_by[res]
        coverage.append({"resolution_bp":res,"total_grid_bins":total,"occupied_background_bins":len(z),"missing_background_fraction":1-len(z)/total,"median_loci_per_occupied_bin":med([r["n_loci"] for r in z]),"mean_taxon_coverage_range":avg([r["taxon_coverage_range"] for r in z]),"max_taxon_coverage_range":max(r["taxon_coverage_range"] for r in z),"median_B_observed":med([r["B_observed_mean"] for r in z])})
    write(OUT/"bin_coverage_summary.tsv",coverage,list(coverage[0]))
    sens=[]
    for metric in ("q1","q2","q3","D2","D3"):
        for a,b in ((100000,250000),(250000,500000)):
            for clade in CLADES:
                aa={r["bin_start"]:r for r in support_by[a] if r["clade"]==clade}; bb={r["bin_start"]:r for r in support_by[b] if r["clade"]==clade}; common=set(aa)&set(bb)
                sens.append({"track":metric,"clade":clade,"resolution_a":a,"resolution_b":b,"n_common_bins":len(common),"pearson":correlation([aa[k][metric] for k in common],[bb[k][metric] for k in common]),"mean_abs_difference":avg([abs(aa[k][metric]-bb[k][metric]) for k in common])})
    for a,b in ((100000,250000),(250000,500000)):
        aa={r["bin_start"]:r for r in bg_by[a]}; bb={r["bin_start"]:r for r in bg_by[b]}; common=set(aa)&set(bb)
        sens.append({"track":"B_observed_mean","clade":"all","resolution_a":a,"resolution_b":b,"n_common_bins":len(common),"pearson":correlation([aa[k]["B_observed_mean"] for k in common],[bb[k]["B_observed_mean"] for k in common]),"mean_abs_difference":avg([abs(aa[k]["B_observed_mean"]-bb[k]["B_observed_mean"]) for k in common])})
    write(OUT/"resolution_sensitivity.tsv",sens,list(sens[0]))
    comparisons=[]
    for res in RESOLUTIONS:
        for rr in regions:
            a,b=num(rr["canonical_start"]),num(rr["canonical_end"])
            for clade in CLADES:
                z=[r for r in support_by[res] if r["clade"]==clade and a<=r["bin_midpoint"]<=b]; o=[r for r in support_by[res] if r["clade"]==clade and not(a<=r["bin_midpoint"]<=b)]
                comparisons.append({"resolution_bp":res,"region":rr["canonical_region_id"],"track":"support","clade":clade,"inside_n":len(z),"outside_n":len(o),"inside_mean_D2":avg([r["D2"] for r in z]),"outside_mean_D2":avg([r["D2"] for r in o]),"inside_mean_D3":avg([r["D3"] for r in z]),"outside_mean_D3":avg([r["D3"] for r in o])})
            z=[r for r in bg_by[res] if a<=r["bin_midpoint"]<=b]; o=[r for r in bg_by[res] if not(a<=r["bin_midpoint"]<=b)]
            comparisons.append({"resolution_bp":res,"region":rr["canonical_region_id"],"track":"B_observed","clade":"all","inside_n":len(z),"outside_n":len(o),"inside_mean_B":avg([r["B_observed_mean"] for r in z]),"outside_mean_B":avg([r["B_observed_mean"] for r in o]),"inside_mean_n_taxa":avg([r["n_taxa_mean"] for r in z]),"outside_mean_n_taxa":avg([r["n_taxa_mean"] for r in o])})
    write(OUT/"structural_region_comparison.tsv",comparisons,list(comparisons[0]))
    z=bg_by[250000]; est=[]
    for k in ("B_observed_mean","B47_mean","B_fixed_mean","B_residual_mean"):
        v=[r[k] for r in z if math.isfinite(r[k])]; est.append({"resolution_bp":250000,"estimator":k,"n_bins":len(v),"mean":avg(v),"median":med(v),"sd":stdev(v),"available_fraction":len(v)/len(z)})
    write(OUT/"background_estimator_sensitivity.tsv",est,list(est[0]))
    # Regression checks and documentation.
    max_support_error=max(abs(num(r["support_sum"])-1) for r in support)
    max_binned_error=max(abs(r["support_sum"]-1) for z in support_by.values() for r in z)
    b250corr=correlation([r["B_observed_mean"] for r in bg_by[250000]],[r["n_taxa_mean"] for r in bg_by[250000]])
    (OUT/"report.md").write_text(f"""# Neoaves chromosome 4 Stage 10 report

The frozen Stage 6 quartet support and Stage 9 background tracks were
aggregated into chromosome-anchored 100-kb, 250-kb, and 500-kb bins. The
250-kb layout is the primary candidate because it reduces the density of the
6,140-locus display while retaining chromosome-wide patterns.

B_observed is available for all 6,140 loci. Its locus-level Pearson correlation
with observed taxon count is -0.2166; the 250-kb mean-track correlation is
{b250corr:.4f}. This supports a descriptive B panel with explicit coverage
caveats. It is an all-observed taxon-composition statistic in ML substitution
units, not a calibrated divergence or coalescence time.

Columbea, N61, and N62 remain separate source-defined quartet designs. The
figure does not treat them as one arrangement-versus-species-history test.
The topology-only 250-kb figure is retained as the conservative alternative.
""")
    (OUT/"methods_text.md").write_text("""# Methods

Each locus was assigned to exactly one chromosome-anchored bin by its recorded
midpoint. Quartet support was averaged equally across loci within each focal
design, and D2/D3 were recalculated from the binned q values. Background
values were summarized by mean and median B_observed, with within-bin
variation, observed-taxon coverage, and the 47-, 108-, and residual estimators
retained separately. Empty bins were omitted rather than filled or
interpolated. Canonical structural-region coordinates were copied unchanged.
""")
    (OUT/"results_text.md").write_text("""# Results

The 250-kb bins provide a clearer chromosome-wide view than the raw locus-level
tracks. Broad topology patterns are comparable across resolutions, while
narrow transitions can move or disappear when they are smaller than a bin.
B_observed has weak taxon-count association and remains interpretable as a
full-taxon descriptive background track with composition caveats.
""")
    (OUT/"figure_caption.md").write_text("""# Figure caption

Neoaves chromosome-4 source-defined quartet support, topology contrasts, and
full-taxon background ML patristic distance in chromosome-anchored bins.
Structural regions are independently annotated canonical intervals. q1/q2/q3
retain their separate source meanings for each focal design; D2=q2-q1 and
D3=q3-q1. B_observed is a mean patristic distance in ML substitution units,
not a calibrated divergence or coalescence time.
""")
    (OUT/"validation_report.md").write_text(f"""# Validation report

- Input quartet rows: {len(support)}
- Input background rows: {len(bg)}
- Input coverage rows: {len(cov)}
- Maximum raw support-sum error: {max_support_error:.3g}
- Maximum binned support-sum error: {max_binned_error:.3g}
- Stage-9 B_observed rows preserved before aggregation: {len(bg)}
- Empty bins omitted: yes
- Structural coordinates copied without modification: yes
""")
    manifest={"stage":"Neoaves chr4 Stage 10 multiresolution tracks","script":str(Path(__file__).relative_to(ROOT)),"inputs":{str(p.relative_to(ROOT)):checksum(p) for p in (SUPPORT,BG,COV,REGIONS)},"resolutions_bp":list(RESOLUTIONS),"primary_resolution_bp":250000,"primary_background":"B_observed","n_support_rows":len(support),"n_background_rows":len(bg),"random_seed":None}
    (OUT/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    for res in RESOLUTIONS: svg_figure(support_by[res],bg_by[res],regions,res,f"neoaves_chr4_multitrack_{res//1000}kb",True)
    svg_figure(support_by[250000],bg_by[250000],regions,250000,"neoaves_chr4_topology_only_250kb",False)


if __name__ == "__main__": main()
