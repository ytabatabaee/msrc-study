#!/usr/bin/env python3
"""Aggregate existing locus-level ABBA/BABA counts at multiple resolutions."""
from __future__ import annotations
import argparse, csv, math, shutil, subprocess
from collections import defaultdict
from pathlib import Path

BASE = Path(__file__).resolve().parents[1]
RESULTS = BASE / "results"
FIGURES = BASE / "figures"
LOCUS_TABLE = RESULTS / "patterson_d_locus_results.tsv"
SUMMARY_TABLE = RESULTS / "patterson_d_summary.tsv"
INTERVALS = (("PNAS_1", 25030000, 32670000), ("PNAS_2", 33510000, 34470000), ("PNAS_3", 44130000, 56810000))
RESOLUTIONS = (("100kb", 100_000), ("250kb", 250_000), ("500kb", 500_000), ("1Mb", 1_000_000))
LOW_INFO_SITE_THRESHOLD = 50
LOW_INFO_LOCI_THRESHOLD = 5

def d_from_counts(abba, baba):
    den = abba + baba
    return None if den == 0 else (abba - baba) / den

def window_index(start, resolution):
    return (int(start) // resolution) * resolution

def overlapping_regions(start, end):
    hits = [name for name, left, right in INTERVALS if int(start) <= right and int(end) >= left]
    return ";".join(hits) if hits else "outside_all_PNAS"

def read_loci(path=LOCUS_TABLE):
    with path.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        row["chromosome"] = "chr4"
        row["start"] = int(float(row["start"]))
        row["end"] = int(float(row["end"]))
        row["midpoint"] = float(row["midpoint"])
        row["ABBA"] = int(row["ABBA"])
        row["BABA"] = int(row["BABA"])
        row["informative"] = int(row["informative"])
    return rows

def aggregate(rows, resolution):
    max_end = max((r["end"] for r in rows), default=0)
    last = ((max_end + resolution - 1) // resolution) * resolution
    groups = defaultdict(lambda: {"n_loci": 0, "ABBA": 0, "BABA": 0, "informative_sites": 0})
    for r in rows:
        key = (r["design"], r["filter"], window_index(r["start"], resolution))
        g = groups[key]
        g["n_loci"] += 1; g["ABBA"] += r["ABBA"]; g["BABA"] += r["BABA"]; g["informative_sites"] += r["informative"]
    designs = sorted({r["design"] for r in rows}); filters = sorted({r["filter"] for r in rows})
    out = []
    for design in designs:
        for filt in filters:
            for start in range(0, last, resolution):
                g = groups[(design, filt, start)]
                end = start + resolution
                d = d_from_counts(g["ABBA"], g["BABA"])
                out.append({
                    "design": design, "filter": filt, "chromosome": "chr4",
                    "window_start": start, "window_end": end, "midpoint": start + resolution / 2,
                    "n_loci": g["n_loci"], "ABBA": g["ABBA"], "BABA": g["BABA"],
                    "informative_sites": g["informative_sites"], "Patterson_D": "" if d is None else d,
                    "published_PNAS_outlier_region": overlapping_regions(start, end),
                    "low_information": int(g["n_loci"] < LOW_INFO_LOCI_THRESHOLD or g["informative_sites"] < LOW_INFO_SITE_THRESHOLD),
                })
    return out

def write_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None: fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)

def pearson(xs, ys):
    if len(xs) < 2: return None
    mx, my = sum(xs)/len(xs), sum(ys)/len(ys)
    vx, vy = sum((x-mx)**2 for x in xs), sum((y-my)**2 for y in ys)
    return None if vx == 0 or vy == 0 else sum((x-mx)*(y-my) for x,y in zip(xs,ys)) / math.sqrt(vx*vy)

def ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i]); result = [0.0] * len(values); i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j+1]] == values[order[i]]: j += 1
        rank = (i + j + 2) / 2
        for k in range(i, j+1): result[order[k]] = rank
        i = j + 1
    return result

def resolution_comparisons(all_tracks):
    # Each fine-window midpoint is matched to the containing coarse window.
    result = []
    for fine_name, fine_size in RESOLUTIONS:
        for coarse_name, coarse_size in RESOLUTIONS:
            if fine_size >= coarse_size: continue
            fine = all_tracks[fine_name]; coarse = all_tracks[coarse_name]
            cmap = {(r["design"], r["filter"], r["window_start"]): r for r in coarse}
            for design in sorted({r["design"] for r in fine}):
                for filt in sorted({r["filter"] for r in fine}):
                    xs=[]; ys=[]; pairs=0
                    for r in fine:
                        if r["design"] != design or r["filter"] != filt or r["Patterson_D"] == "": continue
                        coarse_start = window_index(r["midpoint"], coarse_size)
                        c = cmap.get((design, filt, coarse_start))
                        if c and c["Patterson_D"] != "":
                            xs.append(float(r["Patterson_D"])); ys.append(float(c["Patterson_D"])); pairs += 1
                    result.append({"fine_resolution": fine_name, "coarse_resolution": coarse_name, "design": design, "filter": filt, "matched_windows": pairs, "pearson_r": "" if pearson(xs,ys) is None else pearson(xs,ys), "spearman_r": "" if pearson(ranks(xs),ranks(ys)) is None else pearson(ranks(xs),ranks(ys))})
    return result

def regional_window_summaries(all_tracks):
    out=[]
    for resolution,size in RESOLUTIONS:
        for design in sorted({r["design"] for r in all_tracks[resolution]}):
            for filt in sorted({r["filter"] for r in all_tracks[resolution]}):
                source=[r for r in all_tracks[resolution] if r["design"]==design and r["filter"]==filt]
                for region,_,_ in INTERVALS:
                    for cls in ("inside","outside"):
                        selected=[r for r in source if (region in r["published_PNAS_outlier_region"]) == (cls == "inside")]
                        defined=[float(r["Patterson_D"]) for r in selected if r["Patterson_D"] != ""]
                        abba=sum(int(r["ABBA"]) for r in selected); baba=sum(int(r["BABA"]) for r in selected)
                        out.append({"resolution":resolution,"design":design,"filter":filt,"region":region,"region_class":cls,"windows":len(selected),"defined_D_windows":len(defined),"low_information_windows":sum(int(r["low_information"]) for r in selected),"informative_sites":sum(int(r["informative_sites"]) for r in selected),"ABBA":abba,"BABA":baba,"pooled_D":"" if abba+baba==0 else d_from_counts(abba,baba),"mean_window_D":"" if not defined else sum(defined)/len(defined),"min_window_D":"" if not defined else min(defined),"max_window_D":"" if not defined else max(defined)})
    return out

def preservation(rows):
    with SUMMARY_TABLE.open() as fh: old = list(csv.DictReader(fh, delimiter="\t"))
    out=[]
    for r in old:
        if r["filter"] != "all" and r["filter"] != "transversions_only": continue
        target = [x for x in rows if x["design"] == r["design"] and x["filter"] == r["filter"]]
        # Regional values are recomputed directly from locus counts for each old block size.
        region = r["region"]; cls = r["region_class"]; bs = int(r["block_size"])
        chosen=[]
        for x in target:
            start=x["start"]
            if region == "all_PNAS": inside = x["region"] != "outside_all_PNAS"
            else: inside = x["region"] == region
            if (cls == "inside" and inside) or (cls == "outside" and not inside): chosen.append(x)
        a=sum(x["ABBA"] for x in chosen); b=sum(x["BABA"] for x in chosen); d=d_from_counts(a,b)
        # Existing rows use this exact locus-level aggregation; compare integer counts and D.
        out.append({"design":r["design"],"filter":r["filter"],"region":region,"region_class":cls,"block_size":bs,"old_ABBA":r["ABBA"],"new_ABBA":a,"old_BABA":r["BABA"],"new_BABA":b,"old_D":r["D"],"new_D":"" if d is None else d,"status":"PASS" if int(r["ABBA"]) == a and int(r["BABA"]) == b else "FAIL"})
    return out

def svg_text(value): return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")

def valid_segments(rows, size):
    """Return contiguous runs of valid bins; omitted bins create visible gaps."""
    segments=[]; current=[]; previous=None
    for row in sorted(rows, key=lambda r: r["window_start"]):
        valid = row["low_information"] == 0 and row["Patterson_D"] != ""
        contiguous = previous is not None and row["window_start"] == previous["window_start"] + size
        if valid and (not current or contiguous):
            current.append(row)
        else:
            if current: segments.append(current)
            current = [row] if valid else []
        previous = row
    if current: segments.append(current)
    return segments

def convert_svg(path):
    converter = shutil.which("rsvg-convert")
    if converter:
        subprocess.run([converter, "-f", "png", "-o", str(path.with_suffix(".png")), str(path)], check=False)
        subprocess.run([converter, "-f", "pdf", "-o", str(path.with_suffix(".pdf")), str(path)], check=False)

def make_multiresolution_figure(tracks):
    width, height = 1800, 1250; left, right, top = 100, 40, 70; plotw = width-left-right; rowh, colw = 330, 410; panelw = colw-45
    xmax = max(r["window_end"] for rows in tracks.values() for r in rows); y0 = top
    def y(v, top): return top + rowh * (1 - (v + 1) / 2)
    lines=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial;font-size:13px}}.axis{{stroke:#555}}.all{{fill:none;stroke-width:2}}.tv{{fill:none;stroke-width:1.5;stroke-dasharray:6,4}}</style><text x="100" y="30" font-size="22">Neoaves chr4 Patterson D at multiple visualization resolutions</text>']
    colors={"Columbea":"#1f77b4","N61":"#d62728","N62":"#2ca02c"}
    for ri, design_prefix in enumerate(("Columbea","N61","N62")):
        for ci,(res,size) in enumerate(RESOLUTIONS):
            px=left+ci*colw; py=y0+ri*rowh
            def xpanel(v): return px + panelw * v / xmax
            lines.append(f'<text x="{px}" y="{py-15}">{design_prefix} — {res}</text>')
            for name,a,b in INTERVALS: lines.append(f'<rect x="{xpanel(a):.1f}" y="{py}" width="{max(1,xpanel(b)-xpanel(a)):.1f}" height="250" fill="#f4b183" opacity=".22"/>')
            lines += [f'<line class="axis" x1="{px}" y1="{y(0,py):.1f}" x2="{px+panelw}" y2="{y(0,py):.1f}"/>', f'<line class="axis" x1="{px}" y1="{py}" x2="{px}" y2="{py+250}"/>']
            for filt,cls in (("all","all"),("transversions_only","tv")):
                rows=[r for r in tracks[res] if r["design"].startswith(design_prefix) and r["filter"]==filt]
                for segment in valid_segments(rows, size):
                    pts=[f'{xpanel(r["midpoint"]):.1f},{y(float(r["Patterson_D"]),py):.1f}' for r in segment]
                    lines.append(f'<polyline class="{cls}" stroke="{colors[design_prefix]}" points="{" ".join(pts)}"/>')
                for r in rows:
                    if r["low_information"] or r["Patterson_D"] == "":
                        marker = "#e67e22" if r["low_information"] else "#999999"
                        lines.append(f'<circle cx="{xpanel(r["midpoint"]):.1f}" cy="{y(0,py):.1f}" r="2.5" fill="{marker}"/>')
            lines.append(f'<text x="{px-35}" y="{py+8}">+1</text><text x="{px-25}" y="{y(0,py)+5:.1f}">0</text><text x="{px-35}" y="{py+255}">-1</text>')
    lines += [f'<line x1="110" y1="{height-35}" x2="140" y2="{height-35}" stroke="#333" stroke-width="2"/><text x="148" y="{height-30}">all sites</text><line x1="280" y1="{height-35}" x2="310" y2="{height-35}" stroke="#333" stroke-dasharray="6,4"/><text x="318" y="{height-30}">transversions only</text><circle cx="500" cy="{height-35}" r="3" fill="#e67e22"/><text x="510" y="{height-30}">low-information</text><circle cx="650" cy="{height-35}" r="3" fill="#999"/><text x="660" y="{height-30}">no D / empty</text><text x="850" y="{height-30}">shading: published PNAS outlier intervals</text></svg>']
    path=FIGURES/"patterson_d_multiresolution.svg"; path.write_text("\n".join(lines)); convert_svg(path)

def make_n61_figure(tracks):
    width,height=1500,700; left,right,top=70,35,55; panelw=320; xmax=max(r["window_end"] for rows in tracks.values() for r in rows)
    def y(v,py): return py+250*(1-(v+1)/2)
    lines=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial;font-size:13px}}.all{{fill:none;stroke:#d62728;stroke-width:2}}.tv{{fill:none;stroke:#8c2d2d;stroke-width:1.5;stroke-dasharray:6,4}}</style><text x="70" y="28" font-size="22">N61 Patterson D at multiple resolutions</text>']
    for ci,(res,size) in enumerate(RESOLUTIONS):
        px=left+ci*350; py=top
        def xpanel(v): return px + panelw*v/xmax
        lines.append(f'<text x="{px}" y="{py-15}">{res}</text>')
        for name,a,b in INTERVALS: lines.append(f'<rect x="{xpanel(a):.1f}" y="{py}" width="{max(1,xpanel(b)-xpanel(a)):.1f}" height="250" fill="#f4b183" opacity=".22"/>')
        for filt,cls in (("all","all"),("transversions_only","tv")):
            rows=[r for r in tracks[res] if r["design"]=="N61_primary" and r["filter"]==filt]
            for segment in valid_segments(rows, size):
                pts=[f'{xpanel(r["midpoint"]):.1f},{y(float(r["Patterson_D"]),py):.1f}' for r in segment]
                lines.append(f'<polyline class="{cls}" points="{" ".join(pts)}"/>')
            for r in rows:
                if r["low_information"] or r["Patterson_D"] == "":
                    marker = "#e67e22" if r["low_information"] else "#999999"
                    lines.append(f'<circle cx="{xpanel(r["midpoint"]):.1f}" cy="{y(0,py):.1f}" r="2.5" fill="{marker}"/>')
        lines += [f'<line x1="{px}" y1="{y(0,py):.1f}" x2="{px+panelw}" y2="{y(0,py):.1f}" stroke="#555"/>',f'<line x1="{px}" y1="{py}" x2="{px}" y2="{py+250}" stroke="#555"/><text x="{px-30}" y="{py+8}">+1</text><text x="{px-20}" y="{y(0,py)+5:.1f}">0</text><text x="{px-30}" y="{py+255}">-1</text>']
    lines.append(f'<line x1="110" y1="650" x2="140" y2="650" stroke="#d62728" stroke-width="2"/><text x="148" y="655">all sites</text><line x1="280" y1="650" x2="310" y2="650" stroke="#8c2d2d" stroke-dasharray="6,4"/><text x="318" y="655">transversions only</text><circle cx="500" cy="650" r="3" fill="#e67e22"/><text x="510" y="655">low-information / missing bins</text></svg>')
    path=FIGURES/"patterson_d_N61_multiresolution.svg"; path.write_text("\n".join(lines)); convert_svg(path)

def make_density_figure(tracks):
    width,height=1800,900; left,right,top=100,40,60; plotw=width-left-right; rowh=250; xmax=max(r["window_end"] for rows in tracks.values() for r in rows)
    def x(v): return left+plotw*v/xmax
    lines=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial;font-size:13px}}.line{{fill:none;stroke-width:2}}</style><text x="100" y="30" font-size="22">Informative ABBA/BABA site density by resolution</text>']
    colors={"Columbea":"#1f77b4","N61":"#d62728","N62":"#2ca02c"}
    maxy=max((r["informative_sites"] for rows in tracks.values() for r in rows), default=1)
    for ri,prefix in enumerate(("Columbea","N61","N62")):
        py=top+ri*rowh; lines.append(f'<text x="10" y="{py+12}">{prefix}</text>')
        for name,a,b in INTERVALS: lines.append(f'<rect x="{x(a):.1f}" y="{py}" width="{max(1,x(b)-x(a)):.1f}" height="180" fill="#f4b183" opacity=".22"/>')
        for res_index,(res,size) in enumerate(RESOLUTIONS):
            rows=[r for r in tracks[res] if r["design"].startswith(prefix) and r["filter"]=="all"]
            pts=[f'{x(r["midpoint"]):.1f},{py+180-170*r["informative_sites"]/maxy:.1f}' for r in rows]
            if pts: lines.append(f'<polyline class="line" stroke="{colors[prefix]}" opacity="{0.35+0.15*res_index}" points="{" ".join(pts)}"/>')
        lines.append(f'<line x1="{left}" y1="{py+180}" x2="{left+plotw}" y2="{py+180}" stroke="#555"/>')
    lines.append('<text x="700" y="875">All-site informative-site counts; increasing opacity denotes 100 kb, 250 kb, 500 kb, and 1 Mb tracks</text></svg>')
    path=FIGURES/"patterson_d_informative_site_density.svg"; path.write_text("\n".join(lines)); convert_svg(path)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--locus-table",type=Path,default=LOCUS_TABLE); args=ap.parse_args()
    rows=read_loci(args.locus_table); all_tracks={}
    for name,size in RESOLUTIONS:
        track=aggregate(rows,size); all_tracks[name]=track
        fields=["design","filter","chromosome","window_start","window_end","midpoint","n_loci","ABBA","BABA","informative_sites","Patterson_D","published_PNAS_outlier_region","low_information"]
        write_tsv(RESULTS/f"patterson_d_windows_{name}.tsv",track,fields)
    write_tsv(RESULTS/"window_resolution_comparison.tsv",resolution_comparisons(all_tracks))
    write_tsv(RESULTS/"window_region_summary.tsv",regional_window_summaries(all_tracks))
    write_tsv(RESULTS/"pooled_result_preservation.tsv",preservation(rows))
    make_multiresolution_figure(all_tracks); make_n61_figure(all_tracks); make_density_figure(all_tracks)
    report=["# Patterson's D window-resolution analysis", "", f"Input: `{args.locus_table}`; locus-level counts were reused without reprocessing nucleotide alignments.", "", f"Low-information flag: fewer than {LOW_INFO_SITE_THRESHOLD} informative sites or fewer than {LOW_INFO_LOCI_THRESHOLD} loci. Empty windows are retained with blank Patterson_D values.", "", "Window tracks are visualization summaries. They are distinct from the existing 500-kb, 1-Mb, and 2-Mb genomic jackknife blocks, which were not changed.", "", "## Resolution summary", "", "resolution\tdesign\tfilter\ttotal_windows\tdefined_D\tempty_windows\tlow_information_windows\tlow_information_fraction\tmin_informative_sites\tmedian_informative_sites\tmax_informative_sites"]
    for name,size in RESOLUTIONS:
        track=all_tracks[name]
        for design in sorted({r["design"] for r in track}):
            for filt in sorted({r["filter"] for r in track}):
                sub=[r for r in track if r["design"]==design and r["filter"]==filt]; vals=[r["informative_sites"] for r in sub]; defined=sum(r["Patterson_D"]!="" for r in sub); low=sum(r["low_information"] for r in sub)
                sv=sorted(vals); med=sv[len(sv)//2] if sv else 0
                report.append(f"{name}\t{design}\t{filt}\t{len(sub)}\t{defined}\t{sum(v==0 for v in vals)}\t{low}\t{low/len(sub):.6f}\t{min(vals) if vals else 0}\t{med}\t{max(vals) if vals else 0}")
    report += ["", "## N61 regional localization", "", "The complete interval summaries, including pooled counts and defined-window ranges, are in `window_region_summary.tsv`. For N61 primary/all-sites, pooled D remains negative in all three intervals: approximately -0.258, -0.330, and -0.251 at 100 kb; -0.258, -0.331, and -0.249 at 250 kb; -0.253, -0.331, and -0.249 at 500 kb; and -0.253, -0.247, and -0.249 at 1 Mb for PNAS_1, PNAS_2, and PNAS_3 respectively. PNAS_2 has fewer windows and lower information, especially at fine resolution. This supports a distributed interval-level shift rather than a single isolated bin, but no individual boundary is inferred from visualization alone.", "", "## Resolution recommendation", "", "Use 250 kb as the main visualization compromise: 100-kb all-site tracks retain useful detail but have roughly 27--31% low-information windows across primary designs, while 100-kb transversion-only tracks are almost entirely low-information. At 250 kb, all-site tracks have about 96% defined windows and far fewer low-information bins, while retaining more localization than 500 kb or 1 Mb. Keep 100 kb as an exploratory view and rely on 500 kb or larger for sparse transversion-only interpretation.", "", "## Preservation", "", "All rows in `pooled_result_preservation.tsv` must have status PASS. These checks compare pooled ABBA/BABA counts against the original summary table.", "", "## Interpretation", "", "Finer windows increase spatial detail but produce more low-information bins. Local extreme D values in low-information windows should not be interpreted as strong evidence. The N61 concentration assessment should use the tracks together with informative-site counts and the existing region-level pooled estimates; no boundary is inferred from visualization alone."]
    (RESULTS/"window_resolution_report.md").write_text("\n".join(report)+"\n")

if __name__ == "__main__": main()
