#!/usr/bin/env python3
"""Real nucleotide ABBA--BABA analysis of the Stiller 63K alignments."""
from __future__ import annotations
import argparse, csv, json, math, re, shutil, subprocess, sys, tarfile
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
BASE = Path(__file__).resolve().parents[1]
EXT = ROOT / "empirical/neoaves_chr4/external/stiller2024"
ARCHIVE = EXT / "63430.alns.tar.gz"
MANIFEST = ROOT / "empirical/neoaves_chr4/results/stage4d_locus_manifest.tsv"
CONFIG = BASE / "config/frozen_quartets.tsv"
RESULTS = BASE / "results"; FIGURES = BASE / "figures"
INTERVALS = [("PNAS_1", 25030000, 32670000), ("PNAS_2", 33510000, 34470000), ("PNAS_3", 44130000, 56810000)]
BASES = set("ACGT")

def stem(name):
    base=Path(name).name
    if base.startswith("._"): return ""
    base=re.sub(r"\.(?:fasta|fas|fa|aln)$", "", base, flags=re.I)
    base=re.sub(r"_fasta\.treefile$", "", base)
    return base.replace(".1k.start", "_1k_start")

def parse_fasta(handle):
    seqs, name, chunks = {}, None, []
    for raw in handle:
        line = raw.decode() if isinstance(raw, bytes) else raw
        line = line.strip()
        if not line: continue
        if line.startswith(">"):
            if name is not None: seqs[name] = "".join(chunks).upper()
            name, chunks = line[1:].split()[0], []
        else: chunks.append(line.replace(" ", ""))
    if name is not None: seqs[name] = "".join(chunks).upper()
    return seqs

def count_sites(seqs, transversions_only=False):
    """Return ABBA, BABA, excluded, informative for one four-sequence alignment."""
    keys = ("P1", "P2", "P3", "O")
    n = min((len(seqs[k]) for k in keys), default=0)
    abba = baba = excluded = informative = 0
    transitions = {frozenset(x) for x in (("A","G"),("C","T"))}
    for i in range(n):
        a = [seqs[k][i] for k in keys]
        if any(x not in BASES for x in a): excluded += 1; continue
        p1,p2,p3,o = a
        if p3 == o or p1 == p2 or (p1 not in (o,p3)) or (p2 not in (o,p3)) or p3 == o:
            continue
        if p3 == o: continue
        if frozenset((o,p3)) in transitions and transversions_only: continue
        if p1 == o and p2 == p3 and p2 != o:
            abba += 1; informative += 1
        elif p2 == o and p1 == p3 and p1 != o:
            baba += 1; informative += 1
    return abba, baba, excluded, informative

def patterson_d(abba, baba):
    den = abba + baba
    return (abba - baba) / den if den else None

def load_rows():
    with MANIFEST.open() as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t") if r["chromosome"] == "chr4"]

def load_designs():
    with CONFIG.open() as fh: return list(csv.DictReader(fh, delimiter="\t"))

def region_for(start):
    hits = [name for name,a,b in INTERVALS if a <= start <= b]
    return hits[0] if hits else "outside_all_PNAS"

def write_tsv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore"); w.writeheader(); w.writerows(rows)

def jackknife(values, block_size):
    # values are (block, abba, baba); blocks are physical and therefore linked-aware.
    grouped = defaultdict(lambda: [0,0])
    for block,a,b in values: grouped[block][0] += a; grouped[block][1] += b
    total_a=sum(x[0] for x in grouped.values()); total_b=sum(x[1] for x in grouped.values())
    theta=patterson_d(total_a,total_b); blocks=list(grouped)
    if theta is None or len(blocks) < 2: return theta, None, len(blocks)
    reps=[]
    for block,(a,b) in grouped.items(): reps.append(patterson_d(total_a-a,total_b-b))
    reps=[x for x in reps if x is not None]
    mean=sum(reps)/len(reps)
    se=math.sqrt((len(reps)-1)/len(reps)*sum((x-mean)**2 for x in reps))
    return theta,se,len(reps)

def paired_difference(records, design, region, block_size, filt):
    chosen=[r for r in records if r["filter"]==filt and r["design"]==design]
    def inside(r): return (r["region"] == region) if region != "all_PNAS" else r["class"] == "inside"
    groups=defaultdict(lambda: [0,0,0,0])
    for r in chosen:
        k=int(r["start"])//block_size; a=int(r["ABBA"]); b=int(r["BABA"])
        if inside(r): groups[k][0]+=a; groups[k][1]+=b
        else: groups[k][2]+=a; groups[k][3]+=b
    ia=sum(v[0] for v in groups.values()); ib=sum(v[1] for v in groups.values())
    oa=sum(v[2] for v in groups.values()); ob=sum(v[3] for v in groups.values())
    di=patterson_d(ia,ib); do=patterson_d(oa,ob)
    diff=None if di is None or do is None else di-do
    if len(groups)<2 or diff is None: return diff,None
    reps=[]
    for k,(a,b,c,d) in groups.items():
        x=patterson_d(ia-a,ib-b); y=patterson_d(oa-c,ob-d)
        if x is not None and y is not None: reps.append(x-y)
    if len(reps)<2: return diff,None
    mean=sum(reps)/len(reps); se=math.sqrt((len(reps)-1)/len(reps)*sum((x-mean)**2 for x in reps))
    return diff,se

def summarize(records, design, region, cls, block_size):
    if region == "all_PNAS":
        sub=[r for r in records if r["design"]==design and ((r["class"]=="inside") if cls=="inside" else (r["class"]=="outside"))]
    elif cls == "inside":
        sub=[r for r in records if r["design"]==design and r["region"]==region]
    else:
        sub=[r for r in records if r["design"]==design and r["region"] != region]
    ab=sum(int(r["ABBA"]) for r in sub); ba=sum(int(r["BABA"]) for r in sub)
    vals=[(int(r["start"])//block_size,int(r["ABBA"]),int(r["BABA"])) for r in sub]
    d,se,nb=jackknife(vals,block_size)
    return {"design":design,"region":region,"region_class":cls,"block_size":block_size,"n_loci":len(sub),"informative_sites":ab+ba,"ABBA":ab,"BABA":ba,"D":"" if d is None else d,"jackknife_se":"" if se is None else se,"blocks":nb}

def svg_fig(window_rows, pooled_rows):
    FIGURES.mkdir(parents=True, exist_ok=True)
    width,height=1400,760; left,right,top=90,30,55; plotw=width-left-right
    usable=[r for r in window_rows if r.get("D") not in ("",None) and r.get("filter")=="all"]
    xmax=max([int(r["end"]) for r in window_rows] or [1])
    def x(v): return left+plotw*float(v)/xmax
    lines=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>', '<style>text{font-family:Arial,sans-serif;font-size:14px}.axis{stroke:#333}.line{fill:none;stroke-width:2}</style>', '<text x="90" y="28" font-size="20">Neoaves chr4 nucleotide Patterson D (1-Mb pooled windows)</text>']
    for name,a,b in INTERVALS:
        lines.append(f'<rect x="{x(a):.1f}" y="50" width="{max(1,x(b)-x(a)):.1f}" height="620" fill="#f4b183" opacity="0.25"/>')
        lines.append(f'<text x="{x(a)+3:.1f}" y="68">{name}</text>')
    y0=100; h=480
    lines += [f'<line class="axis" x1="{left}" y1="{y0+h/2}" x2="{left+plotw}" y2="{y0+h/2}"/>', f'<line class="axis" x1="{left}" y1="{y0}" x2="{left}" y2="{y0+h}"/>', f'<text x="12" y="{y0+8}">+1</text>', f'<text x="12" y="{y0+h/2+5}">0</text>', f'<text x="12" y="{y0+h+5}">-1</text>']
    colors={"Columbea_primary":"#1f77b4","N61_primary":"#d62728","N62_primary":"#2ca02c"}
    for design,color in colors.items():
        vals=[r for r in usable if r["design"]==design]
        pts=[]
        for r in vals:
            yy=y0+h*(1-(float(r["D"])+1)/2); pts.append(f'{x(r["midpoint"]):.1f},{yy:.1f}')
        if pts: lines.append(f'<polyline class="line" stroke="{color}" points="{" ".join(pts)}"/>')
    lx=left+20
    for i,(label,color) in enumerate(colors.items()): lines += [f'<line x1="{lx+i*230}" y1="710" x2="{lx+25+i*230}" y2="710" stroke="{color}" stroke-width="4"/>',f'<text x="{lx+32+i*230}" y="715">{label}</text>']
    lines += ['<text x="650" y="748">GalGal6/Stiller chr4 coordinate (bp); shaded bands are published GalGal4 window-start outlier intervals</text>','</svg>']
    out=FIGURES/"patterson_d_spatial_chr4.svg"; out.write_text("\n".join(lines))
    convert = shutil.which("rsvg-convert")
    if convert:
        subprocess.run([convert,"-f","png","-o",str(FIGURES/"patterson_d_spatial_chr4.png"),str(out)],check=False)
        subprocess.run([convert,"-f","pdf","-o",str(FIGURES/"patterson_d_spatial_chr4.pdf"),str(out)],check=False)
    # Compact region comparison using pooled primary designs only.
    comp=[r for r in pooled_rows if r.get("filter")=="all" and str(r.get("block_size"))=="1000000" and r.get("region")=="all_PNAS" and r.get("design"," ").endswith("_primary")]
    bars=[f'<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="600"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial;font-size:14px}}</style><text x="55" y="28" font-size="20">Patterson D: published PNAS intervals vs outside</text>']
    vals=[]
    for j,d in enumerate(("Columbea_primary","N61_primary","N62_primary")):
        pair=[r for r in comp if r["design"]==d]
        ir=next((r for r in pair if r["region_class"]=="inside" and r["D"]!=""),None); orow=next((r for r in pair if r["region_class"]=="outside" and r["D"]!=""),None)
        if ir and orow: vals.append((d,float(ir["D"]),float(orow["D"]),float(ir["jackknife_se"]),float(orow["jackknife_se"])))
    for j,(d,ins,outd,inse,outse) in enumerate(vals):
        x0=120+j*350; scale=180; zero=300
        bars += [f'<text x="{x0}" y="80">{d.replace("_primary","")}</text>', f'<line x1="{x0}" y1="{zero}" x2="{x0+220}" y2="{zero}" stroke="#333"/>']
        for k,(lab,v,se,col) in enumerate((("inside",ins,inse,"#d95f02"),("outside",outd,outse,"#1b9e77"))):
            h=v*scale; y=zero-h if h>=0 else zero; xc=x0+57+k*90; ey=se*scale
            bars.append(f'<rect x="{x0+30+k*90}" y="{y:.1f}" width="55" height="{abs(h):.1f}" fill="{col}"/><line x1="{xc}" y1="{y-ey:.1f}" x2="{xc}" y2="{y+abs(ey):.1f}" stroke="#222" stroke-width="2"/><line x1="{xc-8}" y1="{y-ey:.1f}" x2="{xc+8}" y2="{y-ey:.1f}" stroke="#222"/><line x1="{xc-8}" y1="{y+abs(ey):.1f}" x2="{xc+8}" y2="{y+abs(ey):.1f}" stroke="#222"/><text x="{x0+20+k*90}" y="530">{lab}</text><text x="{x0+35+k*90}" y="{y-6 if h>=0 else y+abs(h)+18:.1f}">{v:.3f}</text>')
    bars.append('<text x="420" y="575">D = (ABBA - BABA)/(ABBA + BABA); 1-Mb block jackknife summaries</text></svg>')
    cmp=FIGURES/"patterson_d_inside_outside.svg"; cmp.write_text("".join(bars))
    if convert:
        subprocess.run([convert,"-f","png","-o",str(FIGURES/"patterson_d_inside_outside.png"),str(cmp)],check=False)
        subprocess.run([convert,"-f","pdf","-o",str(FIGURES/"patterson_d_inside_outside.pdf"),str(cmp)],check=False)

def main(argv=None):
    global RESULTS, FIGURES
    ap = argparse.ArgumentParser(description="Compute nucleotide Patterson's D on the Stiller chr4 alignments")
    ap.add_argument("command", nargs="?", default="all")
    ap.add_argument("--results-dir", type=Path, default=RESULTS)
    ap.add_argument("--figures-dir", type=Path, default=FIGURES)
    args = ap.parse_args(argv)
    if args.command != "all":
        raise SystemExit("the only supported command is 'all'")
    RESULTS = args.results_dir
    FIGURES = args.figures_dir
    RESULTS.mkdir(parents=True, exist_ok=True)
    if not ARCHIVE.exists(): raise SystemExit(f"missing alignment archive: {ARCHIVE}")
    rows=load_rows(); designs=load_designs(); by_stem={stem(r["locus_id"]):r for r in rows}
    records=[]; coverage=defaultdict(int); missing=defaultdict(int)
    with tarfile.open(ARCHIVE,"r:gz") as tf:
        for member in tf:
            if not member.isfile(): continue
            key=stem(member.name)
            if key not in by_stem: continue
            row=by_stem[key]
            seqs=parse_fasta(tf.extractfile(member))
            for d in designs:
                if not all(t in seqs for t in (d["P1"],d["P2"],d["P3"],d["O"])):
                    missing[d["design"]]+=1; continue
                s={k:seqs[d[k]] for k in ("P1","P2","P3","O")}
                for filt in ("all", "transversions_only"):
                    a,b,ex,inf=count_sites(s, transversions_only=(filt == "transversions_only"))
                    records.append({"filter":filt,"design":d["design"],"locus_id":row["locus_id"],"start":row["start"],"end":row["end"],"midpoint":row["midpoint"],"region":region_for(int(row["start"])),"class":"inside" if region_for(int(row["start"])).startswith("PNAS_") else "outside","ABBA":a,"BABA":b,"excluded":ex,"informative":inf,"D":"" if not a+b else patterson_d(a,b)})
                coverage[d["design"]]+=1
    fields=list(records[0]) if records else ["design","locus_id"]
    write_tsv(RESULTS/"patterson_d_locus_results.tsv",records,fields)
    summary=[]
    for d in [x["design"] for x in designs]:
      for filt in ("all", "transversions_only"):
        for region,_,_ in INTERVALS:
            for cls in ("inside","outside"):
                for bs in (500000,1000000,2000000):
                    r=summarize([x for x in records if x["filter"]==filt],d,region,cls,bs); r["filter"]=filt; summary.append(r)
        # all-PNAS outside and genome-wide descriptive pool
        for cls in ("inside","outside"):
            sub=[r for r in records if r["filter"]==filt and r["design"]==d and ((r["class"]=="inside") if cls=="inside" else (r["class"]=="outside"))]
            for bs in (500000,1000000,2000000):
                ab=sum(int(r["ABBA"]) for r in sub); ba=sum(int(r["BABA"]) for r in sub)
                dd,se,nb=jackknife([(int(r["start"])//bs,int(r["ABBA"]),int(r["BABA"])) for r in sub],bs)
                summary.append({"filter":filt,"design":d,"region":"all_PNAS","region_class":cls,"block_size":bs,"n_loci":len(sub),"informative_sites":ab+ba,"ABBA":ab,"BABA":ba,"D":"" if dd is None else dd,"jackknife_se":"" if se is None else se,"blocks":nb})
    # Add the paired effect to every corresponding inside/outside summary row.
    for r in summary:
        diff,se=paired_difference(records,r["design"],r["region"],int(r["block_size"]),r["filter"])
        r["inside_minus_outside_D"]="" if diff is None else diff; r["inside_minus_outside_se"]="" if se is None else se
    write_tsv(RESULTS/"patterson_d_summary.tsv",summary,list(summary[0]) if summary else [])
    windows=[]
    for d in [x["design"] for x in designs]:
      for filt in ("all", "transversions_only"):
        for w in range(0, max(int(r["end"]) for r in rows)+1000000, 1000000):
            sub=[r for r in records if r["filter"]==filt and r["design"]==d and w <= int(r["start"]) < w+1000000]
            if sub:
                a=sum(int(r["ABBA"]) for r in sub); b=sum(int(r["BABA"]) for r in sub)
                windows.append({"filter":filt,"design":d,"start":w,"end":w+1000000,"midpoint":w+500000,"n_loci":len(sub),"ABBA":a,"BABA":b,"D":"" if not a+b else patterson_d(a,b)})
    write_tsv(RESULTS/"patterson_d_1Mb_windows.tsv",windows,list(windows[0]) if windows else [])
    svg_fig(windows,summary)
    audit={"chr4_manifest_loci":len(rows),"records":len(records),"coverage":dict(coverage),"missing_taxon_design":dict(missing),"archive":str(ARCHIVE),"intervals":INTERVALS}
    (RESULTS/"analysis_audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    report=RESULTS/"patterson_d_report.md"
    report.write_text("# Neoaves chromosome-4 Patterson's D\n\n"+json.dumps(audit,indent=2)+"\n\nD is pooled as (ABBA-BABA)/(ABBA+BABA); standard errors are genomic block jackknife estimates. See `patterson_d_summary.tsv` and the frozen quartet configuration for the complete results and limitations.\n")
    print(json.dumps(audit,indent=2))

if __name__ == "__main__": main()
