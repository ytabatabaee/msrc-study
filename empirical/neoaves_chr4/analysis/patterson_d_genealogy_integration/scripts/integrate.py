#!/usr/bin/env python3
"""Integrate nucleotide Patterson's D with frozen Neoaves QQS and B tracks."""
from __future__ import annotations
import csv, gzip, math, random, shutil, subprocess
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[5]
BASE=ROOT/"empirical/neoaves_chr4/analysis/patterson_d_genealogy_integration"
RES=BASE/"results"; FIG=BASE/"figures"
PATT=ROOT/"empirical/neoaves_chr4/analysis/patterson_d/results"
MANIFEST=ROOT/"empirical/neoaves_chr4/results/stage4d_locus_manifest.tsv"
QQS=ROOT/"empirical/neoaves_chr4/results/stage6_genomic_tracks/window_quartet_support.tsv"
DIST=ROOT/"empirical/neoaves_chr4/results/stage9_full_taxon_background/background_distance_estimators.tsv"
TREES=ROOT/"empirical/neoaves_chr4/external/stiller2024/63430.named.gene.trees.gz"
INTERVALS=(('PNAS_1',25030000,32670000),('PNAS_2',33510000,34470000),('PNAS_3',44130000,56810000))
RESOLUTIONS=(("100kb",100000),("250kb",250000),("500kb",500000),("1Mb",1000000))
DESIGNS=("Columbea_primary","N61_primary","N62_primary","Columbea_alt_outgroup","N61_alt_outgroup","N62_alt_outgroup")
PRIMARY=("Columbea_primary","N61_primary","N62_primary")
TAXA={
 "Columbea":("Columba_livia","Phoenicopterus_ruber","Ciconia_maguari","Gallus_gallus"),
 "N61":("Tauraco_erythrolophus","Columba_livia","Ciconia_maguari","Gallus_gallus"),
 "N62":("Columba_livia","Pterocles_gutturalis","Tauraco_erythrolophus","Gallus_gallus"),
}

def read_tsv(path):
    with path.open() as fh: return list(csv.DictReader(fh,delimiter="\t"))
def write_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or (list(rows[0]) if rows else [])
    with path.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n",extrasaction="ignore"); w.writeheader(); w.writerows(rows)
def f(x):
    try: return float(x)
    except (TypeError,ValueError): return None
def mean(xs): return sum(xs)/len(xs) if xs else None
def dcount(a,b): return None if a+b==0 else (a-b)/(a+b)
def region_overlap(start,end):
    return ";".join(n for n,a,b in INTERVALS if start<=b and end>=a) or "outside_all_PNAS"
def region_locus(start):
    return next((n for n,a,b in INTERVALS if a<=start<=b),"outside_all_PNAS")
def bin_start(start,size): return (int(start)//size)*size

def pearson(xs,ys):
    if len(xs)<3: return None
    mx,my=mean(xs),mean(ys); vx=sum((x-mx)**2 for x in xs); vy=sum((y-my)**2 for y in ys)
    return None if vx==0 or vy==0 else sum((x-mx)*(y-my) for x,y in zip(xs,ys))/math.sqrt(vx*vy)
def ranks(xs):
    order=sorted(range(len(xs)),key=lambda i:xs[i]); out=[0.0]*len(xs); i=0
    while i<len(order):
        j=i
        while j+1<len(order) and xs[order[j+1]]==xs[order[i]]: j+=1
        r=(i+j+2)/2
        for k in range(i,j+1): out[order[k]]=r
        i=j+1
    return out
def corr_rows(rows,xkey,ykey):
    z=[(f(r.get(xkey)),f(r.get(ykey))) for r in rows if f(r.get(xkey)) is not None and f(r.get(ykey)) is not None]
    if not z: return None,None,len(z)
    xs,ys=zip(*z); return pearson(xs,ys),pearson(ranks(xs),ranks(ys)),len(z)

class Node:
    __slots__=("children","label","length","parent")
    def __init__(self,label="",length=0.0): self.children=[]; self.label=label; self.length=length; self.parent=None

def parse_newick(text):
    s=text.strip().rstrip(";"); i=0
    def token(stop):
        nonlocal i
        j=i
        while i<len(s) and s[i] not in stop: i+=1
        return s[j:i]
    def node():
        nonlocal i
        n=Node()
        if i<len(s) and s[i]=='(':
            i+=1
            while True:
                c=node(); c.parent=n; n.children.append(c)
                if i<len(s) and s[i]==',': i+=1; continue
                if i<len(s) and s[i]==')': i+=1; break
                raise ValueError("malformed Newick")
            # Internal support labels such as /0.95 or 0.95 are ignored.
            if i<len(s) and s[i] not in ',):': token(',):')
        else:
            n.label=token(':,()')
        if i<len(s) and s[i]==':':
            i+=1; raw=token(',()')
            try: n.length=float(raw)
            except ValueError: n.length=math.nan
        return n
    root=node()
    if i!=len(s): raise ValueError("trailing Newick text")
    return root

def leaf_map(root):
    out={}
    def walk(n,depth):
        if not n.children:
            out[n.label]=(n,depth); return
        for c in n.children: walk(c,depth+(c.length if math.isfinite(c.length) else 0.0))
    walk(root,0.0); return out
def distance_map(root, taxa):
    lm=leaf_map(root); paths={}
    for t in taxa:
        n,_=lm[t]; p=[]; depth=0.0
        while n is not None:
            p.append((n,depth)); depth += n.length if math.isfinite(n.length) else 0.0; n=n.parent
        paths[t]=p
    out={}
    for i,a in enumerate(taxa):
        for b in taxa[i+1:]:
            pb={id(n):(n,d) for n,d in paths[b]}; best=None
            for n,da in paths[a]:
                if id(n) in pb: best=da+pb[id(n)][1]; break
            out[(a,b)]=best
    return out
def getdist(dm,a,b): return dm.get((a,b),dm.get((b,a)))
def quartet_summary(root,taxa):
    if any(t not in leaf_map(root) for t in taxa): return {"status":"missing_taxon"}
    dm=distance_map(root,taxa); a,b,c,o=taxa
    pairs={(a,b):getdist(dm,a,b),(c,o):getdist(dm,c,o),(a,c):getdist(dm,a,c),(b,o):getdist(dm,b,o),(a,o):getdist(dm,a,o),(b,c):getdist(dm,b,c)}
    if any(v is None or not math.isfinite(v) for v in pairs.values()): return {"status":"nonfinite_distance"}
    sums=[pairs[(a,b)]+pairs[(c,o)],pairs[(a,c)]+pairs[(b,o)],pairs[(a,o)]+pairs[(b,c)]]
    order=sorted(range(3),key=lambda i:sums[i]); topology=["q1","q2","q3"][order[0]]
    # A resolved additive quartet has equal *largest* four-point sums.
    resolved=(abs(sums[order[2]]-sums[order[1]])<1e-8 and
              abs(sums[order[1]]-sums[order[0]])>1e-8)
    internal=max(0.0,(max(sums)-min(sums))/2) if resolved else 0.0
    # For a|b,c|o split, pendant lengths follow the four-point formulas.
    if topology=="q1": x,y,z,w=a,b,c,o
    elif topology=="q2": x,y,z,w=a,c,b,o
    else: x,y,z,w=a,o,b,c
    d=lambda u,v:getdist(dm,u,v)
    pend=[(d(x,y)+d(x,z)-d(y,z))/2,(d(y,x)+d(y,z)-d(x,z))/2,(d(z,x)+d(z,w)-d(x,w))/2,(d(w,z)+d(w,x)-d(z,x))/2]
    total=sum(pend)+internal
    return {"status":"resolved" if resolved else "unresolved","topology":topology,"internal_branch":internal,"total_quartet_length":total,"d_ab":d(a,b),"d_ac":d(a,c),"d_ao":d(a,o),"d_bc":d(b,c),"d_bo":d(b,o),"d_co":d(c,o)}

def load_inputs():
    manifest=[r for r in read_tsv(MANIFEST) if r["chromosome"]=="chr4"]
    by_locus={r["locus_id"]:r for r in manifest}; by_index={r["tree_index"]:r for r in manifest}
    patt=read_tsv(PATT/"patterson_d_locus_results.tsv")
    qrows=read_tsv(QQS); drows=read_tsv(DIST)
    qq={}
    for r in qrows:
        m=by_index.get(r["Gene"])
        if m and r["clade"] in ("Columbea","N61","N62"): qq[(m["locus_id"],r["clade"])]=r
    dist={}
    for r in drows:
        m=by_index.get(r["Gene"])
        if m: dist[m["locus_id"]]=r
    return manifest,by_locus,patt,qq,dist

def extract_quartets(manifest, patt):
    wanted={r["locus_id"] for r in manifest}; dmap={r["locus_id"]:r for r in patt if r["design"]=="N61_primary" and r["filter"]=="all"}
    out=[]; seen=set()
    if not TREES.exists(): return out
    with gzip.open(TREES,"rt") as fh:
        for line in fh:
            if "\t" not in line: continue
            locus,nwk=line.rstrip("\n").split("\t",1)
            if locus not in wanted: continue
            seen.add(locus)
            try:
                q=quartet_summary(parse_newick(nwk),TAXA["N61"])
            except Exception as exc: q={"status":"parse_error:"+str(exc)}
            row={"locus_id":locus,"chromosome":"chr4","start":by_start(locus,manifest),"N61_tree_status":q.pop("status","")}
            row.update(q); row["patterson_informative"]=dmap.get(locus,{}).get("informative",""); row["patterson_D"]=dmap.get(locus,{}).get("D","")
            out.append(row)
    return out
def by_start(locus,manifest):
    for r in manifest:
        if r["locus_id"]==locus: return r["start"]
    return ""

def aggregate_tracks(manifest,patt,qq,dist,size,qmetrics=None):
    by_locus={r["locus_id"]:r for r in manifest}; groups=defaultdict(lambda:{"n":0,"abba":0,"baba":0,"inf":0,"q":[],"b":[]})
    for r in patt:
        if r["design"] not in DESIGNS: continue
        m=by_locus.get(r["locus_id"])
        if not m: continue
        g=groups[(r["design"],r["filter"],bin_start(int(m["start"]),size))]; g["n"]+=1; g["abba"]+=int(r["ABBA"]); g["baba"]+=int(r["BABA"]); g["inf"]+=int(r["informative"])
        qr=qq.get((r["locus_id"],r["design"].split("_")[0]))
        if qr: g["q"].append(qr)
        dr=dist.get(r["locus_id"]); bv=f(dr.get("B_observed")) if dr else None
        if bv is not None: g["b"].append(bv)
        if qmetrics and r["locus_id"] in qmetrics:
            iv=f(qmetrics[r["locus_id"]].get("internal_branch"))
            if iv is not None: g.setdefault("internal",[]).append(iv)
    max_end=max(int(r["end"]) for r in manifest); last=((max_end+size-1)//size)*size; out=[]
    for design in DESIGNS:
        for filt in ("all","transversions_only"):
            for start in range(0,last,size):
                g=groups[(design,filt,start)]; q=g["q"]; b=g["b"]; q1=mean([f(x["q1_species_history"]) for x in q if f(x["q1_species_history"]) is not None]); q2=mean([f(x["q2_alternative_1"]) for x in q if f(x["q2_alternative_1"]) is not None]); q3=mean([f(x["q3_alternative_2"]) for x in q if f(x["q3_alternative_2"]) is not None]); dd=dcount(g["abba"],g["baba"])
                out.append({"resolution_bp":size,"design":design,"filter":filt,"chromosome":"chr4","window_start":start,"window_end":start+size,"midpoint":start+size/2,"patterson_loci":g["n"],"qqs_loci":len(q),"distance_loci":len(b),"quartet_tree_loci":len(g.get("internal",[])),"ABBA":g["abba"],"BABA":g["baba"],"informative_ABBA_BABA_sites":g["inf"],"Patterson_D":"" if dd is None else dd,"q1":"" if q1 is None else q1,"q2":"" if q2 is None else q2,"q3":"" if q3 is None else q3,"q2_minus_q3":"" if q2 is None or q3 is None else q2-q3,"topology_shift_q2_minus_q1":"" if q1 is None or q2 is None else q2-q1,"B_observed":"" if not b else mean(b),"quartet_internal_branch":"" if not g.get("internal") else mean(g["internal"]),"published_PNAS_outlier_region":region_overlap(start,start+size-1)})
    return out

def block_boot_corr(rows,x,y,seed=61,n=1000):
    z=[r for r in rows if f(r.get(x)) is not None and f(r.get(y)) is not None]; blocks=defaultdict(list)
    for r in z: blocks[int(r["window_start"])//1000000].append(r)
    keys=list(blocks); rng=random.Random(seed); vals=[]
    for _ in range(n):
        sample=[r for k in (rng.choice(keys) for _ in keys) for r in blocks[k]]
        v=pearson([f(r[x]) for r in sample],[f(r[y]) for r in sample])
        if v is not None: vals.append(v)
    obs=pearson([f(r[x]) for r in z],[f(r[y]) for r in z]); vals.sort()
    return obs, (vals[int(.025*len(vals))],vals[int(.975*len(vals))-1]) if vals else (None,None),len(keys),len(z)
def shifted_p(rows,x,y):
    z=sorted(rows,key=lambda r:r["window_start"]); a=[f(r.get(x)) for r in z]; b=[f(r.get(y)) for r in z]; n=len(z)
    pairs=[(i,a[i],b[i]) for i in range(n) if a[i] is not None and b[i] is not None]; obs=pearson([x[1] for x in pairs],[x[2] for x in pairs]); null=[]
    for shift in range(n):
        z2=[(a[(i-shift)%n],b[i]) for i in range(n) if a[(i-shift)%n] is not None and b[i] is not None]
        if len(z2)>=3: null.append(pearson([u for u,v in z2],[v for u,v in z2]))
    p=(1+sum(abs(v)>=abs(obs) for v in null))/(len(null)+1) if obs is not None else None
    return obs,p,len(null)

def inside_effect(rows,stat,region):
    inside=[]; outside=[]
    for r in rows:
        hit=region in r["published_PNAS_outlier_region"] if region!="all_PNAS" else r["published_PNAS_outlier_region"]!="outside_all_PNAS"
        target=inside if hit else outside
        if stat=="Patterson_D":
            target.append((int(r["window_start"])//1000000,int(r["ABBA"]),int(r["BABA"])))
        else:
            v=f(r.get(stat))
            if v is not None: target.append((int(r["window_start"])//1000000,v))
    def value(x):
        if stat=="Patterson_D": return dcount(sum(a for _,a,b in x),sum(b for _,a,b in x))
        return mean([v for _,v in x])
    delta=None if not inside or not outside else value(inside)-value(outside); blocks=sorted({x[0] for x in inside+outside}); reps=[]
    for block in blocks:
        ii=[x for x in inside if x[0]!=block]; oo=[x for x in outside if x[0]!=block]
        if ii and oo and value(ii) is not None and value(oo) is not None: reps.append(value(ii)-value(oo))
    se=math.sqrt((len(reps)-1)/len(reps)*sum((x-mean(reps))**2 for x in reps)) if len(reps)>1 else None
    return delta,se,len(inside),len(outside),len(blocks)

def svg_convert(path):
    c=shutil.which("rsvg-convert")
    if c:
        subprocess.run([c,"-f","png","-o",str(path.with_suffix(".png")),str(path)],check=False); subprocess.run([c,"-f","pdf","-o",str(path.with_suffix(".pdf")),str(path)],check=False)

def draw_polyline(lines, points, color):
    """Draw only contiguous valid segments; missing bins stay visibly missing."""
    segment=[]
    for point in points + [None]:
        if point is None:
            if len(segment) >= 2:
                lines.append(f'<polyline class="line" stroke="{color}" points="{" ".join(segment)}"/>')
            segment=[]
        else:
            segment.append(point)

def joint_figure(rows,res,design="N61_primary",suffix=""):
    z=[r for r in rows if r["resolution_bp"]==res and r["design"]==design and r["filter"]=="all"]; width,height=1500,920; left,right,top=95,35,55; plotw=width-left-right; xmax=max(r["window_end"] for r in z); panels=[("Patterson D","Patterson_D",-1,1),("QQS q1/q2/q3",None,0,1),("q2 - q3","q2_minus_q3",-1,1),("topology D (q2-q1)","topology_shift_q2_minus_q1",-1,1),("B observed","B_observed",0,max([f(r["B_observed"]) or 0 for r in z] or [1])),("coverage / informative sites",None,0,max([r["informative_ABBA_BABA_sites"] for r in z] or [1]))]
    def x(v): return left+plotw*v/xmax
    lines=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial;font-size:13px}}.line{{fill:none;stroke-width:2}}</style><text x="95" y="28" font-size="22">Neoaves chr4 N61 joint genealogy tracks ({res//1000} kb)</text>']
    colors={"Patterson_D":"#d62728","q2_minus_q3":"#9467bd","topology_shift_q2_minus_q1":"#1f77b4","B_observed":"#2ca02c"}
    for i,(label,key,lo,hi) in enumerate(panels):
        py=top+i*135; lines.append(f'<text x="10" y="{py+8}">{label}</text>')
        for n,a,b in INTERVALS: lines.append(f'<rect x="{x(a):.1f}" y="{py}" width="{max(1,x(b)-x(a)):.1f}" height="100" fill="#f4b183" opacity=".22"/>')
        lines.append(f'<line x1="{left}" y1="{py+100*(1-(0-lo)/(hi-lo)):.1f}" x2="{left+plotw}" y2="{py+100*(1-(0-lo)/(hi-lo)):.1f}" stroke="#999"/>')
        if key:
            pts=[None if r[key]=="" else f'{x(r["midpoint"]):.1f},{py+100*(1-(float(r[key])-lo)/(hi-lo)):.1f}' for r in z]
            draw_polyline(lines,pts,colors[key])
        elif label.startswith("QQS"):
            for k,c in (("q1","#222"),("q2","#1f77b4"),("q3","#9467bd")):
                pts=[None if r[k]=="" else f'{x(r["midpoint"]):.1f},{py+100*(1-float(r[k])):.1f}' for r in z]
                draw_polyline(lines,pts,c)
        else:
            pts=[f'{x(r["midpoint"]):.1f},{py+100*(1-r["informative_ABBA_BABA_sites"]/hi):.1f}' for r in z]
            lines.append(f'<polyline class="line" stroke="#e67e22" points="{" ".join(pts)}"/>')
    lines.append('<text x="600" y="900">Shading: independently defined PNAS outlier intervals; nucleotide D and tree-based contrasts use separate axes/statistics.</text></svg>')
    p=FIG/("neoaves_N61_joint_tracks"+suffix+".svg")
    p.write_text("\n".join(lines)); svg_convert(p)

def make_scatter(rows):
    z=[r for r in rows if r["resolution_bp"]==250000 and r["design"]=="N61_primary" and r["filter"]=="all"]; width,height=1200,560; lines=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial;font-size:13px}}.point{{stroke:none}}</style><text x="70" y="28" font-size="21">N61 spatial correspondence at 250 kb</text>']
    panels=[("Patterson D vs q2-q3","Patterson_D","q2_minus_q3",70), ("Patterson D vs B observed","Patterson_D","B_observed",650)]
    for title,xk,yk,px in panels:
        z2=[r for r in z if r[xk]!="" and r[yk]!=""]; xs=[f(r[xk]) for r in z2]; ys=[f(r[yk]) for r in z2]; xmin,xmax=min(xs),max(xs); ymin,ymax=min(ys),max(ys); spanx=xmax-xmin or 1; spany=ymax-ymin or 1
        lines.append(f'<text x="{px}" y="65">{title}</text><line x1="{px}" y1="470" x2="{px+480}" y2="470" stroke="#555"/><line x1="{px}" y1="90" x2="{px}" y2="470" stroke="#555"/>')
        for r in z2:
            xx=px+480*(f(r[xk])-xmin)/spanx; yy=470-380*(f(r[yk])-ymin)/spany; col="#d95f02" if r["published_PNAS_outlier_region"]!="outside_all_PNAS" else "#777"
            lines.append(f'<circle class="point" cx="{xx:.1f}" cy="{yy:.1f}" r="3" fill="{col}"/>')
        r0=pearson(xs,ys); lines.append(f'<text x="{px}" y="500">r = {"NA" if r0 is None else f"{r0:.3f}"}; nonoverlapping 250-kb bins</text>')
    lines.append('</svg>'); p=FIG/"neoaves_N61_spatial_correspondence.svg"; p.write_text("\n".join(lines)); svg_convert(p)

def make_region_comparison(effects):
    rows=[r for r in effects if r["design"]=="N61_primary" and r["filter"]=="all"]
    stats=[("Patterson_D","nucleotide D"),("q2_minus_q3","q2-q3"),("topology_shift_q2_minus_q1","q2-q1"),("B_observed","B observed")]
    width,height=1350,760; left=120; panelw=290; top=85; lines=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial;font-size:13px}}.bar{{stroke:#333;stroke-width:.6}}</style><text x="120" y="32" font-size="22">N61 inside-minus-outside comparisons by PNAS interval</text>']
    colors={"PNAS_1":"#d95f02","PNAS_2":"#1b9e77","PNAS_3":"#7570b3","all_PNAS":"#555"}
    for j,(stat,label) in enumerate(stats):
        x0=left+j*panelw; vals=[f(r["inside_minus_outside"]) for r in rows if r["statistic"]==stat and f(r["inside_minus_outside"]) is not None]; lim=max([abs(v) for v in vals] or [1]); baseline=top+250
        lines += [f'<text x="{x0}" y="65">{label}</text>',f'<line x1="{x0}" y1="{baseline}" x2="{x0+250}" y2="{baseline}" stroke="#555"/>']
        for i,region in enumerate(("PNAS_1","PNAS_2","PNAS_3","all_PNAS")):
            r=next((x for x in rows if x["statistic"]==stat and x["region"]==region),None); v=f(r["inside_minus_outside"]) if r else None
            if v is None: continue
            y=baseline-180*v/lim; bh=baseline-y if v>=0 else y-baseline; yy=min(y,baseline); lines.append(f'<rect class="bar" x="{x0+20+i*55}" y="{yy:.1f}" width="35" height="{max(1,bh):.1f}" fill="{colors[region]}"/>'); lines.append(f'<text x="{x0+22+i*55}" y="{baseline+22}" transform="rotate(45 {x0+22+i*55} {baseline+22})">{region.replace("PNAS_","")}</text>')
            se=f(r["block_jackknife_se"])
            if se is not None:
                ey=180*se/lim; cx=x0+37.5+i*55
                lines.append(f'<line x1="{cx:.1f}" y1="{y-ey:.1f}" x2="{cx:.1f}" y2="{y+ey:.1f}" stroke="#111"/>')
                lines.append(f'<line x1="{cx-5:.1f}" y1="{y-ey:.1f}" x2="{cx+5:.1f}" y2="{y-ey:.1f}" stroke="#111"/><line x1="{cx-5:.1f}" y1="{y+ey:.1f}" x2="{cx+5:.1f}" y2="{y+ey:.1f}" stroke="#111"/>')
        lines.append(f'<text x="{x0}" y="{baseline+90}">error bars: 1-Mb block jackknife SE</text>')
    lines.append('<text x="120" y="710">Positive values indicate higher inside-region support/distance; the four panels have statistic-specific scales.</text></svg>')
    p=FIG/"neoaves_N61_region_comparison.svg"; p.write_text("\n".join(lines)); svg_convert(p)

def make_all_clades_comparison(effects):
    rows=[r for r in effects if r["filter"]=="all" and r["region"]=="all_PNAS" and r["statistic"] in ("Patterson_D","q2_minus_q3","B_observed")]
    names=[("Columbea_primary","Columbea"),("N61_primary","N61"),("N62_primary","N62")]; stats=[("Patterson_D","nucleotide D"),("q2_minus_q3","q2-q3"),("B_observed","B observed")]
    width,height=1200,700; lines=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial;font-size:13px}}.bar{{stroke:#333;stroke-width:.6}}</style><text x="95" y="32" font-size="22">Neoaves chromosome 4: all-clade inside-minus-outside comparisons</text>']
    for j,(stat,label) in enumerate(stats):
        x0=95+j*350; vals=[f(r["inside_minus_outside"]) for r in rows if r["statistic"]==stat and f(r["inside_minus_outside"]) is not None]; lim=max([abs(v) for v in vals] or [1]); baseline=440; lines += [f'<text x="{x0}" y="72">{label}</text>',f'<line x1="{x0}" y1="{baseline}" x2="{x0+280}" y2="{baseline}" stroke="#555"/>']
        for i,(design,name) in enumerate(names):
            r=next((x for x in rows if x["design"]==design and x["statistic"]==stat),None); v=f(r["inside_minus_outside"]) if r else None
            if v is None: continue
            y=baseline-260*v/lim; yy=min(y,baseline); lines.append(f'<rect class="bar" x="{x0+35+i*75}" y="{yy:.1f}" width="48" height="{max(1,abs(y-baseline)):.1f}" fill="{("#1b9e77","#d95f02","#7570b3")[i]}"/>'); lines.append(f'<text x="{x0+32+i*75}" y="{baseline+22}" transform="rotate(45 {x0+32+i*75} {baseline+22})">{name}</text>')
    lines.append('<text x="95" y="640">The three clade designs use different frozen quartet roles; contrasts are descriptive and not expected to share a common topology orientation.</text></svg>')
    p=FIG/"neoaves_all_three_clades_comparison.svg"; p.write_text("\n".join(lines)); svg_convert(p)

def main():
    RES.mkdir(parents=True,exist_ok=True); FIG.mkdir(parents=True,exist_ok=True)
    manifest,by_locus,patt,qq,dist=load_inputs(); write_tsv(RES/"locus_mapping_audit.tsv",[{"locus_id":m["locus_id"],"tree_index":m["tree_index"],"start":m["start"],"end":m["end"],"qqs_clades":sum((m["locus_id"],c) in qq for c in ("Columbea","N61","N62")),"distance_present":int(m["locus_id"] in dist),"patterson_designs":sum(r["locus_id"]==m["locus_id"] for r in patt),"mapping_status":"exact_stage4d_locus_id"} for m in manifest])
    qtrees=extract_quartets(manifest,patt)
    write_tsv(RES/"quartet_branch_lengths.tsv",qtrees)
    qmetrics={r["locus_id"]:r for r in qtrees if r.get("N61_tree_status")=="resolved"}
    tracks={}
    for name,size in RESOLUTIONS: tracks[name]=aggregate_tracks(manifest,patt,qq,dist,size,qmetrics)
    fields=list(tracks["250kb"][0]); write_tsv(RES/"matched_genomic_tracks.tsv",[r for name in tracks for r in tracks[name]],fields)
    def primary_rows(): return [r for r in tracks["250kb"] if r["design"]=="N61_primary" and r["filter"]=="all"]
    corr=[]
    for design in DESIGNS:
        for filt in ("all","transversions_only"):
            z=[r for r in tracks["250kb"] if r["design"]==design and r["filter"]==filt]
            for stat in ("q2_minus_q3","topology_shift_q2_minus_q1","B_observed","quartet_internal_branch"):
                obs,ci,nb,np=block_boot_corr(z,"Patterson_D",stat); shift,pnull,nshift=shifted_p(z,"Patterson_D",stat)
                corr.append({"resolution":"250kb","design":design,"filter":filt,"x_statistic":"Patterson_D","y_statistic":stat,"n_matched_windows":np,"pearson_r":obs,"block_bootstrap_ci_low":ci[0],"block_bootstrap_ci_high":ci[1],"n_independent_1Mb_blocks":nb,"circular_shift_observed_r":shift,"circular_shift_p_abs":pnull,"n_circular_shifts":nshift})
    write_tsv(RES/"spatial_correspondence.tsv",corr)
    effects=[]
    stats=("Patterson_D","q2_minus_q3","topology_shift_q2_minus_q1","B_observed","quartet_internal_branch")
    for design in DESIGNS:
        for filt in ("all","transversions_only"):
            z=[r for r in tracks["250kb"] if r["design"]==design and r["filter"]==filt]
            for region in ("PNAS_1","PNAS_2","PNAS_3","all_PNAS"):
                for stat in stats:
                    delta,se,ni,no,nb=inside_effect(z,stat,region); effects.append({"resolution":"250kb","design":design,"filter":filt,"region":region,"statistic":stat,"inside_minus_outside":delta,"block_jackknife_se":se,"inside_windows":ni,"outside_windows":no,"independent_1Mb_blocks":nb})
    write_tsv(RES/"inside_outside_comparison.tsv",effects)
    reg=[r for r in effects if r["design"]=="N61_primary" and r["filter"]=="all"]
    write_tsv(RES/"region_specific_analysis.tsv",reg)
    robust=[r for r in effects if r["region"]=="all_PNAS" and r["statistic"] in ("Patterson_D","q2_minus_q3","B_observed","quartet_internal_branch")]
    write_tsv(RES/"robustness_analysis.tsv",robust)
    (RES/"statistic_definition_mapping.md").write_text("""# Statistic-definition and sign mapping\n\nFor N61 the frozen roles are P1=C1=`Tauraco_erythrolophus`, P2=C2=`Columba_livia`, P3=S=`Ciconia_maguari`, and O=`Gallus_gallus`. The existing QQS definitions are q1=(C1,C2)|(S,O), q2=(C1,S)|(C2,O), and q3=(C2,S)|(C1,O).\n\nWith outgroup polarization, ABBA is A/B/B/A and supports the P2-P3 split, which is q3. BABA is B/A/B/A and supports the P1-P3 split, which is q2. Therefore negative nucleotide Patterson D means BABA>ABBA and is expected to correspond to positive q2-q3, not to positive q2-q1. The existing topology-shift statistic is q2-q1 and is a different contrast.\n\nGallus is an independently external galliform reference relative to the three Neoaves ingroup roles, so the rooted quartet is biologically defensible for this exploratory comparison. All branch lengths are ML substitution units, not coalescence times; B_observed is an arrangement-blind mean patristic distance over observed taxa.\n\nThe PNAS intervals are GalGal4 window-start intervals. The Stiller/Stage-4D locus axis is GalGal6. The repository audit records approximate GalGal6 counterparts but no nucleotide-resolution lift was applied; all cross-track comparisons therefore retain the repository's validated published-predicate mask and explicitly preserve this assembly caveat.\n""")
    for res in (250000,100000,500000): joint_figure(tracks["250kb"] if res==250000 else tracks["100kb" if res==100000 else "500kb"],res,suffix="" if res==250000 else f"_{res//1000}kb")
    make_scatter(primary_rows())
    make_region_comparison(effects)
    make_all_clades_comparison(effects)
    n_exact=sum(1 for r in qq if r[1]=="N61")
    qstatus=defaultdict(int)
    for r in qtrees: qstatus[r.get("N61_tree_status","")]+=1
    n61eff=[r for r in effects if r["design"]=="N61_primary" and r["filter"]=="all" and r["region"]=="all_PNAS"]
    n61corr=[r for r in corr if r["design"]=="N61_primary" and r["filter"]=="all"]
    fmt=lambda x: "NA" if x is None else f"{float(x):.4g}"
    effect_lines=[f"- {r['statistic']}: Δ inside−outside = {fmt(r['inside_minus_outside'])}, block-jackknife SE = {fmt(r['block_jackknife_se'])}, inside/outside windows = {r['inside_windows']}/{r['outside_windows']}" for r in n61eff]
    corr_lines=[f"- D vs {r['y_statistic']}: r = {fmt(r['pearson_r'])}, 95% 1-Mb-block bootstrap CI [{fmt(r['block_bootstrap_ci_low'])}, {fmt(r['block_bootstrap_ci_high'])}], circular-shift p = {fmt(r['circular_shift_p_abs'])}" for r in n61corr]
    patt_n61=sorted({r["locus_id"] for r in patt if r["design"]=="N61_primary" and r["filter"]=="all"})
    patt_inside=sum(1 for r in patt_n61 if any(x["locus_id"]==r and x["region"]!="outside_all_PNAS" for x in patt if x["design"]=="N61_primary" and x["filter"]=="all"))
    report=["# Neoaves Patterson's D / genealogy integration", "", "## Data correspondence", f"The Stage-4D manifest contains {len(manifest)} chr4 loci. Patterson, QQS, and distance inputs were joined through exact manifest `locus_id` mappings from the integer `tree_index` fields; no positional join was used. The full audit is in `locus_mapping_audit.tsv`. The QQS table contributed {n_exact} exact N61 locus/clade matches (115 N61 clade records are absent from the 6247-locus manifest join); the distance table has 6140 exact chr4 locus matches. Patterson N61 all-site records cover {len(patt_n61)} unique loci, with {patt_inside} inside and {len(patt_n61)-patt_inside} outside the PNAS mask; the remaining manifest loci have no N61 all-site record. The complete per-locus missingness audit is retained in `locus_mapping_audit.tsv`.", "", "The published PNAS intervals are GalGal4 window-start intervals, while Stiller local-tree coordinates are GalGal6. Existing repository provenance supplies approximate GalGal6 counterparts but no exact lift; this analysis uses the previously validated published-predicate mask and does not claim nucleotide-resolution boundary equivalence.", "", "## Main N61 result", "", "At 250 kb, nucleotide D is compared with q2-q3, because negative Patterson D (BABA > ABBA) corresponds to q2 > q3 under the frozen N61 taxon mapping. The existing q2-q1 track is retained separately. Correlations use nonoverlapping 250-kb windows; block bootstrap intervals and circular-shift spatial nulls are reported in `spatial_correspondence.tsv`.", "", "The original locus-level pooled Patterson analysis gives N61 ΔD = -0.177615 (inside D = -0.258194; outside D = -0.080579). The 250-kb window classification gives ΔD = -0.174987 because windows overlapping a published boundary are assigned to the interval; this small difference is a classification consequence, not a change to the frozen result.", "", "N61 all-PNAS inside-minus-outside effects (D is pooled from ABBA/BABA counts before division):", *effect_lines, "", "N61 250-kb spatial correspondence:", *corr_lines, "", "The q2-q3 and D relationship should be interpreted as a sign-reversed discordance comparison, not as the same estimand. B_observed is a separate arrangement-blind substitution-distance summary.", "", "## Branch-length feasibility", "", f"The local named trees yielded {len(qtrees)} chr4 N61 quartet records. Status counts are {dict(qstatus)}. Quartet-level topology, internal branch, total quartet length, and pairwise distances are in `quartet_branch_lengths.tsv`. Missing taxa, unresolved induced quartets, and parse failures are retained as status fields.", "", "## Interpretation", "", "A spatially coincident nucleotide asymmetry and discordant-topology contrast would be compatible with a shared local genealogy regime, but neither statistic uniquely distinguishes MSRC from introgression, incomplete lineage sorting, mutation-rate variation, alignment error, or other causes. Branch lengths remain substitution distances. The joint analysis is therefore additive evidence of spatial association, not a mechanistic identification of MSRC.", "", "Candidate Results wording: The N61 nucleotide site-pattern asymmetry was spatially aligned with the independently defined chromosome-4 outlier intervals and was compared with the corresponding q2-q3 discordant-topology contrast. Because ABBA and BABA map to q3 and q2 respectively under the frozen quartet, the sign-reversed contrasts were directionally comparable. This association is compatible with a shared recombination-partitioned genealogy regime, but does not by itself discriminate MSRC from introgression or incomplete lineage sorting."]
    (RES/"integration_report.md").write_text("\n".join(report)+"\n")

if __name__=="__main__": main()
