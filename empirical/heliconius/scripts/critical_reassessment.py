#!/usr/bin/env python3
"""Critical, source-separated reassessment of Heliconius all-quartet results.

This stage works on the normalized/assigned local trees produced by
heliconius_pipeline.py.  It uses non-overlapping 50-kb physical bins for
region/run summaries, while retaining the original local-tree assignments for
window-level descriptive counts.  No bin or window is treated as an
independent replicate in the report.
"""
from __future__ import annotations
import argparse, csv, gzip, math, statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PROC = ROOT / "data/heliconius/processed"
EMP = ROOT / "empirical/heliconius"
RES = EMP / "results"
FIG = EMP / "figures"

def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def read_gz_tsv(path):
    with gzip.open(path, "rt", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def f(row, key):
    try: return float(row[key])
    except (KeyError, TypeError, ValueError): return math.nan

def median(values):
    values=[x for x in values if math.isfinite(x)]
    return statistics.median(values) if values else math.nan

def mad(values, center=None):
    values=[x for x in values if math.isfinite(x)]
    if not values: return math.nan
    c=median(values) if center is None else center
    return median([abs(x-c) for x in values])

def group(rows, keys):
    out=defaultdict(list)
    for r in rows: out[tuple(r[k] for k in keys)].append(r)
    return out

def write_tsv(path, rows, fields=None):
    if not rows: return
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields=[]
        for row in rows:
            for key in row:
                if key not in fields: fields.append(key)
    with path.open("w", newline="") as out:
        w=csv.DictWriter(out, fieldnames=fields, delimiter="\t")
        w.writeheader(); w.writerows(rows)

def load():
    assignments=read_gz_tsv(PROC/"all_quartet_assignments.tsv.gz")
    windows=read_tsv(PROC/"all_quartet_window_summary.tsv")
    return assignments, windows

def fixed_quartets(assignments):
    return sorted({r["quartet_id"] for r in assignments})

def binned_windows(windows):
    """Make descriptive non-overlapping bins; preserve source/assembly."""
    bins=defaultdict(list)
    for r in windows:
        try: start=int(float(r["start"])); bs=(start//50000)*50000
        except (KeyError,TypeError,ValueError): continue
        bins[(r["source_id"],r["assembly"],r["chromosome"],bs)].append(r)
    rows=[]
    for (source,assembly,chrom,bs),rs in sorted(bins.items()):
        vals={k:[f(r,k) for r in rs] for k in ("q1","q2","q3","D","mean_pairwise_patristic","mean_quartet_internal_branch","recomb_rate_cM_Mb")}
        row={"source_id":source,"assembly":assembly,"chromosome":chrom,"bin_start":bs,"bin_end":bs+49999,"n_trees":len(rs)}
        for k,v in vals.items(): row[k]=median(v) if k in ("mean_pairwise_patristic","mean_quartet_internal_branch") else (sum(x for x in v if math.isfinite(x))/sum(math.isfinite(x) for x in v) if any(math.isfinite(x) for x in v) else math.nan)
        row["alternative_support"]=row["q2"]+row["q3"] if math.isfinite(row["q2"]) and math.isfinite(row["q3"]) else math.nan
        row["alternative_minus_baseline"]=row["alternative_support"]-row["q1"] if math.isfinite(row["alternative_support"]) and math.isfinite(row["q1"]) else math.nan
        rows.append(row)
    for key,rs in group(rows,["source_id","chromosome"]).items():
        center=median([f(r,"alternative_minus_baseline") for r in rs])
        for r in rs:
            r["chromosome_median_effect"]=center
            r["effect_deviation"]=f(r,"alternative_minus_baseline")-center if math.isfinite(f(r,"alternative_minus_baseline")) and math.isfinite(center) else math.nan
            r["regime"]="alternative_shift" if f(r,"effect_deviation")>0.10 else ("baseline_shift" if f(r,"effect_deviation")<-0.10 else "stable")
    return rows

def contiguous_runs(rows, effect_key, threshold=0.10, min_bins=2):
    """Return signed sustained runs, requiring adjacent 50-kb bins."""
    rows=sorted(rows,key=lambda r:int(r["bin_start"]))
    runs=[]; current=[]; sign=None; prev=None
    for r in rows:
        x=f(r,effect_key)
        good=math.isfinite(x) and abs(x)>=threshold
        s=1 if x>=threshold else (-1 if x<=-threshold else 0)
        adjacent=prev is None or int(r["bin_start"])<=prev+50000
        if good and s and (not current or (s==sign and adjacent)):
            current.append(r); sign=s
        else:
            if len(current)>=min_bins: runs.append(current)
            current=[r] if good else []; sign=s if good else None
        prev=int(r["bin_start"])
    if len(current)>=min_bins: runs.append(current)
    out=[]
    for run in runs:
        effect=[f(r,effect_key) for r in run]
        out.append({"start":run[0]["bin_start"],"end":run[-1]["bin_end"],"span_bp":int(run[-1]["bin_end"])-int(run[0]["bin_start"])+1,"n_bins":len(run),"n_trees":sum(int(r["n_trees"]) for r in run),"mean_effect":sum(effect)/len(effect),"min_effect":min(effect),"max_effect":max(effect),"direction":"alternative" if sum(effect)/len(effect)>0 else "baseline","source_id":run[0]["source_id"],"assembly":run[0]["assembly"],"chromosome":run[0]["chromosome"]})
    return out

def dominance(assignments, windows):
    rows=[]
    for (source,chrom),rs in sorted(group(windows,["source_id","chromosome"]).items()):
        n=len(rs); q2=sum(f(r,"q2")>f(r,"q1") for r in rs); q3=sum(f(r,"q3")>f(r,"q1") for r in rs); alt=sum(f(r,"q2")+f(r,"q3")>f(r,"q1") for r in rs)
        bins=binned_windows(rs); runs=contiguous_runs(bins,"effect_deviation",0.10,2)
        longest=max(runs,key=lambda r:r["span_bp"],default={})
        rows.append({"scope":"aggregate","source_id":source,"chromosome":chrom,"n_windows":n,"fraction_q2_gt_q1":q2/n,"fraction_q3_gt_q1":q3/n,"fraction_alt_gt_q1":alt/n,"longest_run_start":longest.get("start",""),"longest_run_end":longest.get("end",""),"longest_run_span_bp":longest.get("span_bp",""),"n_runs":len(runs)})
    for (source,chrom,quartet),rs in sorted(group(assignments,["source_id","chromosome","quartet_id"]).items()):
        by=defaultdict(list)
        for r in rs: by[r["window_id"]].append(r)
        fixed=[]
        for wid,ar in by.items():
            r=ar[0]; cls=r["topology_class"]; fixed.append({"source_id":source,"assembly":r["assembly"],"chromosome":chrom,"bin_start":(int(float(r["start"]))//50000)*50000,"bin_end":(int(float(r["start"]))//50000)*50000+49999,"n_trees":1,"q1":float(cls=="q1"),"q2":float(cls=="q2"),"q3":float(cls=="q3")})
        n=len(fixed); q2=sum(r["q2"]>r["q1"] for r in fixed); q3=sum(r["q3"]>r["q1"] for r in fixed); alt=sum(r["q2"]+r["q3"]>r["q1"] for r in fixed)
        bins=defaultdict(list)
        for r in fixed: bins[r["bin_start"]].append(r)
        br=[]
        for bs,x in sorted(bins.items()):
            q1=sum(r["q1"] for r in x)/len(x); q2b=sum(r["q2"] for r in x)/len(x); q3b=sum(r["q3"] for r in x)/len(x)
            br.append({"source_id":source,"assembly":x[0]["assembly"],"chromosome":chrom,"bin_start":bs,"bin_end":bs+49999,"n_trees":len(x),"alternative_minus_baseline":q2b+q3b-q1})
        runs=contiguous_runs(br,"alternative_minus_baseline",0.10,2); longest=max(runs,key=lambda r:r["span_bp"],default={})
        rows.append({"scope":"quartet","source_id":source,"chromosome":chrom,"quartet_id":quartet,"n_windows":n,"fraction_q2_gt_q1":q2/n,"fraction_q3_gt_q1":q3/n,"fraction_alt_gt_q1":alt/n,"longest_run_start":longest.get("start",""),"longest_run_end":longest.get("end",""),"longest_run_span_bp":longest.get("span_bp",""),"n_runs":len(runs)})
    return rows

def quartet_bins(assignments):
    bins=defaultdict(list)
    for r in assignments:
        try: bs=(int(float(r["start"]))//50000)*50000
        except (TypeError,ValueError): continue
        bins[(r["source_id"],r["assembly"],r["chromosome"],r["quartet_id"],bs)].append(r)
    out=[]
    for (source,assembly,chrom,quartet,bs),rs in sorted(bins.items()):
        resolved=[r for r in rs if r["topology_class"]!="unresolved"]
        def frac(cls): return sum(r["topology_class"]==cls for r in resolved)/len(resolved) if resolved else math.nan
        out.append({"source_id":source,"assembly":assembly,"chromosome":chrom,"quartet_id":quartet,"bin_start":bs,"bin_end":bs+49999,"n_trees":len(rs),"n_resolved":len(resolved),"q1":frac("q1"),"q2":frac("q2"),"q3":frac("q3"),"alternative_support":frac("q2")+frac("q3") if resolved else math.nan,"alternative_minus_baseline":frac("q2")+frac("q3")-frac("q1") if resolved else math.nan,"mean_internal_branch":median([f(r,"internal_branch_length") for r in resolved]),"mean_pairwise":median([f(r,"mean_pairwise_patristic") for r in resolved])})
    return out

def recombination_summary(windows):
    rows=[]
    for source in sorted(set(r["source_id"] for r in windows)):
        rs=[r for r in windows if r["source_id"]==source and math.isfinite(f(r,"recomb_rate_cM_Mb"))]
        if not rs: continue
        values=[f(r,"recomb_rate_cM_Mb") for r in rs]; q25=statistics.quantiles(values,n=4,method="inclusive")[0] if len(values)>3 else min(values); q75=statistics.quantiles(values,n=4,method="inclusive")[2] if len(values)>3 else max(values)
        for chrom in sorted(set(r["chromosome"] for r in rs),key=lambda x:int(x) if x.isdigit() else 999):
            cr=[r for r in rs if r["chromosome"]==chrom]; low=[r for r in cr if f(r,"recomb_rate_cM_Mb")<=q25]; high=[r for r in cr if f(r,"recomb_rate_cM_Mb")>=q75]
            for label,x in [("low_quantile",low),("high_quantile",high)]: rows.append({"source_id":source,"chromosome":chrom,"group":label,"n_windows":len(x),"mean_recomb_cM_Mb":median([f(r,"recomb_rate_cM_Mb") for r in x]),"mean_q1":median([f(r,"q1") for r in x]),"mean_q2":median([f(r,"q2") for r in x]),"mean_q3":median([f(r,"q3") for r in x]),"mean_D":median([f(r,"D") for r in x]),"mean_B":median([f(r,"mean_pairwise_patristic") for r in x]),"mean_internal":median([f(r,"mean_quartet_internal_branch") for r in x]),"quantile_cut_low":q25,"quantile_cut_high":q75})
    return rows

def region_compare(windows):
    rows=[]
    for source,chrom in [("erato_10kb","2"),("erato_50kb_hera","2"),("erato_50kb_hmel","2")]:
        rs=[r for r in windows if r["source_id"]==source and r["chromosome"]==chrom]
        if not rs: continue
        for region in ("inside","overlap","outside","unclassified"):
            x=[r for r in rs if r["region_class"]==region]
            if not x: continue
            def mean(key):
                vals=[f(r,key) for r in x]; vals=[v for v in vals if math.isfinite(v)]
                return sum(vals)/len(vals) if vals else math.nan
            rows.append({"source_id":source,"chromosome":chrom,"region_class":region,"n_windows":len(x),"mean_q1":mean("q1"),"mean_q2":mean("q2"),"mean_q3":mean("q3"),"mean_D":mean("D"),"mean_B":mean("mean_pairwise_patristic"),"mean_internal":mean("mean_quartet_internal_branch"),"mean_recomb_cM_Mb":mean("recomb_rate_cM_Mb")})
    return rows

def quartet_region_effects(assignments):
    rows=[]
    for (source,chrom,quartet),rs in sorted(group(assignments,["source_id","chromosome","quartet_id"]).items()):
        if chrom!="2": continue
        by={region:[r for r in rs if r["region_class"]==region] for region in ("inside","overlap","outside")}
        if not by["inside"] or not by["outside"]: continue
        def mean(x,key):
            if key in ("q1","q2","q3"):
                vals=[float(r["topology_class"]==key) for r in x]
            else:
                vals=[f(r,key) for r in x]; vals=[v for v in vals if math.isfinite(v)]
            return sum(vals)/len(vals) if vals else math.nan
        inside={k:mean(by["inside"],k) for k in ("q1","q2","q3","internal_branch_length","mean_pairwise_patristic")}
        outside={k:mean(by["outside"],k) for k in ("q1","q2","q3","internal_branch_length","mean_pairwise_patristic")}
        rows.append({"source_id":source,"chromosome":chrom,"quartet_id":quartet,"taxa":rs[0]["taxa"],"n_inside":len(by["inside"]),"n_outside":len(by["outside"]),"inside_q1":inside["q1"],"outside_q1":outside["q1"],"delta_q1":inside["q1"]-outside["q1"],"inside_q2":inside["q2"],"outside_q2":outside["q2"],"delta_q2":inside["q2"]-outside["q2"],"inside_q3":inside["q3"],"outside_q3":outside["q3"],"delta_q3":inside["q3"]-outside["q3"],"inside_alternative":inside["q2"]+inside["q3"],"outside_alternative":outside["q2"]+outside["q3"],"delta_alternative":inside["q2"]+inside["q3"]-outside["q2"]-outside["q3"],"inside_internal":inside["internal_branch_length"],"outside_internal":outside["internal_branch_length"],"delta_internal":inside["internal_branch_length"]-outside["internal_branch_length"],"inside_B":inside["mean_pairwise_patristic"],"outside_B":outside["mean_pairwise_patristic"],"delta_B":inside["mean_pairwise_patristic"]-outside["mean_pairwise_patristic"]})
    return rows

def make_figures(bins, qbins):
    try:
        import os,tempfile
        os.environ.setdefault("MPLCONFIGDIR",str(Path(tempfile.gettempdir())/"msrc-paper-mplconfig"))
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize
    except ImportError: return
    FIG.mkdir(parents=True,exist_ok=True)
    # Heatmaps: alternative support for every quartet, source and focal chromosome.
    for source in sorted(set(r["source_id"] for r in qbins)):
        for chrom in ("2","15"):
            rs=[r for r in qbins if r["source_id"]==source and r["chromosome"]==chrom]
            if not rs: continue
            qs=sorted(set(r["quartet_id"] for r in rs)); xs=sorted(set(int(r["bin_start"]) for r in rs)); lookup={(r["quartet_id"],int(r["bin_start"])):f(r,"alternative_support") for r in rs}
            z=[[lookup.get((q,x),math.nan) for x in xs] for q in qs]
            fig,ax=plt.subplots(figsize=(13,7),constrained_layout=True); im=ax.imshow(z,aspect="auto",interpolation="nearest",cmap="magma",norm=Normalize(0,1)); ax.set_yticks(range(len(qs))); ax.set_yticklabels(qs,fontsize=6); ticks=list(range(0,len(xs),max(1,len(xs)//10))); ax.set_xticks(ticks); ax.set_xticklabels([f"{xs[i]/1e6:.2f}" for i in ticks]); ax.set_xlabel(f"chromosome {chrom} position (Mb)"); ax.set_ylabel("quartet"); ax.set_title(f"All-quartet alternative support: {source}, chromosome {chrom}"); fig.colorbar(im,ax=ax,label="q2 + q3");
            for ext in ("pdf","png","svg"): fig.savefig(FIG/f"heliconius_{source}_chr{chrom}_quartet_heatmap.{ext}",dpi=220 if ext=="png" else None)
            plt.close(fig)
    # Zoom plots use only source-matched, independently documented 2b intervals.
    zooms={
        "erato_10kb":(9114112,9292966),
        "erato_50kb_hera":(9114112,9292966),
        "erato_50kb_hmel":(2147755,2230606),
    }
    for source,(lo,hi) in zooms.items():
        rs=[r for r in bins if r["source_id"]==source and r["chromosome"]=="2" and int(r["bin_end"])>=lo-500000 and int(r["bin_start"])<=hi+500000]
        if not rs: continue
        rs=sorted(rs,key=lambda r:int(r["bin_start"])); x=[(int(r["bin_start"])+int(r["bin_end"]))/(2e6) for r in rs]
        fig,ax=plt.subplots(5,1,figsize=(12,10),sharex=True,constrained_layout=True,gridspec_kw={"height_ratios":[2,1,1,0.8,1.2]})
        for j,k,c in [(0,"q1","#4c72b0"),(0,"q2","#dd8452"),(0,"q3","#55a868"),(1,"D","#c44e52"),(2,"mean_pairwise_patristic","#8172b3")]:
            ax[j].plot(x,[f(r,k) for r in rs],".-",ms=4,label=k,color=c)
        ax[3].plot(x,[f(r,"effect_deviation") for r in rs],".-",ms=4,color="#7f7f7f"); ax[3].axhline(0,color="black",lw=.6); ax[3].set_ylabel("topology\ndeviation")
        ax[4].plot(x,[f(r,"recomb_rate_cM_Mb") for r in rs],".-",ms=4,color="#1b9e77"); ax[4].set_ylabel("recombination\ncM/Mb"); ax[4].set_xlabel("chromosome 2 position (Mb)")
        for a in ax: a.axvspan(lo/1e6,hi/1e6,color="#f0c419",alpha=.25,zorder=0)
        ax[0].legend(ncol=3); ax[0].set_ylabel("support"); ax[1].set_ylabel("D = q2 − q1"); ax[2].set_ylabel("mean patristic"); ax[0].set_title(f"Heliconius {source}: source-matched 2b zoom")
        for ext in ("pdf","png","svg"): fig.savefig(FIG/f"heliconius_{source}_chr2b_zoom.{ext}",dpi=220 if ext=="png" else None)
        plt.close(fig)
    # Genome-wide sustained-regime overview for the 10-kb source.
    rs=[r for r in bins if r["source_id"]=="erato_10kb"]
    if rs:
        fig,ax=plt.subplots(2,1,figsize=(15,7),sharex=True,constrained_layout=True)
        chroms=sorted(set(r["chromosome"] for r in rs),key=lambda x:int(x) if x.isdigit() else 999); offsets={}; off=0; tx=[]
        xx=[]; dd=[]; bb=[]
        for chrom in chroms:
            rr=sorted([r for r in rs if r["chromosome"]==chrom],key=lambda r:int(r["bin_start"])); length=max(int(r["bin_end"]) for r in rr); offsets[chrom]=off; tx.append((off+length/2e6,chrom)); off+=length+500000
            for r in rr: xx.append((offsets[chrom]+(int(r["bin_start"])+int(r["bin_end"]))/2)/1e6); dd.append(f(r,"alternative_minus_baseline")); bb.append(f(r,"mean_pairwise_patristic"))
        ax[0].plot(xx,dd,".",ms=2,color="#c44e52"); ax[0].axhline(0,color="black",lw=.6); ax[0].set_ylabel("(q2 + q3) − q1"); ax[0].set_title("Genome-wide sustained-regime scan: 10-kb Hera source")
        ax[1].plot(xx,bb,".",ms=2,color="#4c72b0"); ax[1].set_ylabel("mean patristic"); ax[1].set_xlabel("cumulative position (Mb)"); ax[1].set_xticks([x for x,_ in tx]); ax[1].set_xticklabels([f"chr{c}" for _,c in tx])
        for ext in ("pdf","png","svg"): fig.savefig(FIG/f"heliconius_genomewide_sustained_regimes.{ext}",dpi=220 if ext=="png" else None)
        plt.close(fig)

def run():
    assignments,windows=load(); bins=binned_windows(windows); qbins=quartet_bins(assignments)
    write_tsv(PROC/"critical_binned_regimes.tsv",bins)
    write_tsv(PROC/"critical_quartet_binned_statistics.tsv",qbins)
    write_tsv(RES/"all_quartet_dominance_summary.tsv",dominance(assignments,windows))
    write_tsv(RES/"recombination_quantile_summary.tsv",recombination_summary(windows))
    write_tsv(RES/"validated_region_comparison.tsv",region_compare(windows))
    write_tsv(RES/"quartet_region_effects_chr2.tsv",quartet_region_effects(assignments))
    runs=[]
    for key,rs in group(bins,["source_id","chromosome"]).items(): runs.extend(contiguous_runs(rs,"effect_deviation",0.10,2))
    write_tsv(RES/"sustained_regimes.tsv",sorted(runs,key=lambda r:(r["source_id"],int(r["chromosome"]) if r["chromosome"].isdigit() else 999,int(r["start"]))))
    branch=[]
    for (source,chrom),rs in group(bins,["source_id","chromosome"]).items():
        center=median([f(r,"mean_pairwise_patristic") for r in rs]); scale=mad([f(r,"mean_pairwise_patristic") for r in rs],center) or 1e-12
        for r in rs:
            z=(f(r,"mean_pairwise_patristic")-center)/(1.4826*scale) if math.isfinite(f(r,"mean_pairwise_patristic")) else math.nan
            if math.isfinite(z): branch.append({"source_id":source,"assembly":r["assembly"],"chromosome":chrom,"bin_start":r["bin_start"],"bin_end":r["bin_end"],"n_trees":r["n_trees"],"B":r["mean_pairwise_patristic"],"B_median_chrom":center,"B_robust_z":z,"topology_effect":r["alternative_minus_baseline"],"branch_only_candidate":str(abs(z)>=2 and abs(f(r,"alternative_minus_baseline"))<0.10).lower()})
    write_tsv(RES/"branch_length_shift_summary.tsv",sorted(branch,key=lambda r:abs(float(r["B_robust_z"])),reverse=True)[:200])
    make_figures(bins,qbins)
    report(assignments,windows,bins)

def report(assignments,windows,bins):
    dom=read_tsv(RES/"all_quartet_dominance_summary.tsv"); runs=read_tsv(RES/"sustained_regimes.tsv"); branch=read_tsv(RES/"branch_length_shift_summary.tsv"); recomb=read_tsv(RES/"recombination_quantile_summary.tsv"); compare=read_tsv(RES/"validated_region_comparison.tsv"); qeffects=read_tsv(RES/"quartet_region_effects_chr2.tsv")
    ag=[r for r in dom if r.get("scope")=="aggregate"]
    top_runs=sorted(runs,key=lambda r:int(r["span_bp"]),reverse=True)[:10]
    top_branch=[r for r in branch if r.get("branch_only_candidate")=="true"][:10]
    lines=["# Critical Heliconius all-quartet reassessment","", "This report supplements `initial_analysis_report.md`. It uses the downloaded Edelman et al. local trees, published species-history resolutions, and source-specific physical coordinates. All run/segment summaries use non-overlapping 50-kb bins; sliding-window trees remain descriptive linked observations.","", "## Coordinate resolution", "", "Thawornwattana et al. Supplementary Table 2 contains two distinct 2b entries. The `etales-8spp` entry pairs Hmel202001o:2,147,755–2,230,606 with Herato0206:880,283–1,059,137, while the `hmelv25-res` entry uses Hmel202001o:4,398,947–6,357,046 and has no H. erato coordinate. They are therefore retained as different dataset-specific definitions, not equated by label. The available Hmel-coordinate erato-clade trees cover both positions, but their source is not the `hmelv25-res` dataset; only the etales-8spp interval is used for their region classification. H. erato 15b remains unresolved in the independent table.","", "## Published topology comparison", "", "Per-quartet q1 now uses the fixed published species-history topology from Edelman et al. Fig. 1A/2 Tree 2, represented in `config/published_species_trees.tsv`; q2 and q3 are the other two deterministic resolutions. The prior background-majority topology is retained in the baseline mapping with an agreement flag. This makes disagreements visible rather than allowing region-specific relabeling.","", "## Aggregate dominance", "", "| source | chromosome | n windows | q2>q1 | q3>q1 | q2+q3>q1 | longest run (bp) |", "|---|---:|---:|---:|---:|---:|---:|"]
    lines += [f"| {r['source_id']} | {r['chromosome']} | {r['n_windows']} | {float(r['fraction_q2_gt_q1']):.3f} | {float(r['fraction_q3_gt_q1']):.3f} | {float(r['fraction_alt_gt_q1']):.3f} | {r['longest_run_span_bp']} |" for r in ag]
    lines += ["", "The complete per-quartet dominance table is `all_quartet_dominance_summary.tsv`; aggregate dominance is not required for a localized MSRC-compatible signal because an effect can involve only a subset of quartets.","", "## Sustained regimes", "", "| source | chromosome | start | end | span | bins | trees | mean alternative-minus-baseline | direction |", "|---|---:|---:|---:|---:|---:|---:|---:|---|"]
    lines += [f"| {r['source_id']} | {r['chromosome']} | {r['start']} | {r['end']} | {r['span_bp']} | {r['n_bins']} | {r['n_trees']} | {float(r['mean_effect']):.3f} | {r['direction']} |" for r in top_runs]
    lines += ["", "The strongest runs are not automatically inversion effects; they are genome-wide candidates defined before looking at the annotated intervals. Missing physical bins are not bridged.","", "## Quartet-specific chr2 effects", "", "The complete table is `quartet_region_effects_chr2.tsv`. It contains every quartet with both inside and outside observations for each source; the ranking below identifies contributors by absolute change in alternative support without discarding the remaining quartets.", "", "| source | quartet | n inside | n outside | delta alternative | delta internal | delta B |", "|---|---|---:|---:|---:|---:|---:|"]
    lines += [f"| {r['source_id']} | {r['quartet_id']} | {r['n_inside']} | {r['n_outside']} | {float(r['delta_alternative']):.3f} | {float(r['delta_internal']):.4f} | {float(r['delta_B']):.4f} |" for r in sorted(qeffects,key=lambda r:abs(float(r['delta_alternative'])),reverse=True)[:12]]
    lines += ["", "## Branch-length-only candidates", "", "| source | chromosome | bin | B robust z | topology effect | candidate |", "|---|---:|---:|---:|---:|---|"]
    lines += [f"| {r['source_id']} | {r['chromosome']} | {r['bin_start']}–{r['bin_end']} | {float(r['B_robust_z']):.2f} | {float(r['topology_effect']):.3f} | {r['branch_only_candidate']} |" for r in top_branch]
    lines += ["", "These are descriptive robust deviations from the chromosome-specific median; patristic distances remain source-tree substitution units, not coalescence times. A relative-rate or alignment-quality analysis would be needed to distinguish mutation-rate effects from genealogical-time effects.","", "## Recombination and validated regions", "", "The H. erato LDhelmet track is compared continuously by chromosome-specific recombination quantiles in `recombination_quantile_summary.tsv`. The current data do not provide a direct crossover measurement in inversion heterozygotes. `validated_region_comparison.tsv` contains topology, B, internal branch, and recombination summaries for the source-matched chr2 intervals. No H. erato 15b comparison is reported because the independent mapping is n/a.","", "## Decision", "", "Recommendation: retain Heliconius as exploratory and do not add it to the MSRC manuscript empirical section. The dataset contains extensive genome-wide discordance and candidate spatial regimes, but the available coordinate-matched comparisons do not yet demonstrate a coherent, independently supported topology-plus-branch-length regime coincident with a verified suppressed-recombination interval. Introgression, ancestral polymorphism/ILS, mutation-rate variation, and tree-estimation error remain strong alternatives.","", "Figures: `heliconius_*_quartet_heatmap.*` and `heliconius_genomewide_sustained_regimes.*`."]
    (RES/"critical_reassessment_report.md").write_text("\n".join(lines)+"\n")
    initial=RES/"initial_analysis_report.md"
    base=initial.read_text() if initial.exists() else ""
    marker="\n## Critical reassessment (published topology, coordinate-matched regions, sustained regimes)\n"
    if marker in base: base=base.split(marker)[0]
    initial.write_text(base+marker+"\n"+"\n".join(lines[2:])+"\n")

if __name__=="__main__":
    argparse.ArgumentParser().parse_args()
    run()
