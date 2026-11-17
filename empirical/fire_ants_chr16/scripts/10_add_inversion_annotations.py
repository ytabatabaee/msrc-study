#!/usr/bin/env python3
"""Create inversion-context tables and versioned fire-ant chr16 figures.

The Stolle supplementary workflow maps nearby breakpoint genes onto its chr16
axis. Its inversion ranges span multiple scaffold placements; the source Rmd
does not include the component-placement CSV. We therefore draw hatched outer
coordinate envelopes and label assembled mapped length, never a solid claim
that every base in each envelope is inverted.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc-paper-mplconfig")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical" / "fire_ants_chr16"
RESULTS = EMP / "results"
DATA = ROOT / "data" / "fire_ants_chr16" / "processed"
PD = EMP / "analysis" / "patterson_d"
PD_RESULTS = PD / "results"
FIG = EMP / "figures"
PD_FIG = PD / "figures"
SUPERGENE = (11_680_438, 27_917_498)

# Exact outer mapped extents from the upstream Stolle Rmd. The source's
# linkage_map_with_inversion.csv (which enumerates component placements) is not
# in the public archive; `assembled_bp_approx` is the reported summed mapped
# component length, rounded to the precision in Supplementary Note 7.
INVERSIONS = [
    {"inversion_id": "In(16)1", "chromosome": "chr16", "start_bp": 17_972_919, "end_bp": 28_149_888,
     "assembled_bp_approx": 9_590_000, "source_study": "Yan et al. (2020); Stolle et al. (2022) mapped annotation",
     "source_assembly": "Yan et al.: S. invicta SB reference; original inversion size/breakpoints",
     "target_assembly": "Stolle Si_gnGA/gng20170922wFex.fa chr16 analysis axis",
     "coordinate_mapping_method": "Stolle upstream twisst_interpretation.rmd formulas; breakpoint-adjacent coding-gene mapping described for Supplementary Fig. 13",
     "mapping_status": "mapped outer coordinate envelope; component placement intervals unavailable",
     "breakpoint_certainty": "outer envelope endpoints exact as encoded in upstream Rmd; internal component gaps unresolved; assembled length rounded",
     "source_figure_or_table": "Stolle Supplementary Fig. 13; Supplementary Note 7; upstream twisst_interpretation.rmd",
     "source_url": "https://wurmlab.com/publications/2022-recurring-adaptive-introgression/rai-supplementary-Information.pdf",
     "notes": "Rmd: start=17621722+351197; end=27455557+694331. Figure labels this In(16)1. Approximate assembled mapped length 9.59 Mb; envelope is not a contiguous inversion-size representation."},
    {"inversion_id": "In(16)2", "chromosome": "chr16", "start_bp": 16_509_882, "end_bp": 17_962_688,
     "assembled_bp_approx": 840_000, "source_study": "Yan et al. (2020); Stolle et al. (2022) mapped annotation",
     "source_assembly": "Yan et al.: S. invicta SB reference; original inversion size/breakpoints",
     "target_assembly": "Stolle Si_gnGA/gng20170922wFex.fa chr16 analysis axis",
     "coordinate_mapping_method": "Stolle upstream twisst_interpretation.rmd formulas; breakpoint-adjacent coding-gene mapping described for Supplementary Fig. 13",
     "mapping_status": "mapped outer coordinate envelope; component placement intervals unavailable",
     "breakpoint_certainty": "outer envelope endpoints exact as encoded in upstream Rmd; internal component gaps unresolved; assembled length rounded",
     "source_figure_or_table": "Stolle Supplementary Fig. 13; Supplementary Note 7; upstream twisst_interpretation.rmd",
     "source_url": "https://wurmlab.com/publications/2022-recurring-adaptive-introgression/rai-supplementary-Information.pdf",
     "notes": "Rmd: start=16509882; end=17621722+340966. Approximate assembled mapped length 0.84 Mb; envelope is not a contiguous inversion-size representation."},
    {"inversion_id": "In(16)3", "chromosome": "chr16", "start_bp": 11_732_910, "end_bp": 16_478_110,
     "assembled_bp_approx": 1_200_000, "source_study": "Yan et al. (2020); Stolle et al. (2022) mapped annotation",
     "source_assembly": "Yan et al.: S. invicta SB reference; original inversion size/breakpoints",
     "target_assembly": "Stolle Si_gnGA/gng20170922wFex.fa chr16 analysis axis",
     "coordinate_mapping_method": "Stolle upstream twisst_interpretation.rmd formulas; breakpoint-adjacent coding-gene mapping described for Supplementary Fig. 13",
     "mapping_status": "mapped outer coordinate envelope; component placement intervals unavailable",
     "breakpoint_certainty": "outer envelope endpoints exact as encoded in upstream Rmd; internal component gaps unresolved; assembled length rounded",
     "source_figure_or_table": "Stolle Supplementary Fig. 13; Supplementary Note 7; upstream twisst_interpretation.rmd",
     "source_url": "https://wurmlab.com/publications/2022-recurring-adaptive-introgression/rai-supplementary-Information.pdf",
     "notes": "Rmd: start=11732910; end=15842782+635328. Stolle window-assignment code instead caps midpoint membership at 16339881; this discrepancy is disclosed. Approximate assembled mapped length 1.2 Mb; envelope is not contiguous inversion size."},
]
COLORS = {"In(16)1": "#8f5da2", "In(16)2": "#56a6c2", "In(16)3": "#3155a4"}
FROZEN_HASHES = {
    "data/fire_ants_chr16/processed/stage1_window_index.tsv": "96ed3fcca3d351a96b2434f17aa7510349d94138d11d60ccb56d45d1b57e237d",
    "data/fire_ants_chr16/processed/stage4a_window_quartet_support.tsv": "b8ff3ea95ba31ae5553f56f0aeb45a393bdb0ba2b19081fc1aca822dfdd70758",
    "empirical/fire_ants_chr16/results/stage4a_region_summary.tsv": "b902d6cdc472fbe2b5749d56e3edb44a093129d8dc418fae05078f55f5bbabca",
    "empirical/fire_ants_chr16/results/stage4b_primary_test.tsv": "310fe14c9cdb0e4d4ba6b1b9355498cfa1e845e56de119061f395c7011d9ec8b",
    "empirical/fire_ants_chr16/results/stage6_cu_branch_length_comparison.tsv": "ab13f0a0af2cf2d7eabd436bb51bb1408f2b04ac426a686cdafe3f6f7d6d337d",
    "empirical/fire_ants_chr16/analysis/patterson_d/results/window_patterson_d.tsv": "c69d4f640d351bf1f4868e131a173b750a579e714f9f88811784eee59a9e374a",
}

def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({k: row.get(k, "") for k in fields} for row in rows)

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

def verify_frozen_inputs():
    for relpath, expected in FROZEN_HASHES.items():
        path=ROOT/relpath
        observed=sha256(path)
        if observed!=expected:
            raise RuntimeError("Frozen input changed: %s (%s != %s)"%(relpath,observed,expected))

def overlap_len(a: int, b: int, c: int, d: int) -> int:
    return max(0, min(b, d) - max(a, c) + 1)

def source_midpoint_membership(mid: float, name: str) -> bool:
    # These reproduce the explicit Stolle upstream Rmd conditions, not a
    # redefinition of the frozen analysis-region membership.
    if name == "In(16)1": return mid > 17_972_920 and mid < 28_149_888
    if name == "In(16)2": return mid > 16_509_882 and mid < 17_962_688
    if name == "In(16)3": return mid > 11_732_910 and mid < 16_339_881
    return False

def add_context(ax, ylim=(0, 1), labels=True):
    ax.set_ylim(*ylim)
    ax.axvspan(SUPERGENE[0]/1e6, SUPERGENE[1]/1e6, color="#777777", alpha=.16,
               label="frozen supergene analysis interval" if labels else None, zorder=0)

def plot_structural_panel(ax):
    """Show the frozen span and each mapped inversion on separate rows."""
    ax.set_xlim(0,30); ax.set_ylim(0,4)
    ax.set_yticks([3.5,2.5,1.5,.5])
    ax.set_yticklabels(["Supergene analysis span","In(16)1 · 9.59 Mb assembled","In(16)2 · 0.84 Mb assembled","In(16)3 · 1.20 Mb assembled"],fontsize=7)
    ax.set_title("Structural context (mapped outer envelopes; not continuous inversion sequence)",loc="left",fontsize=9)
    ax.add_patch(Rectangle((SUPERGENE[0]/1e6,3.18),(SUPERGENE[1]-SUPERGENE[0])/1e6,.64,
                 facecolor="#777777",edgecolor="#666666",alpha=.26,lw=.6,zorder=2))
    for yi,inv in zip((2.5,1.5,.5),INVERSIONS):
        x0,x1=inv["start_bp"]/1e6,inv["end_bp"]/1e6
        ax.add_patch(Rectangle((x0,yi-.28),x1-x0,.56,facecolor=COLORS[inv["inversion_id"]],
                     edgecolor=COLORS[inv["inversion_id"]],alpha=.18,lw=1.1,ls="--",hatch="///",zorder=3))

def save_figure(fig, stem: Path):
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        out = stem.with_suffix("."+ext)
        if ext == "png": fig.savefig(out, dpi=240, bbox_inches="tight")
        elif ext == "pdf":
            with PdfPages(out, metadata={"CreationDate":None,"ModDate":None}) as pdf: pdf.savefig(fig, bbox_inches="tight")
        else: fig.savefig(out, bbox_inches="tight")
    plt.close(fig)

def load_data():
    windows = read_tsv(DATA / "stage1_window_index.tsv")
    support = {r["window_index"]: r for r in read_tsv(DATA / "stage4a_window_quartet_support.tsv")}
    distances = {r["window_index"]: r for r in read_tsv(RESULTS / "stage9_distance" / "stage9_window_patristic_distance.tsv")}
    drows = read_tsv(PD_RESULTS / "window_patterson_d.tsv")
    topology = read_tsv(PD_RESULTS / "induced_quartet_topologies.tsv")
    primary = [r for r in drows if r["quartet_id"] in ("Q1_primary_Argentina", "Q2_primary_Argentina")]
    d_by = {(r["quartet_id"], r["window_index"]): r for r in primary}
    topo_by = {(r["quartet_id"], r["window_index"]): r for r in topology if r["quartet_id"] in ("Q1_primary_Argentina", "Q2_primary_Argentina")}
    return windows, support, distances, d_by, topo_by

def write_tables(windows, d_by, support):
    base = RESULTS / "chr16_inversion_intervals.tsv"
    fields = ["inversion_id","chromosome","start_bp","end_bp","length_bp","assembled_mapped_bp_approx","source_study","source_assembly","target_assembly","coordinate_mapping_method","mapping_status","breakpoint_certainty","source_figure_or_table","source_url","notes"]
    inv_rows=[]
    for inv in INVERSIONS:
        row={k:inv[k] for k in fields if k in inv}
        row["length_bp"] = inv["end_bp"]-inv["start_bp"]+1
        row["assembled_mapped_bp_approx"] = inv["assembled_bp_approx"]
        inv_rows.append(row)
    write_tsv(base,inv_rows,fields)

    overlap_rows=[]; summary=[]; descriptive=[]
    for inv in INVERSIONS:
        lo,hi=inv["start_bp"],inv["end_bp"]
        full=partial=0; assigned=[]
        for w in windows:
            if w["chrom"]!="chr16": continue
            ws,we,wm=int(w["start"]),int(w["end"]),float(w["mid"])
            ov=overlap_len(ws,we,lo,hi)
            cls="none" if not ov else ("full" if ws>=lo and we<=hi else "partial")
            full += cls=="full"; partial += cls=="partial"
            mid_member=source_midpoint_membership(wm,inv["inversion_id"])
            # window-specific D is reported as assigned by the source midpoint
            # rule only; SNPs have no continuous chromosome lift.
            qrows=[d_by.get((q,w["window_index"])) for q in ("Q1_primary_Argentina","Q2_primary_Argentina")]
            assigned.append((w,mid_member,qrows))
            for q,r in zip(("Q1_primary_Argentina","Q2_primary_Argentina"),qrows):
                overlap_rows.append({"inversion_id":inv["inversion_id"],"window_index":w["window_index"],"window_start_bp":ws,"window_end_bp":we,"window_mid_bp":wm,"bounding_envelope_overlap_bp":ov,"window_envelope_overlap_class":cls,"stolle_source_midpoint_assignment":int(mid_member),"region":w["region"],"quartet_id":q,"patterson_D":r["D"] if r else "NA","informative_snps_in_window":r["informative"] if r else "NA","gene_body_span_bp_window":r["target_gene_body_span_bp"] if r else "NA","gene_level_coordinate_overlap":"NA_no_per_gene_target_axis_coordinates","note":"Envelope overlap can include gaps; exact inverted-segment/BUSCO/SNP overlap unavailable"})
        midwins=[(w,qs) for w,m,qs in assigned if m]
        for qidx,q in enumerate(("Q1_primary_Argentina","Q2_primary_Argentina")):
            rs=[d_by.get((q,w["window_index"])) for w,m,qs in assigned if m]
            rs=[r for r in rs if r]
            informative=sum(int(r["informative"]) for r in rs)
            d_inform=[r for r in rs if int(r["informative"])>0]
            summary.append({"inversion_id":inv["inversion_id"],"quartet_id":q,"target_assembly":"Stolle Si_gnGA chr16","mapped_outer_start_bp":lo,"mapped_outer_end_bp":hi,"mapped_envelope_bp":hi-lo+1,"approx_assembled_segment_bp":inv["assembled_bp_approx"],"fully_contained_window_envelopes":full,"partially_overlapping_window_envelopes":partial,"nonoverlapping_window_envelopes":sum(1 for w in windows if w["chrom"]=="chr16")-full-partial,"source_midpoint_assigned_tree_windows":sum(m for w,m,qs in assigned),"source_midpoint_assigned_D_windows_with_informative_sites":sum(1 for w,m,qs in assigned if m and qs[qidx] and int(qs[qidx]["informative"])>0),"informative_snps_in_source_midpoint_assigned_windows":sum(int(qs[qidx]["informative"]) for w,m,qs in assigned if m and qs[qidx]),"informative_snps_in_envelope_overlapping_windows":informative,"callable_gene_bp_overlap":"NA","fraction_inversion_length_covered_by_analyzed_gene_bases":"NA","uncovered_inversion_subintervals":"NA_component_segment_coordinates_unavailable","coverage_note":"Window-envelope overlap is a bounding proxy, not exact inversion overlap; SNPs are assigned only by Stolle's window-midpoint rule. Individual mapped BUSCO gene coordinates on target chr16 are unavailable."})
        for qidx,q in enumerate(("Q1_primary_Argentina","Q2_primary_Argentina")):
            rs=[d_by.get((q,w["window_index"])) for w,m,qs in assigned if m]
            rs=[r for r in rs if r]
            inf=sum(int(r["informative"]) for r in rs)
            informative=[r for r in rs if int(r["informative"])>0]
            abba=sum(float(r["abba"]) for r in informative); baba=sum(float(r["baba"]) for r in informative)
            descriptive.append({"inversion_id":inv["inversion_id"],"quartet_id":q,"tree_windows_by_source_midpoint":sum(m for w,m,qs in assigned),"mean_q_species":"NA" if not midwins else sum(float(support[w["window_index"]]["q_S"]) for w,m,qs in assigned if m)/len(midwins),"mean_q_haplotype":"NA" if not midwins else sum(float(support[w["window_index"]]["q_H"]) for w,m,qs in assigned if m)/len(midwins),"mean_q_third":"NA" if not midwins else sum(float(support[w["window_index"]]["q_3"]) for w,m,qs in assigned if m)/len(midwins),"mean_D_tree_qH_minus_qS":"NA" if not midwins else sum(float(support[w["window_index"]]["D"]) for w,m,qs in assigned if m)/len(midwins),"patterson_D_pooled":"NA" if not inf else (abba-baba)/(abba+baba),"ABBA":"NA" if not inf else abba,"BABA":"NA" if not inf else baba,"informative_snps":inf,"window_assignment":"Stolle source midpoint rule; descriptive, not exact gene/SNP inversion overlap","status":"no assigned windows; NA" if not midwins else ("no informative sites; NA" if not inf else "descriptive only")})
    write_tsv(RESULTS/"chr16_inversion_window_overlap.tsv",overlap_rows,["inversion_id","window_index","window_start_bp","window_end_bp","window_mid_bp","bounding_envelope_overlap_bp","window_envelope_overlap_class","stolle_source_midpoint_assignment","region","quartet_id","patterson_D","informative_snps_in_window","gene_body_span_bp_window","gene_level_coordinate_overlap","note"])
    write_tsv(RESULTS/"chr16_inversion_coverage_summary.tsv",summary,list(summary[0]))
    write_tsv(RESULTS/"chr16_inversion_specific_descriptive_summary.tsv",descriptive,list(descriptive[0]))

def plots(windows,support,distances,d_by,topo_by):
    chr16=sorted([w for w in windows if w["chrom"]=="chr16"],key=lambda r:float(r["mid"]))
    x=[float(w["mid"])/1e6 for w in chr16]
    # Annotated main figure, preserving focal topology and the existing spatial null.
    fig,axs=plt.subplots(3,1,figsize=(10.5,9),gridspec_kw={"height_ratios":[1,1.5,1]})
    ax=axs[0]; plot_structural_panel(ax); ax.set_title("A. Structural context on the Stolle chr16 axis",loc="left",fontsize=10)
    ax.set_xticklabels([])
    ax=axs[1]; add_context(ax,(0,1),False)
    for field,color,label in [("q_S","#345995","q_species"),("q_H","#b23a48","q_haplotype"),("q_3","#2f7f4f","q_third")]:
        ax.plot(x,[float(support[w["window_index"]][field]) for w in chr16],marker="o",ms=2.6,lw=.7,color=color,label=label)
    ax.set_title("B. Frozen focal-quartet topology support",loc="left",fontsize=10);ax.set_ylabel("TWISST weight");ax.set_xlabel("chr16 coordinate (Mb)");ax.legend(frameon=False,ncol=3,fontsize=8,loc="lower left")
    circ=sorted(read_tsv(RESULTS/"stage4b_circular_shifts.tsv"),key=lambda r:float(r.get("delta_D",r.get("Delta_D","0"))))
    ax=axs[2]
    vals=[]
    for r in circ:
        v=r.get("delta_D",r.get("Delta_D",""))
        if v: vals.append(float(v))
    if vals:
        vals.sort();ax.scatter(range(1,len(vals)+1),vals,s=12,color="#345995");ax.axhline(float(read_tsv(RESULTS/"stage4b_primary_test.tsv")[0].get("observed_delta_D",0)),color="#b23a48",lw=1,label="observed")
    ax.set_title("C. Existing circular-shift null (values unchanged)",loc="left",fontsize=10);ax.set_ylabel("D_tree contrast");ax.set_xlabel("ranked circular alignment");ax.legend(frameon=False,fontsize=8)
    fig.tight_layout(h_pad=1.1)
    save_figure(fig,FIG/"fire_ants_main_with_inversions")

    # Aligned topology and distance track.
    fig,axs=plt.subplots(4,1,figsize=(11,9),gridspec_kw={"height_ratios":[1,.95,.75,.95]},sharex=True)
    plot_structural_panel(axs[0]);axs[0].set_title("A. Structural context",loc="left",fontsize=10)
    for ax in axs[1:]: add_context(ax,(0,1),False)
    for field,color,label in [("q_S","#345995","q_species"),("q_H","#b23a48","q_haplotype"),("q_3","#2f7f4f","q_third")]:
        axs[1].plot(x,[float(support[w["window_index"]][field]) for w in chr16],"o-",ms=2.4,lw=.65,color=color,label=label)
    axs[1].set_title("B. Frozen topology support",loc="left",fontsize=10);axs[1].set_ylabel("TWISST weight");axs[1].legend(frameon=False,ncol=3,fontsize=8)
    axs[2].plot(x,[float(support[w["window_index"]]["D"]) for w in chr16],"o-",ms=2.5,lw=.65,color="#8464a0");axs[2].axhline(0,color="black",lw=.5);axs[2].set_title("C. D_tree = q_haplotype − q_species",loc="left",fontsize=10);axs[2].set_ylabel("D_tree")
    axs[3].plot(x,[float(distances[w["window_index"]]["A_relative"]) for w in chr16],"o-",ms=2.4,lw=.65,color="#587b58");axs[3].axhline(0,color="black",lw=.5);axs[3].set_title("D. Patristic-distance contrast (substitution units)",loc="left",fontsize=10);axs[3].set_ylabel("A_relative");axs[3].set_xlabel("chr16 coordinate (Mb)")
    fig.tight_layout(h_pad=.8)
    save_figure(fig,FIG/"fire_ants_chr16_topology_distance_inversions")

    # Matched nucleotide D and exact induced topology D for primary quartets.
    fig,axs=plt.subplots(3,1,figsize=(11,7.8),gridspec_kw={"height_ratios":[1,1.2,1.2]},sharex=True)
    plot_structural_panel(axs[0]);axs[0].set_title("A. Structural context",loc="left",fontsize=10)
    for ax in axs[1:]:add_context(ax,(-1,1),False);ax.axhline(0,color="black",lw=.5)
    colors={"Q1_primary_Argentina":"#2878b5","Q2_primary_Argentina":"#d95f02"}
    for q in colors:
        qx=[];dy=[];ty=[]
        for w in chr16:
            wi=w["window_index"];r=d_by.get((q,wi));t=topo_by.get((q,wi))
            if r and r["D"]!="NA":qx.append(float(w["mid"])/1e6);dy.append(float(r["D"]))
            if t and t["D_topology"]!="NA":ty.append((float(w["mid"])/1e6,float(t["D_topology"])))
        axs[1].scatter(qx,dy,s=16,color=colors[q],label=q.split("_")[0]);axs[2].scatter([a for a,b in ty],[b for a,b in ty],s=16,color=colors[q],label=q.split("_")[0])
    axs[1].set_title("B. Nucleotide Patterson's D (window-pool estimate)",loc="left",fontsize=10);axs[1].set_ylabel("D_site");axs[1].legend(frameon=False,ncol=2,fontsize=8)
    axs[2].set_title("C. Exact induced quartet topology D = q3 − q2",loc="left",fontsize=10);axs[2].set_ylabel("D_topology");axs[2].set_xlabel("chr16 coordinate (Mb)");axs[2].legend(frameon=False,ncol=2,fontsize=8)
    fig.tight_layout(h_pad=.8)
    save_figure(fig,PD_FIG/"figure_B_D_and_TWISST_with_inversions")

def write_docs_and_manifest():
    md=RESULTS/"chr16_inversion_provenance.md"
    md.write_text("""# Chromosome-16 inversion coordinate provenance\n\nThe target coordinate axis is the Stolle et al. chr16 analysis axis (`Si_gnGA` / `gng20170922wFex.fa`). Stolle Supplementary Fig. 13 says Yan et al.'s inversion breakpoint regions were mapped by aligning the closest protein-coding genes to the Stolle assembly. The exact outer envelopes below were recovered from Stolle's public [`Topology weighting/twisst_interpretation.rmd`](https://github.com/wurmlab/2021-fire-ant-social-supergene-introgression/blob/bdc7823941680fa560e61b6467af9f77950f64b0/Topology%20weighting/twisst_interpretation.rmd), in the `inversions` data frame (SHA256 `7192b5819ecef9257ef584b68379e5fac2f45cca34e64d4976876a6d5b63440b`). The supplementary PDF is [here](https://wurmlab.com/publications/2022-recurring-adaptive-introgression/rai-supplementary-Information.pdf) (retrieved PDF SHA256 `f8357d89af425c198d53f3ce8f4de8ef8ead2ad8e49cc48db4f0388d878f797c`); Stolle paper DOI: [10.1038/s41467-022-28806-7](https://doi.org/10.1038/s41467-022-28806-7); Yan paper DOI: [10.1038/s41559-019-1081-1](https://doi.org/10.1038/s41559-019-1081-1).\n\n| Name | Stolle target outer extent (1-based inclusive) | Assembled mapped length in SI |\n|---|---:|---:|\n| In(16)1 | chr16:17,972,919–28,149,888 | about 9.59 Mb |\n| In(16)2 | chr16:16,509,882–17,962,688 | about 0.84 Mb |\n| In(16)3 | chr16:11,732,910–16,478,110 | about 1.2 Mb |\n\nThese coordinates are mapped **outer envelopes** across scaffold placements. Stolle's source calculates assembled lengths by summing component `chr_end-chr_start`; that placement CSV is not in the public source archive. Therefore the envelope span is not the true inverted-base length and may include unplaced/gap space. Figures use hatched outlines and show both the envelope and the reported assembled length. Internal component segments and exact base-level overlap cannot be reconstructed from the available files.\n\nThe upstream Rmd assigns window midpoints using separate strict inequalities: In(16)1 `17,972,920 < mid < 28,149,888`; In(16)2 `16,509,882 < mid < 17,962,688`; In(16)3 `11,732,910 < mid < 16,339,881`. The In(16)3 midpoint cutoff differs from its outer envelope endpoint (16,478,110); this is retained and flagged, not silently harmonized. These rules reproduce Stolle's published window labels. They are separate from our frozen `chr16_supergene` interval, chr16:11,680,438–27,917,498.\n\nYan's source-axis intervals in the existing audit are In(16)1 14,549,064–24,031,576 (9.483 Mb), In(16)2 13,705,210–24,030,990 (10.326 Mb), and In(16)3 12,612,565–13,683,100 (1.071 Mb). These are not overlaid because they are on the Yan SB reference axis. The names are consistent in Stolle Fig. 13 and its Rmd region definitions; the SI prose has an apparent typo calling the largest 9.59-Mb inversion In(16)3, although the figure, code, and later 1.2-Mb description identify the 9.59-Mb region as In(16)1. The mapped assembled lengths differ from Yan's full source inversion sizes, especially for In(16)2, because the Stolle assembly contains only partial mapped sequence; they are not evidence of a renumbering.\n\nCoordinates in the source are GRanges-style one-based inclusive starts/ends. No chromosome-scale individual BUSCO gene coordinates or individual VCF SNP lifts are present in the frozen inputs. Consequently exact analyzed-gene-base coverage, SNP-base overlap, and uncovered inversion component subintervals are reported as unavailable. Window-envelope overlap is tabulated separately as a bounding proxy. Nucleotide informative counts are assigned to windows under Stolle's midpoint rule and are not exact inversion SNP counts.\n""")
    cap="""# Inversion-annotated fire-ant figure captions\n\nThree published inversions, In(16)1–3, distinguish Sb from SB. The analysis continues to treat the frozen chr16:11,680,438–27,917,498 span as one recombination-suppressed supergene analysis region. Hatched outlines mark Stolle-mapped outer coordinate envelopes; labels give the approximate amount of assembled mapped sequence, because exact component intervals are not available in the public source deposit. Window coverage varies: Stolle's source midpoint assignments place 47 windows in In(16)1, four in In(16)2, and none in In(16)3. Our frozen tree and nucleotide tracks are therefore not separately resolved within every inversion. Inversion annotations are structural context, not direct local cM/Mb or a quantitative recombination-rate track. They do not imply that inversions alone cause all recombination suppression, nor that each inversion independently produces a significant topology switch. The published introgression interpretation and alternatives including ancestral polymorphism remain distinct hypotheses.\n"""
    (RESULTS/"chr16_inversion_figure_caption.md").write_text(cap)
    inp=[DATA/"stage1_window_index.tsv",DATA/"stage4a_window_quartet_support.tsv",RESULTS/"stage9_distance"/"stage9_window_patristic_distance.tsv",PD_RESULTS/"window_patterson_d.tsv",PD_RESULTS/"induced_quartet_topologies.tsv"]
    manifest={"task":"inversion structural annotation","analysis_rerun":False,"source_script":"empirical/fire_ants_chr16/scripts/10_add_inversion_annotations.py","inputs_sha256":{str(p.relative_to(ROOT)):sha256(p) for p in inp},"frozen_result_sha256_before_and_after":FROZEN_HASHES,"coordinate_axis":"Stolle Si_gnGA/gng20170922wFex chr16","inversion_intervals":INVERSIONS,"midpoint_assignment_counts":{"In(16)1":47,"In(16)2":4,"In(16)3":0},"limitations":["mapped outer envelopes are not continuous component-interval unions","gene-level target chr16 coordinates are unavailable","individual nucleotide SNPs cannot be lifted to target chr16 in these inputs"],"figure_outputs":["empirical/fire_ants_chr16/figures/fire_ants_main_with_inversions.{png,pdf,svg}","empirical/fire_ants_chr16/figures/fire_ants_chr16_topology_distance_inversions.{png,pdf,svg}","empirical/fire_ants_chr16/analysis/patterson_d/figures/figure_B_D_and_TWISST_with_inversions.{png,pdf,svg}"]}
    (RESULTS/"chr16_inversion_annotation_manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")

def main():
    verify_frozen_inputs()
    windows,support,distances,d_by,topo_by=load_data()
    write_tables(windows,d_by,support)
    plots(windows,support,distances,d_by,topo_by)
    write_docs_and_manifest()
    verify_frozen_inputs()
    print("wrote inversion provenance/coverage tables and 27 versioned figure files")

if __name__=="__main__": main()
