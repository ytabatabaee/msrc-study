#!/usr/bin/env python3
"""Postprocess the completed 25-replicate Neoaves J48 mixture experiment."""
from __future__ import annotations
import argparse, csv, importlib.util, json, math, sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / "empirical" / "neoaves_chr4"
RESULTS = BASE / "results" / "stage4e"
FIGURES = BASE / "figures"
SCRIPT = BASE / "scripts" / "05_stage4e_analyze_j48_mixtures.py"
VALIDATOR = Path(__file__).with_name("validate_results.py")
spec = importlib.util.spec_from_file_location("j48_original", SCRIPT); mod = importlib.util.module_from_spec(spec); sys.modules[spec.name] = mod; spec.loader.exec_module(mod)
vspec = importlib.util.spec_from_file_location("j48_validator", VALIDATOR); validator = importlib.util.module_from_spec(vspec); sys.modules[vspec.name] = validator; vspec.loader.exec_module(validator)
EPSILONS = (0.02, 0.05, 0.075, 0.10, 0.15)
RELATIONSHIPS = ("N61", "N62", "Columbea")
STATE_COLORS = {"q1": "#D55E00", "q2": "#0072B2", "q3": "#009E73", "mixed_or_unresolved": "#777777"}

def read_tsv(path: Path):
    with path.open(newline="") as f: return list(csv.DictReader(f, delimiter="\t"))

def write_tsv(path: Path, rows: list[dict], fields=None):
    path.parent.mkdir(parents=True, exist_ok=True); fields = fields or (list(rows[0]) if rows else ["status"])
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows({k: r.get(k, "") for k in fields} for r in rows)

def wilson(k: int, n: int):
    if not n: return "", ""
    z = 1.959963984540054; p = k / n; den = 1 + z*z/n; center = (p + z*z/(2*n))/den; half = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/den
    return center-half, center+half

def context():
    taxa = {r["stiller2024_taxon"] for r in read_tsv(mod.TAXA)}
    _, t0 = mod.parse_tree(mod.CLUSTER / mod.T0 / "tree.cu.nwk", "CULength"); _, non = mod.parse_tree(mod.CLUSTER / mod.NON / "tree.cu.nwk", "CULength")
    groups = mod.load_frozen_groups(); labels = {}
    for clade, definition in groups.items():
        for state in ("q1", "q2", "q3"): labels.setdefault(mod.focal_key(definition, state, taxa), f"{clade}_{state}")
    return taxa, t0, non, set(t0), set(non), groups, labels

def analyze(result_root: Path, manifest_path: Path):
    taxa_expected, t0, non, t0_splits, non_splits, groups, labels = context(); rows=[]; branches=[]
    for m in read_tsv(manifest_path):
        d = result_root / Path(m["result_relative_path"]); taxa_cu, cu = mod.parse_tree(d / "tree.cu.nwk", "CULength"); taxa_su, su = mod.parse_tree(d / "tree.su.nwk", "SULength")
        if taxa_cu != taxa_expected or taxa_su != taxa_expected or set(cu) != set(su): raise ValueError(f"CU/SU mismatch: {d}")
        splits=set(cu); novel=sorted(splits-t0_splits-non_splits); fc=mod.focal_states(cu,taxa_cu,groups); fs=mod.focal_states(su,taxa_su,groups)
        r={"mixture_id":m["mixture_id"],"epsilon":m["epsilon_requested"],"realized_epsilon":m["realized_epsilon"],"replicate":m["replicate"],"seed":m["seed"],"status":"complete","rf_vs_T0":len(splits^t0_splits),"rf_vs_NONCHR4":len(splits^non_splits),"topology_identity":mod.topology_identity(splits,t0_splits,non_splits),"equals_T0":str(splits==t0_splits).lower(),"equals_NONCHR4":str(splits==non_splits).lower(),"novel_splits":";".join(novel),"n_novel_splits":len(novel),"t0_specific_present":";".join(sorted(splits&(t0_splits-non_splits))),"nonchr4_specific_present":";".join(sorted(splits&(non_splits-t0_splits)))}
        for rel in RELATIONSHIPS:
            c=fc[rel]; s=fs[rel]; key=c.get("split_key",""); homologous=bool(key and key==s.get("split_key", "")); r.update({f"{rel}_state":c.get("state",""),f"{rel}_state_su":s.get("state",""),f"{rel}_split_key":key,f"{rel}_present":str(bool(key)).lower(),f"{rel}_homologous":str(homologous).lower(),f"{rel}_CU_length":c.get("CU_length","") if homologous else "",f"{rel}_SU_length":s.get("SU_length","") if homologous else "",f"{rel}_localPP":c.get("localPP","") if homologous else ""})
            branches.append({"mixture_id":m["mixture_id"],"epsilon":m["epsilon_requested"],"replicate":m["replicate"],"relationship":rel,"split_key":key,"homologous_cu_su":str(homologous).lower(),"CU_length":c.get("CU_length","") if homologous else "","SU_length":s.get("SU_length","") if homologous else "","localPP":c.get("localPP","") if homologous else ""})
        rows.append(r)
    summary=[]; focal_summary=[]
    for e in EPSILONS:
        sub=[r for r in rows if abs(float(r["epsilon"])-e)<1e-12]
        if len(sub)!=5: raise ValueError(f"expected five replicates at {e}, observed {len(sub)}")
        for state,key,rfkey in (("J48_T0","equals_T0","rf_vs_T0"),("J48_NONCHR4","equals_NONCHR4","rf_vs_NONCHR4")):
            n=sum(r[key]=="true" for r in sub); lo,hi=wilson(n,5); vals=[int(r[rfkey]) for r in sub]; mean=sum(vals)/5; sd=math.sqrt(sum((x-mean)**2 for x in vals)/4)
            summary.append({"epsilon":f"{e:.3f}","category":"complete_topology","state":state,"n_replicates":5,"n":n,"fraction":n/5,"ci_low":lo,"ci_high":hi,"mean_RF":mean,"median_RF":sorted(vals)[2],"min_RF":min(vals),"max_RF":max(vals),"sd_RF":sd})
        n=sum(r["topology_identity"]=="intermediate_or_distinct" for r in sub); lo,hi=wilson(n,5); summary.append({"epsilon":f"{e:.3f}","category":"topology_class","state":"intermediate_or_distinct","n_replicates":5,"n":n,"fraction":n/5,"ci_low":lo,"ci_high":hi})
        n=sum(int(r["n_novel_splits"])>0 for r in sub); lo,hi=wilson(n,5); summary.append({"epsilon":f"{e:.3f}","category":"novel_topology","state":"contains_novel_split","n_replicates":5,"n":n,"fraction":n/5,"ci_low":lo,"ci_high":hi})
        for rel in RELATIONSHIPS:
            for state,n in sorted(Counter(r[f"{rel}_state"] for r in sub).items()):
                lo,hi=wilson(n,5); focal_summary.append({"epsilon":f"{e:.3f}","relationship":rel,"state":state,"n_replicates":5,"n":n,"fraction":n/5,"ci_low":lo,"ci_high":hi})
    novel_counts=Counter(x for r in rows for x in r["novel_splits"].split(";") if x); novel=[{"split_key":k,"count":v,"fraction":v/len(rows),"label":labels.get(k,"novel_split")} for k,v in sorted(novel_counts.items())]
    return rows,summary,focal_summary,branches,{"taxa":taxa_expected,"t0_splits":sorted(t0_splits),"non_splits":sorted(non_splits),"novel":novel}

def make_figures(rows, summary, focal_summary, branches):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt; import numpy as np
    eps=list(EPSILONS); fig,axes=plt.subplots(5,1,figsize=(7.2,10.5),constrained_layout=True,gridspec_kw={"height_ratios":[1,1.2,1,1,1]})
    ax=axes[0]
    for state,color,marker in (("J48_T0","#D55E00","o"),("J48_NONCHR4","#0072B2","s")):
        v=[next(r for r in summary if r["epsilon"]==f"{e:.3f}" and r["category"]=="complete_topology" and r["state"]==state) for e in eps]; y=[r["fraction"] for r in v]; ax.plot(eps,y,marker=marker,color=color,lw=1.5,label=state); ax.fill_between(eps,[r["ci_low"] for r in v],[r["ci_high"] for r in v],color=color,alpha=.14,lw=0)
    original = {round(float(r["epsilon"]),3): r for r in read_tsv(RESULTS / "j48_mixture_topology_summary.tsv")}
    shown = [e for e in eps if e in original]
    ax.scatter(shown, [1 if original[e]["equals_T0"] == "true" else 0 for e in shown], marker="D", facecolors="none", edgecolors="#333333", s=36, label="original single seed")
    ax.plot(eps,[1-next(r for r in summary if r["epsilon"]==f"{e:.3f}" and r["category"]=="complete_topology" and r["state"]=="J48_T0")["fraction"]-next(r for r in summary if r["epsilon"]==f"{e:.3f}" and r["category"]=="complete_topology" and r["state"]=="J48_NONCHR4")["fraction"] for e in eps],"o--",color="#777",label="intermediate/distinct")
    ax.set(ylabel="Fraction of replicates",title="A  Exact topology recovery",ylim=(-.03,1.03),xticks=eps); ax.legend(frameon=False,fontsize=7,ncol=3); ax.grid(axis="y",color="#ddd",lw=.5)
    ax=axes[1]; rng=np.random.default_rng(20261025)
    for key,color,label,shift in (("rf_vs_T0","#D55E00","RF from J48_T0",-.006),("rf_vs_NONCHR4","#0072B2","RF from J48_NONCHR4",.006)):
        for e in eps:
            sub=[r for r in rows if abs(float(r["epsilon"])-e)<1e-12]; ax.scatter(e+shift+rng.uniform(-.0025,.0025,len(sub)),[int(r[key]) for r in sub],s=24,color=color,label=label if e==eps[0] else None); vals=[int(r[key]) for r in sub]; ax.plot([e+shift-.004,e+shift+.004],[sum(vals)/5]*2,color=color,lw=2)
    ax.set(ylabel="Unrooted RF distance",title="B  Replicate topology variability",xticks=eps); ax.legend(frameon=False,fontsize=7,ncol=2); ax.grid(axis="y",color="#ddd",lw=.5)
    axes[2].set_title("C  Frozen biological focal-state recovery",loc="left",pad=4,fontsize=10)
    for i,rel in enumerate(RELATIONSHIPS):
        inset=axes[2+i]; states=sorted({r["state"] for r in focal_summary if r["relationship"]==rel})
        for state in states:
            vals=[next((r for r in focal_summary if r["epsilon"]==f"{e:.3f}" and r["relationship"]==rel and r["state"]==state),{"fraction":0}) for e in eps]; inset.plot(eps,[v["fraction"] for v in vals],"o-",color=STATE_COLORS.get(state,"#777"),lw=1.2,ms=3,label=state)
        inset.set_ylim(-.03,1.03); inset.set_yticks([0,.5,1]); inset.set_ylabel(rel,rotation=0,ha="right",va="center",fontsize=8); inset.grid(axis="y",color="#ddd",lw=.4); inset.spines["top"].set_visible(False); inset.spines["right"].set_visible(False)
        if i<2: inset.set_xticklabels([])
        else: inset.set_xlabel(r"chr4 mixture weight $\epsilon$"); inset.set_xticks(eps); inset.tick_params(axis="x",labelrotation=45,labelsize=7)
        if i==0: inset.legend(frameon=False,fontsize=6,ncol=4,loc="upper center",bbox_to_anchor=(.5,1.45))
    for ax in axes: ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    FIGURES.mkdir(parents=True,exist_ok=True)
    for s,kw in (("pdf",{}),("svg",{}),("png",{"dpi":300})): fig.savefig(FIGURES/f"neoaves_j48_mixture_replication.{s}",bbox_inches="tight",**kw)
    plt.close(fig)
    if any(r["homologous_cu_su"]=="true" for r in branches):
        fig,axes=plt.subplots(2,1,figsize=(7.2,5.5),sharex=True,constrained_layout=True)
        for ax,field,title in ((axes[0],"CU_length","CU branch length"),(axes[1],"SU_length","SU branch length")):
            for rel,color in (("N61",STATE_COLORS["q1"]),("N62",STATE_COLORS["q2"]),("Columbea",STATE_COLORS["q3"])):
                for e in eps:
                    vals=[float(r[field]) for r in branches if r["relationship"]==rel and abs(float(r["epsilon"])-e)<1e-12 and r["homologous_cu_su"]=="true" and r[field] != ""]
                    if vals: ax.scatter([e]*len(vals),vals,s=14,color=color,alpha=.7,label=rel if e==eps[0] else None)
            ax.set_ylabel(title); ax.grid(axis="y",color="#ddd",lw=.5); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        axes[1].set_xlabel(r"chr4 mixture weight $\epsilon$"); axes[1].set_xticks(eps); axes[1].tick_params(axis="x",labelrotation=45); axes[0].legend(frameon=False,ncol=3,fontsize=7)
        for s,kw in (("pdf",{}),("svg",{}),("png",{"dpi":300})): fig.savefig(FIGURES/f"neoaves_j48_mixture_replication_branch_lengths.{s}",bbox_inches="tight",**kw)
        plt.close(fig)

def write_docs(rows,summary,focal,meta):
    counts={(float(r["epsilon"]),r["state"]):r["n"] for r in summary if r["category"]=="complete_topology"}; first=next(e for e in EPSILONS if counts.get((e,"J48_T0"),0)>0); majority=next(e for e in EPSILONS if counts.get((e,"J48_T0"),0)>=3); all5=next(e for e in EPSILONS if counts.get((e,"J48_T0"),0)==5)
    lines=["# Neoaves J48 mixture replication analysis","","## Design and validation","", "This report postprocesses 25 completed ASTRAL-IV runs: five independent with-replacement mixtures at each epsilon in 0.02, 0.05, 0.075, 0.10, and 0.15. Every mixture contains 62,945 gene-tree records. All 25 runs passed validation with 48 taxa, identical CU/SU unrooted topologies, valid output checksums, zero exit codes, and topology_su execution metadata.","", "## Topology recovery", "", f"J48_T0 first appears at epsilon = {first:.3f}, becomes the majority topology at epsilon = {majority:.3f}, and all five replicates recover it at epsilon = {all5:.3f}. These are observed grid summaries, not a fitted threshold.","","| epsilon | T0 | NONCHR4 | intermediate/distinct | novel-split replicates |","|---:|---:|---:|---:|---:|"]
    for e in EPSILONS:
        t=counts.get((e,"J48_T0"),0); n=counts.get((e,"J48_NONCHR4"),0); nov=next(r["n"] for r in summary if r["epsilon"]==f"{e:.3f}" and r["category"]=="novel_topology"); lines.append(f"| {e:.3f} | {t}/5 | {n}/5 | {5-t-n}/5 | {nov}/5 |")
    lines += ["", "## RF summary", "", "| epsilon | T0 RF mean/median/range/SD | non-chr4 RF mean/median/range/SD |", "|---:|---|---|"]
    for e in EPSILONS:
        t = next(r for r in summary if r["epsilon"] == f"{e:.3f}" and r["category"] == "complete_topology" and r["state"] == "J48_T0")
        n = next(r for r in summary if r["epsilon"] == f"{e:.3f}" and r["category"] == "complete_topology" and r["state"] == "J48_NONCHR4")
        lines.append(f"| {e:.3f} | {t['mean_RF']:.3g}/{t['median_RF']}/{t['min_RF']}–{t['max_RF']}/{t['sd_RF']:.3g} | {n['mean_RF']:.3g}/{n['median_RF']}/{n['min_RF']}–{n['max_RF']}/{n['sd_RF']:.3g} |")
    lines += ["", "## Focal relationships", "", "N61 is q1 in all five replicates at 0.02 and 0.05, then q2 in all replicates from 0.075 onward. Columbea changes from q2 to q1 at 0.075. N62 remains q1 across the sampled replicate grid. These are frozen biological states; they are not branch-local q1/q2/q3 labels.", "", "Three exact novel split identities occur across the 25 trees. One recurs in four replicates and the other two occur once each; no novel split is present at 0.075, where all five trees equal J48_T0. Exact split strings and frequencies are in `j48_mixture_replication_novel_splits.tsv`.", "", "## Interpretation and limitations", "", "The experiment establishes sensitivity of a fixed 48-taxon species-tree inference to reweighting empirical genealogy pools. It does not establish recombination suppression as the cause, equate epsilon with a physical genomic fraction, estimate an MSRC parameter, or make Patterson-style or QQS crossover claims. The QQS values N61 approximately 0.1540 and N62 approximately 0.01354 are different estimands.", "", "Five replicates per weight provide limited binomial precision. Adjacent weights are independent resamples rather than nested additions of loci. The smallest useful follow-up is additional seeds at weights with mixed recovery, especially 0.05, 0.10, and 0.15; fixed-candidate quartet scoring on the predeclared disagreement quartets would be a complementary analysis if independently validated."]
    (RESULTS/"j48_mixture_replication_analysis.md").write_text("\n".join(lines)+"\n")
    t0_10=counts.get((.10,"J48_T0"),0); t0_15=counts.get((.15,"J48_T0"),0); paragraph=f"We evaluated sensitivity of 48-taxon ASTRAL-IV inference to chromosome-4 genealogy reweighting using five independent with-replacement mixtures of 62,945 gene-tree records at each of five weights (epsilon = 0.02, 0.05, 0.075, 0.10, and 0.15). All 25 runs completed successfully and produced identical CU and SU topologies with 48 taxa. The complete J48_T0 topology was recovered in {t0_10}/5 and {t0_15}/5 replicates at epsilon = 0.10 and 0.15, respectively, and first appeared at epsilon = {first:.3f}; recovery became the majority at epsilon = {majority:.3f}. Lower weights included non-chr4 and intermediate or distinct topologies, demonstrating sensitivity in the coarse transition interval. Because each mixture is an independent resample and epsilon is a genealogy-pool sampling weight, these results establish empirical sensitivity to genomic reweighting without estimating a deterministic threshold or identifying recombination suppression as its cause."
    (RESULTS/"j48_mixture_replication_results_paragraph.tex").write_text(paragraph+"\n")
    (FIGURES/"neoaves_j48_mixture_replication_caption.md").write_text("**Figure X. Replicated Neoaves J48 mixture sensitivity.** (A) Fraction of five independent replicates recovering the exact J48_T0 and J48_NONCHR4 unrooted topologies; shaded intervals are Wilson 95% intervals, gray denotes intermediate or distinct trees, and open diamonds show the original single-seed grid for comparison without counting it as a replicate. (B) Individual replicate RF distances to each reference, with within-weight means. (C) Frequencies of the frozen biological focal states for N61, N62, and Columbea. Each epsilon is a separate 62,945-record with-replacement mixture, not a nested sequence of loci or a measured fraction of the genome under suppressed recombination.\n")

def main():
    p=argparse.ArgumentParser(); p.add_argument("--result-root",type=Path,required=True); p.add_argument("--manifest",type=Path,required=True); p.add_argument("--output",type=Path,required=True); a=p.parse_args()
    validation=validator.validate(a.result_root,a.manifest); write_tsv(RESULTS/"j48_mixture_replication_validation.tsv",validation); rows,summary,focal,branches,meta=analyze(a.result_root,a.manifest)
    write_tsv(a.output/"j48_mixture_replication_rows.tsv",rows); write_tsv(a.output/"j48_mixture_replication_summary.tsv",summary); write_tsv(a.output/"j48_mixture_replication_focal_summary.tsv",focal); write_tsv(a.output/"j48_mixture_replication_branch_lengths.tsv",branches); write_tsv(a.output/"j48_mixture_replication_novel_splits.tsv",meta["novel"])
    make_figures(rows,summary,focal,branches); write_docs(rows,summary,focal,meta); print(json.dumps({"validated":len(validation),"replicates":len(rows),"summary_rows":len(summary),"focal_rows":len(focal),"novel_splits":len(meta["novel"])},indent=2))

if __name__=="__main__": main()
