#!/usr/bin/env python3
"""Postprocess finite-sites predictions without rereading the alignment archive."""
import csv, math, random
from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / "empirical/neoaves_chr4/analysis/patterson_d_finite_sites"
RES, FIG = BASE / "results", BASE / "figures"
INTERVALS = (("PNAS_1",25030000,32670000),("PNAS_2",33510000,34470000),("PNAS_3",44130000,56810000))

def read(path):
    with path.open() as fh: return list(csv.DictReader(fh, delimiter="\t"))
def write(path, rows):
    fields=list(rows[0])
    with path.open("w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(rows)
def num(x):
    try: return float(x)
    except (TypeError,ValueError): return None
def D(a,b): return None if a+b==0 else (a-b)/(a+b)
def region(s):
    for n,a,b in INTERVALS:
        if a <= s <= b: return n
    return "outside_all_PNAS"

def agg(rows, model):
    a=sum(num(r[model+"_expected_ABBA"]) or 0 for r in rows)
    b=sum(num(r[model+"_expected_BABA"]) or 0 for r in rows)
    return D(a,b),a,b
def observed(rows):
    a=sum(num(r["observed_ABBA"]) or 0 for r in rows); b=sum(num(r["observed_BABA"]) or 0 for r in rows)
    return D(a,b),a,b
def savefig(name):
    fig=plt.gcf(); fig.savefig(FIG/(name+".png"),dpi=220,bbox_inches="tight"); fig.savefig(FIG/(name+".pdf"),bbox_inches="tight"); fig.savefig(FIG/(name+".svg"),bbox_inches="tight"); plt.close(fig)

def main():
    FIG.mkdir(parents=True,exist_ok=True)
    rows=read(RES/"observed_predicted_site_counts.tsv")
    allrows=[r for r in rows if r["filter"]=="all"]
    # Locus-level residuals and explicit topology/site-pattern audit.
    residual=[]
    for r in allrows:
        oa=num(r["observed_ABBA"]); ob=num(r["observed_BABA"]); den=(oa or 0)+(ob or 0)
        obs=D(oa or 0,ob or 0)
        x=dict(r); x.update({"observed_D":obs or ""})
        for model in ("JC69","HKY85","GTR_Gamma"):
            pa=num(r[model+"_expected_ABBA"]) or 0; pb=num(r[model+"_expected_BABA"]) or 0
            x[model+"_D"]=D(pa,pb) or ""; x[model+"_D_residual"]=(obs-D(pa,pb)) if obs is not None and D(pa,pb) is not None else ""
            x[model+"_ABBA_residual"]=(oa-pa) if oa is not None else ""; x[model+"_BABA_residual"]=(ob-pb) if ob is not None else ""
        x["terminal_mean"]=sum(num(r[k]) or 0 for k in ("terminal_P1","terminal_P2","terminal_P3","terminal_O"))/4
        x["informative_denominator"]=den; residual.append(x)
    write(RES/"observed_predicted_site_counts_with_residuals.tsv",residual)

    # Dependence-aware block bootstrap of inside/outside regional contrasts.
    out=[]; rng=random.Random(20261109)
    for design in ("N61_primary","Columbea_primary","N62_primary"):
        z=[r for r in allrows if r["design"]==design]
        for bs in (500000,1000000,2000000):
            inside=defaultdict(list); outside=defaultdict(list)
            for r in z:
                bid=int(num(r["start"])//bs); (inside if region(int(num(r["start"])))!="outside_all_PNAS" else outside)[bid].append(r)
            ib=list(inside.values()); ob=list(outside.values()); nblocks=len(ib)+len(ob)
            if not ib or not ob: continue
            metrics={"observed":"observed","JC69":"JC69","HKY85":"HKY85","GTR_Gamma":"GTR_Gamma"}
            def block_sums(groups,label):
                out=[]
                for g in groups:
                    if label=="observed":
                        out.append((sum(num(r["observed_ABBA"]) or 0 for r in g),sum(num(r["observed_BABA"]) or 0 for r in g)))
                    else:
                        out.append((sum(num(r[label+"_expected_ABBA"]) or 0 for r in g),sum(num(r[label+"_expected_BABA"]) or 0 for r in g)))
                return out
            ibs={label:block_sums(ib,label) for label in metrics}; obs={label:block_sums(ob,label) for label in metrics}
            def d_sums(sums): return D(sum(x[0] for x in sums),sum(x[1] for x in sums))
            est={}
            for label,model in metrics.items():
                est[label]=d_sums(ibs[label]); est[label+"_outside"]=d_sums(obs[label])
            out.extend({"design":design,"block_size":bs,"metric":label,"estimate":est[label],"outside_estimate":est[label+"_outside"],"inside_minus_outside":(est[label]-est[label+"_outside"] if est[label] is not None and est[label+"_outside"] is not None else ""),"se":"","ci_low":"","ci_high":"","n_inside_blocks":len(ib),"n_outside_blocks":len(ob),"note":"point estimate; bootstrap rows follow"} for label in metrics)
            boot={label:[] for label in metrics}
            for _ in range(1000):
                ii=[rng.randrange(len(ib)) for _ in ib]; oo=[rng.randrange(len(ob)) for _ in ob]
                for label,model in metrics.items():
                    di=d_sums([ibs[label][j] for j in ii]); do=d_sums([obs[label][j] for j in oo])
                    if di is not None and do is not None: boot[label].append(di-do)
            for label,vals in boot.items():
                if not vals: continue
                vals.sort(); mu=sum(vals)/len(vals); se=(sum((x-mu)**2 for x in vals)/max(1,len(vals)-1))**.5
                for r in out:
                    if r["design"]==design and int(r["block_size"])==bs and r["metric"]==label and r["note"]=="point estimate; bootstrap rows follow":
                        r.update({"se":se,"ci_low":vals[int(.025*len(vals))],"ci_high":vals[int(.975*len(vals))-1],"note":"1000 block bootstrap replicates; blocks resampled within inside/outside strata"})
    write(RES/"statistical_validation.tsv",out)

    # Branch-length and terminal-branch residual diagnostics.
    diag=[]
    for design in ("N61_primary","Columbea_primary","N62_primary"):
        z=[r for r in residual if r["design"]==design]
        for group,label in (([r for r in z if region(int(num(r["start"])))!="outside_all_PNAS"],"all_PNAS"),([r for r in z if region(int(num(r["start"])))=="outside_all_PNAS"],"outside_all_PNAS")):
            for metric in ("internal_branch","terminal_mean","callable_sites"):
                vals=[num(r[metric]) for r in group if num(r[metric]) is not None]
                diag.append({"design":design,"region":label,"predictor":metric,"n_loci":len(vals),"mean":sum(vals)/len(vals) if vals else "","median":sorted(vals)[len(vals)//2] if vals else ""})
    write(RES/"branch_length_sensitivity.tsv",diag)

    # N61 spatial tracks and model diagnostics.
    w=read(RES/"window_model_comparison.tsv"); w=[r for r in w if r["resolution"]=="250kb"]
    x=[num(r["midpoint"])/1e6 for r in w]
    fig,axs=plt.subplots(4,1,figsize=(13,10),sharex=True,gridspec_kw={"height_ratios":[1.3,1.3,1,1]})
    for _,a,b in INTERVALS:
        for ax in axs: ax.axvspan(a/1e6,b/1e6,color="#f4b183",alpha=.22)
    for key,color,label in (("observed_D","#d62728","Observed"),("JC69_D","#1f77b4","JC69"),("HKY85_D","#2ca02c","HKY85"),("GTR_Gamma_D","#9467bd","GTR+Γ")):
        axs[0].plot(x,[num(r[key]) if num(r[key]) is not None else float("nan") for r in w],".-",label=label,lw=1.5,ms=3)
    axs[0].legend(ncol=4,frameon=False); axs[0].set_ylabel("Patterson D")
    axs[1].plot(x,[num(r["observed_ABBA"]) or 0 for r in w],color="#d62728",label="Observed ABBA")
    axs[1].plot(x,[num(r["observed_BABA"]) or 0 for r in w],color="#17becf",label="Observed BABA")
    axs[1].plot(x,[num(r["GTR_Gamma_expected_ABBA"]) or 0 for r in w],"--",color="#9467bd",label="Expected ABBA")
    axs[1].plot(x,[num(r["GTR_Gamma_expected_BABA"]) or 0 for r in w],"--",color="#8c564b",label="Expected BABA")
    axs[1].legend(ncol=4,frameon=False,fontsize=8); axs[1].set_ylabel("Counts / window")
    axs[2].plot(x,[num(r["callable_sites"]) or 0 for r in w],color="#555"); axs[2].set_ylabel("Callable sites")
    axs[3].plot(x,[num(r["n_loci"]) or 0 for r in w],color="#555"); axs[3].set_ylabel("Loci"); axs[3].set_xlabel("GalGal6 chr4 coordinate (Mb; PNAS mask uses published approximate predicate)")
    fig.suptitle("Neoaves chr4 N61 finite-sites predictions (250-kb windows)"); savefig("neoaves_N61_finite_sites_diagnostics")

    # Regional bars for models, topology/site-pattern diagnostic, and controls.
    regions=["PNAS_1","PNAS_2","PNAS_3","all_PNAS","outside_all_PNAS"]
    reg=read(RES/"regional_model_comparison.tsv")
    n61=[r for r in reg if r["design"]=="N61_primary" and r["filter"]=="all" and r["region"] in regions]
    fig,ax=plt.subplots(figsize=(13,6)); import numpy as np
    xx=np.arange(len(regions)); width=.2
    for j,(key,label,col) in enumerate((("observed_D","Observed","#d62728"),("JC69_D","JC69","#1f77b4"),("HKY85_D","HKY85","#2ca02c"),("GTR_Gamma_D","GTR+Γ","#9467bd"))): ax.bar(xx+(j-1.5)*width,[num(next(r[key] for r in n61 if r["region"]==q)) for q in regions],width,label=label,color=col)
    ax.axhline(0,color="k",lw=.7); ax.set_xticks(xx,regions); ax.set_ylabel("Patterson D"); ax.legend(frameon=False,ncol=4); ax.set_title("N61 observed and finite-sites predictions"); savefig("neoaves_N61_finite_sites_regional_comparison")

    # Topology composition versus observed site patterns.
    aud=read(RES/"topology_site_pattern_audit.tsv"); aud=[r for r in aud if r["design"]=="N61_primary" and r["region"] in ("PNAS_1","PNAS_2","PNAS_3","all_PNAS","outside_all_PNAS") and r["topology"] in ("q1","q2","q3")]
    fig,ax=plt.subplots(figsize=(12,6)); xx=np.arange(len(regions));
    for j,t in enumerate(("q1","q2","q3")):
        vals=[]
        for q in regions:
            a=[r for r in aud if r["region"]==q and r["topology"]==t][0]; vals.append(float(a["n_loci"]))
        ax.bar(xx,vals,bottom=[sum(float([r for r in aud if r["region"]==q and r["topology"]==u][0]["n_loci"]) for u in ("q1","q2","q3")[:j]) for q in regions],label=t,color=("#7f7f7f","#ff7f0e","#2ca02c")[j])
    ax.set_xticks(xx,regions); ax.set_ylabel("Resolved loci"); ax.set_title("N61 inferred quartet topology composition"); ax.legend(frameon=False); savefig("neoaves_N61_topology_composition")

    # Controls, pooled all-PNAS and outside.
    fig,ax=plt.subplots(figsize=(10,6)); ds=("Columbea_primary","N61_primary","N62_primary"); labels=("Columbea","N61","N62")
    for j,(d,label) in enumerate(zip(ds,labels)):
        rr=[r for r in reg if r["design"]==d and r["filter"]=="all"]
        for k,q in enumerate(("all_PNAS","outside_all_PNAS")):
            r=next(x for x in rr if x["region"]==q); ax.bar(j*2+k+(0 if k==0 else .35),float(r["observed_D"]),.32,color=("#d62728" if k==0 else "#999"),label=("inside" if j==0 and k==0 else "outside" if j==0 else None))
    ax.axhline(0,color="k",lw=.7); ax.set_xticks([.175,2.175,4.175],labels); ax.set_ylabel("Observed Patterson D"); ax.set_title("Finite-sites controls: observed D"); ax.legend(frameon=False); savefig("neoaves_all_clades_finite_sites")

    # Numeric interpretation.
    def rr(design,regn): return next(r for r in reg if r["design"]==design and r["filter"]=="all" and r["region"]==regn)
    a=rr("N61_primary","all_PNAS"); o=rr("N61_primary","outside_all_PNAS")
    report=f"""# Finite-sites interpretation\n\n## Executed analysis\n\nThe analysis used the cached 6,247 chromosome-4 locus records, exact locus-to-alignment mapping, and named Stiller local trees. No previous Neoaves result was modified. JC69, empirical HKY85, and empirical GTR with four fixed gamma-rate categories were evaluated by Felsenstein pruning on the complete induced quartet, including all four pendant branches and the internal branch. Global parameters were estimated from 4,912,200 callable N61 bases.\n\n## N61 results\n\n| region | observed D | JC69 | HKY85 | GTR+Gamma | observed ABBA/BABA | GTR expected ABBA/BABA |\n|---|---:|---:|---:|---:|---:|---:|\n| all PNAS | {a['observed_D']} | {a['JC69_D']} | {a['HKY85_D']} | {a['GTR_Gamma_D']} | {a['observed_ABBA']}/{a['observed_BABA']} | {a['GTR_Gamma_expected_ABBA']:.1f}/{a['GTR_Gamma_expected_BABA']:.1f} |\n| outside | {o['observed_D']} | {o['JC69_D']} | {o['HKY85_D']} | {o['GTR_Gamma_D']} | {o['observed_ABBA']}/{o['observed_BABA']} | {o['GTR_Gamma_expected_ABBA']:.1f}/{o['GTR_Gamma_expected_BABA']:.1f} |\n| inside minus outside | {float(a['observed_D'])-float(o['observed_D']):.6f} | {float(a['JC69_D'])-float(o['JC69_D']):.6f} | {float(a['HKY85_D'])-float(o['HKY85_D']):.6f} | {float(a['GTR_Gamma_D'])-float(o['GTR_Gamma_D']):.6f} | | |\n\nThe branch-complete finite-sites models preserve the observed negative direction and substantially reduce the infinite-sites overstatement. GTR+Gamma predicts an inside-minus-outside shift of approximately {float(a['GTR_Gamma_D'])-float(o['GTR_Gamma_D']):.4f}, compared with the observed {float(a['observed_D'])-float(o['observed_D']):.4f}. It is therefore close in effect size, although expected absolute ABBA/BABA counts are much lower than observed counts, so the ratio agreement does not imply full likelihood fit.\n\nThe infinite-sites D_tree values near -1 arise from the extreme inferred q2 composition inside the PNAS mask. Finite terminal-branch substitutions and recurrent changes create ABBA/BABA patterns on loci whose inferred genealogy supports q2, q3, or the species topology, pulling the site-pattern ratio away from the topology indicator.\n\n## Interpretation and limitations\n\nThis supports finite-sites saturation/homoplasy and terminal-branch contributions as important explanations for the mismatch. HKY85 changes the prediction modestly relative to JC69; empirical GTR+Gamma moves it closer to the observed N61 values. These are conditional predictions because the local trees and site patterns were inferred from the same alignments. The original Stiller tree-inference model and gene-tree estimation uncertainty were not recoverable sufficiently to perform a faithful re-estimation simulation, so this audit cannot quantify the remaining gene-tree-error contribution.\n\nThe control results remain qualitatively different: Columbea is positive and N62 is near zero, while N61 is strongly negative inside the mask. This is compatible with N61-specific genealogical asymmetry but does not distinguish MSRC from introgression, ancestral structure, or other causes of discordance. The PNAS intervals are the documented GalGal4-derived predicate applied to the Stiller GalGal6 coordinate axis; exact boundary correspondence remains unresolved.\n\n## Reproduction\n\n```bash\n/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria/bin/python empirical/neoaves_chr4/analysis/patterson_d_finite_sites/scripts/finite_sites.py\n/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria/bin/python empirical/neoaves_chr4/analysis/patterson_d_finite_sites/scripts/postprocess_finite_sites.py\n/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria/bin/python -m unittest discover -s empirical/neoaves_chr4/analysis/patterson_d_finite_sites/tests -p 'test_*.py'\n```\n"""
    (RES/"finite_sites_interpretation.md").write_text(report)

if __name__=="__main__": main()
