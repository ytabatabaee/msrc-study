#!/usr/bin/env python3
"""Audit Martin & Van Belleghem's deposited Neurospora windows and quantify shifts."""
from __future__ import annotations

import csv
import gzip
import io
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / "data/neurospora/raw/neurospora"
OUT = ROOT / "empirical/neurospora/results"
FIG = ROOT / "empirical/neurospora/figures"
DATA = RAW / "trees/cras_sit_hist_tet_L10.min20MAC1.phyml_GTR_bionj_w50sites.data.tsv.gz"
TREES = RAW / "trees/cras_sit_hist_tet_L10.min20MAC1.phyml_GTR_bionj_w50sites.trees.gz"
WEIGHTS = RAW / "cras_his_L10.min20MAC1.phyml_GTR_bionj_w50sites.weights.tsv.gz"
BOUNDARY = (1_500_000, 8_000_000)  # L10 mat A/a SR boundaries in 2509a reference coordinates
QUARTETS = [("CJ03", "CJ04"), ("CJ15", "CJ16")]
OUTGROUPS = ("Ncra@999", "8817@999")


def read_tsv_gz(path):
    with gzip.open(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def normalize_weights(vals):
    if not all(math.isfinite(v) for v in vals):
        return [float("nan")] * len(vals)
    total = sum(vals)
    if total <= 0:
        raise ValueError("Nonpositive TWISST count total")
    return [v / total for v in vals]


def read_weights(path):
    # First three lines document the topology labels; fourth is the column header.
    with gzip.open(path, "rt") as fh:
        lines = fh.readlines()
    rows = list(csv.DictReader(lines[3:], delimiter="\t"))
    parsed = []
    for r in rows:
        vals = [float(r[k]) if r[k] != "NA" else float("nan") for k in ("topo1", "topo2", "topo3")]
        vals = normalize_weights(vals)
        parsed.append(dict(zip(("q1", "q2", "q3"), vals)))
    return parsed


def classify_split(tree, a, b, c, d):
    """Return quartet split index: q1=(a,b|c,d), q2=(a,c|b,d), q3=(a,d|b,c)."""
    taxa = {t.name for t in tree.get_terminals()}
    if not {a, b, c, d} <= taxa:
        return None
    dm = quartet_distances(tree, (a, b, c, d))
    sums = [dm[a,b] + dm[c,d], dm[a,c] + dm[b,d], dm[a,d] + dm[b,c]]
    return min(range(3), key=sums.__getitem__), max(0.0, (sorted(sums)[1] - sorted(sums)[0]) / 2.0)


def quartet_distances(tree, names):
    """Compute all focal pairwise distances after only one root-to-tip walk per tip."""
    wanted = set(names); paths = {}
    def visit(node, path, depth):
        current = path + [(node, depth)]
        if node.is_terminal() and node.name in wanted:
            paths[node.name] = current
        for child in node.clades:
            visit(child, current, depth + (child.branch_length or 0.0))
    visit(tree.root, [], 0.0)
    result = {}
    for i, x in enumerate(names):
        for y in names[i+1:]:
            px = paths[x]; py = paths[y]; j = 0
            while j < min(len(px), len(py)) and px[j][0] is py[j][0]: j += 1
            lca_depth = px[j-1][1] if j else 0.0
            result[x,y] = result[y,x] = px[-1][1] + py[-1][1] - 2*lca_depth
    return result


def summarise(rows, label):
    rows = [r for r in rows if all(math.isfinite(r[k]) for k in ("q1", "q2", "q3"))]
    if not rows:
        return None
    q = np.array([[r["q1"], r["q2"], r["q3"]] for r in rows], dtype=float)
    # Window-weighted concordance vector; genomic-span estimate weights each window by its span.
    span = np.array([r["end"] - r["start"] + 1 for r in rows], dtype=float)
    mean = q.mean(axis=0)
    spanmean = np.average(q, axis=0, weights=span)
    discordant = sorted((mean[1], mean[2]))
    cu = -math.log(max(1e-12, 1.5 * sum(discordant))) if mean[0] >= max(mean[1], mean[2]) and abs(mean[1] - mean[2]) < .02 and 1.5 * sum(discordant) < 1 else float("nan")
    return {"subset": label, "n_windows": len(rows), "span_bp": int(span.sum()),
            "q1_species": mean[0], "q2_other": mean[1], "q3_hisp_A": mean[2],
            "D_q3_minus_q1": mean[2] - mean[0], "q1_span_weighted": spanmean[0],
            "q2_span_weighted": spanmean[1], "q3_span_weighted": spanmean[2],
            "dominant_topology": f"q{int(np.argmax(mean)) + 1}", "MSC_CU_if_symmetric": cu,
            "discordant_asymmetry_abs": abs(mean[1] - mean[2])}


def write_tsv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k:(format(v,".8g") if isinstance(v,(float,np.floating)) else v) for k,v in row.items()})


def main():
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    windows = read_tsv_gz(DATA)
    weights = read_weights(WEIGHTS)
    if len(windows) != len(weights):
        raise ValueError(f"Metadata/weight row mismatch: {len(windows)} vs {len(weights)}")
    for i, (r, q) in enumerate(zip(windows, weights), 1):
        r["window"] = i
        r.update({k: int(r[k]) for k in ("start", "end", "mid", "sites")})
        r.update(q); r["D"] = q["q3"] - q["q1"]
    inventory = []
    summaries = []
    inside = lambda r: r["scaffold"] == "LGI" and BOUNDARY[0] <= r["mid"] <= BOUNDARY[1]
    masks = {
        "all": lambda r: True,
        "introgressed_SR_LGI_1.5_8.0Mb": inside,
        "outside_SR": lambda r: not inside(r),
        "LGI_flank_250kb": lambda r: r["scaffold"] == "LGI" and ((1_250_000 <= r["mid"] < BOUNDARY[0]) or (BOUNDARY[1] < r["mid"] <= 8_250_000)),
        "LGI_other": lambda r: r["scaffold"] == "LGI" and not inside(r),
        "other_linkage_groups": lambda r: r["scaffold"] != "LGI",
    }
    for name, predicate in masks.items():
        summaries.append(summarise([r for r in windows if predicate(r)], name))
    write_tsv(OUT / "region_comparison.tsv", [x for x in summaries if x])

    # Spatial inventory contains all metadata and the deposited TWISST frequency vector.
    inventory = [{"window": i + 1, "source_filename": TREES.name, "linkage_group": r["scaffold"], "start": r["start"], "end": r["end"],
                  "mid": r["mid"], "sites": r["sites"], "q1_species": r["q1"], "q2_other": r["q2"],
                  "q3_hisp_A": r["q3"], "D_q3_minus_q1": r["D"],
                  "inside_LGI_1p5_8Mb": int(inside(r))} for i, r in enumerate(windows)]
    write_tsv(OUT / "gene_tree_inventory.tsv", inventory)

    # One independent focal quartet per deposited window for branch-length contrasts.
    # CJ03/CJ04 are an independently documented matched L10 A/a homokaryon pair;
    # CJ15/CJ16 replicate that choice. Outgroups are N. crassa and N. hispaniola.
    branch_rows = []; tree_quality = []; observed_tipsets = {}
    with gzip.open(TREES, "rt") as fh:
        for i, nwk in enumerate(fh):
            if i >= len(windows):
                raise ValueError("More trees than metadata rows")
            r = windows[i]
            nwk = nwk.strip()
            if nwk in ("NA", ""):
                tree_quality.append({"window": i + 1, "linkage_group": r["scaffold"], "start": r["start"], "end": r["end"], "mid": r["mid"],
                                     "tip_set_id": "NA", "tip_count": 0, "missing_nonroot_branch_lengths": "NA",
                                     "polytomy_nodes": "NA", "tree_type": "NA_no_tree"})
                branch_rows.append({"window": i + 1, "scaffold": r["scaffold"], "linkage_group": r["scaffold"], "start": r["start"], "end": r["end"], "mid": r["mid"],
                                   "inside_LGI_1p5_8Mb": int(inside(r)), "published_q1":"NA", "published_q2":"NA", "published_q3":"NA",
                                   "pair1_available": 0, "pair1_A_hisp_minus_a_hisp":"NA", "pair1_B_mean_six_pairwise":"NA", "pair1_internal_branch":"NA", "pair1_tree_split":"NA",
                                   "pair2_available": 0, "pair2_A_hisp_minus_a_hisp":"NA", "pair2_B_mean_six_pairwise":"NA", "pair2_internal_branch":"NA", "pair2_tree_split":"NA"})
                continue
            tree = Phylo.read(io.StringIO(nwk), "newick")
            tips = [t.name for t in tree.get_terminals()]
            tipset=tuple(sorted(tips)); observed_tipsets.setdefault(tipset,len(observed_tipsets)+1); tipset_id=observed_tipsets[tipset]
            internal = [c for c in tree.get_nonterminals() if c is not tree.root]
            tree_quality.append({"window": i + 1, "linkage_group": r["scaffold"], "start": r["start"], "end": r["end"], "mid": r["mid"],
                                 "tip_set_id": tipset_id, "tip_count": len(tips),
                                 "missing_nonroot_branch_lengths": sum(c.branch_length is None for c in tree.find_clades() if c is not tree.root),
                                 "polytomy_nodes": sum(len(c.clades) > 2 for c in internal),
                                 "tree_type": "unrooted_BIONJ"})
            row = {"window": i + 1, "scaffold": r["scaffold"], "linkage_group": r["scaffold"], "start": r["start"], "end": r["end"], "mid": r["mid"],
                   "inside_LGI_1p5_8Mb": int(inside(r)), "published_q1": r["q1"], "published_q2": r["q2"], "published_q3": r["q3"]}
            for pair_i, (a, b) in enumerate(QUARTETS, 1):
                tips = {t.name for t in tree.get_terminals()}
                if not {a, b, *OUTGROUPS} <= tips:
                    row[f"pair{pair_i}_available"] = 0
                    continue
                dm = quartet_distances(tree, (a,b,*OUTGROUPS))
                d_ah = dm[a,OUTGROUPS[1]]; d_bh = dm[b,OUTGROUPS[1]]
                ds = [dm[x,y] for x,y in ((a,b),(a,OUTGROUPS[0]),(a,OUTGROUPS[1]),(b,OUTGROUPS[0]),(b,OUTGROUPS[1]),OUTGROUPS)]
                split = classify_split(tree, a, b, *OUTGROUPS)
                row[f"pair{pair_i}_available"] = 1
                row[f"pair{pair_i}_A_hisp_minus_a_hisp"] = d_ah - d_bh
                row[f"pair{pair_i}_B_mean_six_pairwise"] = float(np.mean(ds))
                row[f"pair{pair_i}_internal_branch"] = split[1] if split else float("nan")
                row[f"pair{pair_i}_tree_split"] = f"q{split[0]+1}" if split else "NA"
            branch_rows.append(row)
    if len(branch_rows) != len(windows):
        raise ValueError(f"Tree/metadata row mismatch: {len(branch_rows)} vs {len(windows)}")
    write_tsv(OUT / "branch_length_inventory.tsv", branch_rows)
    write_tsv(OUT / "tree_quality_inventory.tsv", tree_quality)
    write_tsv(OUT / "tree_taxon_sets.tsv", [{"tip_set_id":i,"tip_count":len(tips),"tip_labels":";".join(tips)} for tips,i in observed_tipsets.items()])

    # Paired quartet branch summaries by genomic subset.
    branch_summary = []
    for pair_i in (1, 2):
        for name, predicate in masks.items():
            subset = [r for r in branch_rows if predicate(r) and r.get(f"pair{pair_i}_available") == 1]
            if not subset: continue
            for metric in (f"pair{pair_i}_A_hisp_minus_a_hisp", f"pair{pair_i}_B_mean_six_pairwise", f"pair{pair_i}_internal_branch"):
                vals = np.array([float(r[metric]) for r in subset])
                branch_summary.append({"pair": f"CJ{3 if pair_i==1 else 15:02d}_CJ{4 if pair_i==1 else 16:02d}", "subset": name, "metric": metric,
                                      "n_windows": len(vals), "mean": float(np.mean(vals)), "median": float(np.median(vals)),
                                      "sd": float(np.std(vals, ddof=1)) if len(vals)>1 else 0.0})
            for split in ("q1", "q2", "q3"):
                vals = [r.get(f"pair{pair_i}_tree_split") == split for r in subset]
                branch_summary.append({"pair": f"CJ{3 if pair_i==1 else 15:02d}_CJ{4 if pair_i==1 else 16:02d}", "subset": name,
                                       "metric": f"tree_topology_{split}_frequency", "n_windows": len(vals),
                                       "mean": float(np.mean(vals)), "median": float(np.mean(vals)), "sd": float("nan")})
    write_tsv(OUT / "branch_length_summary.tsv", branch_summary)

    # One 100-kb grid value per block for the LGI autocorrelation and spatial transition summaries.
    lgi = [r for r in windows if r["scaffold"] == "LGI"]
    bins = {}
    for r in lgi:
        bins.setdefault(r["mid"] // 100_000, []).append(r)
    safe_mean = lambda vals: float(np.mean([v for v in vals if math.isfinite(v)])) if any(math.isfinite(v) for v in vals) else float("nan")
    block = [{"start": k*100_000, "end": (k+1)*100_000-1, "n_windows": len(v),
              "q1": safe_mean([r["q1"] for r in v]), "q2": safe_mean([r["q2"] for r in v]),
              "q3": safe_mean([r["q3"] for r in v]), "D": safe_mean([r["D"] for r in v])} for k,v in sorted(bins.items())]
    write_tsv(OUT / "lgi_100kb_blocks.tsv", block)
    adjacent = [(a["D"],b["D"]) for a,b in zip(block[:-1],block[1:]) if b["start"]-a["start"]==100_000 and math.isfinite(a["D"]) and math.isfinite(b["D"])]
    ac = float(np.corrcoef(np.array(adjacent).T)[0,1]) if len(adjacent)>2 else float("nan")

    # Block bootstrap contrasts and equal-count random matching as descriptive uncertainty controls.
    rng = np.random.default_rng(20261115)
    validq = lambda r: all(math.isfinite(r[k]) for k in ("q1", "q2", "q3"))
    inrows = [r for r in windows if inside(r) and validq(r)]; outrows = [r for r in windows if r["scaffold"] == "LGI" and not inside(r) and validq(r)]
    nmatch = min(len(inrows), len(outrows))
    match = []
    for rep in range(2000):
        ii = rng.choice(len(inrows), nmatch, replace=False); oo = rng.choice(len(outrows), nmatch, replace=False)
        match.append(float(np.mean([inrows[j]["D"] for j in ii]) - np.mean([outrows[j]["D"] for j in oo])))
    bootstrap = []
    # Resample contiguous 100-kb blocks within the two groups, avoiding a false iid-window CI.
    inblocks = [[x for x in v if validq(x)] for k,v in bins.items() if 1_500_000 <= k*100_000 <= 8_000_000]
    outblocks = [[x for x in v if validq(x)] for k,v in bins.items() if not (1_500_000 <= k*100_000 <= 8_000_000)]
    inblocks = [x for x in inblocks if x]; outblocks = [x for x in outblocks if x]
    for _ in range(2000):
        i_s = rng.choice(len(inblocks), len(inblocks), replace=True); o_s = rng.choice(len(outblocks), len(outblocks), replace=True)
        iv = np.concatenate([[x["D"] for x in inblocks[i]] for i in i_s]); ov = np.concatenate([[x["D"] for x in outblocks[i]] for i in o_s])
        bootstrap.append(float(iv.mean()-ov.mean()))
    with open(OUT / "sensitivity_summary.tsv", "w") as fh:
        fh.write("analysis\testimate\tci_low\tci_high\tn\n")
        delta = np.mean([r["D"] for r in inrows]) - np.mean([r["D"] for r in outrows])
        fh.write(f"inside_minus_outside_D\t{delta}\t{np.quantile(bootstrap,.025)}\t{np.quantile(bootstrap,.975)}\t2000_block_resamples\n")
        fh.write(f"equal_count_matched_D_difference\t{np.mean(match)}\t{np.quantile(match,.025)}\t{np.quantile(match,.975)}\t2000_random_matches_n={nmatch}\n")
        fh.write(f"LGI_adjacent_100kb_D_autocorrelation\t{ac}\tNA\tNA\t{len(block)}_blocks\n")
        fh.write(f"LGI_taxon_set_variants\t{len(observed_tipsets)}\tNA\tNA\t{len(branch_rows)}_trees\n")

    # Linked-window thinning: one randomly selected informative window per 100 kb bin.
    thin_rows=[]
    for name,predicate in (("inside",inside),("outside_LGI",lambda r:r["scaffold"]=="LGI" and not inside(r))):
        bybin={}
        for r in windows:
            if predicate(r) and validq(r): bybin.setdefault(r["mid"]//100_000,[]).append(r)
        thinned=[v[int(rng.integers(len(v)))] for v in bybin.values()]
        q=np.array([[r["q1"],r["q2"],r["q3"]] for r in thinned])
        thin_rows.append({"subset":name,"n_100kb_blocks":len(thinned),"q1":q[:,0].mean(),"q2":q[:,1].mean(),"q3":q[:,2].mean(),"D":(q[:,2]-q[:,0]).mean()})
    write_tsv(OUT/"linked_window_thinning.tsv",thin_rows)

    # Dominant topology transition count and longest q3-majority run on LGI.
    informative=sorted([r for r in lgi if validq(r)],key=lambda r:r["mid"])
    labels=[int(np.argmax([r["q1"],r["q2"],r["q3"]]))+1 for r in informative]
    transition=sum(int(np.argmax([a["q1"],a["q2"],a["q3"]]))!=int(np.argmax([b["q1"],b["q2"],b["q3"]]))
                   for a,b in zip(lgi[:-1],lgi[1:]) if b["window"]==a["window"]+1 and validq(a) and validq(b))
    q3blocks=[b for b in block if math.isfinite(b["q3"]) and b["q3"]>b["q1"] and b["q3"]>b["q2"]]
    btrans=sum(int(np.argmax([a["q1"],a["q2"],a["q3"]]))!=int(np.argmax([b["q1"],b["q2"],b["q3"]]))
               for a,b in zip(block[:-1],block[1:]) if b["start"]-a["start"]==100_000 and math.isfinite(a["q1"]) and math.isfinite(b["q1"]))
    bruns=[]; bcurrent=[]
    for b in block:
        if b in q3blocks and (not bcurrent or b["start"]==bcurrent[-1]["end"]+1): bcurrent.append(b)
        else:
            if bcurrent: bruns.append(bcurrent)
            bcurrent=[b] if b in q3blocks else []
    if bcurrent: bruns.append(bcurrent)
    blongest=max(bruns,key=lambda x:x[-1]["end"]-x[0]["start"]) if bruns else []
    wruns=[]; wcurrent=[]; previous_window=None
    for r in lgi:
        label=int(np.argmax([r["q1"],r["q2"],r["q3"]]))+1 if validq(r) else 0
        if label==3 and (not wcurrent or r["window"]==previous_window+1): wcurrent.append(r)
        else:
            if wcurrent: wruns.append(wcurrent)
            wcurrent=[r] if label==3 else []
        previous_window=r["window"]
    if wcurrent: wruns.append(wcurrent)
    wlongest=max(wruns,key=lambda x:x[-1]["end"]-x[0]["start"]) if wruns else []
    with open(OUT/"spatial_signal.tsv","w") as fh:
        fh.write("metric\tvalue\n")
        fh.write(f"LGI_informative_windows\t{len(informative)}\nLGI_dominant_topology_transitions\t{transition}\n")
        fh.write(f"LGI_transitions_per_Mb\t{transition/((max(r['end'] for r in lgi)-min(r['start'] for r in lgi))/1e6):.8g}\n")
        fh.write(f"q3_dominant_runs_in_original_windows\t{len(wruns)}\n")
        fh.write(f"longest_q3_run_original_windows_bp\t{wlongest[-1]['end']-wlongest[0]['start']+1 if wlongest else 0}\n")
        fh.write(f"100kb_blocks_q3_dominant\t{len(q3blocks)}\n100kb_block_transitions\t{btrans}\n")
        fh.write(f"q3_dominant_100kb_runs\t{len(bruns)}\n")
        fh.write(f"longest_q3_run_start\t{blongest[0]['start'] if blongest else 'NA'}\nlongest_q3_run_end\t{blongest[-1]['end'] if blongest else 'NA'}\n")
        fh.write(f"longest_q3_run_span_bp\t{blongest[-1]['end']-blongest[0]['start']+1 if blongest else 0}\n")

    # Figures: topology frequencies, branch metrics, region distributions, and quartet inference sensitivity.
    fig, ax = plt.subplots(figsize=(12, 4.5))
    x = np.array([r["mid"] / 1e6 for r in lgi]); order = np.argsort(x)
    for key, label, color in (("q1","q1 species topology","#3366aa"),("q2","q2 other resolution","#e69f00"),("q3","q3 N. hispaniola + mat A","#c43c39")):
        y = np.array([r[key] for r in lgi])[order]
        ax.scatter(x[order], y, s=2, alpha=.23, color=color)
        if len(y) > 500:
            valid = np.isfinite(y)
            smooth = np.convolve(np.nan_to_num(y), np.ones(101), mode="same") / np.maximum(1, np.convolve(valid.astype(float), np.ones(101), mode="same"))
            ax.plot(x[order], smooth, color=color, lw=1.4, label=label)
    ax.axvspan(BOUNDARY[0]/1e6, BOUNDARY[1]/1e6, color="grey", alpha=.18, label="published L10 SR interval")
    ax.set(xlabel="LGI position (Mb)", ylabel="TWISST topology weight", ylim=(-.03,1.03)); ax.legend(ncol=2, frameon=False); fig.tight_layout()
    for ext in ("png","pdf","svg"): fig.savefig(FIG/f"figure1_spatial_topology.{ext}", dpi=300)
    plt.close(fig)

    fig, axs = plt.subplots(3,1,figsize=(12,8),sharex=True)
    for key,label in (("pair1_A_hisp_minus_a_hisp","A: d(CJ03,hisp) − d(CJ04,hisp)"),("pair1_B_mean_six_pairwise","B: mean six quartet distances"),("pair1_internal_branch","focal quartet internal edge")):
        valid_rows=[r for r in branch_rows if r["linkage_group"]=="LGI" and r.get(key) not in (None,"NA")]
        vals=np.array([float(r[key]) for r in valid_rows]); xx=np.array([r["mid"]/1e6 for r in valid_rows])
        axs[0 if key.startswith("pair1_A") else 1 if key.startswith("pair1_B") else 2].scatter(xx, vals, s=2, alpha=.3)
    for ax in axs: ax.axvspan(1.5,8.0,color="grey",alpha=.18); ax.grid(alpha=.2)
    axs[0].set_ylabel("Patristic contrast\n(substitutions/site)"); axs[1].set_ylabel("Mean patristic distance"); axs[2].set_ylabel("Internal edge"); axs[2].set_xlabel("LGI position (Mb)")
    fig.tight_layout()
    for ext in ("png","pdf","svg"): fig.savefig(FIG/f"figure2_spatial_branch_lengths.{ext}",dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9,5)); names=["outside_SR","LGI_flank_250kb","introgressed_SR_LGI_1.5_8.0Mb","other_linkage_groups"]
    colors=["#3366aa","#e69f00","#c43c39"]; ticklabels=[]; xx=np.arange(len(names))
    for j,key in enumerate(("q1","q2","q3")):
        means=[]
        for name in names:
            rows=[r for r in windows if masks[name](r) and math.isfinite(r[key])]
            means.append(float(np.mean([r[key] for r in rows])))
        ax.bar(xx+(j-1)*.24,means,width=.22,color=colors[j],label=("q1 species","q2 other","q3 hisp + mat A")[j])
    for name in names:
        n=sum(1 for r in windows if masks[name](r) and all(math.isfinite(r[k]) for k in ("q1","q2","q3")))
        ticklabels.append(name.replace("_","\n")+f"\n(n={n:,})")
    ax.set_xticks(xx,ticklabels); ax.set_ylabel("Mean TWISST topology weight"); ax.set_ylim(0,1)
    ax.legend(frameon=False,ncol=3,loc="upper center")
    fig.tight_layout()
    for ext in ("png","pdf","svg"): fig.savefig(FIG/f"figure3_region_comparison.{ext}",dpi=300)
    plt.close(fig)

    fig, axs = plt.subplots(1,2,figsize=(10,4.5)); inf_names=("all","introgressed_SR_LGI_1.5_8.0Mb","outside_SR")
    ss=[next(x for x in summaries if x and x["subset"]==name) for name in inf_names]
    xx=np.arange(3); colors=["#3366aa","#e69f00","#c43c39"]
    for j,key in enumerate(("q1_species","q2_other","q3_hisp_A")):
        axs[0].bar(xx+(j-1)*.24,[s[key] for s in ss],width=.22,color=colors[j],label=("q1","q2","q3")[j])
    axs[0].set_xticks(xx,["all","inside","outside"]); axs[0].set_ylabel("Quartet concordance"); axs[0].legend(frameon=False,title="Resolution")
    for i,s in enumerate(ss):
        cu=s["MSC_CU_if_symmetric"]
        if math.isfinite(cu): axs[1].bar(i,cu,color="#555555",width=.55)
        else: axs[1].text(i,.5,"not\nidentifiable",ha="center",va="center",fontsize=9)
    axs[1].set_xticks(range(3),["all","inside","outside"]); axs[1].set_ylabel("MSC-equivalent internal CU\n(only if discordants symmetric)"); axs[1].set_ylim(0,1.05)
    fig.tight_layout()
    for ext in ("png","pdf","svg"): fig.savefig(FIG/f"figure4_quartet_inference_sensitivity.{ext}",dpi=300)
    plt.close(fig)

    # Matplotlib's SVG writer leaves insignificant trailing spaces on path data;
    # strip them so generated vector figures pass repository whitespace checks.
    for svg in FIG.glob("*.svg"):
        svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")

    # Small file inventory for checksum/provenance report.
    with open(OUT / "source_manifest.tsv", "w") as fh:
        fh.write("relative_path\tbytes\tmd5\trole\tsource_url\tretrieval_date\n")
        import hashlib
        for p,role in ((ROOT/"data/neurospora/raw/neurospora.tar.gz","Zenodo archive"),(DATA,"window coordinates"),(TREES,"local Newick trees"),(WEIGHTS,"published TWISST weights")):
            h=hashlib.md5(p.read_bytes()).hexdigest()
            fh.write(f"{p.relative_to(ROOT)}\t{p.stat().st_size}\t{h}\t{role}\thttps://zenodo.org/records/4946041\t2026-11-15\n")
    valid_tree_count = sum(bool(r["tip_count"]) for r in tree_quality)
    print(f"window_records={len(windows)} deposited_Newick_trees={valid_tree_count} tip_set_variants={len(observed_tipsets)} LGI_records={sum(r['scaffold']=='LGI' for r in windows)} inside_valid_weights={len(inrows)} outside_LGI_valid_weights={len(outrows)} autocorrelation100kb={ac:.4f}")
    for s in summaries:
        if s: print(s)


if __name__ == "__main__":
    main()
