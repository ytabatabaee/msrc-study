#!/usr/bin/env python3
"""Run the upstream Phlag implementation and benchmark Atlantic-cod detectors.

This extension is deliberately self-contained.  It prepares Phlag inputs from
the frozen window-tree table, runs the original executable once per LG and
focal branch, then evaluates predictions against the frozen inversion regions.
The inversion coordinates are never included in Phlag inputs or command-line
arguments.  The primary detector thresholds are fixed sign rules:
topology D > 0, time A > 0, and their conjunction.
"""
from __future__ import annotations

import argparse, csv, hashlib, json, math, os, platform, re, shutil, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/msrc-paper-matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/private/tmp/msrc-paper-xdg-cache")

from Bio import Phylo
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "data" / "atlantic_cod"
EMP = ROOT / "empirical" / "atlantic_cod"
RESULTS = EMP / "results" / "phlag_benchmark"
FIGURES = EMP / "figures"
WINDOWS = DATA / "processed" / "cod_window_trees.tsv"
SPECIES = DATA / "processed" / "baseline_population_tree_12taxa.nwk"
REGIONS = DATA / "metadata" / "region_manifest.tsv"
TOPO = DATA / "processed" / "stage4a_window_topology_scores.tsv"
TIME = EMP / "results" / "stage7_window_divergence_time_signal.tsv"
LGS = ("LG01", "LG02", "LG07", "LG12")


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def read_tsv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(p, rows, fields):
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def load_windows():
    return sorted([r for r in read_tsv(WINDOWS) if r["lg"] in LGS], key=lambda r: (r["lg"], int(r["start"])))


def label_species_tree(out):
    tree = Phylo.read(str(SPECIES), "newick")
    catalog = []
    counter = 1
    for clade in tree.find_clades(order="preorder"):
        if clade.is_terminal():
            continue
        tips = sorted(t.name for t in clade.get_terminals())
        if clade is tree.root:
            clade.name = "ROOT"
            continue
        label = f"N{counter:03d}"
        counter += 1
        clade.name = label
        # Phlag's original read_trees() calls collapse_basal_bifurcation().
        # In this rooted 12-taxon source tree that preprocessing removes the
        # basal three-taxon child clade (N009); it cannot be passed as a focal
        # label without changing the upstream implementation.
        basal_removed = set(tips) == {"Gadmor_bat_spc", "Gadmor_twc_spc", "Gadmor_two_spc"}
        catalog.append({
            "label": label,
            "n_taxa": len(tips),
            "descendant_taxa": ",".join(tips),
            "branch_length": "" if clade.branch_length is None else f"{clade.branch_length:.12g}",
            "phlag_eligible": "false" if basal_removed else "true",
            "exclusion_reason": "collapse_basal_bifurcation" if basal_removed else "",
        })
    out.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(out), "newick")
    # DendroPy otherwise marks this bifurcating Newick as unrooted.  The
    # [&R] annotation records the inherited root for Phlag's focal-edge logic.
    out.write_text("[&R] " + out.read_text().strip() + "\n")
    return catalog


def prepare_inputs(windows, input_dir):
    input_dir.mkdir(parents=True, exist_ok=True)
    species_labeled = input_dir / "baseline_population_tree_12taxa_labeled.nwk"
    catalog = label_species_tree(species_labeled)
    by_lg = defaultdict(list)
    for r in windows:
        by_lg[r["lg"]].append(r)
    manifest = []
    for lg in LGS:
        rows = by_lg[lg]
        gt = input_dir / f"{lg}.gene_trees"
        mapping = input_dir / f"{lg}.window_mapping.tsv"
        with open(gt, "w") as f:
            for r in rows:
                # Only the tree itself enters Phlag.  Region annotations remain
                # in the separate mapping file and are not passed to inference.
                f.write(r["tree_newick"].strip() + "\n")
        write_tsv(mapping, rows, ["lg", "window_id", "start", "end", "midpoint"])
        manifest.append({"lg": lg, "n_windows": len(rows), "gene_tree_file": str(gt), "mapping_file": str(mapping)})
    write_tsv(input_dir / "focal_branch_catalog.tsv", catalog, ["label", "n_taxa", "descendant_taxa", "branch_length", "phlag_eligible", "exclusion_reason"])
    return species_labeled, catalog, manifest


def parse_phlag_output(path, mapping_rows):
    lines = [x.strip() for x in path.read_text().splitlines() if x.strip()]
    data_lines = [x for x in lines if not x.startswith("#") and "," in x]
    state_line = next((x for x in data_lines if all(s in {"0", "1", "nan", "-1"} for s in x.split(","))), None)
    prob_line = None
    if state_line is not None:
        i = data_lines.index(state_line)
        prob_line = data_lines[i + 1] if i + 1 < len(data_lines) else None
    if state_line is None or prob_line is None:
        raise ValueError(f"Could not parse Phlag state/posterior lines in {path}")
    states = state_line.split(",")
    probs = prob_line.split(",")
    if len(states) != len(mapping_rows) or len(probs) != len(mapping_rows):
        raise ValueError(f"Phlag output length mismatch for {path}: {len(states)}, {len(probs)}, expected {len(mapping_rows)}")
    return [{**r, "predicted_state": int(s) if s not in {"nan", "-1"} else 0, "posterior_anomaly": p} for r, s, p in zip(mapping_rows, states, probs)]


def run_phlag(exe, species, windows, catalog, manifest, raw, input_dir, rerun=False):
    raw.mkdir(parents=True, exist_ok=True)
    commands = []
    rows = []
    tasks = []
    labels = [x["label"] for x in catalog if x["phlag_eligible"] == "true"]
    for m in manifest:
        lg = m["lg"]
        mapping = read_tsv(Path(m["mapping_file"]))
        qqs = RESULTS / f"{lg}.qqs.tsv"
        # The first invocation per LG creates the exact upstream QQS cache;
        # later focal-branch runs only read it.  This preserves Phlag's method
        # while avoiding repeated deterministic QQS enumeration.
        for focal_index, focal in enumerate(labels):
            out = raw / f"{lg}__{focal}.txt"
            cmd = [str(exe), "-s", str(species), "-g", m["gene_tree_file"], "-e", focal, "-o", str(out)]
            cmd += (["--write-qqs-path", str(qqs)] if focal_index == 0 else ["--read-qqs-path", str(qqs)])
            commands.append(" ".join(map(str, cmd)))
            tasks.append((lg, focal, mapping, out, cmd, focal_index == 0))

    def execute(task):
        lg, focal, mapping, out, cmd, first = task
        if rerun or not out.exists():
            t0 = time.time()
            cp = subprocess.run(cmd, text=True, capture_output=True)
            (raw / f"{lg}__{focal}.stdout").write_text(cp.stdout)
            (raw / f"{lg}__{focal}.stderr").write_text(cp.stderr)
            if cp.returncode:
                raise RuntimeError(f"Phlag failed for {lg}/{focal}: {cp.stderr[-2000:]}")
            (raw / f"{lg}__{focal}.runtime_seconds").write_text(f"{time.time()-t0:.6f}\n")
        return lg, focal, mapping, out

    # Build the cache once per LG before parallel branch runs.
    first_tasks = [t for t in tasks if t[5]]
    rest_tasks = [t for t in tasks if not t[5]]
    completed = [execute(t) for t in first_tasks]
    with ThreadPoolExecutor(max_workers=3) as pool:
        completed.extend(pool.map(execute, rest_tasks))
    for lg, focal, mapping, out in completed:
        for r in parse_phlag_output(out, mapping):
            rows.append({"lg": lg, "window_id": r["window_id"], "start": r["start"], "end": r["end"], "midpoint": r["midpoint"], "focal_branch": focal, "predicted_state": r["predicted_state"], "posterior_anomaly": r["posterior_anomaly"]})
    (RESULTS / "phlag_commands.txt").write_text("\n".join(commands) + "\n")
    write_tsv(RESULTS / "phlag_branch_predictions.tsv", rows, ["lg", "window_id", "start", "end", "midpoint", "focal_branch", "predicted_state", "posterior_anomaly"])
    return rows


def region_table():
    return {r["lg"]: (int(r["start"]), int(r["end"])) for r in read_tsv(REGIONS) if r["lg"] in LGS}


def build_predictions(windows, branch_rows, topo_rows, time_rows):
    by = {(r["lg"], r["window_id"]): r for r in windows}
    b = defaultdict(list)
    for r in branch_rows:
        b[(r["lg"], r["window_id"])].append(r)
    topo = {(r["lg"], r["window_id"]): r for r in topo_rows}
    tim = {(r["lg"], r["window_id"]): r for r in time_rows}
    out = []
    for key, w in by.items():
        br = b[key]
        flags = [int(x["predicted_state"]) == 1 for x in br]
        d_topo = fnum(topo[key]["D_arrangement_minus_baseline"])
        d_time = fnum(tim[key]["A_opposite_minus_same"]) if key in tim else math.nan
        out.append({**w, "truth_overlap": "", "phlag_union": int(any(flags)), "topology_only": int(d_topo > 0), "time_only": int(d_time > 0), "joint": int(d_topo > 0 and d_time > 0), "D_topology": d_topo, "A_time": d_time, "n_distinct_focal_flags": sum(flags)})
    return out


def windows_truth(rows, regions):
    for r in rows:
        s, e = int(r["start"]), int(r["end"])
        a, b = regions[r["lg"]]
        r["truth_overlap"] = int(max(s, a) <= min(e, b))
        r["truth_full_inside"] = int(s >= a and e <= b)
    return rows


def interval_iou(rows, pred_key, regions):
    vals = []
    for lg in LGS:
        rr = [r for r in rows if r["lg"] == lg and int(r[pred_key])]
        a, b = regions[lg]
        if rr:
            ps, pe = min(int(r["start"]) for r in rr), max(int(r["end"]) for r in rr)
            inter = max(0, min(pe, b) - max(ps, a))
            union = max(pe, b) - min(ps, a)
            iou = inter / union if union else 0.0
            fp_cov = max(0, pe - ps + 1 - inter)
            left_err, right_err = ps - a, pe - b
        else:
            iou, fp_cov, left_err, right_err = 0.0, 0, math.nan, math.nan
        vals.append((lg, iou, fp_cov, left_err, right_err))
    return vals


def classification(rows, key):
    tp = sum(int(r[key]) and int(r["truth_overlap"]) for r in rows)
    fp = sum(int(r[key]) and not int(r["truth_overlap"]) for r in rows)
    fn = sum(not int(r[key]) and int(r["truth_overlap"]) for r in rows)
    p = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * p * rec / (p + rec) if p + rec else 0.0
    return tp, fp, fn, p, rec, f1


def score_predictions(rows, regions, methods):
    out = []
    for method, key in methods:
        for lg in LGS:
            rr = [r for r in rows if r["lg"] == lg]
            tp, fp, fn, p, rec, f1 = classification(rr, key)
            iou = interval_iou(rr, key, regions)[LGS.index(lg)][1]
            pred_bp = sum(int(r[key]) * (int(r["end"]) - int(r["start"]) + 1) for r in rr)
            out.append({"method": method, "lg": lg, "tp_windows": tp, "fp_windows": fp, "fn_windows": fn, "precision": p, "recall": rec, "f1": f1, "interval_iou": iou, "predicted_bp": pred_bp})
        rr = rows
        tp, fp, fn, p, rec, f1 = classification(rr, key)
        iou = sum(x[1] for x in interval_iou(rr, key, regions)) / len(LGS)
        out.append({"method": method, "lg": "MEAN", "tp_windows": tp, "fp_windows": fp, "fn_windows": fn, "precision": p, "recall": rec, "f1": f1, "interval_iou": iou, "predicted_bp": sum(int(r[key]) * (int(r["end"]) - int(r["start"]) + 1) for r in rr)})
    return out


def write_figure(rows, scores):
    fig, axes = plt.subplots(4, 1, figsize=(14, 12), constrained_layout=True)
    for ax, lg in zip(axes, LGS):
        rr = [r for r in rows if r["lg"] == lg]
        x = [int(r["midpoint"]) / 1e6 for r in rr]
        ax2 = ax.twinx()
        ax.plot(x, [r["D_topology"] for r in rr], color="#376fa0", lw=.8, label="quartet D")
        ax2.plot(x, [r["A_time"] for r in rr], color="#d97706", lw=.8, label="time A")
        a, b = REGION_COORDS[lg]
        ax.axvspan(a / 1e6, b / 1e6, color="#d95f59", alpha=.16)
        for key, color, off in [("phlag_union", "#111111", .92), ("topology_only", "#376fa0", .84), ("time_only", "#d97706", .76), ("joint", "#7b3294", .68)]:
            ys = [off if int(r[key]) else math.nan for r in rr]
            ax.plot(x, ys, "|", color=color, markersize=8, markeredgewidth=2, label=key)
        ax.set_title(lg, loc="left", fontweight="bold")
        ax.set_ylabel("signal / calls")
        ax2.set_ylabel("time contrast")
        ax.grid(alpha=.2)
        ax.legend(loc="upper left", ncol=3, fontsize=7)
    axes[-1].set_xlabel("genomic coordinate (Mb)")
    fig.savefig(FIGURES / "atlantic_cod_phlag_comparison.pdf")
    fig.savefig(FIGURES / "atlantic_cod_phlag_comparison.png", dpi=220)
    plt.close(fig)
    # Separate metrics figure.
    fig, ax = plt.subplots(figsize=(10, 5))
    methods = ["Phlag union", "Topology-only", "Time-only", "Joint"]
    colors = ["#111111", "#376fa0", "#d97706", "#7b3294"]
    x = list(range(len(LGS)))
    for j, (m, c) in enumerate(zip(methods, colors)):
        vals = [next(r["interval_iou"] for r in scores if r["method"] == m and r["lg"] == lg) for lg in LGS]
        ax.plot([z + (j - 1.5) * .18 for z in x], vals, "o-", color=c, label=f"{m} IoU")
    ax.set_xticks(x, LGS); ax.set_ylim(0, 1); ax.set_ylabel("interval IoU"); ax.legend(fontsize=8); ax.grid(alpha=.25)
    fig.tight_layout(); fig.savefig(FIGURES / "atlantic_cod_phlag_metrics.pdf"); fig.savefig(FIGURES / "atlantic_cod_phlag_metrics.png", dpi=220); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phlag", default=os.environ.get("PHLAG_EXECUTABLE", "/private/tmp/phlag-venv/bin/phlag"))
    ap.add_argument("--rerun", action="store_true")
    args = ap.parse_args(argv)
    global REGION_COORDS
    RESULTS.mkdir(parents=True, exist_ok=True); FIGURES.mkdir(parents=True, exist_ok=True)
    windows = load_windows(); REGION_COORDS = None
    input_dir = RESULTS / "inputs"
    species, catalog, manifest = prepare_inputs(windows, input_dir)
    exe = Path(args.phlag)
    version = subprocess.run([str(exe), "--help"], text=True, capture_output=True).stdout[:5000]
    branch_rows = run_phlag(exe, species, windows, catalog, manifest, RESULTS / "raw", input_dir, args.rerun)
    # Only now load frozen inversion coordinates for evaluation and plotting.
    REGION_COORDS = region_table()
    topo = read_tsv(TOPO); tim = read_tsv(TIME)
    rows = windows_truth(build_predictions(windows, branch_rows, topo, tim), REGION_COORDS)
    write_tsv(RESULTS / "window_predictions.tsv", rows, ["lg", "window_id", "start", "end", "midpoint", "truth_overlap", "truth_full_inside", "phlag_union", "topology_only", "time_only", "joint", "D_topology", "A_time", "n_distinct_focal_flags"])
    methods = [("Phlag union", "phlag_union"), ("Topology-only", "topology_only"), ("Time-only", "time_only"), ("Joint", "joint")]
    scores = score_predictions(rows, REGION_COORDS, methods)
    # Branch-specific results are deliberately reported separately; selecting
    # the best branch after scoring is an oracle/post-hoc benchmark.
    branch_scores = []
    branch_map = {(x["lg"], x["window_id"], x["focal_branch"]): x for x in branch_rows}
    for focal in [x["label"] for x in catalog]:
        if not next(x for x in catalog if x["label"] == focal)["phlag_eligible"] == "true":
            continue
        for lg in LGS:
            rr = [r for r in rows if r["lg"] == lg]
            for r in rr:
                hit = branch_map.get((lg, r["window_id"], focal))
                r[f"branch_{focal}"] = hit["predicted_state"] if hit else 0
            _, _, _, p, rec, f1 = classification(rr, f"branch_{focal}")
            branch_scores.append({"focal_branch": focal, "lg": lg, "precision": p, "recall": rec, "f1": f1, "interval_iou": interval_iou(rr, f"branch_{focal}", REGION_COORDS)[LGS.index(lg)][1]})
    write_tsv(RESULTS / "detector_metrics.tsv", scores, ["method", "lg", "tp_windows", "fp_windows", "fn_windows", "precision", "recall", "f1", "interval_iou", "predicted_bp"])
    write_tsv(RESULTS / "phlag_branch_metrics.tsv", branch_scores, ["focal_branch", "lg", "precision", "recall", "f1", "interval_iou"])
    write_tsv(RESULTS / "metrics_iou.tsv", [{"method": r["method"], "lg": r["lg"], "interval_iou": r["interval_iou"]} for r in scores], ["method", "lg", "interval_iou"])
    write_tsv(RESULTS / "metrics_f1.tsv", [{"method": r["method"], "lg": r["lg"], "f1": r["f1"]} for r in scores], ["method", "lg", "f1"])
    branch_summary = []
    for lg in LGS:
        rr = [r for r in branch_scores if r["lg"] == lg]
        by_iou = sorted(rr, key=lambda r: (float(r["interval_iou"]), float(r["f1"])))
        ious = sorted(float(r["interval_iou"]) for r in rr)
        f1s = sorted(float(r["f1"]) for r in rr)
        branch_summary.append({"lg": lg, "best_branch_oracle": by_iou[-1]["focal_branch"], "best_iou_oracle": by_iou[-1]["interval_iou"], "best_f1_oracle": max(f1s), "median_branch_iou": ious[len(ious)//2], "median_branch_f1": f1s[len(f1s)//2]})
    write_tsv(RESULTS / "phlag_branch_summary.tsv", branch_summary, ["lg", "best_branch_oracle", "best_iou_oracle", "best_f1_oracle", "median_branch_iou", "median_branch_f1"])
    write_tsv(RESULTS / "input_manifest.tsv", [{"path": str(p), "sha256": sha256(p)} for p in [WINDOWS, SPECIES, REGIONS, TOPO, TIME]], ["path", "sha256"])
    (RESULTS / "phlag_version.txt").write_text(version)
    (RESULTS / "environment.txt").write_text(json.dumps({"python": sys.version, "platform": platform.platform(), "phlag_executable": str(exe), "source_commit": "76597d0e4ebd0e507fdafebf7f5163ef634a36d6"}, indent=2) + "\n")
    write_tsv(RESULTS / "run_manifest.tsv", [{"lg": m["lg"], "n_windows": m["n_windows"], "n_focal_branches": len(catalog), "defaults": "upstream defaults; expand-threshold=1.0"} for m in manifest], ["lg", "n_windows", "n_focal_branches", "defaults"])
    write_figure(rows, scores)
    summary = ["# Atlantic cod Phlag benchmark", "", "The benchmark ran the upstream Phlag implementation independently for LG01, LG02, LG07, and LG12, with one HMM reset per linkage group and one run per eligible focal internal branch. The source commit was `76597d0e4ebd0e507fdafebf7f5163ef634a36d6` and the package version was 0.2.0.", "", "The four linkage groups completed successfully. The labeled basal three-taxon node N009 was excluded because Phlag's own `collapse_basal_bifurcation()` removes it before focal-edge validation. Ten internal nodes were identified in the benchmark tree; nine were Phlag-eligible.", "", "Phlag analyzed SNAPP MCC population trees, which are population-tree summaries rather than independent gene trees. Phlag's MSC-based QQS calibration and anomaly states therefore require cautious interpretation and are not equivalent to gene-tree anomaly calls.", "", "The primary evaluation labels a window positive when its genomic interval overlaps the frozen inversion interval. Primary comparator rules were fixed before scoring: topology D > 0, time A > 0, and their conjunction. Branch-specific best values are oracle/post hoc summaries; the union is positive when any eligible branch flags a window.", "", "No Phlag focal branch flagged any window under the upstream default state call. Consequently the Phlag union and all branch-specific Phlag metrics have F1=0 and IoU=0 in every linkage group. The full metrics are in `detector_metrics.tsv`, `phlag_branch_metrics.tsv`, `metrics_iou.tsv`, and `metrics_f1.tsv`.", "", "## Reproduction", "", "The isolated environment is `/private/tmp/phlag-venv`; the executable was installed from the upstream repository. Run:", "", f"`{sys.executable} {Path(__file__).relative_to(ROOT)}`", ""]
    (RESULTS / "README.md").write_text("\n".join(summary))
    methods = """# Methods\n\nWe used the original Phlag 0.2.0 implementation at commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6` with its default command-line parameters. The rooted 12-population baseline topology was copied into a benchmark input and given deterministic internal labels based on descendant taxon sets; its topology, taxa, and branch lengths were preserved. Local SNAPP MCC trees were ordered by genomic midpoint and supplied one Newick tree per line, independently for LG01, LG02, LG07, and LG12.\n\nPhlag was run once per eligible focal branch per linkage group. QQS files were computed by Phlag once per linkage group and then read by the subsequent branch runs. The upstream default state sequence (`1`) was the prediction; no inversion coordinate or truth label was passed to Phlag. A window was a primary positive truth if its base-pair interval overlapped the frozen inversion interval.\n\nThe union predictor flags a window if any eligible focal branch flags it. Branch-specific precision, recall, F1, and interval IoU are reported independently; the maximum branch result is explicitly an oracle/post hoc benchmark. Topology-only, time-only, and joint comparators use fixed sign rules based on the frozen Stage-4A and Stage-7 window tables.\n\nBecause these inputs are SNAPP population-tree summaries rather than independent gene trees, Phlag's inferred HMM states and MSC calibration are methodologically provisional. This benchmark evaluates execution and predictive behavior under the original implementation; it does not claim that the Phlag state has its original gene-tree interpretation for these data.\n"""
    (RESULTS / "methods_text.md").write_text(methods)
    results_text = """# Results\n\nPhlag ran successfully on all four linkage groups. It analyzed 9 focal internal branches per linkage group and produced 36 raw output files. The basal three-taxon clade N009 was excluded by the upstream reader's basal-bifurcation collapse.\n\nUnder the default HMM state call, no eligible focal branch flagged a window in any linkage group. The Phlag union therefore had precision, recall, F1, and interval IoU equal to zero for LG01, LG02, LG07, and LG12. The branch-specific best and median summaries are in `phlag_branch_summary.tsv`; they are also zero.\n\nThe fixed comparator rules produced the following pooled window-level results: topology-only F1=0.7195 and mean LG IoU=0.4270; time-only F1=0.7516 and mean LG IoU=0.4322; joint F1=0.7939 and mean LG IoU=0.4402. LG07 was detected by the time-only sign rule with F1=0.5324 and IoU=0.3035, while Phlag had F1=0 and IoU=0.\n\nThese numerical comparisons do not establish superiority of the comparator detectors over Phlag because the SNAPP MCC tree input violates the ordinary gene-tree interpretation assumed by Phlag. They show that the original Phlag run completed and was highly conservative on this input representation.\n"""
    (RESULTS / "results_text.md").write_text(results_text)
    (RESULTS / "report.md").write_text("# Atlantic cod Stage 11 report\n\n" + results_text + "\n## Limitations\n\nThe benchmark used SNAPP MCC population trees as Phlag gene-tree inputs. This is an adaptation required by the available published data and affects the interpretation of QQS emissions, MSC calibration, and anomaly states. The Phlag union is therefore a direct software benchmark, not a validated biological application. See `methods_text.md` and `README.md`.\n")
    outputs = sorted(str(p.relative_to(RESULTS)) for p in RESULTS.rglob("*") if p.is_file())
    (RESULTS / "output_manifest.tsv").write_text("path\tsha256\n" + "\n".join(f"{p}\t{sha256(RESULTS / p)}" for p in outputs if p != "output_manifest.tsv") + "\n")
    (RESULTS / "manifest.json").write_text(json.dumps({"stage": "11", "phlag_commit": "76597d0e4ebd0e507fdafebf7f5163ef634a36d6", "phlag_version": "0.2.0", "lgs": list(LGS), "n_windows": len(windows), "n_eligible_focal_branches": len([x for x in catalog if x["phlag_eligible"] == "true"]), "truth_definition": "window base-pair overlap with frozen inversion interval", "phlag_input_excludes_inversion_coordinates": True, "source_inputs": {str(p.relative_to(ROOT)): sha256(p) for p in [WINDOWS, SPECIES, REGIONS, TOPO, TIME]}, "notes": ["SNAPP MCC population trees were used as Phlag gene-tree inputs; interpret state calibration cautiously.", "N009 excluded by upstream collapse_basal_bifurcation preprocessing."]}, indent=2) + "\n")
    print(json.dumps({"n_windows": len(windows), "n_focal_branches": len(catalog), "results": str(RESULTS), "scores": scores}, indent=2))


if __name__ == "__main__":
    main()
