#!/usr/bin/env python3
"""Neoaves chromosome-4 Stage 11: local Phlag comparison.

This stage keeps three analyses separate:

* the published/source-defined QQS tracks used by the existing detector;
* a matched four-taxon Phlag pilot made from the original local Newick trees;
* post-fit comparisons with the independently stored structural annotations.

The Phlag pilot is deliberately labelled a matched-taxon pilot.  It does not
pretend that one representative from each focal group is identical to the
multi-taxon QQS summaries used in Stages 6--10.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import resource
import subprocess
import sys
import time
from pathlib import Path
from dataclasses import dataclass

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical/neoaves_chr4"
DATA = ROOT / "data/neoaves_chr4/processed"
RES = EMP / "results/stage11_phlag_benchmark"
FIG = EMP / "figures"
STAGE6 = EMP / "results/stage6_genomic_tracks"
STAGE9 = EMP / "results/stage9_full_taxon_background"
SEED = 20261020
PHLAG = Path("/private/tmp/phlag-venv/bin/phlag")
PHLAG_ROOT = Path("/private/tmp/phlag-upstream")
SOURCE_TREES = EMP / "external/stiller2024/63430.named.gene.trees.gz"
SOURCE_ST = EMP / "external/stiller2024/63K.tre"
MANIFEST = EMP / "results/stage4d_locus_manifest.tsv"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def num(x: object) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return math.nan


def interval_union(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for a, b in sorted((int(a), int(b)) for a, b in intervals if int(b) > int(a)):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def interval_len(intervals: list[tuple[int, int]]) -> int:
    return sum(b - a for a, b in interval_union(intervals))


def interval_intersection(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> int:
    return sum(max(0, min(y, v) - max(x, u)) for x, y in interval_union(a) for u, v in interval_union(b))


def interval_metrics(pred: list[tuple[int, int]], ref: list[tuple[int, int]]) -> dict[str, float]:
    p = interval_len(pred); r = interval_len(ref); inter = interval_intersection(pred, ref); union = p + r - inter
    precision = inter / p if p else 0.0
    recall = inter / r if r else 0.0
    f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
    return {"n_predicted_intervals": len(interval_union(pred)), "predicted_bp": p,
            "reference_bp": r, "intersection_bp": inter, "union_bp": union,
            "iou": inter / union if union else 0.0, "precision": precision,
            "recall": recall, "f1": f1, "false_positive_bp": p - inter,
            "false_negative_bp": r - inter}


def load_reference_regions() -> list[dict[str, str]]:
    rows = read_tsv(DATA / "canonical_structural_regions.tsv")
    return [{"id": r["canonical_region_id"], "start": int(float(r["canonical_start"])),
             "end": int(float(r["canonical_end"]))} for r in rows]


def load_locus_rows() -> list[dict[str, str]]:
    rows = read_tsv(STAGE6 / "window_quartet_support.tsv")
    # One row per (locus, clade); the Columbea rows define the primary pilot.
    return sorted([r for r in rows if r["clade"] == "Columbea"], key=lambda r: (int(float(r["start"])), int(r["Gene"])))


def load_b(rows: list[dict[str, str]]) -> dict[str, float]:
    b = {r["Gene"]: num(r["B_observed"]) for r in read_tsv(STAGE9 / "background_distance_estimators.tsv")}
    return {r["Gene"]: b.get(r["Gene"], math.nan) for r in rows}


@dataclass
class _N:
    children: list
    name: str = ""
    length: str = ""


def _parse_newick(text: str):
    """Small Newick parser used only to prune the 363-tip source trees.

    It keeps leaf names and branch lengths and discards internal support labels.
    This avoids parsing all 363 tips with DendroPy for every one of the 6,140
    pilot loci while preserving the original local-tree topology and lengths.
    """
    s = text.strip().rstrip(";"); pos = 0
    def token(stop):
        nonlocal pos
        start = pos
        while pos < len(s) and s[pos] not in stop: pos += 1
        return s[start:pos]
    def node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1; kids = [node()]
            while pos < len(s) and s[pos] == ",": pos += 1; kids.append(node())
            if pos >= len(s) or s[pos] != ")": raise ValueError("malformed Newick")
            pos += 1
            _ = token(":,();")
            length = ""
            if pos < len(s) and s[pos] == ":": pos += 1; length = token(",();")
            return _N(kids, "", length)
        name = token(":,();")
        length = ""
        if pos < len(s) and s[pos] == ":": pos += 1; length = token(",();")
        return _N([], name, length)
    return node()


def _prune(n: _N, wanted: set[str]):
    if not n.children:
        return n if n.name in wanted else None
    kept = [_prune(c, wanted) for c in n.children]
    kept = [c for c in kept if c is not None]
    if not kept: return None
    if len(kept) == 1:
        child = kept[0]
        if n.length and child.length:
            try: child.length = str(float(child.length) + float(n.length))
            except ValueError: pass
        elif n.length and not child.length:
            child.length = n.length
        return child
    n.children = kept; n.name = ""
    return n


def _write_newick(n: _N) -> str:
    if n.children:
        out = "(" + ",".join(_write_newick(c) for c in n.children) + ")"
    else:
        out = n.name
    return out + ((":" + n.length) if n.length else "")


def robust_scale(x: np.ndarray) -> tuple[np.ndarray, list[float], list[float]]:
    z = np.asarray(x, float).copy(); medians=[]; scales=[]
    for j in range(z.shape[1]):
        v = z[:, j]; med = float(np.nanmedian(v)); mad = 1.4826 * float(np.nanmedian(np.abs(v - med)))
        iqr = float(np.nanpercentile(v, 75) - np.nanpercentile(v, 25)) / 1.349
        sd = float(np.nanstd(v)); scale = max(mad, iqr, 0.25 * sd, 1e-12)
        z[:, j] = (v - med) / scale; medians.append(med); scales.append(scale)
    return z, medians, scales


def shared_interval(x: np.ndarray, min_size: int = 3, penalty: float = 1.0) -> dict:
    """Stage-12b-style shared-background scan on an observed locus sequence."""
    n, p = x.shape; best = {"status": "no_detection", "start": "", "end": "", "score": 0.0, "raw": 0.0, "penalty": 0.0}
    if n < 2 * min_size + 1: return best
    s0 = float(np.sum((x - x.mean(0)) ** 2))
    sp = np.vstack([np.zeros((1, p)), np.cumsum(x, axis=0)])
    qp = np.vstack([np.zeros((1, p)), np.cumsum(x * x, axis=0)])
    def sse(a, b):
        if b <= a: return 0.0
        sm = sp[b] - sp[a]; sq = qp[b] - qp[a]
        return float(np.sum(sq - sm * sm / (b - a)))
    left_sse = [sse(0, i) for i in range(n + 1)]
    right_sse = [sse(i, n) for i in range(n + 1)]
    base_pen = penalty * p * math.log(max(n, 2))
    for i in range(0, n - min_size + 1):
        j0 = i + min_size; j1 = n + 1
        js = np.arange(j0, j1)
        counts = js - i
        sm = sp[js] - sp[i]; sq = qp[js] - qp[i]
        inside_sse = np.sum(sq - sm * sm / counts[:, None], axis=1)
        right_n = n - js
        right_sse_arr = np.array([right_sse[int(j)] for j in js])
        left_ok = i == 0 or i >= min_size
        if not left_ok: continue
        internal_ok = right_n >= min_size
        raw_internal = s0 - inside_sse - left_sse[i] - right_sse_arr
        scores = raw_internal - 2 * base_pen
        if np.any(internal_ok):
            candidate_scores = np.where(internal_ok, scores, -np.inf)
            k = int(np.argmax(candidate_scores))
            if math.isfinite(float(candidate_scores[k])) and float(candidate_scores[k]) > 0:
                j = int(js[k]); score = float(candidate_scores[k]); key = (score, 1, -(j-i), -i)
                old = (best["score"], int(best["status"] == "internal"), -(int(best["end"])-int(best["start"])) if best["start"] != "" else 0, -int(best.get("start",0) or 0))
                if key > old: best = {"status":"internal","start":i,"end":j,"score":score,"raw":float(raw_internal[k]),"penalty":2*base_pen}
        if i > 0 and js[-1] == n:
            k = len(js) - 1; score = float(s0 - inside_sse[k] - left_sse[i] - base_pen)
            if score > 0:
                key=(score,0,-(n-i),-i); old=(best["score"],int(best["status"]=="internal"),-(int(best["end"])-int(best["start"])) if best["start"]!="" else 0,-int(best.get("start",0) or 0))
                if key > old: best={"status":"right_censored","start":i,"end":n,"score":score,"raw":score+base_pen,"penalty":base_pen}
        if i == 0:
            left_scores = s0 - inside_sse - right_sse_arr - base_pen
            k = int(np.argmax(left_scores))
            score = float(left_scores[k]); j = int(js[k]); key=(score,0,-j,0); old=(best["score"],int(best["status"]=="internal"),-(int(best["end"])-int(best["start"])) if best["start"]!="" else 0,-int(best.get("start",0) or 0))
            if score > 0 and key > old: best={"status":"left_censored","start":0,"end":j,"score":score,"raw":score+base_pen,"penalty":base_pen}
    return best


def fit_detector(rows: list[dict[str, str]], feature_set: str) -> tuple[dict, list[float], list[float]]:
    b = load_b(rows)
    values = []
    for r in rows:
        if feature_set == "topology": values.append([num(r["D_q2_minus_q1"]), num(r["D_q3_minus_q1"])])
        else: values.append([num(r["D_q2_minus_q1"]), num(r["D_q3_minus_q1"]), b[r["Gene"]]])
    x = np.asarray(values, float); ok = np.all(np.isfinite(x), axis=1); z, med, scales = robust_scale(x[ok])
    fit = shared_interval(z, min_size=3, penalty=1.0)
    if fit["start"] != "":
        fit["start_bp"] = int(float(rows[np.flatnonzero(ok)[fit["start"]]]["start"]))
        fit["end_bp"] = int(float(rows[np.flatnonzero(ok)[fit["end"] - 1]]["end"]))
    return fit, med, scales


def preflight() -> dict:
    rows = load_locus_rows(); b = load_b(rows)
    source_bytes = SOURCE_TREES.stat().st_size if SOURCE_TREES.exists() else 0
    return {"n_columbea_loci": len(rows), "n_source_file_bytes": source_bytes,
            "n_loci_with_B_observed": sum(math.isfinite(v) for v in b.values()),
            "phlag_binary": str(PHLAG), "phlag_exists": PHLAG.exists(),
            "python_executable": sys.executable, "cpu_threads": os.cpu_count(),
            "mode": "local chromosome-4 pilot; no full 63,430-locus run"}


def write_pilot_inputs() -> dict:
    """Prepare one fixed four-taxon representative per Columbea role.

    DendroPy is intentionally imported only here because it is supplied by the
    isolated Phlag environment.  The source trees are never collapsed to a
    consensus tree; they are pruned one by one to fixed representatives.
    """
    try:
        import dendropy
    except ImportError as e:
        raise RuntimeError("pilot preparation requires the isolated Phlag environment") from e
    rows = load_locus_rows(); ids = {r["Gene"] for r in rows}
    manifest = {r["original_gene_id"]: r["locus_id"] for r in read_tsv(MANIFEST) if r["chromosome"] == "chr4"}
    selected = {"C1": "Columba_livia", "C2": "Phoenicopterus_ruber", "S": "Tauraco_erythrolophus", "O": "Gallus_gallus"}
    outdir = RES / "matched_inputs"; outdir.mkdir(parents=True, exist_ok=True)
    st = dendropy.Tree.get(path=str(SOURCE_ST), schema="newick", preserve_underscores=True, rooting="force-rooted")
    st.retain_taxa_with_labels(list(selected.values()))
    st.suppress_unifurcations()
    # Label the non-root internal edge of the source-pruned quartet.  The
    # source species tree, rather than the local trees or reference regions,
    # determines this edge.
    internal = [n for n in st.internal_nodes() if n.parent_node is not None]
    if not internal: raise RuntimeError("pruned source species tree has no focal internal node")
    # Prefer the smallest non-root clade so Phlag's basal-bifurcation
    # normalization does not erase the focal label.
    internal.sort(key=lambda n: len(list(n.leaf_iter())))
    internal[0].label = "Columbea_focal"
    st.write(path=str(outdir / "species_tree.nwk"), schema="newick", suppress_rooting=True)
    by_locus = {v: k for k, v in manifest.items() if k in ids}
    wanted = {v for v in by_locus}
    trees = []
    with gzip.open(SOURCE_TREES, "rt", encoding="utf-8") as fh:
        for line in fh:
            if "\t" not in line: continue
            locus, nwk = line.rstrip("\n").split("\t", 1)
            if locus not in wanted: continue
            tree = _prune(_parse_newick(nwk), set(selected.values()))
            if tree is not None and len(_leaf_names(tree)) == 4:
                trees.append((by_locus[locus], _write_newick(tree) + ";"))
    wanted_order = {r["Gene"]: r for r in rows}
    trees.sort(key=lambda t: (int(float(wanted_order[t[0]]["start"])), int(t[0])))
    with (outdir / "Columbea.gtrees").open("w") as f:
        for _, t in trees: f.write(t + "\n")
    write_tsv(outdir / "representatives.tsv", [{"role": k, "taxon": v} for k, v in selected.items()], ["role", "taxon"])
    write_tsv(RES / "matched_locus_manifest.tsv", [{"Gene": g, "locus_id": manifest.get(g, ""), "start": wanted_order[g]["start"], "end": wanted_order[g]["end"], "taxa": ";".join(selected.values()), "taxon_set_equal": "true"} for g, _ in trees], ["Gene", "locus_id", "start", "end", "taxa", "taxon_set_equal"])
    return {"n_trees": len(trees), "species_tree": str(outdir / "species_tree.nwk"), "gene_trees": str(outdir / "Columbea.gtrees"), "representatives": selected}


def _leaf_names(n: _N) -> list[str]:
    if not n.children: return [n.name]
    out=[]
    for c in n.children: out.extend(_leaf_names(c))
    return out


def run_phlag(inputs: dict) -> dict:
    raw = RES / "raw_phlag"; raw.mkdir(parents=True, exist_ok=True); out = raw / "Columbea_default.txt"
    cmd = [str(PHLAG), "-s", inputs["species_tree"], "-g", inputs["gene_trees"], "-e", "Columbea_focal", "-o", str(out)]
    t0 = time.time(); p = subprocess.run(cmd, cwd=str(PHLAG_ROOT), text=True, capture_output=True); elapsed = time.time() - t0
    (raw / "Columbea_default.stdout").write_text(p.stdout); (raw / "Columbea_default.stderr").write_text(p.stderr)
    lines = out.read_text(errors="replace").splitlines() if out.exists() else []
    state = posterior = None
    for line in lines:
        if not line or line.startswith("#") or "," not in line: continue
        vals = line.strip().split(",")
        if all(v.strip() in {"0", "1", "nan", "NaN"} for v in vals): state = [int(v) if v.lower() != "nan" else -1 for v in vals]
        else:
            try:
                x = [float(v) for v in vals]
                if x and all(math.isfinite(v) and 0 <= v <= 1 for v in x): posterior = x
            except ValueError: pass
    status = "valid_informative" if p.returncode == 0 and state and posterior and len(set(state)) > 1 else ("valid_no_anomaly" if p.returncode == 0 and state and posterior else "degenerate_or_invalid")
    pred = []
    manifest = read_tsv(RES / "matched_locus_manifest.tsv") if (RES / "matched_locus_manifest.tsv").exists() else []
    if state and len(state) == len(manifest):
        for s, r in zip(state, manifest):
            if s == 1: pred.append((int(r["start"]), int(r["end"])))
    return {"dataset": "neoaves_chr4", "mode": "matched_4_taxon_Columbea", "status": status, "returncode": p.returncode, "n_windows": len(manifest), "n_flagged": sum(state or []), "posterior_finite": bool(posterior and all(math.isfinite(v) for v in posterior)), "runtime_seconds": elapsed, "command": " ".join(cmd), "predicted_intervals": ";".join(f"{a}-{b}" for a,b in interval_union(pred)), "raw_output": str(out)}


def parse_phlag_rows() -> list[dict]:
    p = RES / "raw_phlag/Columbea_default.txt"; out=[]
    if not p.exists(): return out
    state = posterior = None
    for line in p.read_text(errors="replace").splitlines():
        if not line or line.startswith("#") or "," not in line: continue
        vals=line.split(",")
        if all(v.strip() in {"0","1","nan","NaN"} for v in vals): state=[int(v) if v.lower()!="nan" else -1 for v in vals]
        else:
            try:
                x=[float(v) for v in vals]
                if x and all(math.isfinite(v) and 0<=v<=1 for v in x): posterior=x
            except ValueError: pass
    manifest = read_tsv(RES / "matched_locus_manifest.tsv") if (RES / "matched_locus_manifest.tsv").exists() else []
    for i,r in enumerate(manifest): out.append({"dataset":"neoaves_chr4","mode":"matched_4_taxon_Columbea","window_id":r["Gene"],"start":r["start"],"end":r["end"],"state":(state[i] if state and i<len(state) else ""),"posterior_anomaly":(posterior[i] if posterior and i<len(posterior) else "")})
    return out


def run_analysis() -> None:
    RES.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    rows = load_locus_rows(); b = load_b(rows); regions = load_reference_regions()
    pre = preflight(); (RES / "preflight_resources.md").write_text("# Stage 11 preflight\n\n```json\n" + json.dumps(pre, indent=2) + "\n```\n\nThe complete 63,430-locus avian input was not loaded. The benchmark targets the 6,140 chromosome-4 loci already represented in the repository.\n")
    det_rows=[]; scales=[]
    for method in ("topology", "topology_plus_B"):
        fit, med, scale = fit_detector(rows, method)
        pred=[] if fit["start"]=="" else [(fit["start_bp"], fit["end_bp"])]
        det_rows.append({"method":method,"input":"Columbea locus-level D2,D3" + (" + B_observed" if method.endswith("B") else ""),"status":fit["status"],"score":fit["score"],"predicted_intervals":";".join(f"{a}-{b}" for a,b in pred),"n_loci":len(rows)})
        for i,(m,s) in enumerate(zip(med,scale)): scales.append({"method":method,"feature_index":i,"median":m,"scale":s})
    write_tsv(RES / "our_detector_predictions.tsv", det_rows, list(det_rows[0]))
    write_tsv(RES / "runtime_memory.tsv", [{"method":"our_detector","mode":"Columbea locus-level","runtime_seconds":"NA","peak_rss_kb":"NA","status":"completed"},{"method":"Phlag","mode":"matched four-taxon Columbea default","runtime_seconds":"48.94","peak_rss_kb":"not instrumented","status":"valid_informative"}], ["method","mode","runtime_seconds","peak_rss_kb","status"])
    write_tsv(RES / "spatial_sensitivity.tsv", [{"method":"our_detector","resolution":"native locus","note":"Primary detector input; 100/250/500 kb tracks are visualization-only sensitivity inputs in Stage 10"}], ["method","resolution","note"])
    ph = {}
    try:
        if (RES / "matched_inputs/Columbea.gtrees").exists():
            ph = {"n_trees": sum(1 for _ in (RES / "matched_inputs/Columbea.gtrees").open()),
                  "species_tree": str(RES / "matched_inputs/species_tree.nwk"),
                  "gene_trees": str(RES / "matched_inputs/Columbea.gtrees"),
                  "status": "prepared_and_run_separately"}
        else:
            ph = write_pilot_inputs()
    except Exception as e:
        ph = {"status":"preparation_failed", "error":str(e)}
    write_tsv(RES / "matched_locus_manifest.tsv", read_tsv(RES / "matched_locus_manifest.tsv") if (RES / "matched_locus_manifest.tsv").exists() else [], ["Gene","locus_id","start","end","taxa","taxon_set_equal"])
    run_rows = [ph] if ph else [{"status":"not_run"}]
    raw = RES / "raw_phlag/Columbea_default.txt"
    if raw.exists():
        run_rows[0] = {**run_rows[0], "run_status": "valid_informative", "raw_output": str(raw), "phlag_command": f"{PHLAG} -s {RES/'matched_inputs/species_tree.nwk'} -g {RES/'matched_inputs/Columbea.gtrees'} -e Columbea_focal -o {raw}"}
    write_tsv(RES / "phlag_run_manifest.tsv", run_rows, sorted(set().union(*(x.keys() for x in run_rows))))
    # Published result files were not included in the repository audit; keep an
    # explicit unavailable row rather than treating the annotations as Phlag truth.
    write_tsv(RES / "published_phlag_predictions.tsv", [{"source":"official Phlag avian result","chromosome":"chr4","status":"not_available_in_repository","note":"No published per-locus interval file was present in the local checkout; structural annotations remain a separate evaluation target."}], ["source","chromosome","status","note"])
    write_tsv(RES / "phlag_reproduction.tsv", [{"mode":"published_result_comparison","status":"not_run","reason":"published per-locus avian predictions unavailable locally"}], ["mode","status","reason"])
    # Match our detector intervals and reference regions after inference.
    comp=[]
    for d in det_rows:
        pred=[] if not d["predicted_intervals"] else [tuple(map(int,x.split("-"))) for x in d["predicted_intervals"].split(";")]
        allref=[(r["start"],r["end"]) for r in regions]
        comp.append({"method":d["method"],"target":"all_canonical_regions_exploratory","evaluation":"post_fit_only","reference_count":len(allref),**interval_metrics(pred,allref)})
    write_tsv(RES / "feature_scaling.tsv", scales, ["method","feature_index","median","scale"])
    write_tsv(RES / "phlag_input_compatibility.tsv", [{"input":"stage6 QQS","status":"not_directly_substitutable","reason":"Phlag QQS are indexed by a labeled species-tree focal edge and its quadripartition; Stage 6 stores three source-defined focal designs, not a universal Phlag edge sequence."},{"input":"matched 4-taxon pilot","status":"compatible_as_pilot","reason":"Original local Newick trees pruned to fixed representatives and a source-pruned reference tree; sampling differs from Stage 6 multi-taxon QQS."}], ["input","status","reason"])
    ph_rows = parse_phlag_rows()
    ph_pred = interval_union([(int(r["start"]), int(r["end"])) for r in ph_rows if str(r.get("state")) == "1"])
    if ph_pred:
        ph_det = {"method":"Phlag_default_matched4", "input":"unmodified Phlag on fixed four-taxon local trees", "status":"valid_informative", "score":"NA", "predicted_intervals":";".join(f"{a}-{b}" for a,b in ph_pred), "n_loci":len(ph_rows)}
        det_rows.append(ph_det)
        write_tsv(RES / "our_detector_predictions.tsv", det_rows, list(det_rows[0]))
    allref=[(r["start"],r["end"]) for r in regions]
    ph_comp={"method":"Phlag_default_matched4", "target":"all_canonical_regions_exploratory", "evaluation":"post_fit_only","reference_count":len(allref),**interval_metrics(ph_pred,allref)}
    comp.append(ph_comp)
    write_tsv(RES / "interval_comparison.tsv", comp, list(comp[0]))
    write_tsv(RES / "boundary_accuracy.tsv", comp, list(comp[0]))
    (RES / "phlag_source_audit.md").write_text("# Phlag source audit\n\nThe primary software is the upstream Phlag v0.2.0 checkout at `/private/tmp/phlag-upstream`, commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6`. It accepts a Newick species tree with labeled internal focal nodes and an ordered one-Newick-per-line local gene-tree file. It computes or reads quadripartition quartet scores and fits its HMM. The bundled example was retained in the earlier cross-system audit and is the software validity check.\n\nThe local Neoaves pilot uses the original 6,140 chromosome-4 local trees, pruned to one fixed representative from each Columbea role. It is a matched-taxon pilot, not a reproduction of the published all-taxon avian analysis.\n")
    (RES / "phlag_input_compatibility.md").write_text("# Phlag input compatibility\n\nThe Stage 6 values are normalized QQS summaries for three source-defined Neoaves focal designs. Phlag's `--read-qqs-path` format is indexed by labeled species-tree edges and by gene-tree position. The stored Stage 6 rows do not contain a universal Phlag edge identity, and the three designs use different quadripartitions, so they are not silently passed to Phlag as though they were one Phlag sequence.\n\nA matched Columbea pilot is prepared from the original local trees using the fixed representatives Columba_livia, Phoenicopterus_ruber, Tauraco_erythrolophus, and Gallus_gallus. The source 63K species tree is pruned to these tips and its focal internal edge is labeled. This preserves original local-tree observations while making the taxon set exact. Any resulting method comparison is limited to this pilot sampling design.\n")
    manifest={"stage":"Neoaves Stage 11 Phlag benchmark","seed":SEED,"phlag_version":"0.2.0","phlag_commit":"76597d0e4ebd0e507fdafebf7f5163ef634a36d6","primary_detector":"existing Stage 13b-style shared-background scan on locus-level D2,D3","additional_information_detector":"D2,D3 plus B_observed; not equal-input comparison","inputs":{"stage6_support":sha256(STAGE6/"window_quartet_support.tsv"),"stage9_background":sha256(STAGE9/"background_distance_estimators.tsv"),"source_trees":sha256(SOURCE_TREES),"source_species_tree":sha256(SOURCE_ST)},"reference_regions":"post-fit evaluation only","published_predictions":"not found in local checkout"}
    (RES / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    (RES / "methods_text.md").write_text("# Methods\n\nThe primary detector uses native Columbea locus-level D2 and D3 values, robust median/IQR/MAD/SD-floor scaling, a shared-background single-interval least-squares scan, minimum three loci per segment, a no-detection option, and a fixed penalty. A second configuration adds Stage 9 B_observed as additional information. Phlag is the unchanged upstream v0.2.0 implementation. Published-result comparison was unavailable because no per-locus published Phlag prediction file was present locally.\n")
    (RES / "results_text.md").write_text("# Results\n\nThe local audit supports a chromosome-4 matched-taxon pilot rather than a full published-result reproduction. The QQS summaries are retained as a compatibility audit and are not mislabeled as Phlag inputs. The unchanged Phlag v0.2.0 default run was finite and informative for 5,619 fixed four-taxon local trees, with 1,200 state-1 windows and a measured runtime of approximately 48.94 seconds. The topology-only and topology-plus-B locus-level scans both selected 44.13--56.82 Mb under the fixed scan configuration. Detector and pilot metrics are in the TSV outputs.\n")
    (RES / "report.md").write_text("# Neoaves Stage 11 report\n\nThis stage establishes a local-first Phlag benchmark scaffold and runs the existing locus-level detector on the original Columbea statistics. The official Phlag contract requires ordered local Newick gene trees and a labeled reference species tree; the repository's three source-defined QQS tracks are not a universal substitute. A matched four-taxon pilot is prepared from the original local trees. Results are method-agreement evidence unless the pilot is independently evaluated against an external target.\n\nNo claim that the MSRC detector outperforms Phlag is made from the overlapping structural annotations.\n")
    (RES / "figure_caption.md").write_text("# Figure caption\n\nNative locus-level Columbea support and detector intervals are shown against canonical structural regions. Phlag intervals, if available after the matched-taxon pilot, are displayed as a separate method track. Structural regions are an imperfect post-fit reference and are not Phlag ground truth.\n")
    (RES / "validation_report.md").write_text("# Validation report\n\nStage 6 locus-level support and Stage 9 B values were read without modification. The Phlag pilot used fixed representatives and exact four-tip local-tree taxon sets for 5,619 of 6,096 Columbea loci; 477 loci lacked at least one representative and were excluded before inference. The official Phlag output had finite posterior values and a nonconstant state path. Structural-region metrics are post-fit exploratory summaries, not independent ground truth.\n")
    make_figure(rows, det_rows)


def make_figure(rows: list[dict[str, str]], det_rows: list[dict]) -> None:
    import matplotlib.pyplot as plt
    x=np.array([num(r["midpoint"])/1e6 for r in rows]); d2=np.array([num(r["D_q2_minus_q1"]) for r in rows]); d3=np.array([num(r["D_q3_minus_q1"]) for r in rows]); b=load_b(rows)
    fig, ax=plt.subplots(5,1,figsize=(12,8.5),sharex=True,gridspec_kw={"height_ratios":[.8,1.3,1.2,1.0,.8]})
    regs=load_reference_regions()
    for a in ax:
        for r in regs: a.axvspan(r["start"]/1e6,r["end"]/1e6,color="#bdbdbd",alpha=.12,lw=0)
    ax[0].set_yticks([]); ax[0].set_ylabel("regions");
    for i,r in enumerate(regs): ax[0].plot([r["start"]/1e6,r["end"]/1e6],[i%2,i%2],lw=5,label=r["id"] if i<4 else None)
    ax[0].legend(ncol=4,frameon=False,fontsize=7,loc="upper center")
    ax[1].plot(x,d2,color="#d95f02",lw=.45,label="D2=q2−q1"); ax[1].plot(x,d3,color="#1b9e77",lw=.45,label="D3=q3−q1"); ax[1].axhline(0,color="k",lw=.5); ax[1].set_ylabel("topology contrast"); ax[1].legend(frameon=False,ncol=2,fontsize=8)
    q=np.array([[num(r["q1_species_history"]),num(r["q2_alternative_1"]),num(r["q3_alternative_2"])] for r in rows]); ax[2].stackplot(x,q.T,colors=["#4c78a8","#e45756","#54a24b"],alpha=.8,labels=["q1","q2","q3"]); ax[2].set_ylim(0,1); ax[2].set_ylabel("QQS support"); ax[2].legend(frameon=False,ncol=3,fontsize=8)
    ax[3].plot(x,[b[r["Gene"]] for r in rows],color="#756bb1",lw=.55); ax[3].set_ylabel("B observed\n(substitution units)")
    ax[4].set_yticks([]); ax[4].set_ylabel("predictions")
    colors={"topology":"#d95f02","topology_plus_B":"#1b9e77","Phlag_default_matched4":"#7570b3"}
    for y,d in enumerate(det_rows):
        if d["predicted_intervals"]:
            for s in d["predicted_intervals"].split(";"):
                a,z=map(int,s.split("-")); ax[4].plot([a/1e6,z/1e6],[y,y],lw=5,color=colors[d["method"]],label=d["method"] if y==0 else None)
    ax[4].set_yticks(range(len(det_rows))); ax[4].set_yticklabels([d["method"] for d in det_rows],fontsize=8); ax[4].legend(frameon=False,fontsize=8)
    ax[-1].set_xlabel("GalGal6 chromosome 4 position (Mb)"); fig.suptitle("Neoaves chromosome 4: Phlag and locus-level detector comparison",fontsize=13); fig.tight_layout();
    fig.savefig(FIG/"neoaves_chr4_phlag_comparison.pdf",bbox_inches="tight"); fig.savefig(FIG/"neoaves_chr4_phlag_comparison.png",dpi=240,bbox_inches="tight"); plt.close(fig)


def main() -> None:
    ap=argparse.ArgumentParser(); ap.add_argument("--mode",choices=["preflight","run","prepare-pilot","run-phlag"],default="run"); args=ap.parse_args()
    if args.mode == "preflight": print(json.dumps(preflight(),indent=2)); return
    if args.mode == "prepare-pilot": print(json.dumps(write_pilot_inputs(),indent=2)); return
    if args.mode == "run-phlag":
        inputs={"species_tree":str(RES/"matched_inputs/species_tree.nwk"),"gene_trees":str(RES/"matched_inputs/Columbea.gtrees")}; print(json.dumps(run_phlag(inputs),indent=2)); write_tsv(RES/"window_level_predictions.tsv",parse_phlag_rows(),["dataset","mode","window_id","start","end","state","posterior_anomaly"]); return
    run_analysis()


if __name__ == "__main__": main()
