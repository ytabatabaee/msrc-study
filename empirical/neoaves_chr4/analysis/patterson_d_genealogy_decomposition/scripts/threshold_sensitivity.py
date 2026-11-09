#!/usr/bin/env python3
"""Lightweight informative-site threshold sensitivity using cached predictions."""
from pathlib import Path
import importlib.util

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("decompose", HERE / "decompose.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

rows = m.read_tsv(m.RES / "locus_genealogy_predictions.tsv")
out = []
for design in m.DESIGNS:
    for filt in m.FILTERS:
        for threshold in (0, 10, 50, 100):
            keep = [r for r in rows if r["design"] == design and r["topology_status"] == "resolved" and m.num(r.get(f"observed_{filt}_informative")) is not None and m.num(r[f"observed_{filt}_informative"]) >= threshold]
            inside = m.pred_from_rows([r for r in keep if r["region"] != "outside_all_PNAS"], design, filt)
            outside = m.pred_from_rows([r for r in keep if r["region"] == "outside_all_PNAS"], design, filt)
            out.append({"analysis": "informative_site_threshold", "design": design, "filter": filt, "resolution": "locus", "min_informative_sites": threshold, "region": "all_PNAS", "D_site_inside": inside["D_site"], "D_site_outside": outside["D_site"], "D_site_inside_minus_outside": None if inside["D_site"] is None or outside["D_site"] is None else inside["D_site"] - outside["D_site"], "D_topology_inside": inside["D_topology"], "D_topology_outside": outside["D_topology"], "D_tree_inside": inside["D_tree"], "D_tree_outside": outside["D_tree"], "n_inside": inside["n_prediction_loci"], "n_outside": outside["n_prediction_loci"]})
m.write_tsv(m.RES / "informative_site_threshold_sensitivity.tsv", out, list(out[0]))
