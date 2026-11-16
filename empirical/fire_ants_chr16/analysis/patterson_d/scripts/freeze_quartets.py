#!/usr/bin/env python3
"""Freeze a small Stage-1 set of rooted fire-ant D quartet designs."""
import csv, random
from pathlib import Path
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / "empirical/fire_ants_chr16/analysis/patterson_d"
AUDIT = BASE / "results/sample_mapping_audit.tsv"
MAP = ROOT / "data/fire_ants_chr16/processed/stage6c_individual_to_group.tsv"
TREE = ROOT / "data/fire_ants_chr16/processed/stage3_background_tree.nwk"

def read_tsv(path):
    with path.open() as f: return list(csv.DictReader(f, delimiter="\t"))

def main():
    audit = read_tsv(AUDIT)
    mapping = {r["sample"]: r for r in read_tsv(MAP)}
    tree = Phylo.read(str(TREE), "newick")
    random_ids = {"USP3-1-bigB-m-1", "AR223-1-bigB-m-1", "Copa2-1-bigB-m-1", "Par1-1-bigB-m-1", "AR112-1-bigB-p-1", "SS1-2-bigB-m-1"}
    eligible = {}
    for r in audit:
        sid = r["vcf_sample_id"]
        if (r["local_tree_tip"] == "1" and r["caste"].lower() == "male" and
            float(r["source_sample_missing_fraction"]) < .05 and
            r["random_allele_treatment"] == "0" and r["pooled_worker_sample"] == "0" and
            sid not in random_ids):
            group = mapping[sid]["astral_group"] if sid in mapping else r["twisst_group"]
            eligible.setdefault(group, []).append(r)

    # The two primary tests compare the two arrangements within each major species group,
    # using the other species group's Sb population as the third ingroup.
    designs = [
      ("Q1_primary_Argentina", "invicta/macdonaghi_SB", "invicta/macdonaghi_Sb", "richteri_Sb", "saevissima", "Argentina", "all"),
      ("Q2_primary_Argentina", "richteri_SB", "richteri_Sb", "invicta/macdonaghi_Sb", "saevissima", "Argentina", "all"),
      ("Q1_sinvicta_only", "invicta/macdonaghi_SB", "invicta/macdonaghi_Sb", "richteri_Sb", "saevissima", "Argentina", "S.invicta"),
      ("Q2_sinvicta_only", "richteri_SB", "richteri_Sb", "invicta/macdonaghi_Sb", "saevissima", "Argentina", "S.invicta"),
      ("Q1_all_geography", "invicta/macdonaghi_SB", "invicta/macdonaghi_Sb", "richteri_Sb", "saevissima", "all", "all"),
      ("Q2_all_geography", "richteri_SB", "richteri_Sb", "invicta/macdonaghi_Sb", "saevissima", "all", "all"),
      ("Q1_outgroup_Brazil", "invicta/macdonaghi_SB", "invicta/macdonaghi_Sb", "richteri_Sb", "saevissima", "Argentina", "outgroup:Brazil"),
      ("Q1_outgroup_French_Guiana", "invicta/macdonaghi_SB", "invicta/macdonaghi_Sb", "richteri_Sb", "saevissima", "Argentina", "outgroup:French Guiana"),
    ]
    selected=[]; memberships=[]
    for name,p1,p2,p3,o,country,species_filter in designs:
        groups=[p1,p2,p3,o]; chosen=[]
        for g in groups:
            rows=eligible.get(g, [])
            if country != "all" and g != "saevissima": rows=[r for r in rows if r["country"] == country]
            if species_filter == "S.invicta" and g.startswith("invicta/"):
                rows=[r for r in rows if r["species"] == "S.invicta"]
            if species_filter.startswith("outgroup:") and g == "saevissima":
                rows=[r for r in rows if r["country"] == species_filter.split(":",1)[1]]
            if not rows: raise RuntimeError(f"Empty {name} group {g}")
            chosen.append(rows)
        # Unique exact IDs, preserving sorted order.
        labels=(p1,p2,p3,o)
        ids=[sorted({r["vcf_sample_id"] for r in rs}) for rs in chosen]
        all_ids=sum(ids, [])
        if len(all_ids) != len(set(all_ids)): raise RuntimeError(f"Non-disjoint quartet {name}")
        # Root/topology check against independently frozen chromosome 1-15 ASTRAL tree.
        def lca_tips(xs): return {t.name for t in tree.common_ancestor(xs).get_terminals()}
        p12=lca_tips(ids[0]+ids[1]); ing=lca_tips(ids[0]+ids[1]+ids[2])
        topology_ok=not bool(p12.intersection(ids[2]+ids[3])) and not bool(ing.intersection(ids[3]))
        if not topology_ok: raise RuntimeError(f"Baseline rooted quartet conditions fail for {name}")
        selected.append({"quartet_id":name,"P1":p1,"P2":p2,"P3":p3,"O":o,
          "n_P1":len(ids[0]),"n_P2":len(ids[1]),"n_P3":len(ids[2]),"n_O":len(ids[3]),
          "baseline_sister_check":"PASS","outgroup_root_check":"PASS",
          "baseline_tree":"stage3_background_tree.nwk; independent chr1-15 ASTRAL history",
          "hypothesis":"P1/P2 arrangement-stratified groups; P3 contrasting population; D sign reflects excess P1-P3 vs P2-P3 site patterns",
          "rationale":"predeclared arrangement/ecotype contrast; sample membership restricted by metadata, callability, caste and source filtering",
          "confounders":"geography; ancestry/species mixture; linkage; upstream random-allele treatment excluded; outgroup geographic sampling"})
        for role,(g,sids) in enumerate(zip(labels,ids),1):
            for sid in sids:
                memberships.append({"quartet_id":name,"role":f"P{role}" if role<4 else "O","population_group":g,"sample_id":sid})
    out=BASE/"config"; out.mkdir(parents=True,exist_ok=True)
    for fn,rows in [("frozen_quartets.tsv",selected),("frozen_sample_membership.tsv",memberships)]:
        with (out/fn).open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
    print(f"froze {len(selected)} quartets and {len(memberships)} membership rows")

if __name__ == "__main__": main()
