#!/usr/bin/env python3
"""Stage 7: audit and test independent synteny-state distance contrasts."""
from __future__ import annotations
import csv, gzip, hashlib, json, math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import mean

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
EMP = ROOT / "empirical/neoaves_chr4"
RES = EMP / "results/stage7_arrangement_divergence"
FIG = EMP / "figures"
DATA = ROOT / "data/neoaves_chr4/processed"
TREES = EMP / "external/stiller2024/63430.named.gene.trees.gz"
MANIFEST = EMP / "results/stage4d_locus_manifest.tsv"
LOCUS = DATA / "locus_table_chr4.tsv"
SUPPORT = EMP / "results/stage6_genomic_tracks/window_quartet_support.tsv"
DIST = EMP / "results/stage6_genomic_tracks/window_distance_tracks.tsv"
STATES = DATA / "canonical_arrangement_states.tsv"
REGIONS = DATA / "canonical_structural_regions.tsv"
ROLES = DATA / "stage3b_v2_quartet_role_mapping.tsv"
J48 = EMP / "data/stage4d/jarvis48_taxa.tsv"
CLADES = ["Columbea", "N61", "N62"]
SPECIES = {
    "Stork": "Ciconia_maguari", "Flamingo": "Phoenicopterus_ruber",
    "Dove": "Columba_livia", "Sandgrouse": "Pterocles_gutturalis",
    "Turaco": "Tauraco_erythrolophus", "Cuckoo": "Cuculus_canorus",
}


def read(path):
    with path.open() as h: return list(csv.DictReader(h, delimiter="\t"))


def write(path, data, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        w.writeheader(); w.writerows(data)


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def parse_newick_distances(text, wanted):
    """Return patristic distances for wanted leaves using a small Newick parser."""
    parent, edge, label, children = [], [], [], []
    i = 0
    def new(par=-1):
        k=len(parent); parent.append(par); edge.append(0.0); label.append(""); children.append([]); return k
    def name():
        nonlocal i
        j=i
        while i < len(text) and text[i] not in ",():;": i += 1
        return text[j:i].strip()
    def length(k):
        nonlocal i
        if i < len(text) and text[i] == ':':
            i += 1; j=i
            while i < len(text) and text[i] not in ",();": i += 1
            try: edge[k]=float(text[j:i])
            except ValueError: edge[k]=math.nan
    def parse(par=-1):
        nonlocal i
        k=new(par)
        if i < len(text) and text[i]=='(':
            i += 1; c=parse(k); children[k].append(c)
            while i < len(text) and text[i]==',':
                i += 1; c=parse(k); children[k].append(c)
            if i >= len(text) or text[i]!=')': raise ValueError('unbalanced Newick')
            i += 1
            if i < len(text) and text[i] not in ':,();': label[k]=name()
        else: label[k]=name()
        length(k); return k
    root=parse()
    depth=[0.0]*len(parent)
    for k in range(1,len(parent)):
        depth[k]=depth[parent[k]]+edge[k] if math.isfinite(edge[k]) else math.nan
    leaves={label[k]:k for k in range(len(parent)) if not children[k] and label[k]}
    present=[x for x in wanted if x in leaves]
    paths={}
    for x in present:
        k=leaves[x]; d={}
        while k >= 0: d[k]=depth[k]; k=parent[k]
        paths[x]=d
    out={}
    for ai,a in enumerate(present):
        for b in present[ai+1:]:
            common=set(paths[a]) & set(paths[b])
            if not common: continue
            lca=max(common,key=lambda k: depth[k])
            v=depth[leaves[a]]+depth[leaves[b]]-2*depth[lca]
            if math.isfinite(v): out[(a,b)]=v
    return out, len(leaves)


def main():
    RES.mkdir(parents=True, exist_ok=True); FIG.mkdir(parents=True, exist_ok=True)
    locus=read(LOCUS); manifest={r['original_gene_id']:r for r in read(MANIFEST) if r['chromosome']=='chr4'}
    states=read(STATES); regions=read(REGIONS); roles=read(ROLES)
    state_by_region={}
    for r in states:
        state_by_region.setdefault(r['canonical_region_id'], {})[r['species']]=r
    wanted=list(SPECIES.values())
    jarvis={r['stiller2024_taxon']:r for r in read(J48)}
    mapping=[]
    for common,tip in SPECIES.items():
        mapping.append({'structural_species':common,'local_tree_tip':tip,'jarvis_mapping_status':jarvis.get(tip,{}).get('mapping_status','not_in_47_taxon_mapping'),'local_tree_exact_label':tip,'note':'structural-state label mapped through frozen role definitions'})
    write(RES/'structural_state_taxon_mapping.tsv',mapping,list(mapping[0]))

    eligibility=[]
    for reg in regions:
        rid=reg['canonical_region_id']; s=state_by_region.get(rid,{})
        present=[x for x in SPECIES if x in s]
        counts=Counter(s[x]['normalized_arrangement_state'] for x in present)
        pairs=[(a,b) for ai,a in enumerate(present) for b in present[ai+1:]]
        same=[p for p in pairs if s[p[0]]['normalized_arrangement_state']==s[p[1]]['normalized_arrangement_state']]
        diff=[p for p in pairs if s[p[0]]['normalized_arrangement_state']!=s[p[1]]['normalized_arrangement_state']]
        confident=sum(s[x].get('confidence','')=='confident' for x in present)
        eligible=len(present)>=4 and bool(same) and len(diff)>=2 and len(counts)>=2
        status='exploratory_eligible' if eligible else 'unsupported_insufficient_or_singleton_states'
        if eligible and confident < len(present): status='exploratory_eligible_state_uncertainty'
        eligibility.append({'canonical_region_id':rid,'start':reg['canonical_start'],'end':reg['canonical_end'],'region_type':reg['region_type'],'n_state_taxa':len(present),'state_counts':';'.join(f'{k}:{v}' for k,v in sorted(counts.items())),'n_same_pairs':len(same),'n_different_state_pairs':len(diff),'n_state_groups':len(counts),'n_confident_states':confident,'n_uncertain_states':len(present)-confident,'eligible':eligible,'analysis_status':status,'same_pairs':';'.join(f'{a}-{b}' for a,b in same),'different_pairs':';'.join(f'{a}-{b}' for a,b in diff)})
    write(RES/'region_arrangement_eligibility.tsv',eligibility,list(eligibility[0]))

    # Connect independently defined structural taxa to the frozen quartet roles.
    fmap=[]
    for cl in CLADES:
        rr=[r for r in roles if r['clade']==cl]
        focal=set()
        for r in rr:
            focal |= set(x for x in r.get('structural_species','').split(';') if x)
        for e in eligibility:
            st=state_by_region[e['canonical_region_id']]
            overlap=sorted(focal & set(st))
            fmap.append({'clade':cl,'canonical_region_id':e['canonical_region_id'],'focal_structural_taxa':';'.join(sorted(focal)),'overlap_taxa':';'.join(overlap),'overlap_n':len(overlap),'region_analysis_status':e['analysis_status'],'compatibility':'same independently defined species universe; no qH reassignment' if overlap else 'no focal structural taxon overlap'})
    write(RES/'focal_design_region_mapping.tsv',fmap,list(fmap[0]))

    # Map source tree records to unique frozen genes and calculate six-species pairs.
    unique={r['Gene']:r for r in locus}; tree_metrics={}
    with gzip.open(TREES,'rt') as h:
        for line in h:
            if not line.startswith('chr4_') or '\t' not in line: continue
            loc,nwk=line.rstrip('\n').split('\t',1); key=loc.removesuffix('_fasta.treefile')
            parts=loc.split('_');
            # use the manifest's coordinate/offset identity to recover Gene
            hit=None
            for gene,m in manifest.items():
                expected=f"chr4_{m['start']}_{m['end']}_1k_start{m['locus_start_within_window']}"
                if key==expected: hit=gene; break
            if hit is None or hit not in unique or hit in tree_metrics: continue
            tree_metrics[hit]=parse_newick_distances(nwk,wanted)[0]

    eligible_ids={e['canonical_region_id'] for e in eligibility if e['eligible']}
    distance_rows=[]
    for rid in sorted(eligible_ids):
        st=state_by_region[rid]; pairs=[(a,b) for ai,a in enumerate(SPECIES) for b in list(SPECIES)[ai+1:] if a in st and b in st]
        same=[p for p in pairs if st[p[0]]['normalized_arrangement_state']==st[p[1]]['normalized_arrangement_state']]
        diff=[p for p in pairs if st[p[0]]['normalized_arrangement_state']!=st[p[1]]['normalized_arrangement_state']]
        for gene,r in sorted(unique.items(),key=lambda x:float(x[1]['midpoint'])):
            dm=tree_metrics.get(gene,{})
            sv=[dm[(SPECIES[p[0]],SPECIES[p[1]])] for p in same if (SPECIES[p[0]],SPECIES[p[1]]) in dm]; dv=[dm[(SPECIES[p[0]],SPECIES[p[1]])] for p in diff if (SPECIES[p[0]],SPECIES[p[1]]) in dm]; av=[dm[(SPECIES[p[0]],SPECIES[p[1]])] for p in pairs if (SPECIES[p[0]],SPECIES[p[1]]) in dm]
            if not sv or not dv: continue
            same_bal=mean([mean([dm[(SPECIES[p[0]],SPECIES[p[1]])] for p in same if (SPECIES[p[0]],SPECIES[p[1]]) in dm])])
            # Equal weighting of observed state-pair categories avoids a large
            # state group receiving an extra weight in the sensitivity output.
            diff_groups=defaultdict(list)
            for p in diff:
                if (SPECIES[p[0]],SPECIES[p[1]]) in dm:
                    a,b=p; sa=st[a]['normalized_arrangement_state']; sb=st[b]['normalized_arrangement_state']; diff_groups[tuple(sorted((sa,sb)))].append(dm[(SPECIES[p[0]],SPECIES[p[1]])])
            diff_bal=mean(mean(v) for v in diff_groups.values()) if diff_groups else math.nan
            same_mean=mean(sv); diff_mean=mean(dv); b=mean(av); a=diff_mean-same_mean
            distance_rows.append({'Gene':gene,'canonical_region_id':rid,'start':r['start'],'end':r['end'],'midpoint':r['midpoint'],'n_same_pairs':len(sv),'n_different_state_pairs':len(dv),'n_eligible_pairs':len(av),'mean_same_state':same_mean,'mean_different_state':diff_mean,'A_distance':a,'B_region':b,'A_relative':a/b if b>1e-12 else math.nan,'same_state_balanced':same_bal,'different_state_balanced':diff_bal,'A_balanced':diff_bal-same_bal,'distance_status':'exploratory_multistate_equal_vs_unequal_state'})
    if distance_rows:
        fields=list(distance_rows[0]); write(RES/'window_arrangement_distance.tsv',distance_rows,fields)
        write(RES/'window_matched_background_distance.tsv',distance_rows,['Gene','canonical_region_id','start','end','midpoint','n_eligible_pairs','B_region','distance_status'])
        write(RES/'window_relative_distance.tsv',distance_rows,['Gene','canonical_region_id','start','end','midpoint','A_distance','B_region','A_relative','A_balanced'])
    else:
        for name in ['window_arrangement_distance.tsv','window_matched_background_distance.tsv','window_relative_distance.tsv']:
            (RES/name).write_text('Gene\tcanonical_region_id\tstatus\n')

    sens=[]
    for rid in sorted(eligible_ids):
        z=[r for r in distance_rows if r['canonical_region_id']==rid]
        for field,label in [('A_distance','raw'),('A_balanced','state_pair_balanced')]:
            vals=[float(r[field]) for r in z if math.isfinite(float(r[field]))]
            sens.append({'canonical_region_id':rid,'sensitivity':label,'n_windows':len(vals),'mean':mean(vals) if vals else math.nan,'median':sorted(vals)[len(vals)//2] if vals else math.nan,'min':min(vals) if vals else math.nan,'max':max(vals) if vals else math.nan})
    write(RES/'distance_sensitivity.tsv',sens,list(sens[0]) if sens else ['canonical_region_id','sensitivity','n_windows','mean','median','min','max'])

    # Biological comparison table: merge frozen QQS tracks with region-specific A/B.
    support=read(SUPPORT); comp=[]
    for e in eligibility:
        rid=e['canonical_region_id']; s,eend=float(e['start']),float(e['end']); zd=[r for r in distance_rows if r['canonical_region_id']==rid]
        for cl in CLADES:
            q=[r for r in support if r['clade']==cl]; inside=[r for r in q if s<=float(r['midpoint'])<=eend]; outside=[r for r in q if not(s<=float(r['midpoint'])<=eend)]
            def av(rs,k): return mean(float(x[k]) for x in rs) if rs else math.nan
            zin=[r for r in zd if s<=float(r['midpoint'])<=eend]; zout=[r for r in zd if not (s<=float(r['midpoint'])<=eend)]
            comp.append({'canonical_region_id':rid,'clade':cl,'inside_q1':av(inside,'q1_species_history'),'outside_q1':av(outside,'q1_species_history'),'inside_D2':av(inside,'D_q2_minus_q1'),'outside_D2':av(outside,'D_q2_minus_q1'),'inside_D3':av(inside,'D_q3_minus_q1'),'outside_D3':av(outside,'D_q3_minus_q1'),'n_inside_q':len(inside),'n_outside_q':len(outside),'n_distance_windows':len(zd),'inside_A':av(zin,'A_distance'),'outside_A':av(zout,'A_distance'),'inside_A_relative':av(zin,'A_relative'),'outside_A_relative':av(zout,'A_relative'),'inside_B_region':av(zin,'B_region'),'outside_B_region':av(zout,'B_region'),'analysis_status':e['analysis_status']})
    write(RES/'inside_outside_comparison.tsv',comp,list(comp[0]))

    write_docs(eligibility,mapping,tree_metrics,distance_rows)
    make_figures(support,regions,distance_rows)
    inputs=[LOCUS,STATES,REGIONS,ROLES,TREES,MANIFEST,SUPPORT,DIST]
    outputs=list(RES.glob('*.tsv'))+list(RES.glob('*.md'))
    (RES/'manifest.json').write_text(json.dumps({'stage':'Neoaves Stage 7 arrangement-informed divergence audit','script':str(Path(__file__).relative_to(ROOT)),'reproduction_command':'MPLBACKEND=Agg python empirical/neoaves_chr4/scripts/07_arrangement_divergence.py','inputs':{str(x.relative_to(ROOT)):sha(x) for x in inputs},'n_tree_metrics':len(tree_metrics),'n_arrangement_distance_rows':len(distance_rows),'candidate_regions':sorted(eligible_ids),'arrangement_A_primary_supported':False,'notes':['States are independent synteny/anchor classifications.','A is exploratory only where multistate same-versus-different state pairs are numerically estimable; no binary qH is assigned.'],'output_sha256':{str(x.relative_to(ROOT)):sha(x) for x in outputs}},indent=2)+'\n')


def write_docs(eligibility,mapping,trees,distance):
    counts=Counter(e['analysis_status'] for e in eligibility)
    (RES/'structural_state_provenance.md').write_text(f'''# Structural-state provenance and compatibility audit

The canonical states are region-specific synteny states derived from GalGal6-to-focal-genome alignment anchors and block order/orientation. They are independent of the local ML trees and QQS values used here. They are multistate labels (`A0`, `A1`, ...), not a universal binary arrangement assignment.

There are {len(eligibility)} canonical regions. Eligibility counts are {dict(counts)}. Most regions contain singleton states or low-confidence assignments; only candidate regions with at least four mapped structural species, at least one same-state pair, at least two different-state pairs, and at least two state groups receive exploratory distance rows. State equality is used as the explicit multistate comparison; biologically distinct state labels are not collapsed into a binary orientation.

The six structural species map to exact local-tree labels: {', '.join(f'{r["structural_species"]}={r["local_tree_tip"]}' for r in mapping)}. Ciconia is present in the local trees but is absent from the exact 47-taxon Jarvis mapping used by Stage 6, so Stage 7 records this distinction rather than silently expanding the Stage 6 set. The structural states are species-level and cross-clade; there are no replicate individuals within a state in these trees. Consequently the raw contrasts remain confounded with deep taxonomic divergence and are exploratory rather than a direct analogue of population-level A(w).

The local trees have substitution-unit ML branch lengths. The QQS table contains normalized support but not raw informative quartet counts. Structural-region coordinates and state assignments are inherited from the frozen processed tables; no region was selected using the local topology or distance signal.
''')
    (RES/'validation_report.md').write_text(f'''# Validation report

- Local tree records with six-species pair distances: {len(trees)}.
- Exploratory region-distance rows: {len(distance)}.
- Candidate regions: {', '.join(sorted({e['canonical_region_id'] for e in eligibility if e['eligible']})) or 'none'}.
- Primary arrangement-informed A(w): not accepted as a fully defensible population-level statistic because the states are multistate species-level synteny classes with singleton/low-confidence assignments and no within-state population replication.
- Stage 6 support values and global B(w) were read without modification.
- Pairwise distances use the same eligible species universe for A_r and B_r within each region; incomplete pairs are recorded by counts.
''')
    (RES/'methods_text.md').write_text('''# Methods

Canonical region-specific synteny states were joined to exact local-tree species labels. For each eligible region and locus, `A_distance = mean(different-state pairs) - mean(same-state pairs)` and `B_region` is the mean across the union of those eligible pairs. `A_relative = A_distance/B_region` is missing for near-zero denominators. Because states are multistate, same versus different is an equality comparison; state labels are not collapsed into a binary orientation. A sensitivity estimator gives equal weight to observed state-pair categories. The primary interpretation is exploratory because the comparisons are across deep species-level lineages.
''')
    (RES/'results_text.md').write_text('''# Results

The state tables provide independent, region-specific synteny classifications, but their representation is too sparse and taxonomically deep for a clean population-level arrangement-informed A(w) claim. Candidate region-specific tracks are retained as exploratory sensitivity analyses. The primary figure keeps the three QQS designs separate and retains the Stage 6 global arrangement-blind distance as a supplementary reference.
''')
    (RES/'report.md').write_text(f'''# Neoaves chromosome 4 Stage 7 report

## Scope

This stage audits whether the independently defined structural-state tables support an arrangement-informed divergence contrast analogous to A(w) in the other empirical systems. The local trees are ML trees with substitution-unit branch lengths. They are not calibrated coalescence-time observations.

## Findings

- The structural labels are region-specific, multistate synteny/anchor classes and are independent of the local ML-tree topology used for the QQS tracks.
- Exact structural-species labels were mapped to six local-tree tips. The 47-taxon Jarvis mapping is retained as a separate provenance reference; it does not contain Ciconia maguari.
- Of {len(eligibility)} canonical regions, the numerically eligible exploratory regions are {', '.join(sorted({e['canonical_region_id'] for e in eligibility if e['eligible']})) or 'none'}. The other regions have only singleton or otherwise insufficient state replication.
- The exploratory contrast is `A_distance = mean(different-state pairs) - mean(same-state pairs)` with a matched region-specific background `B_region`; `A_relative = A_distance/B_region`.
- These tracks compare deep species-level states and lack within-state population replication. They therefore do not support a primary population-level A(w) claim and should not be interpreted as direct evidence of recombination suppression.

## Figure recommendation

The topology overview is the primary Stage 7 figure: it keeps Columbea, N61, and N62 as separate source-defined quartet designs and shows both alternative contrasts. Region-specific distance tracks are shown separately as exploratory sensitivity results. The Stage 6 fixed-taxon global B(w) remains the appropriate arrangement-blind reference; no binary qH or universal A(w) was manufactured.
''')
    (RES/'figure_caption.md').write_text('''# Figure caption

Neoaves chromosome-4 topology and structural-state distance audit. Columbea, N61, and N62 are separate source-defined quartet designs; q1 is the species-history resolution and q2/q3 are separate alternatives. Exploratory region-specific A tracks are shown only for numerically eligible multistate synteny regions and are labeled as same-versus-different-state substitution-distance contrasts. They are not calibrated coalescence times or direct evidence of recombination suppression. Structural regions are independently derived synteny regions, not necessarily nucleotide-resolved breakpoints.
''')


def make_figures(support,regions,distance):
    # Primary topology overview: structural annotation, three support rows,
    # and the two contrasts for each source-defined focal design.
    fig,axs=plt.subplots(7,1,figsize=(13,15),sharex=True,gridspec_kw={'height_ratios':[.45,1.0,1.0,1.0,1.0,1.0,1.0]})
    xmax=max(float(r['end']) for r in support)/1e6
    cols={'q1':'#2673B8','q2':'#D65F5F','q3':'#4A9B68'}
    for ax in axs:
        for rr in regions:
            s,e=float(rr['canonical_start'])/1e6,float(rr['canonical_end'])/1e6
            ax.axvspan(s,e,color='#999',alpha=.14,lw=0); ax.axvline(s,color='#666',ls='--',lw=.4); ax.axvline(e,color='#666',ls='--',lw=.4)
    axs[0].set_ylim(0,1); axs[0].set_yticks([]); axs[0].set_ylabel('regions',fontsize=8)
    for rr in regions:
        s,e=float(rr['canonical_start'])/1e6,float(rr['canonical_end'])/1e6; axs[0].bar((s+e)/2,.65,width=e-s,color='#777',alpha=.55,align='center')
    for i,cl in enumerate(CLADES,1):
        z=sorted([r for r in support if r['clade']==cl],key=lambda r:float(r['midpoint']))
        x=[float(r['midpoint'])/1e6 for r in z]
        axs[i].stackplot(x,[float(r['q1_species_history']) for r in z],[float(r['q2_alternative_1']) for r in z],[float(r['q3_alternative_2']) for r in z],colors=[cols['q1'],cols['q2'],cols['q3']],alpha=.85,labels=['q1 species history','q2 alternative 1','q3 alternative 2'])
        axs[i].set_ylim(0,1); axs[i].set_ylabel(f'{cl}\nq support',fontsize=8); axs[i].grid(alpha=.12)
        j=i+3
        axs[j].plot(x,[float(r['D_q2_minus_q1']) for r in z],color=cols['q2'],lw=.5,label='D2=q2−q1')
        axs[j].plot(x,[float(r['D_q3_minus_q1']) for r in z],color=cols['q3'],lw=.5,label='D3=q3−q1')
        axs[j].axhline(0,color='#333',lw=.4); axs[j].set_ylim(-1,1); axs[j].set_ylabel(f'{cl}\nD',fontsize=8); axs[j].grid(alpha=.12)
    axs[-1].set_xlabel('GalGal6 chromosome 4 position (Mb)'); axs[-1].set_xlim(0,xmax)
    axs[1].legend(loc='upper center',bbox_to_anchor=(.5,1.32),ncol=3,frameon=False,fontsize=8)
    axs[4].legend(loc='upper center',bbox_to_anchor=(.5,1.32),ncol=2,frameon=False,fontsize=8)
    fig.suptitle('Neoaves chromosome 4: structural regions, quartet support, and focal contrasts',y=.997,fontsize=13); fig.tight_layout(rect=[0,.02,1,.98])
    fig.savefig(FIG/'neoaves_chr4_stage7_topology_overview.pdf',bbox_inches='tight'); fig.savefig(FIG/'neoaves_chr4_stage7_topology_overview.png',dpi=240,bbox_inches='tight'); plt.close(fig)
    # Exploratory distance figure, only candidate region tracks.
    if distance:
        fig,ax=plt.subplots(figsize=(13,4.5));
        for rid in sorted({r['canonical_region_id'] for r in distance}):
            z=sorted([r for r in distance if r['canonical_region_id']==rid],key=lambda r:float(r['midpoint']))
            x=[float(r['midpoint'])/1e6 for r in z]; ax.plot(x,[float(r['A_distance']) for r in z],lw=.6,label=f'{rid} A distance'); ax.plot(x,[float(r['A_relative']) for r in z],lw=.6,ls='--',label=f'{rid} A relative')
        for rr in regions:
            ax.axvspan(float(rr['canonical_start'])/1e6,float(rr['canonical_end'])/1e6,color='#999',alpha=.1,lw=0)
        ax.axhline(0,color='#333',lw=.5); ax.set_xlabel('GalGal6 chromosome 4 position (Mb)'); ax.set_ylabel('exploratory contrast'); ax.legend(ncol=2,frameon=False,fontsize=8); ax.set_title('Exploratory multistate synteny distance contrasts'); fig.tight_layout(); fig.savefig(FIG/'neoaves_chr4_stage7_exploratory_arrangement_distance.pdf'); fig.savefig(FIG/'neoaves_chr4_stage7_exploratory_arrangement_distance.png',dpi=240); plt.close(fig)
    # Supplementary global background track retained from Stage 6.
    global_distance=read(DIST)
    fig,ax=plt.subplots(figsize=(13,3.8))
    z=sorted([r for r in global_distance if r.get('fixed_mapped_set_complete')=='True' and r.get('B_fixed_mapped') not in ('','nan','NaN')],key=lambda r:float(r['midpoint']))
    if z:
        ax.plot([float(r['midpoint'])/1e6 for r in z],[float(r['B_fixed_mapped']) for r in z],color='#555',lw=.55,label='Stage 6 global B(w)')
    for rr in regions:
        ax.axvspan(float(rr['canonical_start'])/1e6,float(rr['canonical_end'])/1e6,color='#999',alpha=.1,lw=0)
    ax.set_xlabel('GalGal6 chromosome 4 position (Mb)'); ax.set_ylabel('mean patristic distance\n(substitution units)'); ax.set_title('Supplementary arrangement-blind global distance track'); ax.legend(frameon=False,fontsize=8); fig.tight_layout(); fig.savefig(FIG/'neoaves_chr4_stage7_global_B_supplement.pdf'); fig.savefig(FIG/'neoaves_chr4_stage7_global_B_supplement.png',dpi=240); plt.close(fig)


if __name__=='__main__': main()
