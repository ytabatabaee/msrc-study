#!/usr/bin/env python3
"""Run/rebuild coordinate-indexed local-tree treeness analyses in dataset-local folders."""
from __future__ import annotations
import argparse,csv,gzip,hashlib,io,math,re,statistics,sys
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'empirical'/'shared'))
from treeness_core import read_tsv,write_tsv,fmt,calculate,mark_region,compare_inside_outside

DATA=ROOT/'data';EMP=ROOT/'empirical'

def cod_rows():
    p=EMP/'atlantic_cod/analysis/treeness/results/atlantic_cod_treeness_250kb.tsv'
    source=DATA/'atlantic_cod/processed/cod_window_trees.tsv'
    v={r['window_id']:r for r in read_tsv(p)};src=read_tsv(source);out=[]
    for raw in src:
        r=v[raw['window_id']];valid=r['valid']=='true';m=calculate(raw['tree_newick'],rooted=True)
        if valid and (not math.isclose(float(r['treeness']),m['treeness'],rel_tol=1e-8,abs_tol=1e-8) or not math.isclose(float(r['internal_branch_sum']),m['internal_branch_length_sum'],rel_tol=1e-8,abs_tol=1e-8) or not math.isclose(float(r['total_branch_sum']),m['total_tree_length'],rel_tol=1e-8,abs_tol=1e-8)):
            raise ValueError(f"Atlantic-cod stored treeness failed source-tree validation: {raw['window_id']}")
        total=float(r['total_branch_sum']) if r.get('total_branch_sum') not in ('NA','') else math.nan
        inter=float(r['internal_branch_sum']) if r.get('internal_branch_sum') not in ('NA','') else math.nan
        out.append({'dataset':'atlantic_cod','tree_set':'snapp_mcc_250kb','chromosome':r['lg'],'window_id':r['window_id'],'start':int(r['start']),'end':int(r['end']),'midpoint':int(r['midpoint']),'window_size_bp':int(r['window_size_bp']),'region_class':r['region_class'],'region_id':'cod_suppressed_recombination','treeness':float(r['treeness']) if valid else math.nan,'internal_branch_length_sum':inter,'terminal_branch_length_sum':(total-inter if math.isfinite(total) and math.isfinite(inter) else math.nan),'total_tree_length':total,'max_terminal_branch_length':m['max_terminal_branch_length'],'max_terminal_branch_fraction':m['max_terminal_branch_fraction'],'n_taxa':int(r['n_tips']) if r.get('n_tips') not in ('NA','') else math.nan,'n_internal_branches':m['n_internal_branches'],'n_terminal_branches':m['n_terminal_branches'],'n_zero_length_internal_branches':m['n_zero_length_internal_branches'],'n_zero_length_terminal_branches':m['n_zero_length_terminal_branches'],'n_missing_branch_lengths':m['n_missing_branch_lengths'],'n_invalid_branch_lengths':m['n_invalid_branch_lengths'],'n_valid_nonroot_branches':m['n_valid_nonroot_branches'],'root_stem_length_excluded':m['root_stem_length_excluded'],'n_outgroup_tips':0,'treeness_no_outgroups':float(r['treeness']) if valid else math.nan,'delta_treeness_all_minus_no_outgroups':0 if valid else math.nan,'taxa':m['taxa'],'taxon_set_sha256':m['taxon_set_sha256'],'valid':'true' if valid else 'false','validation_issue':'' if valid else 'negative branch recorded in existing validation','branch_length_units':'Ma','tree_method':'SNAPP BEAST MCC time tree','tree_type':'population tree','rooting':'rooted BEAST/SNAPP crown-age MCC; [&R] marker absent','rooted':True,'data_source':'data/atlantic_cod/processed/cod_window_trees.tsv','region_precision':'validated frozen gadMor2 region manifest','existing_treeness_verified':True})
    return out

def load_anopheles():
    grid={r['window_id']:r for r in read_tsv(DATA/'anopheles_2la/processed/stage2r_window_grid.tsv')}
    species_map={r['haplotype_id']:r['species'] for r in read_tsv(DATA/'anopheles_2la/processed/stage4r_haplotype_species_map.tsv')}
    ints=[('2La',20524058,42165532)];out=[]
    sources=[('stage2r_12sample_local_NJ','stage2r_local_nj_trees.nwk.gz'),('stage4r_all72_haplotype_NJ','stage4r_all72_haplotype_local_trees.nwk.gz')]
    for tree_set,filename in sources:
        p=DATA/'anopheles_2la/processed'/filename
        with gzip.open(p,'rt') as f:
            for line in f:
                parts=line.rstrip('\n').split('\t',1)
                if len(parts)!=2:continue
                wid,nwk=parts
                if wid not in grid:continue
                g=grid[wid];st,en=int(g['start']),int(g['end']);cl,rid=mark_region(st,en,ints)
                if tree_set.startswith('stage4r'):
                    outgroups={tip+'|'+sp for tip,sp in species_map.items() if sp=='quadriannulatus'}
                else:
                    outgroups=set(re.findall(r'[^,():;]*\|quadriannulatus\|[^,():;]*',nwk))
                out.append(base_row('anopheles_2la',tree_set,'2L',wid,st,en,cl,rid,nwk,'published NJ substitution-distance scale (uncalibrated)','unrooted Neighbor-Joining local genealogy',f'data/anopheles_2la/processed/{filename}',outgroups))
    return out

def load_fireants():
    src=read_tsv(DATA/'fire_ants_chr16/processed/stage1_window_trees.tsv');ints=[('In16_3',11732910,16478110),('In16_2',16509882,17962688),('In16_1',17972919,28149888)];out=[]
    for r in src:
        chrom=r['chrom'];st=int(float(r['start']));en=int(float(r['end']));mid=float(r['mid']);
        if chrom=='chr16':cl,rid=mark_region(st,en,ints)
        else:cl,rid='unclassified',''
        out.append(base_row('fire_ants_chr16','RAxML_213_local_gene_trees',chrom,str(r['window_index']),st,en,cl,rid,r['tree_newick'],'source RAxML gene-tree branch-length units; substitution-scale, not time calibrated','unrooted local individual genealogy','data/fire_ants_chr16/processed/stage1_window_trees.tsv',()))
    return out

def load_mouse():
    tre=DATA/'house_mouse_t_complex/processed/house_mouse_t_complex_ml_5kb.tre';meta=read_tsv(DATA/'house_mouse_t_complex/processed/house_mouse_t_complex_ml_5kb_metadata.tsv')
    lines=[x.strip() for x in tre.open() if x.strip()];assert len(lines)==len(meta),(len(lines),len(meta))
    ints=[(r['inversion_name'],int(r['start_bp']),int(r['end_bp'])) for r in read_tsv(EMP/'house_mouse_t_complex/results/t_complex_inversion_intervals.tsv') if r.get('comparison_eligible')=='true'];out=[]
    for nwk,r in zip(lines,meta):
        st,en=int(r['start_bp']),int(r['end_bp']);cl,rid=mark_region(st,en,ints)
        # M. spretus is the explicit outgroup to the three Mus musculus subspecies.
        outgroups=set(re.findall(r'[^,():;]*Ms_SPRE[^,():;]*',nwk))
        out.append(base_row('house_mouse_t_complex','published_IQTREE_ML_5kb','chr17',r['locus_id'],st,en,cl,rid,nwk,'published ML substitution-distance units','unrooted IQ-TREE maximum-likelihood local tree','data/house_mouse_t_complex/processed/house_mouse_t_complex_ml_5kb.tre',outgroups))
    return out

def load_neurospora():
    inv=read_tsv(EMP/'neurospora/results/gene_tree_inventory.tsv');p=DATA/'neurospora/raw/neurospora/trees/cras_sit_hist_tet_L10.min20MAC1.phyml_GTR_bionj_w50sites.trees.gz'
    with gzip.open(p,'rt') as f:trees=[x.strip() for x in f if x.strip()]
    if len(inv)!=len(trees):raise ValueError(f'Neurospora local-tree/window mismatch: {len(trees)} trees vs {len(inv)} coordinate rows')
    out=[];ints=[('LGI_1p5_8Mb',1500000,8000000)]
    for nwk,r in zip(trees,inv):
        st,en=int(r['start']),int(r['end']);cl,rid=mark_region(st,en,ints)
        out.append(base_row('neurospora','L10_50SNP_local_BIONJ','LGI',str(r['window']),st,en,cl,rid,nwk,'PhyML GTR substitution/site scale; 50-SNP window estimate','unrooted BIONJ local genealogy','data/neurospora/raw/neurospora/trees/cras_sit_hist_tet_L10.min20MAC1.phyml_GTR_bionj_w50sites.trees.gz',()))
    return out

def load_sunflowers():
    p=DATA/'sunflowers/raw/10.2kbphysical.trees.txt';out=[]
    with p.open() as f:
        for line in f:
            a=line.rstrip('\n').split('\t',1)
            if len(a)!=2:continue
            key,nwk=a
            m=re.match(r'(Ha412HOChr[^-]+)-(\d+)-(\d+)',key)
            if not m:continue
            chrom='chr'+m.group(1).replace('Ha412HOChr','');st,en=int(m.group(2)),int(m.group(3))
            out.append(base_row('sunflowers','StudyB_IQTREE_10kbphysical',chrom,key,st,en,'unclassified','',nwk,'IQ-TREE 2.0.6 substitutions/site','unrooted maximum-likelihood local tree','data/sunflowers/raw/10.2kbphysical.trees.txt',()))
    return out

def load_heliconius():
    rows=read_tsv(DATA/'heliconius/processed/local_trees.tsv');reg=read_tsv(EMP/'heliconius/config/regions.tsv');out=[]
    for r in rows:
        src=r['source_id'];assembly=r['assembly'];chrom=str(r['chromosome']);st=int(r['start']);en=int(r['end'])
        intervals=[]
        for x in reg:
            if x['status']!='verified' or x['start'] in ('NA','') or x['end'] in ('NA',''):continue
            if src not in x['applies_to_source_ids'].split(','):continue
            if str(x['chromosome'])!=chrom:continue
            # The local trees are never projected between distinct assemblies.
            if x['assembly'] != assembly:continue
            intervals.append((x['region_id'],int(x['start']),int(x['end'])))
        cl,rid=mark_region(st,en,intervals)
        nwk=r['tree_newick']
        outgroups={tip for tip in ('HmelRef','Etal') if re.search(r'(?<![A-Za-z0-9_])'+re.escape(tip)+r'(?=[:,)])',nwk)}
        tree_set=src+'_'+assembly
        out.append(base_row('heliconius',tree_set,chrom,r['window_id'],st,en,cl,rid,nwk,'source-tree patristic units; not converted or normalized','published local ML gene genealogy; source-rooting marker absent','data/heliconius/processed/local_trees.tsv',outgroups))
    return out

def load_neoaves():
    """Map named chromosome-4 gene trees to the frozen Stage-4D loci."""
    tree_path=EMP/'neoaves_chr4/external/stiller2024/63430.named.gene.trees.gz'
    manifest=read_tsv(EMP/'neoaves_chr4/results/stage4d_locus_manifest.tsv')
    by_id={r['locus_id']:r for r in manifest if r['chromosome']=='chr4'}
    regions=read_tsv(DATA/'neoaves_chr4/processed/canonical_structural_regions.tsv')
    intervals=[(r['canonical_region_id'],int(r['canonical_start']),int(r['canonical_end'])) for r in regions if r['reference_chrom']=='chr4']
    out=[];seen=set()
    with gzip.open(tree_path,'rt') as f:
        for line in f:
            a=line.rstrip('\n').split('\t',1)
            if len(a)!=2:continue
            locus_id,nwk=a
            r=by_id.get(locus_id)
            if r is None or locus_id in seen:continue
            seen.add(locus_id)
            st,en=int(r['start']),int(r['end']);cl,rid=mark_region(st,en,intervals)
            out.append(base_row('neoaves_chr4','stiller2024_named_local_ML','chr4',locus_id,st,en,cl,rid,nwk,'published ML substitution-distance units; uncalibrated','published local maximum-likelihood gene tree; support labels embedded in internal node labels','empirical/neoaves_chr4/external/stiller2024/63430.named.gene.trees.gz',()))
    if len(out)<6000:raise ValueError(f'Neoaves exact tree/locus-ID join unexpectedly small: {len(out)}')
    return out

def base_row(dataset,tree_set,chrom,wid,st,en,region_class,region_id,nwk,units,method,source,outgroups):
    row={'dataset':dataset,'tree_set':tree_set,'chromosome':str(chrom),'window_id':str(wid),'start':int(st),'end':int(en),'midpoint':(int(st)+int(en))/2,'window_size_bp':int(en)-int(st)+1,'region_class':region_class,'region_id':region_id,'branch_length_units':units,'tree_method':method,'tree_type':'local genealogy','rooting':'source tree unrooted; internal edges partition tips 2+ versus 2+; root stem excluded','rooted':False,'data_source':source,'input_tree_newick':nwk,'explicit_outgroups':'|'.join(sorted(outgroups))}
    return row

LOADERS={'anopheles_2la':load_anopheles,'fire_ants_chr16':load_fireants,'house_mouse_t_complex':load_mouse,'neurospora':load_neurospora,'sunflowers':load_sunflowers,'heliconius':load_heliconius,'neoaves_chr4':load_neoaves,'atlantic_cod':cod_rows}

def analyze(dataset, plot_all=True):
    rows=LOADERS[dataset]()
    validation=[];processed=[]
    for r in rows:
        nwk=r.pop('input_tree_newick',None);outs=set(r.pop('explicit_outgroups','').split('|')) if r.get('explicit_outgroups') else set()
        try:
            if nwk is not None:
                m=calculate(nwk,outs,rooted=bool(r.get('rooted',False)));r.update(m);r['n_taxa']=m['n_tips'];r['valid']='true' if math.isfinite(m['treeness']) else 'false';r['validation_issue']='' if r['valid']=='true' else ('zero_total_branch_length' if m['n_missing_branch_lengths']==0 and m['n_invalid_branch_lengths']==0 else 'invalid_or_missing_branch_length')
            elif r.get('existing_treeness_verified'):
                if 'root_stem_length_excluded' not in r:r['root_stem_length_excluded']=False
                if 'n_valid_nonroot_branches' not in r:r['n_valid_nonroot_branches']=int(r.get('n_internal_branches',0) or 0)+int(r.get('n_terminal_branches',0) or 0)
            else:raise ValueError('missing_tree_input')
        except Exception as exc:
            r.update({'treeness':math.nan,'internal_branch_length_sum':math.nan,'terminal_branch_length_sum':math.nan,'total_tree_length':math.nan,'n_taxa':math.nan,'n_internal_branches':math.nan,'n_terminal_branches':math.nan,'n_zero_length_internal_branches':0,'n_zero_length_terminal_branches':0,'n_missing_branch_lengths':0,'n_invalid_branch_lengths':1,'n_outgroup_tips':0,'treeness_no_outgroups':math.nan,'delta_treeness_all_minus_no_outgroups':math.nan,'taxa':'','taxon_set_sha256':'','n_valid_nonroot_branches':0,'root_stem_length_excluded':False,'valid':'false','validation_issue':str(exc)})
        # Summaries of long terminal branches and tree-completeness sensitivity.
        total=r.get('total_tree_length');term=r.get('terminal_branch_length_sum')
        r['terminal_branch_share']=(float(term)/float(total) if total not in (None,0) and math.isfinite(float(total)) else math.nan)
        processed.append(r)
        validation.append({k:r.get(k) for k in ['dataset','tree_set','chromosome','window_id','start','end','region_class','region_id','n_taxa','taxon_set_sha256','n_missing_branch_lengths','n_invalid_branch_lengths','n_zero_length_internal_branches','n_zero_length_terminal_branches','n_outgroup_tips','treeness','treeness_no_outgroups','delta_treeness_all_minus_no_outgroups','root_stem_length_excluded','valid','validation_issue']})
    base=EMP/dataset/'analysis/treeness';res=base/'results';fig=base/'figures';res.mkdir(parents=True,exist_ok=True);fig.mkdir(parents=True,exist_ok=True)
    fields=['dataset','tree_set','chromosome','window_id','start','end','midpoint','window_size_bp','region_class','region_id','treeness','internal_branch_length_sum','terminal_branch_length_sum','total_tree_length','terminal_branch_share','n_taxa','n_internal_branches','n_terminal_branches','n_zero_length_internal_branches','n_zero_length_terminal_branches','max_terminal_branch_length','max_terminal_branch_fraction','n_missing_branch_lengths','n_invalid_branch_lengths','n_valid_nonroot_branches','n_outgroup_tips','treeness_no_outgroups','delta_treeness_all_minus_no_outgroups','taxa','taxon_set_sha256','root_stem_length_excluded','branch_length_units','tree_method','tree_type','rooting','data_source','valid','validation_issue']
    write_tsv(res/'treeness_windows.tsv',processed,fields);write_tsv(res/'tree_validation.tsv',validation)
    comparisons=compare_inside_outside(processed, reps=2000, seed=20261124)
    write_tsv(res/'inside_outside_summary.tsv',comparisons)
    # Branch component differences mirror the region comparison and clarify mechanism.
    components=[]
    for c in comparisons:
        components.append({k:c.get(k) for k in ['tree_set','chromosome','region_id','n_inside','n_outside','median_internal_inside','median_internal_outside','delta_internal_inside','median_terminal_inside','median_terminal_outside','delta_terminal_inside','median_treeness_inside','median_treeness_outside','delta_median_treeness']})
    write_tsv(res/'branch_component_comparison.tsv',components)
    # Taxon/root audit at tree-set/chromosome granularity.
    audits=[]
    groups=defaultdict(list)
    for r in processed:groups[(r['tree_set'],r['chromosome'])].append(r)
    for (ts,ch),g in groups.items():
        tax=Counter(r['taxon_set_sha256'] for r in g if r['taxon_set_sha256']);modal=max(tax,key=tax.get) if tax else ''
        audits.append({'tree_set':ts,'chromosome':ch,'n_trees':len(g),'n_valid':sum(r['valid']=='true' for r in g),'n_taxon_sets':len(tax),'modal_taxon_set_n':tax.get(modal,0),'modal_taxon_set_fraction':tax.get(modal,0)/len(g) if g else math.nan,'min_taxa':min((int(r['n_taxa']) for r in g if r['n_taxa'] not in (None,'NA') and math.isfinite(float(r['n_taxa']))),default=0),'max_taxa':max((int(r['n_taxa']) for r in g if r['n_taxa'] not in (None,'NA') and math.isfinite(float(r['n_taxa']))),default=0),'has_named_outgroups':any(int(r['n_outgroup_tips'] or 0)>0 for r in g),'median_max_terminal_branch_fraction':statistics.median([float(r['max_terminal_branch_fraction']) for r in g if r.get('max_terminal_branch_fraction') is not None and math.isfinite(float(r['max_terminal_branch_fraction']))]) if any(r.get('max_terminal_branch_fraction') is not None and math.isfinite(float(r['max_terminal_branch_fraction'])) for r in g) else math.nan,'median_treeness_outgroup_pruned':statistics.median([float(r['treeness_no_outgroups']) for r in g if r['treeness_no_outgroups'] is not None and math.isfinite(float(r['treeness_no_outgroups']))]) if any(r['treeness_no_outgroups'] is not None and math.isfinite(float(r['treeness_no_outgroups'])) for r in g) else math.nan,'median_delta_all_minus_no_outgroups':statistics.median([float(r['delta_treeness_all_minus_no_outgroups']) for r in g if r['delta_treeness_all_minus_no_outgroups'] is not None and math.isfinite(float(r['delta_treeness_all_minus_no_outgroups']))]) if any(r['delta_treeness_all_minus_no_outgroups'] is not None and math.isfinite(float(r['delta_treeness_all_minus_no_outgroups'])) for r in g) else math.nan})
    write_tsv(res/'tree_set_audit.tsv',audits)
    if plot_all:plot_tracks(dataset,processed,fig)
    write_readme(dataset,processed,comparisons,audits,base)
    return processed,comparisons,audits

def plot_tracks(dataset,rows,figdir):
    by=defaultdict(list)
    for r in rows:by[(r['tree_set'],r['chromosome'])].append(r)
    palette={'inside':'#245b78','outside':'#778899','boundary':'#c27a2c','unclassified':'#7d3c98'}
    for (ts,chrom),g in by.items():
        g=sorted(g,key=lambda x:int(x['start']));pdf=figdir/f'{ts}_{chrom}_treeness_tracks.pdf'
        with PdfPages(pdf) as pp:
            fig,axs=plt.subplots(3,1,figsize=(11,8),sharex=True)
            for ax,measure,label in zip(axs,['treeness','internal_branch_length_sum','terminal_branch_length_sum'],['Treeness (dimensionless)','Sum of internal branch lengths','Sum of terminal branch lengths']):
                for cl in sorted(set(r['region_class'] for r in g)):
                    valid=[r for r in g if r['region_class']==cl and r.get(measure) is not None and math.isfinite(float(r[measure]))]
                    if not valid:continue
                    # Break lines at missing windows or large coordinate gaps.
                    run=[];last=None
                    widths=[int(r['window_size_bp']) for r in valid];gap=max(statistics.median(widths)*2,1)
                    for r in valid:
                        st=int(r['start'])
                        if run and st-last>gap:
                            ax.plot([float(x['midpoint'])/1e6 for x in run],[float(x[measure]) for x in run],'.-',ms=2,lw=.7,color=palette.get(cl,'#7d3c98'),label=cl if run[0] is valid[0] else None);run=[]
                        run.append(r);last=int(r['end'])
                    if run:ax.plot([float(x['midpoint'])/1e6 for x in run],[float(x[measure]) for x in run],'.-',ms=2,lw=.7,color=palette.get(cl,'#7d3c98'),label=cl if run[0] is valid[0] else None)
                # Amber region shading only when validated interval classes exist; shade plotted inside span.
                inside=[r for r in g if r['region_class']=='inside']
                if inside:
                    ax.axvspan(min(int(r['start']) for r in inside)/1e6,max(int(r['end']) for r in inside)/1e6,color='#e6ad4b',alpha=.16)
                ax.set_ylabel(label);ax.grid(alpha=.18)
            axs[-1].set_xlabel(f'{chrom} genomic position (Mb; source assembly coordinates)')
            fig.suptitle(f'{dataset} | {ts} | {chrom}: local-tree treeness and branch components');fig.tight_layout();pp.savefig(fig);fig.savefig(figdir/f'{ts}_{chrom}_treeness_tracks.png',dpi=170);plt.close(fig)

def write_readme(dataset,rows,comparisons,audits,base):
    valid=sum(r['valid']=='true' for r in rows);bad=sum(r['valid']!='true' for r in rows);sources=sorted(set(r['tree_set'] for r in rows));masks=sorted(set(r['region_id'] for r in rows if r['region_id']))
    text=f'''# {dataset} treeness analysis\n\n## Method\n\nTreeness is internal branch-length sum divided by total tree length, following the Systematic Biology CASTLES/CoalBL definition. For rooted trees, internal edges are non-root edges leading to nonterminal nodes; for unrooted trees, internal edges split at least two retained tips on each side. Terminal edges lead to one tip. Root stems are excluded; no missing edge length is set to zero. Zero-length edges remain zero. Missing/nonfinite/negative branches are counted and trees with any such branch have NA treeness; a zero-total tree also has NA treeness. Internal and terminal branch sums are retained separately.\n\n## Inputs and provenance\n\nTree sets: {', '.join(sources)}. Coordinate-indexed local trees processed: {len(rows)}; valid: {valid}; invalid/zero-total: {bad}. The per-window file records tree set, assembly chromosome/coordinates, branch-length units, tree method/type, rooting convention, taxon-set hash, branch-length counts, zero internal/terminal branches, and explicit outgroup sensitivity where taxon labels supported it. No species-wide reference tree is used. Posterior tree samples were not available as coordinate-matched samples in the included local-tree sources; these are point trees, not posterior summaries.\n\nSuppressed-recombination masks used when validated and source-compatible: {', '.join(masks) if masks else 'none for this source; no inside/outside contrast is calculated'}. Boundaries are not transferred across assemblies or source-specific coordinate axes. Region labels from unvalidated/approximate masks are not treated as validated comparisons.\n\n## Comparisons\n\nSee `results/inside_outside_summary.tsv` and `results/branch_component_comparison.tsv`. Contrasts are descriptive median/mean differences. Mean-difference intervals use a 2-Mb moving-block bootstrap over genomic runs; windows are not treated as independent biological replicates. Inference is not reported where there are fewer than three observations in either class. Multiple chromosomes/source groups remain separate. Trees with nonmodal taxon sets are excluded from inside/outside estimates and counted in `n_excluded_nonmodal_taxa`; full tracks retain every valid tree.\n\nOutgroup-pruned values are provided only where the outgroup was explicitly identifiable from source labels/metadata; otherwise no outgroup is guessed. Rooted/unrooted sensitivity can be assessed by the unrooted edge-split rule; root stem exclusion is invariant. Branch-length rescaling is not applied. See `tree_set_audit.tsv` for taxon consistency and outgroup influence.\n\n## Reproduction\n\nFrom repository root: `python3 empirical/{dataset}/analysis/treeness/scripts/run_treeness.py`. Shared math is in `empirical/shared/treeness_core.py`; per-dataset results and figures remain here. Cross-dataset summary is `empirical/results/treeness_cross_dataset_summary.tsv`. Tests: `python3 -m unittest discover -s empirical/shared/tests`.\n'''
    (base/'README.md').write_text(text)

def main():
    p=argparse.ArgumentParser();p.add_argument('--dataset',choices=['all',*LOADERS],default='all');a=p.parse_args()
    names=list(LOADERS) if a.dataset=='all' else [a.dataset];summaries=[]
    for d in names:
        rows,cmp,aud=analyze(d)
        for c in cmp:
            summaries.append({'dataset':d,**c,'tree_type':next((x['tree_type'] for x in rows if x['tree_set']==c['tree_set']),''),'branch_length_units':next((x['branch_length_units'] for x in rows if x['tree_set']==c['tree_set']),''),'valid_windows_in_group':sum(r['valid']=='true' and r['tree_set']==c['tree_set'] and r['chromosome']==c['chromosome'] for r in rows),'direction_median':'increase_inside' if c['delta_median_treeness']>0 else ('decrease_inside' if c['delta_median_treeness']<0 else 'no_difference'),'important_limitations':'see dataset-local README; source-specific taxon/root/units and region masks'})
        if not cmp and rows:
            for ts in sorted(set(x['tree_set'] for x in rows)):
                grp=[x for x in rows if x['tree_set']==ts]
                summaries.append({'dataset':d,'tree_set':ts,'chromosome':'multiple' if len(set(x['chromosome'] for x in grp))>1 else grp[0]['chromosome'],'region_id':'NA_no_validated_mask','n_inside':0,'n_outside':0,'median_treeness_inside':math.nan,'median_treeness_outside':math.nan,'delta_median_treeness':math.nan,'mean_treeness_inside':math.nan,'mean_treeness_outside':math.nan,'delta_mean_treeness':math.nan,'valid_windows_in_group':sum(x['valid']=='true' for x in grp),'direction_median':'not_estimable','important_limitations':'coordinate-indexed trees available; no source-compatible validated suppressed-recombination boundary'})
        print(d,len(rows),'trees',sum(r['valid']=='true' for r in rows),'valid',len(cmp),'inside/outside groups')
    if a.dataset=='all':
        fields=list(summaries[0]) if summaries else []
        write_tsv(EMP/'results/treeness_cross_dataset_summary.tsv',summaries,fields)
        inventory=[{'dataset':d,'analysis_status':'analyzed' if d in names else 'not analyzed','notes':note} for d,note in [
          ('atlantic_cod','426 published SNAPP MCC population trees, 250-kb windows; reused existing treeness track'),
          ('anopheles_2la','978 Stage-4R all-haplotype local NJ trees on 2L; 50-kb grid'),
          ('fire_ants_chr16','213 RAxML local individual gene trees with coordinates; chr1 and chr16'),
          ('house_mouse_t_complex','4,046 local ML trees, 5-kb chr17 windows'),
          ('neurospora','L10 local BIONJ trees linked by existing window inventory; LGI, 50 SNPs/tree'),
          ('sunflowers','3,242 coordinate-indexed local IQ-TREE trees; no validated precise suppressed-region mask'),
          ('heliconius','57,488 coordinate-indexed local trees in source-specific assemblies; verified 2b masks only where source-compatible'),
          ('neoaves_chr4','not analyzed: deposited tree-support names/coordinate inventory do not include coordinate-indexed branch-length Newicks; existing topology summaries cannot supply branch lengths'),
          ('primula','not analyzed: gene trees exist but source tree-to-coordinate crosswalk is explicitly unresolved'),
          ('microbotryum','not analyzed: gene trees exist but are not coordinate-indexed local genomic windows'),
          ('monkeyflowers','not analyzed: window-weight products exist but no coordinate-indexed local phylogenetic tree file found'),
          ('songbirds_pips','not analyzed: gene-tree code exists but no local coordinate-indexed tree set currently available'),
        ]]
        write_tsv(EMP/'results/treeness_dataset_inventory.tsv',inventory)
        (EMP/'results/treeness_cross_dataset_report.md').write_text(cross_report(summaries,inventory))

def cross_report(summaries,inventory):
    lines=['# Cross-dataset local-tree treeness audit','', 'Treeness is calculated within each tree set, but absolute values are not compared across datasets because taxon sampling, units, tree types, and inference procedures differ. Inferences are local genomic-window contrasts, not independent biological replicates.','', '| Dataset / source / chromosome | valid windows | median inside | median outside | difference | direction | limits |','|---|---:|---:|---:|---:|---|---|']
    for r in summaries:lines.append(f"| {r['dataset']} / {r['tree_set']} / {r['chromosome']} | {r['valid_windows_in_group']} | {fmt(r['median_treeness_inside'])} | {fmt(r['median_treeness_outside'])} | {fmt(r['delta_median_treeness'])} | {r['direction_median']} | {r['important_limitations']} |")
    lines += ['', '## Dataset discovery / unavailable analyses','', '| Dataset | Status | Audit note |','|---|---|---|']
    for r in inventory:lines.append(f"| {r['dataset']} | {r['analysis_status']} | {r['notes']} |")
    lines += ['', 'Details, figures, source provenance, and limitations remain in each dataset’s `analysis/treeness/` folder. The Atlantic-cod analysis reuses and checks the previously generated MCC treeness output. Sunflower has coordinate-indexed trees but the available recombination-suppression boundaries are rounded or unresolved, so tracks are produced without a validated inside/outside estimate. Heliconius uses only validated, source-compatible 2b intervals; 15b on the H. erato assembly is unresolved and excluded.']
    return '\n'.join(lines)+'\n'

if __name__=='__main__':main()
