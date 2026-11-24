#!/usr/bin/env python3
"""PCA, clustering, and time-weighted rooted-triplet analyses of cod SNAPP MCCs."""
from __future__ import annotations
import csv, itertools, math, random, statistics, sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import Phylo
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score, adjusted_rand_score

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
ROOT=BASE.parents[3]
DATA=ROOT/'data/atlantic_cod/processed/cod_window_trees.tsv'
BASELINE=ROOT/'data/atlantic_cod/processed/baseline_population_tree_12taxa.nwk'
OUT=BASE/'results';FIG=BASE/'figures';
sys.path.insert(0,str(HERE))
from time_tree_core import extract_mrca_ages, triplet_features

LGS=('LG01','LG02','LG07','LG12')
BLOCK_BP=2_000_000
RNG_SEED=20261124
COLORS={'inside':'#b43c5b','outside':'#326b9a','boundary':'#999999'}

def read_tsv(path):
    with path.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(path,rows,fields=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    if fields is None: fields=list(rows[0]) if rows else []
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader()
        for r in rows:w.writerow({k:fmt(r.get(k)) for k in fields})

def fmt(x):
    if x is None:return 'NA'
    if isinstance(x,(float,np.floating)):
        return 'NA' if not math.isfinite(float(x)) else f'{float(x):.10g}'
    return x

def boolish(x):return str(x).strip().lower() in {'true','1','yes'}

def local_tree_table():
    src=read_tsv(DATA); expected=None;tree_rows=[]; audits=[]
    for raw in src:
        taxa_here=tuple(sorted(x.name for x in Phylo.read(__import__('io').StringIO(raw['tree_newick']),'newick').get_terminals()))
        if expected is None:expected=taxa_here
        f=extract_mrca_ages(raw['tree_newick'],expected_taxa=expected)
        tree=f.get('tree')
        branch_lengths=[float(c.branch_length) for c in tree.find_clades() if c is not tree.root and c.branch_length is not None] if tree else []
        zeros=sum(x==0 for x in branch_lengths)
        inside=boolish(raw['inside_inversion']);boundary=boolish(raw['overlaps_inversion']) and not inside
        region='inside' if inside else ('boundary' if boundary else 'outside')
        row={k:raw[k] for k in ('lg','window_id','start','end','midpoint','window_size_bp','inversion_id')}
        row.update({'region_class':region,'root_age':f['root_age'],'ultrametric_spread':f['ultrametric_spread'],'n_taxa':len(f['taxa']),'n_zero_nonroot_branches':zeros,'valid':f['valid'],'validation_issue':f['issue'],'source_tree_newick':raw['tree_newick'],'tree':tree,'pair_ages':f['pair_ages']})
        tree_rows.append(row)
        audits.append({k:row[k] for k in ('lg','window_id','start','end','region_class','root_age','ultrametric_spread','n_taxa','n_zero_nonroot_branches','valid','validation_issue')})
    return tree_rows,expected,audits

def pair_names(taxa):return [f'{a}__{b}' for a,b in itertools.combinations(taxa,2)]

def bh_adjust(pvalues):
    vals=np.asarray(pvalues,float);n=len(vals);order=np.argsort(vals);q=np.empty(n);running=1.
    for rank in range(n,0,-1):
        ix=order[rank-1];running=min(running,vals[ix]*n/rank);q[ix]=running
    return q

def shifted_effect(values,mask):
    values=np.asarray(values,float);mask=np.asarray(mask,bool)
    a=values[mask];b=values[~mask]
    if len(a)==0 or len(b)==0:return math.nan
    return float(np.nanmean(a)-np.nanmean(b))

def circular_null(values,mask):
    """Enumerate cyclic shifts of a contiguous-region mask against a fixed track."""
    values=np.asarray(values,float);mask=np.asarray(mask,bool);n=len(mask)
    obs=shifted_effect(values,mask)
    null=np.asarray([shifted_effect(values,np.roll(mask,k)) for k in range(n)],float)
    null=null[np.isfinite(null)]
    if not math.isfinite(obs) or not len(null):return obs,math.nan,math.nan,math.nan
    p=(np.sum(np.abs(null)>=abs(obs)-1e-15))/len(null)
    return obs,float(p),float(np.quantile(null,.025)),float(np.quantile(null,.975))

def stratified_circular_null(values,mask,groups,reps=20000,seed=RNG_SEED):
    """Monte Carlo null with an independent cyclic shift inside each LG."""
    values=np.asarray(values,float);mask=np.asarray(mask,bool);groups=np.asarray(groups)
    ids=[np.where(groups==g)[0] for g in LGS]
    def effect(m):
        es=[shifted_effect(values[ix],m[ix]) for ix in ids]
        return float(np.mean(es)) if all(math.isfinite(x) for x in es) else math.nan
    obs=effect(mask);rng=np.random.default_rng(seed);null=np.empty(reps)
    for b in range(reps):
        shifted=mask.copy()
        for ix in ids:shifted[ix]=np.roll(mask[ix],int(rng.integers(0,len(ix))))
        null[b]=effect(shifted)
    null=null[np.isfinite(null)]
    p=(np.sum(np.abs(null)>=abs(obs)-1e-15)+1)/(len(null)+1) if len(null) and math.isfinite(obs) else math.nan
    return obs,float(p),float(np.quantile(null,.025)),float(np.quantile(null,.975))

def prepare_vectors(trees,taxa):
    pairs=list(itertools.combinations(taxa,2));names=pair_names(taxa);valid=[r for r in trees if r['valid']]
    absolute=np.asarray([[r['pair_ages'][p] for p in pairs] for r in valid],float)
    roots=np.asarray([r['root_age'] for r in valid],float)
    relative=absolute/roots[:,None]
    for i,r in enumerate(valid):
        r['vector_index']=i;r['feature_vector_abs']=absolute[i];r['feature_vector_rel']=relative[i]
    return valid,pairs,names,absolute,relative,roots

def fit_pcas(valid,names,absolute,relative,roots):
    result={};scores=[];loadings=[];summary=[]
    for analysis,X in [('absolute_time',absolute),('relative_time',relative),('relative_feature_zscore_sensitivity',StandardScaler().fit_transform(relative))]:
        pca=PCA(n_components=min(X.shape),svd_solver='full').fit(X);Z=pca.transform(X);result[analysis]=(pca,Z)
        for i,r in enumerate(valid):scores.append({'analysis':analysis,'lg':r['lg'],'window_id':r['window_id'],'start':r['start'],'end':r['end'],'midpoint':r['midpoint'],'region_class':r['region_class'],'root_age':roots[i],'PC1':Z[i,0],'PC2':Z[i,1]})
        for pc in range(min(5,Z.shape[1])):
            for j,name in enumerate(names):loadings.append({'analysis':analysis,'PC':pc+1,'feature':name,'loading':pca.components_[pc,j],'abs_loading':abs(pca.components_[pc,j]),'explained_variance_ratio':pca.explained_variance_ratio_[pc]})
        corr=pearsonr(Z[:,0],roots) if len(roots)>2 and np.std(Z[:,0])>0 else (math.nan,math.nan)
        summary.append({'analysis':analysis,'n_windows':len(valid),'n_features':len(names),'PC1_variance':pca.explained_variance_ratio_[0],'PC2_variance':pca.explained_variance_ratio_[1],'PC1_root_age_pearson_r':corr[0],'PC1_root_age_p':corr[1],'scaling':'no feature scaling' if analysis!='relative_feature_zscore_sensitivity' else 'z-scored features (sensitivity only)'})
    write_tsv(OUT/'pca_scores.tsv',scores);write_tsv(OUT/'pca_loadings.tsv',loadings);write_tsv(OUT/'pca_summary.tsv',summary)
    # Pair features that dominate each component, retained separately for easy reading.
    top=[]
    for a in result:
        pca,_=result[a]
        for pc in (0,1):
            ids=np.argsort(np.abs(pca.components_[pc]))[::-1][:15]
            for rank,j in enumerate(ids,1):top.append({'analysis':a,'PC':pc+1,'rank':rank,'feature':names[j],'loading':pca.components_[pc,j],'explained_variance_ratio':pca.explained_variance_ratio_[pc]})
    write_tsv(OUT/'pca_top_population_pairs.tsv',top)
    return result,scores,summary

def pca_figures(pcas,scores,valid):
    # Absolute and normalized PCA: 2x2 linkage-group panels.
    for analysis,title in [('absolute_time','Absolute MRCA-age PCA'),('relative_time','Root-age-normalized MRCA-age PCA')]:
        subset=[r for r in scores if r['analysis']==analysis]
        fig,axs=plt.subplots(2,2,figsize=(10,8),sharex=False,sharey=False)
        for ax,lg in zip(axs.flat,LGS):
            rows=[r for r in subset if r['lg']==lg]
            for cl in ('outside','inside','boundary'):
                q=[r for r in rows if r['region_class']==cl]
                if q:ax.scatter([r['PC1'] for r in q],[r['PC2'] for r in q],s=18,c=COLORS[cl],label=cl,alpha=.82,edgecolor='none')
            ax.set_title(lg);ax.set_xlabel('PC1');ax.set_ylabel('PC2');ax.grid(alpha=.2)
        axs[0,0].legend(frameon=False,fontsize=8)
        pca=pcas[analysis][0];fig.suptitle(f'{title} | PC1 {pca.explained_variance_ratio_[0]:.1%}, PC2 {pca.explained_variance_ratio_[1]:.1%}')
        fig.tight_layout();fig.savefig(FIG/f'{analysis}_pca_scatter.png',dpi=220);fig.savefig(FIG/f'{analysis}_pca_scatter.pdf');plt.close(fig)
    for analysis in ('absolute_time','relative_time'):
        sub=[r for r in scores if r['analysis']==analysis]
        fig,axs=plt.subplots(4,2,figsize=(12,10),sharex=False)
        for i,lg in enumerate(LGS):
            rows=sorted((r for r in sub if r['lg']==lg),key=lambda r:int(r['start']))
            ints=[r for r in rows if r['region_class']=='inside']
            for j,pc in enumerate(('PC1','PC2')):
                ax=axs[i,j]
                for cl in ('outside','inside','boundary'):
                    q=[r for r in rows if r['region_class']==cl]
                    ax.scatter([float(r['midpoint'])/1e6 for r in q],[r[pc] for r in q],s=13,c=COLORS[cl],label=cl,alpha=.8)
                if ints:ax.axvspan(min(float(r['start']) for r in ints)/1e6,max(float(r['end']) for r in ints)/1e6,color='#d8a12a',alpha=.14)
                ax.set_ylabel(f'{lg} {pc}');ax.grid(alpha=.18)
        axs[-1,0].set_xlabel('Position (Mb)');axs[-1,1].set_xlabel('Position (Mb)');axs[0,0].legend(frameon=False,ncol=3,fontsize=8)
        fig.suptitle(f'{analysis}: genomic PCA tracks (inversion intervals shaded)');fig.tight_layout();fig.savefig(FIG/f'{analysis}_pc_tracks.png',dpi=220);fig.savefig(FIG/f'{analysis}_pc_tracks.pdf');plt.close(fig)

def pca_spatial_effects(scores):
    out=[]
    for analysis in ('absolute_time','relative_time'):
        for lg in LGS:
            rows=sorted((r for r in scores if r['analysis']==analysis and r['lg']==lg),key=lambda r:int(r['start']))
            mask=np.asarray([r['region_class']=='inside' for r in rows])
            for pc in ('PC1','PC2','root_age'):
                vals=np.asarray([r[pc] for r in rows],float);obs,p,lo,hi=circular_null(vals,mask)
                out.append({'analysis':analysis,'lg':lg,'PC':pc,'n_inside':int(mask.sum()),'n_outside':int((~mask).sum()),'inside_minus_outside_mean':obs,'circular_shift_p':p,'null_q025':lo,'null_q975':hi,'spatial_null':'all cyclic shifts of frozen inversion mask over ordered available windows'})
    for analysis in sorted(set(r['analysis'] for r in out)):
        for pc in sorted(set(r['PC'] for r in out if r['analysis']==analysis)):
            ix=[i for i,r in enumerate(out) if r['analysis']==analysis and r['PC']==pc and math.isfinite(float(r['circular_shift_p']))]
            qs=bh_adjust([float(out[i]['circular_shift_p']) for i in ix])
            for j,i in enumerate(ix):out[i]['BH_q']=qs[j]
    write_tsv(OUT/'pca_inside_outside_spatial_tests.tsv',out)

def cluster_analyses(valid,relative):
    cluster_rows=[];sil_rows=[];composition=[];enrich=[];stability=[];chosen={}
    ys=np.asarray(relative,float);all_lg=[r['lg'] for r in valid]
    idx_by_lg={lg:np.asarray([i for i,r in enumerate(valid) if r['lg']==lg]) for lg in LGS}
    for metric in ('euclidean','cityblock'):
        D=squareform(pdist(ys,metric=metric));Z=linkage(squareform(D,checks=False),method='average')
        candidates=[]
        for k in range(2,8):
            lab=fcluster(Z,k,criterion='maxclust')
            if len(set(lab))<2:continue
            sil=silhouette_score(ys,lab,metric=metric,sample_size=min(425,len(ys)),random_state=RNG_SEED)
            candidates.append((sil,k,lab));sil_rows.append({'metric':metric,'k':k,'silhouette':sil})
        best=max(candidates,key=lambda x:x[0]);sil,k,labels=best;chosen[metric]=(k,labels,sil)
        for i,r in enumerate(valid):cluster_rows.append({'metric':metric,'chosen_k':k,'cluster':int(labels[i]),'lg':r['lg'],'window_id':r['window_id'],'start':r['start'],'end':r['end'],'midpoint':r['midpoint'],'region_class':r['region_class']})
        for c in sorted(set(labels)):
            for lg in LGS:
                ix=[i for i,r in enumerate(valid) if labels[i]==c and r['lg']==lg]
                for reg in ('inside','outside'):
                    n=sum(valid[i]['region_class']==reg for i in ix);den=sum(r['region_class']==reg and r['lg']==lg for r in valid)
                    composition.append({'metric':metric,'chosen_k':k,'cluster':int(c),'lg':lg,'region_class':reg,'n_cluster_region':n,'n_region':den,'proportion':n/den if den else math.nan})
            mask=np.asarray([r['region_class']=='inside' for r in valid]);vals=(labels==c).astype(float)
            observed,p,lo,hi=stratified_circular_null(vals,mask,all_lg,seed=RNG_SEED+int(c)+(0 if metric=='euclidean' else 100))
            enrich.append({'metric':metric,'chosen_k':k,'cluster':int(c),'n_inside':int(mask.sum()),'n_outside':int((~mask).sum()),'inside_minus_outside_cluster_proportion':observed,'circular_shift_p':p,'null_q025':lo,'null_q975':hi,'n_LGs_equal_weighted':len(LGS),'spatial_null':'independent cyclic interval-mask shift within each LG; equal-LG mean effect'})
        # 80%-coverage contiguous block subsampling for cluster-label stability.
        rng=random.Random(RNG_SEED+(0 if metric=='euclidean' else 1));ari=[]
        for rep in range(80):
            selected=set()
            for lg,ix in idx_by_lg.items():
                target=max(1,int(.8*len(ix)));picked=set()
                while len(picked)<target:
                    start=rng.randrange(len(ix))
                    for j in range(8):picked.add(int(ix[(start+j)%len(ix)]))
                    if len(picked)>=len(ix):break
                selected.update(picked)
            sel=np.asarray(sorted(selected));d=pdist(ys[sel],metric=metric);zz=linkage(d,method='average');lab=fcluster(zz,k,criterion='maxclust')
            ari.append(adjusted_rand_score(labels[sel],lab))
        stability.append({'metric':metric,'chosen_k':k,'resamples':len(ari),'block_size_windows':8,'subsample_fraction_target':.8,'median_ARI_vs_full':statistics.median(ari),'q025_ARI':float(np.quantile(ari,.025)),'q975_ARI':float(np.quantile(ari,.975))})
    if enrich: q=bh_adjust([x['circular_shift_p'] for x in enrich]);[x.update({'BH_q':q[i]}) for i,x in enumerate(enrich)]
    write_tsv(OUT/'cluster_assignments.tsv',cluster_rows);write_tsv(OUT/'cluster_silhouette.tsv',sil_rows);write_tsv(OUT/'cluster_composition.tsv',composition);write_tsv(OUT/'cluster_spatial_enrichment.tsv',enrich);write_tsv(OUT/'cluster_stability.tsv',stability)
    cluster_figures(valid,cluster_rows,composition,chosen);cluster_distance_heatmaps(valid,ys,chosen)
    representative_trees(valid,chosen['euclidean'][1])
    return chosen,sil_rows,enrich,stability

def cluster_distance_heatmaps(valid,features,chosen):
    fig,axs=plt.subplots(1,2,figsize=(12,5))
    for ax,metric in zip(axs,('euclidean','cityblock')):
        D=squareform(pdist(features,metric=metric));Z=linkage(squareform(D,checks=False),method='average')
        order=__import__('scipy.cluster.hierarchy',fromlist=['leaves_list']).leaves_list(Z)
        ordered=D[np.ix_(order,order)];labels=chosen[metric][1][order]
        im=ax.imshow(ordered,origin='lower',interpolation='nearest',aspect='equal',cmap='magma')
        changes=np.where(labels[1:]!=labels[:-1])[0]+1
        for b in changes:ax.axhline(b-.5,color='cyan',lw=.45,alpha=.7);ax.axvline(b-.5,color='cyan',lw=.45,alpha=.7)
        ax.set_title(f'{metric} MRCA-vector distances | k={chosen[metric][0]}');ax.set_xlabel('Windows ordered by average-linkage clustering');ax.set_ylabel('Windows ordered by clustering')
        fig.colorbar(im,ax=ax,fraction=.046,pad=.04,label='Pairwise distance')
    fig.suptitle('Hierarchical clustering distance heatmaps (cyan lines mark selected clusters)');fig.tight_layout();fig.savefig(FIG/'cluster_distance_heatmaps.png',dpi=220);fig.savefig(FIG/'cluster_distance_heatmaps.pdf');plt.close(fig)

def cluster_figures(valid,clusters,composition,chosen):
    fig,axs=plt.subplots(2,4,figsize=(16,7),sharex=False)
    cmap=plt.get_cmap('tab10')
    for i,metric in enumerate(('euclidean','cityblock')):
        rows=[r for r in clusters if r['metric']==metric]
        for j,lg in enumerate(LGS):
            ax=axs[i,j];q=sorted((r for r in rows if r['lg']==lg),key=lambda r:int(r['start']))
            inside=[r for r in q if r['region_class']=='inside']
            if inside:ax.axvspan(min(int(r['start']) for r in inside)/1e6,max(int(r['end']) for r in inside)/1e6,color='#d8a12a',alpha=.16)
            ax.scatter([float(r['midpoint'])/1e6 for r in q],[int(r['cluster']) for r in q],c=[cmap((int(r['cluster'])-1)%10) for r in q],marker='|',s=100)
            ax.set_title(f'{metric} | {lg}');ax.set_ylabel('Cluster');ax.set_xlabel('Position (Mb)');ax.grid(alpha=.18)
    
    fig.suptitle('Hierarchical clusters of root-age-normalized MRCA vectors');fig.tight_layout();fig.savefig(FIG/'cluster_assignment_tracks.png',dpi=220);fig.savefig(FIG/'cluster_assignment_tracks.pdf');plt.close(fig)
    # Cluster composition inside/outside grouped by LG.
    # One compact composition figure per distance metric, faceted by LG.
    for metric in ('euclidean','cityblock'):
        fig,axs=plt.subplots(2,2,figsize=(10,7))
        for ax,lg in zip(axs.flat,LGS):
            rows=[r for r in composition if r['metric']==metric and r['lg']==lg]
            cs=sorted(set(r['cluster'] for r in rows));x=np.arange(len(cs));width=.38
            for j,reg in enumerate(('inside','outside')):
                v=[next((r['proportion'] for r in rows if r['cluster']==c and r['region_class']==reg),0) for c in cs]
                ax.bar(x+(j-.5)*width,v,width,label=reg,color=COLORS[reg])
            ax.set_title(lg);ax.set_xticks(x,[str(c) for c in cs]);ax.set_xlabel('Cluster');ax.set_ylabel('Proportion of region windows');ax.grid(axis='y',alpha=.2)
        axs[0,0].legend(frameon=False);fig.suptitle(f'{metric} clusters: inside/outside composition');fig.tight_layout();fig.savefig(FIG/f'cluster_composition_{metric}.png',dpi=220);fig.savefig(FIG/f'cluster_composition_{metric}.pdf');plt.close(fig)

def representative_trees(valid,labels):
    fig,axs=plt.subplots(3,2,figsize=(15,16));cids=sorted(set(labels))
    for ax,c in zip(axs.flat,cids):
        idx=np.where(labels==c)[0];mat=np.vstack([valid[i]['feature_vector_rel'] for i in idx]);center=np.median(mat,axis=0);med=idx[np.argmin(np.linalg.norm(mat-center,axis=1))];r=valid[med]
        Phylo.draw(r['tree'],axes=ax,do_show=False,show_confidence=False,branch_labels=None);ax.set_title(f'Euclidean cluster {c}; medoid {r["window_id"]}');ax.set_xlabel('Branch time (Ma)')
    for ax in list(axs.flat)[len(cids):]:ax.axis('off')
    fig.suptitle('Observed SNAPP MCC medoid trees by Euclidean hierarchical cluster');fig.tight_layout();fig.savefig(FIG/'cluster_representative_observed_trees.png',dpi=220);fig.savefig(FIG/'cluster_representative_observed_trees.pdf');plt.close(fig)

def published_topology_map(taxa):
    tree=Phylo.read(str(BASELINE),'newick');present=set(t.name for t in tree.get_terminals());
    if set(taxa)-present:raise ValueError(f'published baseline lacks local taxa: {sorted(set(taxa)-present)}')
    tips={t.name:t for t in tree.get_terminals()};out={}
    for a,b,c in itertools.combinations(taxa,3):
        trip=set((a,b,c));possible=[]
        for x,y,z in ((a,b,c),(a,c,b),(b,c,a)):
            anc=tree.common_ancestor(tips[x],tips[y]);desc={t.name for t in anc.get_terminals()}&trip
            if desc=={x,y}:possible.append(f'{x}+{y}|{z}')
        out[(a,b,c)]=possible[0] if len(possible)==1 else 'unresolved'
    return out

def triplet_analyses(valid,taxa):
    baseline=published_topology_map(taxa);long=[];summaries=[];per_trip_lg={}
    for r in valid:
        feats=triplet_features(taxa,r['pair_ages'],r['root_age'])
        for f in feats:
            a,b,c=f['taxa'];base=baseline[(a,b,c)]; agrees=f['topology']==base and f['topology']!='unresolved'
            rec={'lg':r['lg'],'window_id':r['window_id'],'start':r['start'],'end':r['end'],'midpoint':r['midpoint'],'region_class':r['region_class'],'triplet':f['triplet'],'topology':f['topology'],'baseline_topology':base,'agrees_with_baseline':agrees,'temporal_separation':f['temporal_separation'],'normalized_temporal_separation':f['normalized_temporal_separation'],'unresolved':f['unresolved']}
            for name,value in f['q_contrasts'].items():rec['Q_'+name]=value
            for name,value in f['separation_by_resolution'].items():rec['S_'+name]=value
            for name,value in f['normalized_separation_by_resolution'].items():rec['Snorm_'+name]=value
            for name,value in f['normalized_q_contrasts'].items():rec['Qnorm_'+name]=value
            for name,value in f['resolution_topologies'].items():rec['topology_'+name]=value
            long.append(rec);per_trip_lg.setdefault((r['lg'],f['triplet']),[]).append(rec)
    categories=defaultdict(list)
    for lg in LGS:
        for trip in itertools.combinations(taxa,3):
            tname=','.join(trip);rows=per_trip_lg.get((lg,tname),[])
            if not rows:continue
            ins=[r for r in rows if r['region_class']=='inside'];outs=[r for r in rows if r['region_class']=='outside']
            if not ins or not outs:continue
            row={'lg':lg,'triplet':tname,'baseline_topology':rows[0]['baseline_topology'],'n_inside':len(ins),'n_outside':len(outs),'n_boundary':sum(r['region_class']=='boundary' for r in rows)}
            for resolution,topology in enumerate([f'{trip[0]}+{trip[1]}|{trip[2]}',f'{trip[0]}+{trip[2]}|{trip[1]}',f'{trip[1]}+{trip[2]}|{trip[0]}'],1):
                key=f'r{resolution}'
                for side,group in [('inside',ins),('outside',outs)]:row[f'freq_{key}_{side}']=sum(r['topology']==topology for r in group)/len(group)
                delta=row[f'freq_{key}_inside']-row[f'freq_{key}_outside']
                vals=np.asarray([float(r['topology']==topology) for r in rows]);mask=np.asarray([r['region_class']=='inside' for r in rows]);_,p,lo,hi=circular_null(vals,mask)
                row[f'delta_freq_{key}']=delta;row[f'p_freq_{key}']=p;row[f'null_q025_freq_{key}']=lo;row[f'null_q975_freq_{key}']=hi
                sepkey=f'Snorm_{key}'
                for side,group in [('inside',ins),('outside',outs)]:
                    all_sep=np.asarray([float(r[sepkey]) for r in group]);supported=np.asarray([float(r[sepkey]) for r in group if r['topology']==topology])
                    row[f'mean_{sepkey}_{side}']=float(np.mean(all_sep))
                    row[f'median_{sepkey}_{side}']=float(np.median(all_sep))
                    row[f'mean_supported_{sepkey}_{side}']=float(np.mean(supported)) if len(supported) else math.nan
                row[f'delta_{sepkey}_inside_minus_outside']=row[f'mean_{sepkey}_inside']-row[f'mean_{sepkey}_outside']
                v=np.asarray([float(r[sepkey]) for r in rows]);obs,p,lo,hi=circular_null(v,mask)
                row[f'p_{sepkey}_circular']=p;row[f'null_q025_{sepkey}']=lo;row[f'null_q975_{sepkey}']=hi
            for side,group in [('inside',ins),('outside',outs)]:
                row[f'freq_unresolved_{side}']=sum(r['unresolved'] for r in group)/len(group)
                resolved=[r for r in group if not r['unresolved']]
                row[f'mean_sep_Ma_{side}']=float(np.mean([r['temporal_separation'] for r in resolved])) if resolved else math.nan
                row[f'mean_sep_norm_{side}']=float(np.mean([r['normalized_temporal_separation'] for r in resolved])) if resolved else math.nan
                row[f'median_sep_norm_{side}']=float(np.median([r['normalized_temporal_separation'] for r in resolved])) if resolved else math.nan
                row[f'agreement_baseline_{side}']=sum(r['agrees_with_baseline'] for r in group)/len(group)
            row['delta_unresolved_inside_minus_outside']=row['freq_unresolved_inside']-row['freq_unresolved_outside']
            row['delta_mean_sep_Ma_inside_minus_outside']=row['mean_sep_Ma_inside']-row['mean_sep_Ma_outside']
            row['delta_mean_sep_norm_inside_minus_outside']=row['mean_sep_norm_inside']-row['mean_sep_norm_outside']
            row['delta_agreement_inside_minus_outside']=row['agreement_baseline_inside']-row['agreement_baseline_outside']
            resolved_vals=np.asarray([r['normalized_temporal_separation'] if not r['unresolved'] else np.nan for r in rows],float);mask=np.asarray([r['region_class']=='inside' for r in rows])
            obs,p,lo,hi=circular_null(resolved_vals,mask);row.update({'p_sep_norm_circular':p,'null_q025_sep_norm':lo,'null_q975_sep_norm':hi})
            agree_vals=np.asarray([float(r['agrees_with_baseline']) for r in rows]);obs,p,lo,hi=circular_null(agree_vals,mask);row.update({'p_baseline_agreement_circular':p,'null_q025_baseline_agreement':lo,'null_q975_baseline_agreement':hi})
            summaries.append(row);categories[tname].append(row)
    # FDR within each metric family across all LG x triplet comparisons.
    for prefix in ('p_sep_norm_circular','p_baseline_agreement_circular'):
        ix=[i for i,r in enumerate(summaries) if math.isfinite(float(r[prefix]))]
        qs=bh_adjust([float(summaries[i][prefix]) for i in ix])
        for j,i in enumerate(ix):summaries[i][prefix.replace('_circular','_BH_q')]=qs[j]
    s_families=sorted({k for r in summaries for k in r if k.startswith('p_Snorm_') and k.endswith('_circular')})
    for pkey in s_families:
        ix=[i for i,r in enumerate(summaries) if math.isfinite(float(r.get(pkey,math.nan)))];qs=bh_adjust([float(summaries[i][pkey]) for i in ix]);qkey=pkey.replace('p_','BH_q_',1).replace('_circular','')
        for j,i in enumerate(ix):summaries[i][qkey]=qs[j]
    for topology_key in sorted({k for r in summaries for k in r if k.startswith('p_freq_')}):
        ix=[i for i,r in enumerate(summaries) if math.isfinite(float(r.get(topology_key,math.nan)))];qs=bh_adjust([float(summaries[i][topology_key]) for i in ix])
        qkey=topology_key.replace('p_','BH_q_',1)
        for j,i in enumerate(ix):summaries[i][qkey]=qs[j]
    write_tsv(OUT/'rooted_triplets_per_window.tsv',long);write_tsv(OUT/'rooted_triplet_inside_outside.tsv',summaries)
    # Criteria fixed in advance: rank complete triplet set by max LG-specific
    # normalized separation effect, topology-frequency shift, and baseline mismatch.
    rankings=[];triples=sorted(categories)
    for t in triples:
        rs=categories[t]
        sep=max(abs(r['delta_mean_sep_norm_inside_minus_outside']) for r in rs)
        topo=max(abs(r[k]) for r in rs for k in r if k.startswith('delta_freq_') and not k.startswith('delta_freq_unresolved'))
        mismatch=max(r['agreement_baseline_outside']-r['agreement_baseline_inside'] for r in rs)
        rankings.append({'triplet':t,'max_abs_LG_delta_normalized_separation':sep,'max_abs_LG_topology_frequency_shift':topo,'max_LG_added_baseline_disagreement_inside':mismatch})
    for metric in ('max_abs_LG_delta_normalized_separation','max_abs_LG_topology_frequency_shift','max_LG_added_baseline_disagreement_inside'):
        for rank,r in enumerate(sorted(rankings,key=lambda x:x[metric],reverse=True),1):r['rank_'+metric]=rank
    write_tsv(OUT/'rooted_triplet_rankings.tsv',rankings)
    selected=[]
    for metric in ('max_abs_LG_delta_normalized_separation','max_abs_LG_topology_frequency_shift','max_LG_added_baseline_disagreement_inside'):
        selected.extend(r['triplet'] for r in sorted(rankings,key=lambda x:x[metric],reverse=True)[:2])
    selected=list(dict.fromkeys(selected));write_tsv(OUT/'selected_triplets_for_figures.tsv',[{'rank':i+1,'triplet':t,'selection_criteria':'top two per prespecified ranking: normalized separation shift, topology-frequency shift, baseline disagreement'} for i,t in enumerate(selected)])
    triplet_figures(long,summaries,selected,taxa)
    return long,summaries,rankings,selected

def triplet_figures(long,summaries,selected,taxa):
    # Signed time contrasts, with each linkage group in its own panel.
    fig,axs=plt.subplots(len(selected),4,figsize=(16,max(2.0*len(selected),9)),squeeze=False,sharey=True)
    qcolors=('#0072B2','#D55E00','#009E73')
    short=lambda x:x.replace('Gadmor_','').replace('_spc','')
    for i,t in enumerate(selected):
        trip=tuple(t.split(','));topos=[f'{trip[0]}+{trip[1]}|{trip[2]}',f'{trip[0]}+{trip[2]}|{trip[1]}',f'{trip[1]}+{trip[2]}|{trip[0]}']
        for j,lg in enumerate(LGS):
            ax=axs[i,j];rr=sorted((r for r in long if r['triplet']==t and r['lg']==lg),key=lambda r:int(r['start']))
            inside=[r for r in rr if r['region_class']=='inside']
            if inside:ax.axvspan(min(int(r['start']) for r in inside)/1e6,max(int(r['end']) for r in inside)/1e6,color='#d8a12a',alpha=.15)
            for q,top in enumerate(topos):
                ax.plot([float(r['midpoint'])/1e6 for r in rr],[float(r[f'Qnorm_r{q+1}']) for r in rr],'.-',ms=1.8,lw=.65,color=qcolors[q],alpha=.9,label=f'Q{q+1}' if i==0 and j==0 else None)
            ax.axhline(0,color='black',lw=.6,alpha=.55);ax.grid(alpha=.15);ax.set_ylim(-1.05,1.05)
            if i==0:ax.set_title(lg)
            if j==0:ax.set_ylabel(f'T{i+1}: {short(trip[0])}+{short(trip[1])}|{short(trip[2])}',fontsize=8)
            if i==len(selected)-1:ax.set_xlabel('Position (Mb)')
    fig.suptitle('Signed root-age-normalized Q contrasts for objectively selected triplets',y=.995)
    fig.legend(loc='upper center',bbox_to_anchor=(.55,.972),ncol=3,frameon=False)
    fig.tight_layout(rect=(.01,.01,.99,.94));fig.savefig(FIG/'rooted_triplet_selected_tracks.png',dpi=220);fig.savefig(FIG/'rooted_triplet_selected_tracks.pdf');plt.close(fig)
    # Show all three rooted-resolution frequencies plus unresolved observations.
    fig,axs=plt.subplots(len(selected),4,figsize=(17,max(2.4*len(selected),8)),squeeze=False)
    cmap=plt.get_cmap('Set2')
    from matplotlib.patches import Patch
    legend_handles=[Patch(facecolor=cmap(i),label=f'r{i+1}') for i in range(3)]+[Patch(facecolor=cmap(3),edgecolor='white',label='unresolved'),Patch(facecolor='#999',label='inside'),Patch(facecolor='#999',hatch='//',label='outside')]
    for i,t in enumerate(selected):
        triple=tuple(t.split(','));topos=[f'{triple[0]}+{triple[1]}|{triple[2]}',f'{triple[0]}+{triple[2]}|{triple[1]}',f'{triple[1]}+{triple[2]}|{triple[0]}','unresolved']
        for j,lg in enumerate(LGS):
            row=next(r for r in summaries if r['lg']==lg and r['triplet']==t);ax=axs[i,j];x=np.arange(4);w=.38
            for h,side in enumerate(('inside','outside')):
                if side=='inside':vals=[row[f'freq_r{z}_inside'] for z in range(1,4)]+[row['freq_unresolved_inside']]
                else:vals=[row[f'freq_r{z}_outside'] for z in range(1,4)]+[row['freq_unresolved_outside']]
                ax.bar(x+(h-.5)*w,vals,w,color=[cmap(z) for z in range(4)],alpha=.92 if side=='inside' else .55,hatch='' if side=='inside' else '//',edgecolor='white',label=side if i==0 and j==0 else None)
            ax.set_xticks(x,['r1','r2','r3','unres'],rotation=30,ha='right');ax.set_ylim(0,1);ax.grid(axis='y',alpha=.17)
            if j==0:
                short=lambda x:x.replace('Gadmor_','').replace('_spc','')
                ax.set_ylabel(f"{short(triple[0])}+{short(triple[1])},{short(triple[0])}+{short(triple[2])},{short(triple[1])}+{short(triple[2])}",fontsize=6.5)
            if i==0:ax.set_title(lg)
    fig.suptitle('Inside/outside rooted-triplet topology frequencies; row labels list r1–r3 sister pairs in order',y=.995)
    fig.legend(handles=legend_handles,loc='upper center',bbox_to_anchor=(.63,.975),ncol=6,frameon=False,fontsize=8)
    fig.tight_layout(rect=(.01,.01,.99,.94));fig.savefig(FIG/'selected_triplet_topology_frequencies.png',dpi=220);fig.savefig(FIG/'selected_triplet_topology_frequencies.pdf');plt.close(fig)

def report(trees,valid,taxa,pca_summary,cluster_choice,cluster_sil,cluster_enrich,cluster_stability,trip_summaries,selected):
    invalid=[r for r in trees if not r['valid']]
    lines=['# Atlantic cod SNAPP time-tree structure analysis','',
      'This analysis summarizes the existing published 250-kb SNAPP MCC population time trees. It does not rerun SNAPP or infer a new joint genome-wide tree. A window is the spatial observation; neighboring windows are dependent. No window-matched posterior samples were available, so uncertainty below is spatial-null uncertainty, not posterior uncertainty.','',
      '## Input and tree audit','',f'- Source rows: {len(trees)} across LG01, LG02, LG07, LG12; valid rooted ultrametric trees: {len(valid)}; invalid/excluded: {len(invalid)}.',
      f'- All trees use the same {len(taxa)} cod-population tips; no outgroup tips occur in the local MCC trees. Thus outgroup features do not drive this PCA, and cod-population-only is the primary analysis.',
      '- Branch lengths are the published SNAPP time units (Ma); root stem, if present, is excluded from root age. Zero non-root branches are retained. Missing/negative lengths, non-ultrametric trees, or zero root age are not repaired.',
      '- Per-window posterior files were not available in the source inventory. Supergene-wide posterior files are not matched to individual windows and are not treated as independent loci.',
      '- Invalid tree rows: '+(', '.join(f"{r['lg']}:{r['window_id']} ({r['validation_issue']})" for r in invalid) if invalid else 'none'),'',
      '## Time-tree PCA','']
    for r in pca_summary:
        lines.append(f"- {r['analysis']}: PC1 explains {float(r['PC1_variance']):.1%}, PC2 {float(r['PC2_variance']):.1%}; PC1/root-age correlation r={float(r['PC1_root_age_pearson_r']):.3f}.")
    for analysis in ('absolute_time','relative_time'):
        top=[r for r in read_tsv(OUT/'pca_top_population_pairs.tsv') if r['analysis']==analysis and r['PC']=='1'][:5]
        lines.append(f"- Largest PC1 pair loadings in {analysis}: "+', '.join(r['feature'].replace('__','–') for r in top)+'.')
    for lg in LGS:
        q=[r for r in read_tsv(OUT/'pca_inside_outside_spatial_tests.tsv') if r['analysis']=='absolute_time' and r['lg']==lg and r['PC']=='root_age']
        if q:lines.append(f"- {lg} mean root age inside minus outside = {float(q[0]['inside_minus_outside_mean']):.4g} Ma (cyclic-shift p={float(q[0]['circular_shift_p']):.4g}, BH q={float(q[0]['BH_q']):.4g}).")
    lines += ['', 'Per-linkage-group PC1/PC2 contrasts and cyclic-shift nulls are in `pca_inside_outside_spatial_tests.tsv`. The cyclic-shift test rotates the frozen interval mask over the ordered available windows and leaves each observed score track intact. Relative-time PCA removes overall root-age scaling; the separate feature-z-score PCA is only a scaling sensitivity.', '', '## Hierarchical clustering','']
    for metric,(k,labels,sil) in cluster_choice.items():
        st=next(x for x in cluster_stability if x['metric']==metric)
        lines.append(f'- {metric}: best silhouette among k=2–7 is k={k} (silhouette {sil:.3f}); 80% contiguous-block subsample median ARI={st["median_ARI_vs_full"]:.3f} (2.5–97.5% {st["q025_ARI"]:.3f}–{st["q975_ARI"]:.3f}).')
    lines += ['', 'Cluster count/enrichment tables use independent cyclic shifts of the interval mask within each LG, with equal-LG weighting. Do not read cluster membership as a discrete biological history; inspect the representative observed medoid trees and the MRCA features. Cluster enrichment is descriptive and window-local, with only four chromosomes.', '', 'Spatial cluster enrichment:']
    for e in sorted(cluster_enrich,key=lambda x:float(x['BH_q']))[:8]:
        lines.append(f"- {e['metric']} cluster {e['cluster']}: inside-minus-outside proportion {float(e['inside_minus_outside_cluster_proportion']):+.3f}; cyclic-shift p={float(e['circular_shift_p']):.3g}, BH q={float(e['BH_q']):.3g}.")
    lines += ['', '## Rooted-triplet analysis','',f'- All {math.comb(len(taxa),3)} cod-population triplets were enumerated. Q contrasts and temporal separations are included for all three possible rooted resolutions; triplets with zero/near-zero internal separation are retained as unresolved.',
      '- Local topology frequencies are frequencies among MCC window point trees, not posterior support for a genome-wide species tree. Temporal separation is reported separately from topology frequency, in Ma and divided by root age.',
      '- Published-tree agreement is computed from the induced rooted topology in `baseline_population_tree_12taxa.nwk`; it compares topology only, not times. The source baseline includes an artificial root-stem value that is irrelevant to the induced topology.',
      '- For all triplets and linkage groups, the summary includes inside/outside topology frequencies, unresolved frequencies, mean/median temporal separation, baseline agreement, effects, and circular-shift null p-values (with BH q-values by metric family).',
      '- Track examples were selected objectively: top two triplets each for maximum LG-specific normalized separation shift, topology-frequency shift, and added baseline disagreement inside. Selected triplets: '+', '.join(selected)+'.', '',
      'Per-linkage-group triplet summaries:']
    for lg in LGS:
        q=[r for r in trip_summaries if r['lg']==lg]
        def modal(side,r):
            vals={}
            for k,v in r.items():
                if k.startswith('freq_') and k.endswith('_'+side) and 'unresolved' not in k:vals[k[len('freq_'):-len('_'+side)]]=float(v)
            vals['unresolved']=float(r[f'freq_unresolved_{side}'])
            return max(vals,key=vals.get)
        flips=sum(modal('inside',r)!=modal('outside',r) for r in q)
        top=max(q,key=lambda r:abs(float(r['delta_mean_sep_norm_inside_minus_outside'])))
        agree_i=statistics.mean(float(r['agreement_baseline_inside']) for r in q);agree_o=statistics.mean(float(r['agreement_baseline_outside']) for r in q)
        mean_delta=statistics.mean(float(r['delta_mean_sep_norm_inside_minus_outside']) for r in q)
        lines.append(f"- {lg}: inside/outside modal rooted topology differs for {flips}/{len(q)} triplets; average normalized separation shift across triplets is {mean_delta:+.3f}; agreement with the published baseline is {agree_i:.1%} inside vs {agree_o:.1%} outside. The largest normalized mean-separation shift is {float(top['delta_mean_sep_norm_inside_minus_outside']):+.3f} for {top['triplet']} (p={float(top['p_sep_norm_circular']):.3g}, BH q={float(top.get('p_sep_norm_BH_q',math.nan)):.3g}).")
    lines += ['', '## Interpretation','',
      'These analyses can establish that the observed local MCC time trees differ in absolute divergence ages, relative MRCA structure, triplet ordering frequencies, or temporal separation across the frozen inversion intervals. They cannot establish that inversions caused those changes, that the MCC trees represent posterior uncertainty, or that a local effect biases a genome-wide divergence-time estimate. Comparisons are descriptive, spatially adjusted by cyclic shifts, and not independent replication across 426 windows.',
      '', '## Reproduction','', '`python3 empirical/atlantic_cod/analysis/time_tree_structure/scripts/run_time_tree_structure.py`',
      '', 'Validation tests: `python3 -m unittest discover -s empirical/atlantic_cod/analysis/time_tree_structure/tests -v`.',
      '']
    (OUT/'interpretation_report.md').write_text('\n'.join(lines))

def main():
    OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
    trees,taxa,audits=local_tree_table();write_tsv(OUT/'tree_availability_validation.tsv',audits)
    valid,pairs,names,absolute,relative,roots=prepare_vectors(trees,taxa)
    vector_rows=[]
    for i,r in enumerate(valid):
        v={'lg':r['lg'],'window_id':r['window_id'],'start':r['start'],'end':r['end'],'midpoint':r['midpoint'],'region_class':r['region_class'],'root_age':roots[i]}
        v.update({names[j]:absolute[i,j] for j in range(len(names))});vector_rows.append(v)
    write_tsv(OUT/'window_mrca_age_vectors_absolute.tsv',vector_rows)
    vector_rel=[]
    for i,r in enumerate(valid):
        v={'lg':r['lg'],'window_id':r['window_id'],'start':r['start'],'end':r['end'],'midpoint':r['midpoint'],'region_class':r['region_class'],'root_age':roots[i]}
        v.update({names[j]:relative[i,j] for j in range(len(names))});vector_rel.append(v)
    write_tsv(OUT/'window_mrca_age_vectors_relative.tsv',vector_rel)
    pcas,scores,psummary=fit_pcas(valid,names,absolute,relative,roots);pca_figures(pcas,scores,valid);pca_spatial_effects(scores)
    cluster_choice,sils,enrich,stability=cluster_analyses(valid,relative)
    long,summaries,rankings,selected=triplet_analyses(valid,taxa)
    report(trees,valid,taxa,psummary,cluster_choice,sils,enrich,stability,summaries,selected)
    print(f'{len(trees)} trees; {len(valid)} valid; {len(taxa)} populations; {math.comb(len(taxa),3)} triplets; selected {len(selected)} examples')

if __name__=='__main__':main()
