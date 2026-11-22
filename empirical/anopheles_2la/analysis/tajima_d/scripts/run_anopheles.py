#!/usr/bin/env python3
"""Ag3.10 all-site Tajima D tracks for population/background groups on 2L."""
from pathlib import Path
from collections import defaultdict
import sys,csv,json,time,argparse
import numpy as np,pandas as pd,matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[5]
sys.path.insert(0,str(ROOT/'empirical/analysis/tajimas_d'))
from tajimas import tajima_from_sufficient_stats
BASE=ROOT/'empirical/anopheles_2la/analysis/tajima_d'; DATA=ROOT/'data/anopheles_2la'; JOINT=ROOT/'empirical/anopheles_2la/analysis/joint_genomic_tracks'
MAN=DATA/'metadata/sample_manifest.tsv'; GRID=DATA/'processed/stage2r_window_grid.tsv'; INV=tuple(int(x) for x in pd.read_csv(DATA/'processed/stage0r_2la_interval.tsv',sep='\t').iloc[0][['start','end']])
BUCKET='gs://vo_agam_release_master_us_central1'; SAMPLE_SET='fontaine-2015-rebuild'; CHUNK=2_000_000

def make_groups(rows):
    groups={}
    species=['gambiae','coluzzii']
    for sp in species:
        sr=[r for r in rows if r['species']==sp]
        groups['species_pooled_'+sp]=sr
        for state,label in [('A0_homozygous','standard'),('A1_homozygous','inverted')]:
            rr=[r for r in sr if r['2La_state']==state]
            groups['species_'+label+'_'+sp]=rr
        for loc in sorted(set(r['locality'] for r in sr)):
            lr=[r for r in sr if r['locality']==loc]
            groups['population_pooled_'+sp+'_'+loc]=lr
            for state,label in [('A0_homozygous','standard'),('A1_homozygous','inverted')]:
                ar=[r for r in lr if r['2La_state']==state]
                if len(ar)>=2: groups['population_'+label+'_'+sp+'_'+loc]=ar
    return {k:v for k,v in groups.items() if len(v)>=2}

def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--summary-only',action='store_true'); args=parser.parse_args()
    rows=list(csv.DictReader(MAN.open(),delimiter='\t')); rows=[r for r in rows if r['species'] in ('gambiae','coluzzii')]
    groups=make_groups(rows); sample_ids=list(dict.fromkeys(r['sample_id'] for r in rows))
    if not args.summary_only:
        import malariagen_data
        ag=malariagen_data.Ag3(url=BUCKET,check_location=False,show_progress=False); meta=ag.sample_metadata(sample_sets=SAMPLE_SET)
        mdids=[str(x) for x in meta.sample_id]; meta_ix={s:i for i,s in enumerate(mdids)}
        missing=set(sample_ids)-set(meta_ix)
        if missing: raise ValueError('manifest IDs absent from Ag3 metadata: '+repr(sorted(missing)))
        indices=[meta_ix[s] for s in sample_ids]
    grid=pd.read_csv(GRID,sep='\t'); maxend=int(grid.end.max()); rows50=grid.to_dict('records')
    win_defs={50_000:rows50,100_000:[]}
    for lo in range(1,maxend+1,100_000):
        hi=min(maxend,lo+99_999);win_defs[100_000].append({'window_id':f'2L_{lo}_{hi}','chrom':'2L','start':lo,'end':hi,'region_class':'inside' if (lo+hi)//2>=INV[0] and (lo+hi)//2<=INV[1] else 'outside'})
    stats=defaultdict(lambda:[0,0.,0,0,0,0,0]) # S, pi, callable, multiallelic, missing, mask fail, n positions
    audit={'source':'MalariaGEN Ag3.10 fontaine-2015-rebuild; AgamP4; Ag3.10 all-site genotype API','sample_count_species':{},'group_sample_counts':{k:len(v) for k,v in groups.items()},'filters':'gamb_colu site mask; complete genotype calls for the fixed group at a site; exclude observed multiallelic sites; genotypes diploid','window_resolutions_bp':[50000,100000],'interval_1based_inclusive':INV,'chunks':[]}
    # Keep stable sample indices; chunking bounds memory. Ag3.10 call arrays include all non-N reference positions.
    chunks=[(lo,min(maxend,lo+CHUNK-1)) for lo in range(1,maxend+1,CHUNK)]
    state_path=BASE/'results/anopheles_tajima_d_sufficient_stats.json'
    if args.summary_only:
        state=json.load(state_path.open())
        stats=defaultdict(lambda:[0,0.,0,0,0,0,0])
        for item in state['stats']:
            stats[(int(item[0]),item[1],item[2])]=item[3]
        audit=json.load((BASE/'results/anopheles_tajima_d_checkpoint.json').open())['audit']
    else:
        for ci,(lo,hi) in enumerate(chunks,1):
            region=f'2L:{lo}-{hi}'; t=time.time(); print(f'query {ci}/{len(chunks)} {region}',flush=True)
            ds=ag.snp_calls(region=region,sample_sets=SAMPLE_SET,sample_indices=indices,site_mask=None,inline_array=True,chunks='native')
            returned=[str(x) for x in ds.sample_id.values.tolist()]
            if returned!=sample_ids: raise RuntimeError('Ag3 sample order mismatch')
            pos=np.asarray(ds.variant_position.values,dtype=np.int64); gt=np.asarray(ds.call_genotype.values,dtype=np.int8)
            passed=np.asarray(ds.variant_filter_pass_gamb_colu.values,dtype=bool)
            official=ag.site_filters(region=region,mask='gamb_colu',field='filter_pass'); official=np.asarray(official.compute() if hasattr(official,'compute') else official,dtype=bool)
            if not np.array_equal(passed,official): raise RuntimeError(f'site mask mismatch in {region}')
            ref=np.asarray(ag.genome_sequence(region).compute()).astype(str); nonn=np.char.upper(ref)!='N'; expect=np.flatnonzero(nonn)+lo
            if not np.array_equal(pos,expect): raise RuntimeError(f'Ag3 records not all and only non-N reference positions in {region}')
            if len(pos)!=len(gt) or len(pos)!=len(passed): raise RuntimeError('variant arrays have unequal lengths')
            win_maps={}
            for res,defs in win_defs.items():
                bids=(pos-1)//res
                win_maps[res]=bids
            for gname,rs in groups.items():
                ix=np.array([sample_ids.index(r['sample_id']) for r in rs],dtype=int); n=2*len(ix)
                sub=gt[:,ix,:]; complete=np.all(sub>=0,axis=(1,2)); good=passed&complete
                allele_counts=np.stack([(sub==a).sum(axis=(1,2)) for a in range(4)],axis=1)
                nstates=(allele_counts>0).sum(axis=1); biallelic=nstates<=2
                include=good&biallelic
                variable=include&(nstates==2)
                pi=np.divide(n*n-(allele_counts*allele_counts).sum(axis=1),n*(n-1),out=np.zeros(len(pos)),where=n>1)
                for res,defs in win_defs.items():
                    bid=win_maps[res]; bins=np.unique(bid)
                    for bi in bins:
                        wstart=int(bi*res+1); wend=min(maxend,int((bi+1)*res)); # align with frozen 50kb IDs where possible
                        # Window ID and class from the authoritative 50-kb grid or the physical 100-kb grid.
                        if res==50_000:
                            match=defs[int(bi)] if int(bi)<len(defs) else None
                            if match is None: continue
                            wid=match['window_id']; wstart=int(match['start']); wend=int(match['end']); cls=match.get('region_class','')
                        else:
                            di=int(bi); match=defs[di] if di<len(defs) else None
                            if match is None: continue
                            wid=match['window_id']; wstart=int(match['start']); wend=int(match['end']); cls=match['region_class']
                        m=(bid==bi); key=(res,gname,wid); acc=stats[key]; acc[6]+=int(m.sum()); acc[5]+=int((m&~passed).sum()); acc[4]+=int((m&passed&~complete).sum()); acc[3]+=int((m&good&~biallelic).sum()); acc[2]+=int((m&include).sum()); acc[0]+=int((m&variable).sum()); acc[1]+=float(pi[m&include].sum())
            audit['chunks'].append({'region':region,'positions':len(pos),'non_N_reference_positions':int(nonn.sum()),'site_mask_pass':int(passed.sum()),'seconds':round(time.time()-t,2)})
            if ci%5==0 or ci==len(chunks):
                checkpoint={'completed_chunks':ci,'total_chunks':len(chunks),'completed_at':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'audit':audit}
                (BASE/'results').mkdir(parents=True,exist_ok=True); (BASE/'results/anopheles_tajima_d_checkpoint.json').write_text(json.dumps(checkpoint,indent=2)+'\n')
            del ds,gt
        state={'stats':[[k[0],k[1],k[2],v] for k,v in stats.items()]}
        state_path.write_text(json.dumps(state,separators=(',',':'))+'\n')
    # Emit all windows/groups, preserving undefined values and reasons.
    out=[]
    for res,defs in win_defs.items():
      for d in defs:
       for gname,rs in groups.items():
        x=stats[(res,gname,d['window_id'])]; S,pi,callable_n,multi,missing_n,mask_fail,total=x; n=2*len(rs)
        result=tajima_from_sufficient_stats(n,int(S),float(pi),int(callable_n),int(callable_n) if callable_n else None)
        mid=(int(d['start'])+int(d['end']))//2
        out.append({'resolution_bp':res,'window_id':d['window_id'],'chrom':'2L','start':d['start'],'end':d['end'],'mid':mid,'interval_class':'inside_2La' if INV[0]<=mid<=INV[1] else 'outside_2La','group':gname,'species':rs[0]['species'],'group_definition':'independent metadata; arrangement groups use homozygous karyotypes; pooled groups include heterokaryotypes','n_individuals':len(rs),'n_haplotypes':n,'ploidy':'diploid','S':int(S),'pi_total':pi,'theta_W':result['theta_W'],'tajimas_D':result['tajimas_D'],'D_undefined_reason':result['D_undefined_reason'],'n_complete_biallelic_callable_sites':int(callable_n),'callable_bases':int(callable_n) if callable_n else None,'n_positions_seen':int(total),'n_mask_fail':int(mask_fail),'n_missing_genotype':int(missing_n),'n_multiallelic_excluded':int(multi),'pi_per_callable_base':result['pi_per_callable_base'],'thetaW_per_callable_base':result['thetaW_per_callable_base'],'callability_definition':'Ag3.10 non-N reference positions passing gamb_colu mask with all fixed-group genotypes called; multiallelic positions excluded from numerator and denominator','ascertainment':'Ag3.10 all-site genotype dataset; per chunk verified one record per non-N reference position; segregating sites complete conditional on source callset and gamb_colu mask'})
    out=pd.DataFrame(out); (BASE/'results').mkdir(parents=True,exist_ok=True)
    for res in (50000,100000): out[out.resolution_bp==res].to_csv(BASE/f'results/anopheles_tajima_d_{res//1000}kb.tsv',sep='\t',index=False,na_rep='NA')
    # Regional estimates with 500-kb spatial block bootstrap, independent within each group.
    summary=[]; rng=np.random.default_rng(20261121)
    for (res,g),df in out.groupby(['resolution_bp','group']):
        for metric in ['tajimas_D','pi_per_callable_base','thetaW_per_callable_base']:
            a=df[df.interval_class=='inside_2La'].dropna(subset=[metric]); b=df[df.interval_class=='outside_2La'].dropna(subset=[metric])
            a=a.copy();b=b.copy();a['block']=(a.mid//500000).astype(int);b['block']=(b.mid//500000).astype(int)
            av=a.groupby('block')[metric].mean().values; bv=b.groupby('block')[metric].mean().values
            boots=[]
            if len(av)>1 and len(bv)>1:
                for _ in range(5000): boots.append(rng.choice(av,len(av),replace=True).mean()-rng.choice(bv,len(bv),replace=True).mean())
            summary.append({'resolution_bp':res,'group':g,'statistic':metric,'inside_mean':a[metric].mean(),'outside_mean':b[metric].mean(),'inside_minus_outside':a[metric].mean()-b[metric].mean(),'n_inside_windows':len(a),'n_outside_windows':len(b),'n_inside_blocks':len(av),'n_outside_blocks':len(bv),'block_bootstrap_95_low':np.quantile(boots,.025) if boots else np.nan,'block_bootstrap_95_high':np.quantile(boots,.975) if boots else np.nan,'interpretation':'descriptive spatial contrast; not neutrality-test p-value; arrangement groups confounded with geography'})
    pd.DataFrame(summary).to_csv(BASE/'results/anopheles_tajima_d_inside_outside.tsv',sep='\t',index=False,na_rep='NA')
    if not args.summary_only or not (BASE/'results/anopheles_tajima_d_joined_tracks.tsv').exists():
        # Exact 50-kb joins to existing population-genetic and genealogy tracks.
        track=out[out.resolution_bp==50000].copy(); pi=pd.read_csv(JOINT/'results/anopheles_2la_pi_windows.tsv',sep='\t'); dx=pd.read_csv(JOINT/'results/anopheles_2la_dxy_fst_windows.tsv',sep='\t'); tree=pd.read_csv(JOINT/'results/anopheles_2la_integrated_genomic_tracks.tsv',sep='\t')
        track=track.merge(pi,on='window_id',how='left',suffixes=('','_pi')).merge(dx,on='window_id',how='left',suffixes=('','_dxy')).merge(tree,on='window_id',how='left',suffixes=('','_tree'))
        track.to_csv(BASE/'results/anopheles_tajima_d_joined_tracks.tsv',sep='\t',index=False,na_rep='NA')
        correlations=[]
        for group,df in track.groupby('group'):
            cand=[c for c in df.columns if any(t in c.lower() for t in ['pi_per_callable','dxy_per_callable','hudson_fst','q1','q2','q3','branch','d_topology']) and c!='group']
            for col in cand:
                x=pd.to_numeric(df.tajimas_D,errors='coerce');y=pd.to_numeric(df[col],errors='coerce');ok=x.notna()&y.notna()
                if ok.sum()>=10: correlations.append({'group':group,'track':col,'n_windows':int(ok.sum()),'spearman_rho':x[ok].corr(y[ok],method='spearman'),'note':'descriptive window correlation; spatial dependence and non-independent overlapping estimands'})
        pd.DataFrame(correlations).to_csv(BASE/'results/anopheles_tajima_d_track_correlations.tsv',sep='\t',index=False,na_rep='NA')
    # Plot core groups for each species, preserving missing gaps.
    fig,axs=plt.subplots(4,1,figsize=(13,10),sharex=True,constrained_layout=True)
    core=['species_pooled_gambiae','species_standard_gambiae','species_inverted_gambiae','species_pooled_coluzzii']
    labels=['A. gambiae pooled','A. gambiae standard homozygotes','A. gambiae inverted homozygotes','A. coluzzii pooled']
    for ax,g,label in zip(axs,core,labels):
        d=out[(out.resolution_bp==50000)&(out.group==g)].sort_values('start');ax.axvspan(INV[0]/1e6,INV[1]/1e6,color='#e5a344',alpha=.18,lw=0);ax.plot(d.mid/1e6,d.tajimas_D,color='#235d7a',lw=.7);ax.axhline(0,color='.5',lw=.6);ax.set_ylabel("Tajima's D\n"+label);ax.grid(axis='y',alpha=.2)
    axs[-1].set_xlabel('AgamP4 2L coordinate (Mb)');fig.suptitle("Anopheles 2La population and arrangement SFS tracks\nAg3.10 all-site calls; arrangement contrasts are geographically confounded",fontsize=12)
    for ext in ('pdf','svg','png'): fig.savefig(BASE/'figures'/f'anopheles_tajima_d_tracks.{ext}',dpi=300)
    audit['n_samples']=len(sample_ids);audit['sample_counts_by_species']={sp:sum(r['species']==sp for r in rows) for sp in ('gambiae','coluzzii')};audit['groups_omitted_below_two_individuals']=True;audit['result_windows']={'50kb':len(rows50),'100kb':len(win_defs[100000])};audit['groups']=audit['group_sample_counts'];audit['coordinate_assembly']='AgamP4';audit['analysis_limitations']='Species-wide arrangement groups span geography; locality and arrangement are partly confounded. Local groups included when n>=2 diploid individuals. Classical D is descriptive; no neutral p-values.'
    (BASE/'results/anopheles_tajima_d_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('COMPLETE',len(out),'rows',len(groups),'groups',len(summary),'regional summaries')
    print(pd.DataFrame(summary).query("resolution_bp==50000 and statistic=='tajimas_D'").to_string(index=False))
if __name__=='__main__':main()
