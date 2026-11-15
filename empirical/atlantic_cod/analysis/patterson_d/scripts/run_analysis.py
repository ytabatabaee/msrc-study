#!/usr/bin/env python3
"""Run the frozen Atlantic cod nucleotide Patterson D analysis."""
from __future__ import annotations
import csv, gzip, hashlib, json, math, random, sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from patterson_d_core import (read_designs,read_tree_windows,scan_vcf,add_window_stats,
    pool_windows,block_jackknife,d_from_sums,window_id,triple_topology,d_correlation,WINDOW)

REPO=HERE.parents[4]
AN=REPO/'empirical/atlantic_cod/analysis/patterson_d'
RAW=REPO/'data/atlantic_cod/raw/zenodo/local_downloads'
RESULTS=AN/'results'; FIG=AN/'figures'
TREES=REPO/'data/atlantic_cod/processed/cod_window_trees.tsv'
REGIONS=REPO/'data/atlantic_cod/metadata/region_manifest.tsv'
DESIGNFILE=AN/'config/frozen_quartets.tsv'


def read_tsv(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))

def write_tsv(path, rows, fields=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    if fields is None: fields=list(rows[0]) if rows else []
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader()
        for r in rows:w.writerow({k:fmt(v) for k,v in r.items()})

def fmt(v):
    if v is None:return 'NA'
    if isinstance(v,float):return f'{v:.10g}' if math.isfinite(v) else 'NA'
    return v

def read_regions():
    return {r['lg']:(int(r['start']),int(r['end'])) for r in read_tsv(REGIONS)}

def record_genome(path):
    h=hashlib.md5()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def sum_scenarios(designs, tree_windows, inside_results, outside_results):
    combined={}
    for d in designs:
        key=d['design_id']; lg=d['lg']; win=defaultdict(dict)
        add_window_stats(win,outside_results.get(key,{}).get('windows',{}))
        add_window_stats(win,inside_results.get(key,{}).get('windows',{}))
        combined[key]=dict(win)
    return combined

def classify_windows(lg, wins, interval):
    a,b=interval;out={}
    for wid,r in wins[lg].items():
        s,e=r['start'],r['end']
        if s>=a and e<=b: lab='inside'
        elif e<a or s>b: lab='outside'
        else: lab='boundary_overlap'
        out[wid]=lab
    return out

def pooled(windows, ids):return pool_windows(windows,ids)

def scenario_run(name, designs, tree_windows, regions, outside_path, outgroups=None,
                 sample_mode='frozen', min_call_fraction=1.0, transversions_only=False,
                 inside_source='snapp', keep_sites=False):
    dmap={d['lg']:d for d in designs}
    og=outgroups or {d['lg']:d['O_sample_id'] for d in designs}
    out=scan_vcf(outside_path,dmap,tree_windows,outgroups=og,min_call_fraction=min_call_fraction,
                 sample_mode=sample_mode,transversions_only=transversions_only,keep_sites=keep_sites)
    ins={}
    for d in designs:
        lg=d['lg']
        if inside_source=='snapp':fname=f'gadus_morhua_supergenes_{lg.lower()}_snapp.vcf.gz'
        else:fname=f'gadus_morhua_supergene_{lg.lower()}_dsuite.vcf.gz'
        path=RAW/fname
        ins.update(scan_vcf(path,{lg:d},tree_windows,outgroups={lg:og[lg]},min_call_fraction=min_call_fraction,
                            sample_mode=sample_mode,transversions_only=transversions_only,keep_sites=keep_sites))
    windows=sum_scenarios(designs,tree_windows,ins,out)
    rows={}
    for d in designs:
        key=d['design_id'];lg=d['lg'];w=windows[key];cl=classify_windows(lg,tree_windows,regions[lg])
        by={k:{wid for wid,v in cl.items() if v==k} for k in ('inside','outside','boundary_overlap')}
        vals={k:pooled(w,ids) for k,ids in by.items()}
        a,b=regions[lg]
        left={wid for wid,r in tree_windows[lg].items() if r['start']>=max(1,a-1_000_000) and r['end']<a}
        right={wid for wid,r in tree_windows[lg].items() if r['start']>b and r['end']<=b+1_000_000}
        vals['left_flank_1Mb']=pooled(w,left);vals['right_flank_1Mb']=pooled(w,right)
        rows[key]={'design':d,'windows':w,'class':cl,'regions':vals,'inside_results':ins.get(key,{}),'outside_results':out.get(key,{})}
        if keep_sites:
            rows[key]['sites']=out.get(key,{}).get('sites',[])+ins.get(key,{}).get('sites',[])
    return rows

def make_primary_window_table(primary, designs, tree_windows, regions):
    rows=[];lookup={}
    for d in designs:
        key=d['design_id'];lg=d['lg'];w=primary[key]['windows'];cl=primary[key]['class']
        for wid,meta in sorted(tree_windows[lg].items(),key=lambda x:x[1]['start']):
            s=w.get(wid,{})
            row={'lg':lg,'design_id':key,'window_id':wid,'start':meta['start'],'end':meta['end'],'midpoint':meta['midpoint'],
                 'region_class':cl[wid],'inside_inversion':cl[wid]=='inside','overlaps_inversion':cl[wid]=='boundary_overlap',
                 'n_accessible_biallelic_snps':int(s.get('biallelic_snps',0)),'n_polarized_sites':int(s.get('polarized_sites',0)),
                 'n_informative_sites':int(s.get('informative_sites',0)),'n_outgroup_hom_ref':int(s.get('outgroup_hom_ref',0)),
                 'n_outgroup_hom_alt':int(s.get('outgroup_hom_alt',0)),'n_outgroup_heterozygous':int(s.get('outgroup_heterozygous',0)),
                 'n_outgroup_missing':int(s.get('outgroup_missing',0)),'n_genotype_missingness_exclusions':int(s.get('missingness_exclusions',0)),
                 'sum_ABBA':s.get('sum_abba',0.0),'sum_BABA':s.get('sum_baba',0.0),
                 'D_site':s.get('D'),'n_P1':d['n_P1'],'n_P2':d['n_P2'],'n_P3':d['n_P3'],'n_O':d['n_O'],
                 'outgroup':'Borsai','data_source':'Zenodo SNAPP complete-call VCFs; positions restricted to exact published tree windows'}
            rows.append(row);lookup[(lg,wid)]=row
    return rows,lookup

def enrich_filter_inventory(path, tree_windows, rows):
    """Add raw-record/filter counts without reparsing all genotype fields."""
    idx={(r['lg'],r['window_id']):r for r in rows}
    countkeys=('n_vcf_records','n_nonpass_records','n_multiallelic_exclusions','n_non_snp_exclusions')
    regions=read_regions()
    for r in rows:
        for k in countkeys:r[k]=0
    for source in [path]+[RAW/f'gadus_morhua_supergenes_{lg.lower()}_snapp.vcf.gz' for lg in tree_windows]:
        if not source.exists():continue
        with gzip.open(source,'rt') as f:
            for line in f:
                if line.startswith('#'):continue
                v=line.split('\t',8);lg=v[0]
                if lg not in tree_windows:continue
                pos=int(v[1]);a,b=regions[lg]
                is_inside=source.name!='gadus_morhua_outside_supergenes_snapp.vcf.gz'
                if is_inside != (a<=pos<=b):continue
                wid,_,_=window_id(lg,pos)
                r=idx.get((lg,wid))
                if r is None:continue
                r['n_vcf_records']+=1
                if v[6] not in ('PASS','.'):r['n_nonpass_records']+=1;continue
                alt=v[4];ref=v[3]
                if ',' in alt:r['n_multiallelic_exclusions']+=1;continue
                if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT':r['n_non_snp_exclusions']+=1;continue
    for r in rows:
        r['n_outgroup_polarization_failures']=int(r['n_outgroup_missing'])+int(r['n_outgroup_heterozygous'])
        r['n_callable_polarized_sites']=int(r['n_polarized_sites'])
    return rows

def enrich_existing_filter_counts():
    designs=read_designs(DESIGNFILE);trees=read_tree_windows(TREES)
    path=RESULTS/'window_patterson_d.tsv';rows=read_tsv(path)
    rows=enrich_filter_inventory(RAW/'gadus_morhua_outside_supergenes_snapp.vcf.gz',trees,rows)
    write_tsv(path,rows)

def write_inside_outside(primary, designs, tree_windows):
    rows=[]
    for d in designs:
        key=d['design_id'];lg=d['lg'];regions=primary[key]['regions']
        inside=regions['inside'];outside=regions['outside'];delta=(inside['D']-outside['D']) if inside['D'] is not None and outside['D'] is not None else None
        rr={'lg':lg,'design_id':key,'primary_outgroup':'Borsai','inside_D':inside['D'],'outside_D':outside['D'],'inside_minus_outside_D':delta,
            'inside_ABBA_weight':inside['sum_abba'],'inside_BABA_weight':inside['sum_baba'],'outside_ABBA_weight':outside['sum_abba'],'outside_BABA_weight':outside['sum_baba'],
            'inside_informative_sites':int(inside['informative_sites']),'outside_informative_sites':int(outside['informative_sites']),
            'inside_accessible_snps':int(inside['biallelic_snps']),'outside_accessible_snps':int(outside['biallelic_snps']),
            'left_flank_1Mb_D':regions['left_flank_1Mb']['D'],'left_flank_1Mb_informative_sites':int(regions['left_flank_1Mb']['informative_sites']),
            'right_flank_1Mb_D':regions['right_flank_1Mb']['D'],'right_flank_1Mb_informative_sites':int(regions['right_flank_1Mb']['informative_sites']),
            'inside_250kb_windows':sum(1 for x in primary[key]['class'].values() if x=='inside'),
            'outside_250kb_windows':sum(1 for x in primary[key]['class'].values() if x=='outside'),
            'boundary_overlap_windows_excluded':sum(1 for x in primary[key]['class'].values() if x=='boundary_overlap'),
            'comparison_note':'Pooled ABBA/BABA contributions across complete 250-kb tree windows; D not averaged. Outside is same-linkage-group complement, not a neutral control.'}
        rows.append(rr)
    return rows

def boundary_rows(primary, designs, regions):
    rows=[]
    for d in designs:
        key=d['design_id'];lg=d['lg'];a,b=regions[lg];sites=primary[key]['sites']
        chrom_max=max((int(w.split('_')[2]) for _,w,_,_,_,_,_ in sites),default=b)
        for offset in (-100000,-50000,0,50000,100000):
            lo=max(1,a+offset);hi=min(chrom_max,b+offset);summ=defaultdict(float);counts=defaultdict(int)
            for pos,wid,abba,baba,pol,info,ost in sites:
                region='inside' if lo<=pos<=hi else 'outside'
                counts[region+'_accessible']+=1
                if pol:
                    counts[region+'_polarized']+=1;summ[region+'_abba']+=abba;summ[region+'_baba']+=baba
                    if info:counts[region+'_informative']+=1
            di=d_from_sums(summ['inside_abba'],summ['inside_baba']);do=d_from_sums(summ['outside_abba'],summ['outside_baba'])
            rows.append({'lg':lg,'design_id':key,'boundary_shift_bp':offset,'basis':'exact 1-based site coordinates within tree-window span',
                'inside_start':lo,'inside_end':hi,'inside_D':di,'outside_D':do,'inside_minus_outside_D':di-do if di is not None and do is not None else None,
                'inside_informative_sites':counts['inside_informative'],'outside_informative_sites':counts['outside_informative'],
                'inside_polarized_sites':counts['inside_polarized'],'outside_polarized_sites':counts['outside_polarized'],
                'inside_accessible_snps':counts['inside_accessible'],'outside_accessible_snps':counts['outside_accessible']})
    return rows

def spatial_null_rows(primary, designs, tree_windows, regions, window_rows):
    rows=[]
    bylookup={(r['lg'],r['window_id']):r for r in window_rows}
    for d in designs:
        lg=d['lg'];key=d['design_id'];a,b=regions[lg];maxend=max(r['end'] for r in tree_windows[lg].values())
        starts=list(range(1,maxend+1,WINDOW)); vals=[];labels=[]
        for s in starts:
            wid=f'{lg}_{s:09d}_{s+WINDOW-1:09d}'
            r=bylookup.get((lg,wid),{})
            vals.append((float(r.get('sum_ABBA',0) or 0),float(r.get('sum_BABA',0) or 0)))
            e=s+WINDOW-1
            labels.append('inside' if s>=a and e<=b else ('outside' if e<a or s>b else 'boundary'))
        null=[];n=len(starts)
        for shift in range(n):
            sums={'inside':[0.,0.],'outside':[0.,0.]}
            for i,lab in enumerate(labels):
                if lab not in sums:continue
                ab,ba=vals[(i-shift)%n];sums[lab][0]+=ab;sums[lab][1]+=ba
            di=d_from_sums(*sums['inside']);do=d_from_sums(*sums['outside'])
            if di is not None and do is not None:null.append(di-do)
        obs=primary[key]['regions']['inside']['D']-primary[key]['regions']['outside']['D']
        nullmean=float(np.mean(null)) if null else None
        p=(sum(abs(x-nullmean)>=abs(obs-nullmean) for x in null)+1)/(len(null)+1) if null and nullmean is not None else None
        rows.append({'lg':lg,'design_id':key,'method':'circular shift of 250-kb ABBA/BABA window contributions; inversion window mask fixed','n_shifts_valid':len(null),
            'observed_inside_minus_outside_D':obs,'null_mean':nullmean,'null_min':min(null) if null else None,'null_max':max(null) if null else None,
            'empirical_two_sided_p_centered':p,'note':'Physical 250-kb window order retained; missing-tree windows contribute zero; boundary-overlapping windows excluded.'})
    return rows

def resolution_rows(primary, designs, regions):
    rows=[]
    for d in designs:
        lg=d['lg'];key=d['design_id'];sites=primary[key]['sites']
        for size in (50000,100000,250000,500000):
            bins=defaultdict(lambda:defaultdict(float))
            for pos,wid,abba,baba,pol,info,og in sites:
                s=((pos-1)//size)*size+1;e=s+size-1;x=bins[s];x['accessible']+=1
                x['abba']+=abba;x['baba']+=baba
                if pol:x['polarized']+=1
                if info:x['informative']+=1
            a,b=regions[lg]
            for s,x in sorted(bins.items()):
                e=s+size-1
                cls='inside' if s>=a and e<=b else ('outside' if e<a or s>b else 'boundary_overlap')
                rows.append({'lg':lg,'design_id':key,'resolution_bp':size,'window_start':s,'window_end':e,'region_class':cls,
                    'accessible_biallelic_snps':int(x['accessible']),'polarized_sites':int(x['polarized']),'informative_sites':int(x['informative']),
                    'sum_ABBA':x['abba'],'sum_BABA':x['baba'],'D_site':d_from_sums(x['abba'],x['baba'])})
    return rows

def circular_block_corr(rows, xname, yname, nboot=1000, seed=20261115):
    from scipy.stats import pearsonr,spearmanr
    pairs=[r for r in rows if r.get(xname) is not None and r.get(yname) is not None]
    if len(pairs)<3:return len(pairs),None,None,None,None
    x=np.array([float(r[xname]) for r in pairs]);y=np.array([float(r[yname]) for r in pairs])
    if np.std(x)==0 or np.std(y)==0:return len(pairs),None,None,None,None
    rp=float(pearsonr(x,y)[0]);rs=float(spearmanr(x,y)[0])
    blocks=defaultdict(list)
    for r in pairs:blocks[(int(r['start'])-1)//1_000_000].append(r)
    keys=list(blocks);rng=random.Random(seed);pr=[];sr=[]
    if len(keys)>=3:
        for _ in range(nboot):
            sampled=[rng.choice(keys) for i in keys];sub=[r for k in sampled for r in blocks[k]]
            xx=[float(r[xname]) for r in sub];yy=[float(r[yname]) for r in sub]
            if len(xx)>=3 and np.std(xx)>0 and np.std(yy)>0:
                pr.append(float(pearsonr(xx,yy)[0]));sr.append(float(spearmanr(xx,yy)[0]))
    ci=lambda v:(float(np.quantile(v,.025)),float(np.quantile(v,.975))) if v else (None,None)
    return len(pairs),rp,rs,ci(pr),ci(sr)


def main():
    RESULTS.mkdir(parents=True,exist_ok=True); FIG.mkdir(parents=True,exist_ok=True)
    designs=read_designs(DESIGNFILE); trees=read_tree_windows(TREES); regions=read_regions()
    outside=RAW/'gadus_morhua_outside_supergenes_snapp.vcf.gz'
    primary=scenario_run('primary',designs,trees,regions,outside,keep_sites=True)
    wrows,wlookup=make_primary_window_table(primary,designs,trees,regions)
    wrows=enrich_filter_inventory(outside,trees,wrows)
    write_tsv(RESULTS/'window_patterson_d.tsv',wrows)
    inside_rows=write_inside_outside(primary,designs,trees)
    write_tsv(RESULTS/'inside_outside_summary.tsv',inside_rows)
    jack=[]
    for d in designs:
        x=primary[d['design_id']]
        for size in (250000,500000,1000000):
            for r in block_jackknife(x['windows'],x['class'],size): jack.append({'lg':d['lg'],'design_id':d['design_id'],**r})
    write_tsv(RESULTS/'block_jackknife.tsv',jack)
    write_tsv(RESULTS/'spatial_null_validation.tsv',spatial_null_rows(primary,designs,trees,regions,wrows))
    write_tsv(RESULTS/'boundary_sensitivity.tsv',boundary_rows(primary,designs,regions))
    write_tsv(RESULTS/'resolution_patterson_d.tsv',resolution_rows(primary,designs,regions))

    # SNAPP trees omit Borsai; report a rooted three-population genealogy contrast only.
    topology=[]; corrs=[]
    stage4={r['window_id']:r for r in read_tsv(REPO/'data/atlantic_cod/processed/stage4a_window_topology_scores.tsv')}
    stage7={r['window_id']:r for r in read_tsv(REPO/'empirical/atlantic_cod/results/stage7_window_divergence_time_signal.tsv')}
    stage7b={r['window_id']:r for r in read_tsv(REPO/'empirical/atlantic_cod/results/stage7b_msrc_coalescent_window_signal.tsv')}
    for d in designs:
        lg=d['lg'];key=d['design_id'];pairset={"q1":0,"q2":0,"q3":0}; rows=[]
        for tr in trees[lg].values():
            labels={d['P1']+'_spc':d['P1'],d['P2']+'_spc':d['P2'],d['P3']+'_spc':d['P3']}
            q=triple_topology(tr['tree_newick'],tuple(labels)) if tr.get('tree_newick') else None
            dat=wlookup.get((lg,tr['window_id']),{})
            s4=stage4.get(tr['window_id'],{});s7=stage7.get(tr['window_id'],{});s7b=stage7b.get(tr['window_id'],{})
            r={'lg':lg,'design_id':key,'window_id':tr['window_id'],'start':tr['start'],'end':tr['end'],'region_class':dat.get('region_class'),
               'D_site':dat.get('D_site'),'sum_ABBA':dat.get('sum_ABBA'),'sum_BABA':dat.get('sum_BABA'),
               'q1_triple':int(q=='q1'),'q2_triple':int(q=='q2'),'q3_triple':int(q=='q3'),
               'rooted_triple_topology':q,'rooted_triple_contrast_q3_minus_q2':(1 if q=='q3' else (-1 if q=='q2' else 0)) if q else None,
               'D_topology_exact_four_taxon':None,'exact_match_status':'not estimable: Borsai outgroup absent from all 426 local population trees',
               'D_arrangement_minus_baseline':s4.get('D_arrangement_minus_baseline'),
               'A_opposite_minus_same':s7.get('A_opposite_minus_same'),'A_res':s7b.get('A_res'),
               'tree_branch_length_units':'published SNAPP population divergence-time units; not substitution or coalescent units'}
            topology.append(r);rows.append(r)
            if q:pairset[q]+=1
        for signal in ('rooted_triple_contrast_q3_minus_q2','D_arrangement_minus_baseline','A_opposite_minus_same','A_res'):
            n,rp,rs,ci_p,ci_s=circular_block_corr(rows,'D_site',signal)
            corrs.append({'lg':lg,'design_id':key,'signal':signal,'n_matched_windows':n,'pearson_r':rp,'spearman_rho':rs,
                          'pearson_1Mb_block_bootstrap_ci_low':ci_p[0],'pearson_1Mb_block_bootstrap_ci_high':ci_p[1],
                          'spearman_1Mb_block_bootstrap_ci_low':ci_s[0],'spearman_1Mb_block_bootstrap_ci_high':ci_s[1],
                          'inference_note':'Descriptive spatial association; 1-Mb block bootstrap, not independent-site significance.'})
    write_tsv(RESULTS/'topology_D_comparison.tsv',topology)
    write_tsv(RESULTS/'spatial_correlations.tsv',corrs)

    # Robustness was predeclared: other available single-individual outgroups, sample thinning,
    # transversions, and Dsuite's more-missing genotype call set at 50% population callability.
    robustness_path=RESULTS/'robustness_analysis.tsv'
    robust=read_tsv(robustness_path) if robustness_path.exists() else []
    scenarios=[('primary_fixed_samples_Borsai',{},'frozen',1.0,False,'snapp'),
      ('transversions_Borsai',{},'frozen',1.0,True,'snapp'),
      ('one_individual_per_population_Borsai',{},'first_available',1.0,False,'snapp')]
    # Gadoga is selected from the independently documented baseline outgroup set;
    # the remaining candidates are inventoried for future extended sensitivity runs.
    scenarios.append(('outgroup_Gadoga',{lg:'Gadoga' for lg in ('LG01','LG02','LG07','LG12')},'frozen',1.0,False,'snapp'))
    scenarios.append(('Dsuite_roster_50pct_callability',{},'all_available',0.5,False,'dsuite'))
    for name,ogs,mode,frac,tv,src in scenarios:
        already={r['lg'] for r in robust if r['scenario']==name}
        if already=={d['lg'] for d in designs}:
            print('Robustness scenario already checkpointed:',name,flush=True)
            continue
        robust=[r for r in robust if r['scenario']!=name]
        res=scenario_run(name,designs,trees,regions,outside,outgroups=ogs or None,sample_mode=mode,min_call_fraction=frac,
                         transversions_only=tv,inside_source=src,keep_sites=False)
        for d in designs:
            x=res[d['design_id']]['regions'];i=x['inside'];o=x['outside']
            robust.append({'scenario':name,'lg':d['lg'],'design_id':d['design_id'],'outgroup':(ogs or {}).get(d['lg'],'Borsai'),
              'sample_mode':mode,'minimum_population_call_fraction':frac,'transversions_only':tv,'inside_source':src,
              'inside_D':i['D'],'outside_D':o['D'],'inside_minus_outside_D':i['D']-o['D'] if i['D'] is not None and o['D'] is not None else None,
              'inside_informative_sites':int(i['informative_sites']),'outside_informative_sites':int(o['informative_sites']),
              'caveat':'single-individual alternate outgroups; Dsuite roster replaces Gadmor_lfo2 with Gadmor_lfo3; 50% call rule.'})
        write_tsv(robustness_path,robust)
        print('Robustness scenario complete:',name,flush=True)

    # Determine input digests and record precise command-level provenance.
    inputs=[outside]+[RAW/f'gadus_morhua_supergenes_{d["lg"].lower()}_snapp.vcf.gz' for d in designs]
    manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'analysis':'Patterson D from diploid population allele frequencies; fixed homozygous outgroup polarization',
       'primary_quartet_config':str(DESIGNFILE.relative_to(REPO)),'window_table':str(TREES.relative_to(REPO)),
       'reference':'gadMor2','coordinates':'1-based inclusive','primary_window_size_bp':250000,
       'genotype_policy':'biallelic SNPs; PASS/.; at least one called diploid in each ingroup; outgroup genotype must be homozygous; missing/polymorphic outgroup excluded',
       'inputs':[{'path':str(p.relative_to(REPO)),'bytes':p.stat().st_size,'md5':record_genome(p)} for p in inputs],
       'python':sys.version,'outputs':['window_patterson_d.tsv','inside_outside_summary.tsv','block_jackknife.tsv','spatial_null_validation.tsv','topology_D_comparison.tsv','robustness_analysis.tsv']}
    (RESULTS/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Completed',len(wrows),'matched 250-kb windows')
    for r in inside_rows: print(r['lg'],'inside',fmt(r['inside_D']),'outside',fmt(r['outside_D']),'delta',fmt(r['inside_minus_outside_D']))

if __name__=='__main__':
    if sys.argv[1:]==['--filter-audit-only']:enrich_existing_filter_counts()
    else:main()
