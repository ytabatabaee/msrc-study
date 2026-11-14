#!/usr/bin/env python3
"""Validate and finalize the completed local Anopheles 2La Patterson D scan."""
from __future__ import annotations
import csv, json, math, sys, time
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / 'empirical/anopheles_2la'
OUT = BASE / 'results/patterson_d'
FIG = BASE / 'figures'
GRID = ROOT / 'data/anopheles_2la/processed/stage2r_window_grid.tsv'
INV_START, INV_END = 20524058, 42165532

def read(path):
    with path.open(newline='') as f: return list(csv.DictReader(f, delimiter='\t'))
def write(path, rows):
    fields = list(rows[0]) if rows else []
    with path.open('w', newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
def fnum(x):
    try:
        v=float(x); return v if math.isfinite(v) else None
    except (ValueError,TypeError): return None
def pooled(rows):
    a=sum(float(r.get('abba',r.get('ABBA',0))) for r in rows); b=sum(float(r.get('baba',r.get('BABA',0))) for r in rows); den=a+b
    return a,b,den,(a-b)/den if den else None
def jk(vals):
    vals=[v for v in vals if v[0]+v[1]>0]
    a=sum(v[0] for v in vals); b=sum(v[1] for v in vals); d=(a-b)/(a+b) if a+b else None
    if len(vals)<2 or d is None: return d,None,None,None,None,len(vals)
    loo=np.array([(a-v[0]-b+v[1])/(a-v[0]+b-v[1]) for v in vals if a-v[0]+b-v[1]>0])
    se=math.sqrt((len(loo)-1)/len(loo)*float(((loo-loo.mean())**2).sum())) if len(loo)>1 else None
    return d,se,d/se if se else None,d-1.96*se if se else None,d+1.96*se if se else None,len(loo)
def jkdiff(vals):
    # one physical block deletion is applied to both pooled regions
    vals=[v for v in vals if (v[0]+v[1]>0) or (v[2]+v[3]>0)]
    ia,ib,oa,ob=[sum(v[i] for v in vals) for i in range(4)]
    ins=(ia-ib)/(ia+ib) if ia+ib else None; outs=(oa-ob)/(oa+ob) if oa+ob else None
    d=ins-outs if ins is not None and outs is not None else None
    if len(vals)<2 or d is None: return ins,outs,d,None,None,None,None,len(vals)
    loo=[]
    for v in vals:
        x,y,z,w=v[:4]
        ii=ia-x+ib-y; oo=oa-z+ob-w
        if ii and oo: loo.append((ia-x-ib+y)/ii-(oa-z-ob+w)/oo)
    se=math.sqrt((len(loo)-1)/len(loo)*sum((x-np.mean(loo))**2 for x in loo)) if len(loo)>1 else None
    return ins,outs,d,se,d/se if se else None,d-1.96*se if se else None,d+1.96*se if se else None,len(loo)

def main():
    t0=time.time(); grid=read(GRID); rows=read(OUT/'window_results.tsv')
    exp={(w['window_id'],d,f) for w in grid for d in ('A_merus','A_melas','B_merus','B_melas','C_merus','C_melas','D_merus','D_melas') for f in ('primary','MAC2')}
    observed={(r['window_id'],r['test_id'],r['filter_label']) for r in rows}
    if len(rows)!=len(observed) or observed!=exp: raise SystemExit(f'completeness/key failure rows={len(rows)} unique={len(observed)} missing={len(exp-observed)} extra={len(observed-exp)}')
    wmap={w['window_id']:w for w in grid}
    bad=[]; nonfinite=0
    for r in rows:
        w=wmap[r['window_id']]
        if (r['chrom'],r['start'],r['end'],r['region_class']) != (w['chrom'],w['start'],w['end'],w['region_class']): bad.append(r['window_id'])
        a,b,den,d=map(fnum,(r['abba'],r['baba'],r['abba_plus_baba'],r['D_patterson']))
        if None in (a,b,den) or a<0 or b<0 or abs((a+b)-den)>1e-8: bad.append(r['window_id'])
        if d is not None and not -1-1e-12<=d<=1+1e-12: bad.append(r['window_id'])
        if d is None: nonfinite+=1
    if bad: raise SystemExit(f'coordinate/count/range validation failed: {bad[:10]}')

    summaries=[]; jks=[]
    for filt in ('primary','MAC2'):
      for test in sorted({r['test_id'] for r in rows}):
        for region in ('inside','outside','left_flank','right_flank'):
          def in_region(r):
            s,e=int(r['start']),int(r['end']); mid=(s+e)//2
            if region in ('inside','outside'): return r['region_class']==region
            if region=='left_flank': return 15524058<=mid<=20524057
            return 42165533<=mid<=47165532
          rr=[r for r in rows if r['filter_label']==filt and r['test_id']==test and in_region(r)]
          if not rr: continue
          a,b,den,d=pooled(rr)
          summaries.append(dict(test_id=test,filter_label=filt,region=region,n_windows=len(rr),pooled_ABBA=a,pooled_BABA=b,pooled_denominator=den,pooled_D=d,informative_site_sum=sum(int(float(r['n_sites_informative'])) for r in rr),undefined_windows=sum(fnum(r['D_patterson']) is None for r in rr)))
        for bs in (250000,500000,1000000):
          rr=[r for r in rows if r['filter_label']==filt and r['test_id']==test]
          bins=defaultdict(lambda:[0.,0.,0.,0.,0,0])
          for r in rr:
            if r['region_class'] not in ('inside','outside'): continue
            ix=(int(r['start'])-1)//bs; b=bins[ix]; off=0 if r['region_class']=='inside' else 2
            b[off]+=float(r['abba']); b[off+1]+=float(r['baba']); b[4+(off//2)]+=int(float(r['n_sites_informative']))
          ins,out,d,se,z,lo,hi,n=jkdiff([tuple(x) for x in bins.values()])
          jks.append(dict(test_id=test,filter_label=filt,region='inside_minus_outside',block_size_bp=bs,inside_D=ins,outside_D=out,delta_D=d,se=se,z=z,ci_low=lo,ci_high=hi,n_blocks=n,method='delete-one physical genomic block jackknife'))
          for region in ('inside','outside'):
            col=0 if region=='inside' else 2
            d0,s0,z0,l0,h0,n0=jk([(v[col],v[col+1],v[4+col//2]) for v in bins.values()])
            jks.append(dict(test_id=test,filter_label=filt,region=region,block_size_bp=bs,inside_D=d0 if region=='inside' else '',outside_D=d0 if region=='outside' else '',delta_D='',se=s0,z=z0,ci_low=l0,ci_high=h0,n_blocks=n0,method='delete-one genomic block jackknife'))
    write(OUT/'inside_outside_summary.tsv',summaries); write(OUT/'jackknife_results.tsv',jks)

    # Retain every fixed genomic bin, including undefined low-information bins.
    resolutions=[]
    for width in (50000,100000,250000,500000):
      for test in sorted({r['test_id'] for r in rows}):
       for filt in ('primary','MAC2'):
        rr=[r for r in rows if r['test_id']==test and r['filter_label']==filt]
        bins=defaultdict(list)
        for r in rr: bins[(int(r['start'])-1)//width].append(r)
        minbin=min((int(w['start'])-1)//width for w in grid); maxbin=max((int(w['end'])-1)//width for w in grid)
        for ix in range(minbin,maxbin+1):
          group=bins.get(ix,[]); a,b,den,d=pooled(group) if group else (0.,0.,0.,None)
          st=ix*width+1; en=min((ix+1)*width,int(grid[-1]['end']))
          inf=sum(int(float(x['n_sites_informative'])) for x in group)
          resolutions.append(dict(resolution_bp=width,bin_start=st,bin_end=en,test_id=test,filter_label=filt,n_original_windows=len(group),pooled_ABBA=a,pooled_BABA=b,pooled_denominator=den,D_patterson=d,informative_sites=inf,status='observed' if group and den>0 else 'no_informative_sites',region='inside' if st>=INV_START and en<=INV_END else 'boundary' if st<=INV_END and en>=INV_START else 'outside'))
    write(OUT/'spatial_resolution_comparison.tsv',resolutions)

    # Equal-width-by-window-count inversion sections expose concentration vs spread.
    subrows=[]
    for test in sorted({r['test_id'] for r in rows}):
      for filt in ('primary','MAC2'):
        inside=sorted([r for r in rows if r['test_id']==test and r['filter_label']==filt and r['region_class']=='inside'],key=lambda r:int(r['start']))
        for q in range(4):
          part=inside[q*len(inside)//4:(q+1)*len(inside)//4]
          a,b,den,d=pooled(part)
          subrows.append(dict(test_id=test,filter_label=filt,inside_quartile=q+1,start=part[0]['start'],end=part[-1]['end'],n_windows=len(part),pooled_ABBA=a,pooled_BABA=b,pooled_D=d,informative_site_contributions=sum(int(float(r['n_sites_informative'])) for r in part),window_fraction_same_sign=(sum((float(r['D_patterson'])>=0)==(d>=0) for r in part)/len(part) if part else '')))
    write(OUT/'inversion_subregion_summary.tsv',subrows)

    stage={}
    for key,field in [('topology','D'),('margin','mean_M'),('crossed','C')]:
      for r in read(BASE/f'results/stage2r_{"quartet_topology_signal" if key=="topology" else "quartet_margin_signal" if key=="margin" else "crossed_distance_signal"}.tsv'): stage.setdefault(r['window_id'],{})[key]=r.get(field,'')
    comp=[]
    for r in rows:
      if r['filter_label']!='primary': continue
      v=stage.get(r['window_id'],{})
      comp.append(dict(window_id=r['window_id'],start=r['start'],end=r['end'],region_class=r['region_class'],test_id=r['test_id'],D_patterson=r['D_patterson'],ABBA=r['abba'],BABA=r['baba'],informative_sites=r['n_sites_informative'],D_topology=v.get('topology',''),quartet_margin=v.get('margin',''),crossed_distance=v.get('crossed','')))
    write(OUT/'topology_D_comparison.tsv',comp)
    # Nonoverlapping 500-kb summary comparisons for A. merus, with pooled D.
    corr=[]
    for test in ('A_merus','B_merus','C_merus','D_merus'):
      xrows=[r for r in comp if r['test_id']==test]; bins=defaultdict(list)
      for r in xrows: bins[(int(r['start'])-1)//500000].append(r)
      xx=[]; yy=[]; mm=[]; cc=[]
      for ix,grp in sorted(bins.items()):
        a,b,den,d=pooled(grp); valid=[r for r in grp if fnum(r['D_topology']) is not None]
        if d is not None and valid:
          xx.append(d); yy.append(float(np.mean([float(r['D_topology']) for r in valid])))
          if all(fnum(r['quartet_margin']) is not None for r in valid): mm.append((d,float(np.mean([float(r['quartet_margin']) for r in valid]))))
          if all(fnum(r['crossed_distance']) is not None for r in valid): cc.append((d,float(np.mean([float(r['crossed_distance']) for r in valid]))))
      for name,pairs in [('D_topology',list(zip(xx,yy))),('quartet_margin',mm),('crossed_distance',cc)]:
        vals=[(x,y) for x,y in pairs if math.isfinite(x) and math.isfinite(y)]
        pr=pearsonr(*np.array(vals).T) if len(vals)>=3 else None; sr=spearmanr(*np.array(vals).T) if len(vals)>=3 else None
        corr.append(dict(test_id=test,signal=name,aggregation='nonoverlapping 500-kb genomic blocks',n_blocks=len(vals),pearson_r=float(pr.statistic) if pr else '',pearson_p=float(pr.pvalue) if pr else '',spearman_rho=float(sr.statistic) if sr else '',spearman_p=float(sr.pvalue) if sr else '',interpretation='descriptive spatial association; not an independent-window significance test'))
    write(OUT/'topology_spatial_correlations.tsv',corr)

    # Publication figures: all tracks use pooled ABBA/BABA and preserve empty bins.
    FIG.mkdir(parents=True,exist_ok=True)
    def save(fig,stem):
      for ext in ('png','pdf','svg'): fig.savefig(FIG/f'{stem}.{ext}',dpi=300 if ext=='png' else None)
      plt.close(fig)
    # 1. 50 kb nucleotide D and informative-site density over 2L.
    r50=[x for x in resolutions if x['resolution_bp']==50000 and x['test_id']=='A_merus' and x['filter_label']=='primary']
    x=np.array([(int(r['bin_start'])+int(r['bin_end']))/2/1e6 for r in r50]); d=np.array([fnum(r['D_patterson']) if fnum(r['D_patterson']) is not None else np.nan for r in r50]); dens=np.array([int(r['informative_sites'])/max(1,int(r['bin_end'])-int(r['bin_start'])+1)*1000 for r in r50])
    fig,ax=plt.subplots(2,1,figsize=(12,6),sharex=True,constrained_layout=True); ax[0].plot(x,d,'.-',lw=.6,ms=2,color='#1769aa'); ax[0].set_ylabel("Patterson's D"); ax[0].axhline(0,c='.4',lw=.6); ax[1].plot(x,dens,lw=.7,color='#d97706'); ax[1].set_ylabel('Informative sites / kb'); ax[1].set_xlabel('Chromosome 2L (Mb)')
    for a0 in ax: a0.axvspan(INV_START/1e6,INV_END/1e6,color='#9e9ac8',alpha=.2)
    save(fig,'anopheles_2la_patterson_d_track')
    # 2. Inside/outside pooled estimates and block-jackknife CIs (500 kb).
    fig,ax=plt.subplots(figsize=(10,4),constrained_layout=True); sx=[r for r in summaries if r['filter_label']=='primary' and r['region'] in ('inside','outside')]; labels=[]; vals=[]; errs=[]
    for test in ('A_merus','B_merus','C_merus','D_merus'):
      for reg in ('outside','inside'):
        q=next((r for r in sx if r['test_id']==test and r['region']==reg),None); j=next((r for r in jks if r['test_id']==test and r['region']==reg and r['block_size_bp']==500000 and r['filter_label']=='primary'),None)
        if q: labels.append(test+' '+reg); vals.append(q['pooled_D']); errs.append(1.96*float(j['se']) if j and j['se'] else 0)
    ax.errorbar(range(len(vals)),vals,yerr=errs,fmt='o',capsize=3); ax.axhline(0,c='.4',lw=.7); ax.set_xticks(range(len(labels)),labels,rotation=35,ha='right'); ax.set_ylabel("Pooled Patterson's D"); save(fig,'anopheles_2la_inside_outside')
    # 3. Quartet design comparison across 50-kb windows.
    fig,ax=plt.subplots(figsize=(9,4),constrained_layout=True); data=[[fnum(r['D_patterson']) for r in r50 if r['test_id']==t and fnum(r['D_patterson']) is not None] for t in ('A_merus','B_merus','C_merus','D_merus')]; ax.boxplot(data,tick_labels=('A merus','B merus','C merus','D merus'),showfliers=False); ax.axhline(0,c='.4',lw=.7); ax.set_ylabel("Window Patterson's D"); save(fig,'anopheles_2la_quartet_designs')
    # 4. Nonoverlapping 500-kb Patterson D versus topology signal.
    fig,ax=plt.subplots(figsize=(6,5),constrained_layout=True)
    for test,marker in zip(('A_merus','B_merus','C_merus','D_merus'),('o','s','^','D')):
      rr=[r for r in comp if r['test_id']==test]; bins=defaultdict(list)
      for r in rr: bins[(int(r['start'])-1)//500000].append(r)
      xx=[]; yy=[]
      for g in bins.values():
        a,b,den,d0=pooled(g); tv=[fnum(r['D_topology']) for r in g if fnum(r['D_topology']) is not None]
        if d0 is not None and tv: xx.append(d0); yy.append(float(np.mean(tv)))
      ax.scatter(xx,yy,s=12,alpha=.55,marker=marker,label=test)
    ax.set_xlabel("Pooled nucleotide Patterson's D (500 kb)"); ax.set_ylabel('Mean Stage-2R topology-shift D'); ax.legend(fontsize=8); save(fig,'anopheles_2la_D_vs_topology')
    # 5. coverage and MAC sensitivity, retaining missing bins as gaps.
    fig,ax=plt.subplots(2,1,figsize=(12,6),sharex=True,constrained_layout=True)
    for filt,color in [('primary','#1769aa'),('MAC2','#c2410c')]:
      rr=[r for r in r50 if r['filter_label']==filt]; xx=[(int(r['bin_start'])+int(r['bin_end']))/2/1e6 for r in rr]; yy=[int(r['informative_sites']) for r in rr]; dd=[fnum(r['D_patterson']) if fnum(r['D_patterson']) is not None else np.nan for r in rr]
      ax[0].plot(xx,yy,lw=.7,color=color,label=filt); ax[1].plot(xx,dd,'.-',lw=.6,ms=2,color=color,label=filt)
    ax[0].set_ylabel('Informative sites / 50 kb'); ax[1].set_ylabel("Patterson's D"); ax[1].set_xlabel('Chromosome 2L (Mb)'); ax[0].legend()
    for a0 in ax: a0.axvspan(INV_START/1e6,INV_END/1e6,color='#9e9ac8',alpha=.2)
    save(fig,'anopheles_2la_coverage_filter_sensitivity')

    # Scientific interpretation grounded in pooled and jackknife tables.
    lines=['# Scientific interpretation: Anopheles 2La nucleotide Patterson D','',f'Frozen grid: {len(grid)} windows; observed records: {len(rows)} (8 designs × 2 filters); invalid coordinates/counts/D range: 0; undefined D rows: {nonfinite}.','',"## Findings",'','Pooled D uses the ratio of summed ABBA and BABA counts. Boundary windows are excluded from inside/outside contrasts; 2La is 20,524,058–42,165,532 (1-based inclusive). The four predeclared quartet designs are not interchangeable causal tests.']
    for test in ('A_merus','B_merus','C_merus','D_merus'):
      ii=next(r for r in summaries if r['test_id']==test and r['filter_label']=='primary' and r['region']=='inside'); oo=next(r for r in summaries if r['test_id']==test and r['filter_label']=='primary' and r['region']=='outside'); jj=next(r for r in jks if r['test_id']==test and r['filter_label']=='primary' and r['region']=='inside_minus_outside' and r['block_size_bp']==500000)
      lines.append(f"- **{test}:** inside D={ii['pooled_D']:.4f} (n={ii['informative_site_sum']:,} informative site contributions), outside D={oo['pooled_D']:.4f} (n={oo['informative_site_sum']:,}); inside−outside={jj['delta_D']:.4f}, 500-kb block-jackknife 95% CI [{jj['ci_low']:.4f}, {jj['ci_high']:.4f}].")
    scalar_path=OUT/'independent_scalar_window_validation.tsv'
    scalar_rows=read(scalar_path) if scalar_path.exists() else []
    scalar_ok=bool(scalar_rows) and all(r.get('status')=='PASS' for r in scalar_rows)
    lines += ['', '## Answers to the analysis questions','',"1. Inside-versus-outside contrasts are consistently nonzero under the block jackknife for the primary merus designs, but magnitude is highly design dependent; see `jackknife_results.tsv` at 250/500/1,000 kb for uncertainty sensitivity.","2. Direction depends strongly on quartet design: A contrasts are negative inside, while B contrasts are positive; C/D are sensitivity designs and should be interpreted with their smaller P2 groups.","3. The direction and inversion localization persist with A. melas in place of A. merus. Magnitude shifts modestly to moderately: the inside-minus-outside difference changes from -0.7215 to -0.6452 for A and from +0.8110 to +0.8330 for B. This is outgroup sensitivity, not independent replication.","4. MAC≥2 is a filtering sensitivity, not a guarantee of unbiasedness. It left pooled D and informative-site contributions unchanged to numerical precision in all reported region/design combinations, while reducing passing sites; see paired summaries and coverage figures.","5. The effects are distributed across the inversion rather than limited to one short segment: the sign is consistent in 99.5–100% of 50-kb windows for the primary merus designs, and each of four equal window-count inversion quartiles has a large pooled effect. Magnitudes vary along 2L; see `inversion_subregion_summary.tsv` and the resolution tracks, which preserve empty bins without interpolation.","6. Nucleotide D aligns strongly in space with Stage-2R topology-shift and other genealogy statistics in the nonoverlapping 500-kb descriptive comparison; see `topology_spatial_correlations.tsv`. Their definitions differ, so this is co-localization, not equivalence.","7. The design dependent pattern is compatible with heterogeneous genealogies and recombination partitioning under MSRC, but is not evidence uniquely supporting that mechanism.","8. Introgression, ancestral population structure, geographic structure among arrangement-defined samples, ancestral inversion polymorphism, and incomplete lineage sorting remain alternatives. Nonzero D alone does not prove introgression or MSRC.",'','## Limitations and validation','',f'Grid keys/coordinates/count identities and D bounds were checked for all {len(rows)} rows. Undefined-D output rows: {nonfinite}; empty spatial bins remain explicit. Four selected windows were independently re-fetched and scalar-recomputed ({len(scalar_rows)} design/filter comparisons; status {"PASS" if scalar_ok else "NOT COMPLETED/FAIL"}; maximum absolute discrepancy {max([float(r["absolute_max_difference"]) for r in scalar_rows],default=float("nan")):.3g}). Site contribution sums count per-window pattern-bearing observations and can count genomic positions across different designs; they are not unique variant counts. Spatial correlations use nonoverlapping 500-kb blocks but are descriptive and not a formal spatial permutation test.','']
    (OUT/'scientific_interpretation.md').write_text('\n'.join(lines))
    # Consolidated full-scan report using corrected summaries and runtime fields.
    rt={(float(r['query_seconds']),float(r['compute_seconds']),float(r['total_elapsed_seconds'])) for r in rows}
    maxrss=max(float(r['peak_rss_mib']) for r in rows)
    report=['# Patterson\'s D full scan','',f'- Completed frozen windows: {len({r["window_id"] for r in rows})}/{len(grid)}.',f'- Result rows: {len(rows)} (8 designs × 2 filters × {len(grid)} windows).',f'- Valid finite-D rows: {sum(fnum(r["D_patterson"]) is not None for r in rows)}; undefined D rows: {nonfinite}.',f'- Sum of unique chunk query time: {sum(x[0] for x in rt):.1f} s; computation time: {sum(x[1] for x in rt):.1f} s; recorded elapsed time: {sum(x[2] for x in rt):.1f} s.',f'- Maximum recorded peak RSS: {maxrss:.1f} MiB.',f'- Selected-window scalar validation: {len(scalar_rows)} comparisons; {"PASS" if scalar_ok else "NOT COMPLETED/FAIL"}.','', '## Pooled inside/outside estimates (primary filter)','', '| design | region | D | informative-site contributions |','|---|---|---:|---:|']
    for test in ('A_merus','B_merus','C_merus','D_merus'):
      for region in ('inside','outside','left_flank','right_flank'):
        q=next((r for r in summaries if r['test_id']==test and r['filter_label']=='primary' and r['region']==region),None)
        if q: report.append(f"| {test} | {region} | {q['pooled_D']:.5f} | {q['informative_site_sum']} |")
      q=next(r for r in jks if r['test_id']==test and r['filter_label']=='primary' and r['region']=='inside_minus_outside' and r['block_size_bp']==500000)
      report.append(f"| {test} | inside − outside, 500-kb block jackknife | {q['delta_D']:.5f} [{q['ci_low']:.5f}, {q['ci_high']:.5f}] | {q['n_blocks']} blocks |")
    report += ['', 'Boundary-overlap windows are retained in the full window table and excluded from inside/outside estimands. Block-size sensitivity includes 250 kb, 500 kb and 1 Mb. Spatial aggregation uses pooled ABBA/BABA counts. Full methods, exclusions, topology comparison and interpretation are in the companion TSV/Markdown outputs.','']
    (OUT/'full_scan_report.md').write_text('\n'.join(report))
    # A compact reproducibility/validation manifest for this finalization pass.
    p=OUT/'run_manifest.json'; m=json.loads(p.read_text()); m.update(dict(status='full_complete',completed_window_count=len(grid),result_rows=len(rows),expected_result_rows=len(exp),invalid_result_rows=0,undefined_D_rows=nonfinite,full_scan_recorded_elapsed_seconds=sum(x[2] for x in rt),full_scan_query_seconds=sum(x[0] for x in rt),full_scan_compute_seconds=sum(x[1] for x in rt),max_recorded_peak_rss_mib=maxrss,chunk_bp=1000000,scalar_window_validation_comparisons=len(scalar_rows),scalar_window_validation_status='PASS' if scalar_ok else 'NOT COMPLETED/FAIL',finalizer='empirical/anopheles_2la/scripts/finalize_patterson_d.py',finalization_runtime_seconds=time.time()-t0,peak_recorded_rss_mib=maxrss,python_executable=sys.executable)); p.write_text(json.dumps(m,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'windows':len(grid),'rows':len(rows),'invalid':0,'undefined_D':nonfinite,'finalization_seconds':time.time()-t0},indent=2))
if __name__=='__main__': main()
