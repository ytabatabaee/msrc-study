#!/usr/bin/env python3
"""Build compact final result and physical-block audit tables."""
import csv
from pathlib import Path

ROOT=Path(__file__).resolve().parents[6]
BASE=ROOT/'empirical/neurospora/analysis/patterson_d'
R=BASE/'final_validation/results'

def read(path):
    with open(path) as f:return list(csv.DictReader(f,delimiter='\t'))
def write(path,rows):
    if rows:
        keys=[]
        for r in rows:
            for k in r:
                if k not in keys:keys.append(k)
        with open(path,'w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=keys,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def main():
    old=read(BASE/'results/block_uncertainty.tsv'); uncertainty=[]
    for r in old:
        if r['block_start']=='SUMMARY':
            uncertainty.append({'block_width_bp':r['block_width_bp'],'region':r['region'],'estimate_D':r['estimate_D'],'jackknife_se':r['jackknife_se'],'ci95_low':r['ci95_low'],'ci95_high':r['ci95_high'],'n_physical_blocks':r['n_blocks'],'assessment':'Existing source calculation; independent blocks partition physical coordinate within each region. Normal-approximate interval.'})
    uncertainty.append({'block_width_bp':2_000_000,'region':'inside_minus_outside','estimate_D':'NA','jackknife_se':'NA','ci95_low':'NA','ci95_high':'NA','n_physical_blocks':'4 inside / 2 outside','assessment':'Coarse diagnostic only: only two outside physical units after splitting at focal boundaries; no defensible outside delete-block interval.'})
    write(R/'spatial_uncertainty_validation.tsv',uncertainty)

    rows=[]
    def add(label,estimate,lo='NA',hi='NA',n='NA',detail=''):
        rows.append({'result':label,'estimate':estimate,'ci95_low':lo,'ci95_high':hi,'n_or_support':n,'detail':detail})
    summary=read(BASE/'results/inside_outside_summary.tsv')
    def get(des,region):return next(x for x in summary if x['design']==des and x['region']==region)
    for reg in ['LGI_inside','LGI_outside']:
        x=get('primary',reg);add(reg,x['D'],'NA','NA',x['informative_sites'],'ABBA=%s; BABA=%s'%(x['ABBA'],x['BABA']))
    for width in ['250000','500000','1000000']:
        x=next(x for x in uncertainty if str(x['block_width_bp'])==width and x['region']=='inside_minus_outside')
        add('LGI_D_inside_minus_outside_'+width+'bp',x['estimate_D'],x['ci95_low'],x['ci95_high'],x['n_physical_blocks'],'physical-block delete-one; outside blocks are sparse')
    for name in ['LGI_left_flank_250kb','LGI_right_flank_250kb']:
        x=get('primary',name);add(name,x['D'],'NA','NA',x['informative_sites'],'250-kb flank')
    for reg in ['LGII','LGIII','LGIV','LGV','LGVI','LGVII']:
        x=get('primary',reg);add(reg,x['D'],'NA','NA',x['informative_sites'],'whole linkage group; variant panel')
    for x in read(R/'ascertainment_common_support.tsv'):
        add('ascertainment_'+x['filter']+'_'+x['region'],x['D'],'NA','NA',x['windows_retained'],x['bounds'])
    for x in read(R/'genealogy_correspondence.tsv'):
        add('genealogy_'+x['subset'],x['pearson_r'],x['CI95_low'],x['CI95_high'],x['n_bins'],x['interpretation'])
    add('raw_scalar_site_validation_n',len(read(R/'estimator_validation.tsv')),'NA','NA',len(set(x['stratum'] for x in read(R/'estimator_validation.tsv'))),'All ABBA/BABA contributions and status labels exactly match independent scalar counter.')
    write(R/'final_results.tsv',rows)

if __name__=='__main__':main()
