#!/usr/bin/env python3
"""Cross-check cod f_dM scan against frozen D windows and render summaries."""
import csv,json,hashlib
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[5]
BASE=ROOT/'empirical/atlantic_cod/analysis/joint_genomic_tracks';RES=BASE/'results';FIG=BASE/'figures'
def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 new=read(RES/'atlantic_cod_fdm_windows.tsv');old=read(ROOT/'empirical/atlantic_cod/analysis/patterson_d/results/window_patterson_d.tsv')
 oldmap={(r['design_id'],r['window_id']):r for r in old};cross=[]
 for r in new:
  wid=f"{r['lg']}_{int(r['window_start']):09d}_{int(r['window_end']):09d}";o=oldmap.get((r['design_id'],wid))
  if o is None:continue
  if r['D_Patterson']=='NA' or o['D_site']=='NA':
   cross.append(dict(lg=r['lg'],design_id=r['design_id'],window_id=wid,D_from_fdm_scan=r['D_Patterson'],D_original=o['D_site'],absolute_difference='NA',status='BOTH_UNDEFINED' if r['D_Patterson']==o['D_site'] else 'MISMATCH_UNDEFINED'))
  else:
   x=float(r['D_Patterson']);y=float(o['D_site'])
   cross.append(dict(lg=r['lg'],design_id=r['design_id'],window_id=wid,D_from_fdm_scan=x,D_original=y,absolute_difference=abs(x-y),status='PASS' if abs(x-y)<1e-8 else 'FAIL'))
 write(RES/'atlantic_cod_fdm_D_crosscheck.tsv',cross)
 fails=sum(r['status'] not in ('PASS','BOTH_UNDEFINED') for r in cross)
 summaries=read(RES/'atlantic_cod_fdm_regions.tsv');window_rows=new
 # Separate track figures use one statistic and one y-axis per panel.
 FIG.mkdir(parents=True,exist_ok=True);plt.rcParams.update({'font.size':9,'axes.spines.top':False})
 colors={'D':'#31688e','fdm':'#d1495b'}
 region={r['lg']:(int(r['start']),int(r['end'])) for r in read(ROOT/'data/atlantic_cod/metadata/region_manifest.tsv')}
 for stat_key,stat_label,color,stem in (
     ('D_Patterson',"Patterson's D",colors['D'],'atlantic_cod_patterson_D_tracks'),
     ('f_dM','$f_{dM}$',colors['fdm'],'atlantic_cod_fdm_tracks')):
  fig,axs=plt.subplots(4,1,figsize=(11,9),sharex=False,constrained_layout=True)
  values=[float(r[stat_key]) for r in window_rows if r[stat_key]!='NA']
  ylim=max(abs(v) for v in values)*1.05
  for ax,lg in zip(axs,['LG01','LG02','LG07','LG12']):
   rr=sorted([r for r in window_rows if r['lg']==lg],key=lambda r:int(r['window_start']))
   x=[(int(r['window_start'])+int(r['window_end']))/2/1e6 for r in rr]
   y=[float(r[stat_key]) if r[stat_key]!='NA' else float('nan') for r in rr]
   lo,hi=region[lg]
   ax.axvspan(lo/1e6,hi/1e6,color='#f1c75b',alpha=.22)
   ax.scatter(x,y,s=10,color=color)
   ax.set_ylim(-ylim,ylim)
   ax.axhline(0,color='0.45',lw=.7)
   ax.set_ylabel(f'{lg}\n{stat_label}')
   ax.set_xlabel('Linkage-group coordinate (Mb)')
  fig.suptitle(f'Atlantic cod: {stat_label} across 250-kb windows (shading = frozen inversion interval)')
  for ext,kwargs in (('png',{'dpi':300}),('pdf',{}),('svg',{})):
   fig.savefig(FIG/f'{stem}.{ext}',**kwargs)
  plt.close(fig)
 # Inside/outside region estimates, physical-block jackknife.
 fig,ax=plt.subplots(figsize=(9,4.8),constrained_layout=True)
 lgs=['LG01','LG02','LG07','LG12'];cols=['#31688e','#d1495b','#55a868','#8172b3']
 for i,lg in enumerate(lgs):
  for j,reg in enumerate(['inside','outside']):
   r=next(x for x in summaries if x['lg']==lg and x['region']==reg)
   val=float(r['fdm']);lo=float(r['ci_low']);hi=float(r['ci_high']);y=i+(j-.5)*.28
   ax.errorbar(val,y,xerr=[[val-lo],[hi-val]],fmt='o',capsize=3,color=cols[i],label=lg if j==0 else None)
 ax.axvline(0,color='black',lw=.7);ax.set_yticks(range(4),lgs);ax.invert_yaxis();ax.set_xlabel('$f_{dM}$ (pooled numerator / denominator; 500-kb physical-block jackknife 95% CI)');ax.set_title('Frozen interval (upper) and same-linkage-group outside windows (lower)');ax.legend(frameon=False,ncol=4)
 fig.savefig(FIG/'atlantic_cod_fdm_regions.png',dpi=300);fig.savefig(FIG/'atlantic_cod_fdm_regions.pdf');fig.savefig(FIG/'atlantic_cod_fdm_regions.svg');plt.close(fig)
 # Machine-readable provenance from exact source hashes.
 paths=[ROOT/'data/atlantic_cod/raw/zenodo/local_downloads/gadus_morhua_outside_supergenes_snapp.vcf.gz']
 for lg in lgs:paths.append(ROOT/f'data/atlantic_cod/raw/zenodo/local_downloads/gadus_morhua_supergenes_{lg.lower()}_snapp.vcf.gz')
 provenance={'status':'completed','estimator':'Dsuite Dinvestigate f_dM site terms; fixed-outgroup polarization','window_bp':250000,'jackknife_block_bp':500000,'sample_design':'empirical/atlantic_cod/analysis/patterson_d/config/frozen_quartets.tsv','window_inventory':'data/atlantic_cod/processed/cod_window_trees.tsv','sources':[{'path':str(p.relative_to(ROOT)),'sha256':sha(p)} for p in paths], 'crosscheck_rows':len(cross),'crosscheck_failures':fails,'max_abs_D_difference':max((float(r['absolute_difference']) for r in cross if r['absolute_difference']!='NA'),default=None),'dsuite_executable_validation':'not_run; Dsuite CLI is not installed. Piecewise terms tested against official Dsuite/D.cpp source transcription.'}
 (RES/'atlantic_cod_fdm_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
 if fails:raise SystemExit(f'{fails} Patterson D cross-checks failed')
 print('crosscheck',len(cross),'PASS; max abs difference',provenance['max_abs_D_difference'])
if __name__=='__main__':main()
