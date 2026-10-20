#!/usr/bin/env python3
"""Boundary-aware single- and multi-interval validation for mouse tracks."""
from pathlib import Path
import csv, json, hashlib, random
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]; OUT=ROOT/'empirical/house_mouse_t_complex/results/stage4_distance'; FIG=ROOT/'empirical/house_mouse_t_complex/figures'
def read(p):
 with open(p) as h:return list(csv.DictReader(h,delimiter='\t'))
def scan(rows,f,min_n=5,penalty=.05):
 z=sorted([r for r in rows if r.get(f,'') not in ('','nan')],key=lambda r:float(r['midpoint_bp']));
 # A deterministic 25-kb candidate grid makes the scan reproducible while
 # preserving the 5-kb observations and avoids pretending breakpoints are
 # resolved more finely than the broad local-tree signal supports.
 z=z[::5] if len(z)>1000 else z
 y=np.array([float(r[f]) for r in z]); base=((y-y.mean())**2).sum(); best=(-1,None,None,0)
 sy=np.r_[0,np.cumsum(y)]; sy2=np.r_[0,np.cumsum(y*y)]; n=len(z)
 def sse(a,b):
  nn=b-a; sm=sy[b]-sy[a]; ss=sy2[b]-sy2[a]; return ss-sm*sm/nn if nn else 0
 for i in range(len(z)):
  for j in range(i+min_n-1,len(z)-1):
   s=sse(0,i)+sse(i,j+1)+sse(j+1,n); score=base-s-penalty*base
   if score>best[0]:best=(score,i,j,s)
 if best[1] is None:return {'feature':f,'score':0,'start_bp':'','end_bp':'','n_inside':0,'status':'no_detection','penalty':penalty}
 return {'feature':f,'score':best[0],'start_bp':z[best[1]]['start_bp'],'end_bp':z[best[2]]['end_bp'],'n_inside':best[2]-best[1]+1,'status':'detected' if best[0]>0 else 'no_detection','penalty':penalty}
def iou(a,b,c,d):return max(0,min(b,d)-max(a,c)+1)/(max(b,d)-min(a,c)+1)
def main():
 rows=read(OUT/'stage4_window_distance.tsv'); methods=[scan(rows,f) for f in ('A_distance','A_relative','q_t_alt','mean_pair_distance_no_spretus')]
 intervals=read(ROOT/'empirical/house_mouse_t_complex/results/t_complex_inversion_intervals.tsv'); preds=[]
 for m in methods:
  for iv in intervals: preds.append({**m,'inversion_name':iv['inversion_name'],'IoU':iou(int(m['start_bp']),int(m['end_bp']),int(iv['start_bp']),int(iv['end_bp']))})
 fields=list(preds[0]);
 with (OUT/'stage4_boundary_predictions.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(preds)
 # Fixed one-interval scans are supplemented by an explicitly labelled multi-interval greedy scan.
 multi=[]; remaining=rows[:]
 for k in range(1,5):
  p=scan(remaining,'A_distance');multi.append({**p,'n_intervals':k}); a=float(p['start_bp']);b=float(p['end_bp']);remaining=[r for r in remaining if not(a<=float(r['midpoint_bp'])<=b)]
 with (OUT/'stage4_multi_interval_predictions.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=list(multi[0]),delimiter='\t');w.writeheader();w.writerows(multi)
 # Coordinate-aware interval summaries are descriptive alignment tests. They
 # do not test whether any change point exists.
 null=[]; rng=random.Random(20261019)
 for iv in intervals:
  a,b=int(iv['start_bp']),int(iv['end_bp']); inside=[float(r['A_distance']) for r in rows if a<=float(r['midpoint_bp'])<=b]; outside=[float(r['A_distance']) for r in rows if float(r['end_bp'])<a or float(r['start_bp'])>b]
  obs=np.mean(inside)-np.mean(outside); vals=np.array([float(r['A_distance']) for r in rows]); k=len(inside); sims=[]
  for _ in range(999):
   ix=rng.sample(range(len(vals)),k); mask=np.zeros(len(vals),bool);mask[ix]=1;sims.append(float(np.mean(vals[mask])-np.mean(vals[~mask])))
  null.append({'inversion_name':iv['inversion_name'],'inside_minus_outside_A_distance':obs,'two_sided_spatial_p':(1+sum(abs(v)>=abs(obs) for v in sims))/1000,'n_inside':k,'n_outside':len(outside)})
 with (OUT/'stage4_spatial_null.tsv').open('w',newline='') as h:w=csv.DictWriter(h,fieldnames=list(null[0]),delimiter='\t');w.writeheader();w.writerows(null)
 x=np.array([float(r['midpoint_bp'])/1e6 for r in rows]);fig,ax=plt.subplots(2,1,figsize=(12,7),sharex=True,constrained_layout=True);ax[0].plot(x,[float(r['q_t_alt']) if r['q_t_alt'] not in ('','nan') else np.nan for r in rows],'.',ms=1);ax[0].set_ylabel('topology signal');ax[1].plot(x,[float(r['A_distance']) if r['A_distance'] not in ('','nan') else np.nan for r in rows],'.',ms=1);ax[1].set_ylabel('A_distance');ax[1].set_xlabel('chr17 position (Mb)')
 for iv in intervals:
  for a in ax:a.axvspan(int(iv['start_bp'])/1e6,int(iv['end_bp'])/1e6,alpha=.12)
 fig.savefig(FIG/'house_mouse_stage4_boundary_distance.pdf');fig.savefig(FIG/'house_mouse_stage4_boundary_distance.png',dpi=220);plt.close(fig);(OUT/'boundary_report.md').write_text('# House mouse boundary validation\n\nThe single-interval and explicitly labelled greedy multi-interval scans use fixed least-squares penalties and do not use inversion coordinates during fitting. Frozen intervals are applied only for post-inference IoU summaries.\n'); (OUT/'report.md').write_text((OUT/'boundary_report.md').read_text()); (OUT/'methods_text.md').write_text((OUT/'boundary_report.md').read_text()); (OUT/'results_text.md').write_text((OUT/'boundary_report.md').read_text()); (OUT/'manifest.json').write_text(json.dumps({'seed':20261019,'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ROOT/'data/house_mouse_t_complex/processed/house_mouse_t_complex_ml_5kb.tre',ROOT/'data/house_mouse_t_complex/processed/house_mouse_t_complex_ml_5kb_metadata.tsv',ROOT/'data/house_mouse_t_complex/metadata/tip_mapping.tsv',ROOT/'empirical/house_mouse_t_complex/results/t_complex_inversion_intervals.tsv']},'n_windows':len(rows)},indent=2)+'\n')
if __name__=='__main__':main()
