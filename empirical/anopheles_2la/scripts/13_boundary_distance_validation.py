#!/usr/bin/env python3
"""Boundary-aware validation of Anopheles distance and topology tracks."""
from pathlib import Path
import csv, json, math
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[3]; OUT=ROOT/'empirical/anopheles_2la/results/stage12_distance'; FIG=ROOT/'empirical/anopheles_2la/figures'
def read(p):
 with open(p) as h:return list(csv.DictReader(h,delimiter='\t'))
def best(rows,feature,min_n=3,penalty=0.05):
 z=[r for r in rows if r.get(feature,'') not in ('','nan','NA')]; z=sorted(z,key=lambda r:float(r['midpoint']))
 if len(z)<min_n*2: return {'feature':feature,'score':0,'start':'','end':'','n_inside':0,'status':'no_detection','penalty':penalty}
 y=np.array([float(r[feature]) for r in z]); base=((y-y.mean())**2).sum(); best=(0,0,0,base)
 sy=np.r_[0,np.cumsum(y)]; sy2=np.r_[0,np.cumsum(y*y)]; n=len(z)
 def sse(a,b):
  nn=b-a; sm=sy[b]-sy[a]; ss=sy2[b]-sy2[a]; return ss-sm*sm/nn if nn else 0
 for i in range(len(z)):
  for j in range(i+min_n-1,len(z)-1):
   s=sse(0,i)+sse(i,j+1)+sse(j+1,n); score=base-s-penalty*base
   if score>best[0]:best=(score,i,j,s)
 return {'feature':feature,'score':best[0],'start':z[best[1]]['start'],'end':z[best[2]]['end'],'n_inside':best[2]-best[1]+1,'status':'detected' if best[0]>0 else 'no_detection','penalty':penalty}
def main():
 rows=read(OUT/'stage12_window_distance.tsv'); methods=[]
 for f in ('A_distance','A_relative','D','mean_pair_distance'):
  methods.append(best(rows,f))
 out=OUT/'stage12_boundary_predictions.tsv'; fields=list(methods[0]);
 with out.open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(methods)
 # prespecified frozen interval is loaded only for the post-inference alignment summary.
 iv=(20524058,42165532); acc=[]
 for r in methods:
  a,b=int(r['start']),int(r['end']); inter=max(0,min(b,iv[1])-max(a,iv[0])+1); union=max(b,iv[1])-min(a,iv[0])+1; acc.append({'feature':r['feature'],'IoU':inter/union,'overlap_bp':inter})
 with (OUT/'stage12_boundary_accuracy.tsv').open('w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=list(acc[0]),delimiter='\t');w.writeheader();w.writerows(acc)
 x=np.array([float(r['midpoint'])/1e6 for r in rows]); fig,ax=plt.subplots(2,1,figsize=(12,6),sharex=True,constrained_layout=True); ax[0].plot(x,[float(r['D']) if r['D'] not in ('','nan') else np.nan for r in rows],'.',ms=2);ax[0].set_ylabel('D(w)'); ax[1].plot(x,[float(r['A_distance']) if r['A_distance'] not in ('','nan') else np.nan for r in rows],'.',ms=2);ax[1].set_ylabel('A_distance');ax[1].set_xlabel('2L position (Mb)');
 ax[0].axvspan(iv[0]/1e6,iv[1]/1e6,color='#8c6bb1',alpha=.15);ax[1].axvspan(iv[0]/1e6,iv[1]/1e6,color='#8c6bb1',alpha=.15);fig.savefig(FIG/'anopheles_2la_stage12_boundary_distance.pdf');fig.savefig(FIG/'anopheles_2la_stage12_boundary_distance.png',dpi=220);plt.close(fig)
 (OUT/'boundary_report.md').write_text('# Anopheles boundary validation\n\nThe interval scan uses a fixed least-squares improvement penalty and a minimum of three observed windows. The frozen 2La coordinates are used only after prediction for alignment summaries. The distance features are substitution-distance summaries, not coalescence times.\n')
if __name__=='__main__':main()
