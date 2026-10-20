"""Reusable spatial interval segmentation for Stage 13.

The detector works on the observed window order.  Missing windows are removed
from the fitted sequence, while their physical gaps remain visible in output
coordinates.  For a fixed number of segments, the dynamic program is exact
for the stated least-squares objective.  A one-interval internal model has
three free-mean segments; endpoint models have two.  Multi-interval models use
2*K+1 free-mean segments and a fixed penalty per additional segment.  This is
an explicit piecewise-constant model rather than a biological claim that both
flanks have identical means.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np

@dataclass
class Segmentation:
    status: str
    intervals: list[tuple[int, int]]
    score: float
    raw_improvement: float
    penalty: float
    n_segments: int
    feature_contributions: dict[str, float]
    scale: dict[str, float]

def robust_standardize(matrix, names):
    x=np.asarray(matrix,dtype=float).copy(); scale={}
    for j,name in enumerate(names):
        col=x[:,j]; med=float(np.nanmedian(col)); mad=float(np.nanmedian(np.abs(col-med))); s=1.4826*mad
        if not math.isfinite(s) or s<=1e-12: s=float(np.nanstd(col)) or 1.0
        x[:,j]=(col-med)/s; scale[name]={'median':med,'scale':s}
    return x,scale

def _sse_prefix(x):
    # One prefix pair per feature, so range SSE is O(1).
    return np.vstack([np.zeros((1,x.shape[1])),np.cumsum(x,axis=0)]), np.vstack([np.zeros((1,x.shape[1])),np.cumsum(x*x,axis=0)])
def _cost(sum_p,sq_p,a,b):
    n=b-a; sm=sum_p[b]-sum_p[a]; ss=sq_p[b]-sq_p[a]
    return float(np.sum(ss-sm*sm/n)) if n else float('inf')
def _dp(x,nseg,min_size):
    n=x.shape[0]; sp,qp=_sse_prefix(x); inf=float('inf'); dp=np.full((nseg+1,n+1),inf); prev=np.full((nseg+1,n+1),-1,int); dp[0,0]=0
    for k in range(1,nseg+1):
        for end in range(k*min_size,n+1):
            lo=(k-1)*min_size; hi=end-min_size
            cuts=np.arange(lo,hi+1,dtype=int); nn=end-cuts
            sm=sp[end]-sp[cuts]; ss=qp[end]-qp[cuts]
            vals=dp[k-1,cuts]+np.sum(ss-sm*sm/nn[:,None],axis=1)
            j=int(np.argmin(vals)); dp[k,end]=float(vals[j]);prev[k,end]=int(cuts[j])
    if not math.isfinite(dp[nseg,n]):return float('inf'),[]
    cuts=[];end=n
    for k in range(nseg,0,-1): cuts.append((prev[k,end],end));end=prev[k,end]
    return float(dp[nseg,n]),list(reversed(cuts))
def _intervals_for(cuts, mode):
    if mode=='internal': return [(cuts[i][0],cuts[i][1]) for i in range(1,len(cuts)-1,2)]
    if mode=='left': return [(cuts[0][0],cuts[0][1])]
    if mode=='right': return [(cuts[-1][0],cuts[-1][1])]
    return []
def segment_track(values, names, min_size=3, max_intervals=4, penalty=None):
    x,scale=robust_standardize(values,names); n=len(x)
    base=_cost(*_sse_prefix(x),0,n); pen=float(penalty if penalty is not None else 2*math.log(max(n,2))*x.shape[1]); candidates=[Segmentation('no_detection',[],0,0,0,1,{},scale)]
    for k in range(1,max_intervals+1):
        cost,cuts=_dp(x,2*k+1,min_size)
        if math.isfinite(cost):
            raw=base-cost; candidates.append(Segmentation('internal',_intervals_for(cuts,'internal'),raw-k*pen,raw,k*pen,2*k+1,{},scale))
    if n>=2*min_size:
        for mode in ('left','right'):
            cost,cuts=_dp(x,2,min_size)
            if math.isfinite(cost):
                raw=base-cost;candidates.append(Segmentation(mode,_intervals_for(cuts,mode),raw-pen,raw,pen,2,{},scale))
    best=max(candidates,key=lambda z:z.score)
    if best.score<=0:return Segmentation('no_detection',[],best.score,best.raw_improvement,best.penalty,1,{},scale)
    if best.intervals:
        contrib={}
        for j,name in enumerate(names):
            col=x[:,j]; mu=float(np.mean(col)); inside=np.zeros(n,bool)
            for a,b in best.intervals:inside[a:b]=True
            contrib[name]=float(np.sum((col[~inside]-mu)**2)+np.sum((col[inside]-col[inside].mean())**2))
        best.feature_contributions=contrib
    return best
def intervals_to_coordinates(seg, rows):
    out=[]
    for a,b in seg.intervals:
        out.append({'start':rows[a]['start'],'end':rows[b-1]['end'],'start_is_censored':int(a==0),'end_is_censored':int(b==len(rows))})
    return out
