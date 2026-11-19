"""Frequency-based modified f_d statistic following Dsuite Dinvestigate."""
from __future__ import annotations
from dataclasses import dataclass
import math
from collections import defaultdict
from typing import Iterable, Optional, Sequence

@dataclass(frozen=True)
class FdMComponents:
    numerator: float
    denominator: float
    abba: float
    baba: float
    included: bool

@dataclass(frozen=True)
class FdMEstimate:
    numerator: float
    denominator: float
    fdm: Optional[float]
    n_sites: int
    n_informative_sites: int

def _freq(value: Optional[float], name: str) -> Optional[float]:
    if value is None: return None
    x=float(value)
    if not math.isfinite(x) or x < 0 or x > 1: raise ValueError(f"{name} must be finite and in [0,1]")
    return x

def derived_frequency_from_alt(alt_frequency: Optional[float], outgroup_alt_frequency: Optional[float], *, require_fixed: bool=True) -> Optional[float]:
    """Orient an ALT frequency using the outgroup state; reject polymorphic O by default."""
    af=_freq(alt_frequency,'alt_frequency');of=_freq(outgroup_alt_frequency,'outgroup_alt_frequency')
    if af is None or of is None:return None
    if require_fixed and of not in (0.0,1.0):return None
    # For a nonfixed outgroup and require_fixed=False, use majority-rule ancestry,
    # matching Dsuite's polarization convention; p_O is handled separately.
    ancestral_alt=(of>=0.5)
    return 1.0-af if ancestral_alt else af

def fdm_site(p1: Optional[float], p2: Optional[float], p3: Optional[float], p_outgroup: Optional[float]=0.0) -> FdMComponents:
    """Site terms using derived frequencies; fixed ancestral O has p_outgroup=0."""
    p1=_freq(p1,'p1');p2=_freq(p2,'p2');p3=_freq(p3,'p3');po=_freq(p_outgroup,'p_outgroup')
    if None in (p1,p2,p3,po): return FdMComponents(0.,0.,0.,0.,False)
    if (p1==p2==p3==0.) or (p1==p2==p3==1.): return FdMComponents(0.,0.,0.,0.,False)
    w=1.-po
    abba=(1.-p1)*p2*p3*w; baba=p1*(1.-p2)*p3*w; num=abba-baba
    # Equivalent to the four branches in Dsuite D.cpp / Dinvestigate.
    if p1<=p2:
        if p2>p3: den=((1.-p1)*p2*p2-p1*(1.-p2)*p2)*w
        else: den=((1.-p1)*p3*p3-p1*(1.-p3)*p3)*w
    elif p1>p3:
        den=-(((1.-p1)*p2*p1-p1*(1.-p2)*p1)*w)
    else:
        den=-(((1.-p3)*p2*p3-p3*(1.-p2)*p3)*w)
    return FdMComponents(num,den,abba,baba,True)

def pool_fdm(frequencies: Iterable[Sequence[Optional[float]]]) -> FdMEstimate:
    """Pool site terms before division; rows have p1,p2,p3[,p_outgroup]."""
    n=ni=0;num=den=0.
    for row in frequencies:
        if len(row) not in (3,4): raise ValueError('rows require p1,p2,p3[,p_outgroup]')
        c=fdm_site(*row)
        if not c.included:continue
        n+=1;num+=c.numerator;den+=c.denominator;ni+=int(c.abba+c.baba>0)
    return FdMEstimate(num,den,num/den if den else None,n,ni)

def window_fdm(records: Iterable[tuple[int,float,float,float,float]],window_bp: int,block_bp: int=500_000):
    """Aggregate 1-based (position,p1,p2,p3,pO) rows into physical windows/blocks."""
    if window_bp<=0 or block_bp<=0:raise ValueError('window/block sizes must be positive')
    sums=defaultdict(lambda:[0.,0.,0,0]);blocks=defaultdict(lambda:[0.,0.])
    for pos,p1,p2,p3,po in records:
        pos=int(pos)
        if pos<1:raise ValueError('positions must be positive 1-based coordinates')
        c=fdm_site(p1,p2,p3,po)
        if not c.included:continue
        start=((pos-1)//window_bp)*window_bp+1;q=sums[start]
        q[0]+=c.numerator;q[1]+=c.denominator;q[2]+=1;q[3]+=int(c.abba+c.baba>0)
        b=((pos-1)//block_bp)*block_bp+1;blocks[b][0]+=c.numerator;blocks[b][1]+=c.denominator
    rows=[]
    for start,(num,den,n,ni) in sorted(sums.items()):
        rows.append(dict(start=start,end=start+window_bp-1,numerator=num,denominator=den,fdm=num/den if den else None,n_sites=n,n_informative_sites=ni))
    return rows,dict(blocks)

def jackknife_ratio(blocks: dict[int,Sequence[float]],z: float=1.96):
    """Delete-one physical-block jackknife for ratio of pooled components."""
    vals={k:(float(v[0]),float(v[1])) for k,v in blocks.items() if len(v)>=2 and v[1]!=0}
    n=len(vals);num=sum(x[0] for x in vals.values());den=sum(x[1] for x in vals.values());est=num/den if den else None
    if est is None or n<2:return dict(fdm=est,se=None,ci_low=None,ci_high=None,n_blocks=n)
    leave=[(num-a)/(den-b) for a,b in vals.values() if den-b!=0]
    if len(leave)<2:return dict(fdm=est,se=None,ci_low=None,ci_high=None,n_blocks=n)
    m=sum(leave)/len(leave);se=math.sqrt((len(leave)-1)/len(leave)*sum((x-m)**2 for x in leave))
    return dict(fdm=est,se=se,ci_low=est-z*se,ci_high=est+z*se,n_blocks=n)
