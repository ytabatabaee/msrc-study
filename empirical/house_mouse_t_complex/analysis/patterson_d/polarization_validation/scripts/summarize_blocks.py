#!/usr/bin/env python3
"""Block-bootstrap corrected full-chr17 inversion/segment contrasts."""
import csv,random,statistics,math
from pathlib import Path
from collections import defaultdict
BASE=Path(__file__).resolve().parents[1];OUT=BASE/'results';rows=list(csv.DictReader(open(OUT/'corrected_window_D.tsv'),delimiter='\t'))
DESIGNS=sorted({x['design'] for x in rows});INV=[('In(17)1',5747331,6957919),('In(17)2',7435308,13016510),('In(17)3',14068229,18051507),('In(17)4',19401000,40049964)]
SEGS=INV+[('left_flank',5000000,5747331),('collinear_gap_1',6957919,7435308),('collinear_gap_2',13016510,14068229),('collinear_gap_3',18051507,19401000),('right_flank',40049964,42050000)]
def d(a,b):return (a-b)/(a+b) if a+b else float('nan')
def pooled(rr):return sum(float(x['ABBA']) for x in rr),sum(float(x['BABA']) for x in rr)
def boot(ins,outs,des,reg,bp,n=1000):
 def aggregate(rr):
  g=defaultdict(lambda:[0.,0.])
  for x in rr:
   k=int(x['start'])//bp;g[k][0]+=float(x['ABBA']);g[k][1]+=float(x['BABA'])
  return g
 a,b=aggregate(ins),aggregate(outs);ka,kb=sorted(a),sorted(b)
 if len(ka)<2 or len(kb)<2:return float('nan'),float('nan'),float('nan'),len(ka),len(kb)
 rng=random.Random(sum((i+1)*ord(c) for i,c in enumerate(f'{des}|{reg}|{bp}')));vals=[]
 for _ in range(n):
  sa=[rng.choice(ka) for _ in ka];sb=[rng.choice(kb) for _ in kb]
  da=d(sum(a[k][0] for k in sa),sum(a[k][1] for k in sa));db=d(sum(b[k][0] for k in sb),sum(b[k][1] for k in sb));
  if math.isfinite(da-db):vals.append(da-db)
 vals.sort();point=d(*pooled(ins))-d(*pooled(outs));return point,vals[int(.025*len(vals))],vals[min(len(vals)-1,int(.975*len(vals)))],len(ka),len(kb)
def main():
 out=[]
 for des in DESIGNS:
  allr=[x for x in rows if x['design']==des]; union=[x for x in allr if any(int(x['start'])<b and int(x['end'])>a for _,a,b in INV)];
  for name,a,b in [('combined_union',0,0)]+SEGS:
   inside=union if name=='combined_union' else [x for x in allr if int(x['start'])<b and int(x['end'])>a]
   starts={(x['start']) for x in inside};outside=[x for x in allr if x['start'] not in starts]
   for bp in (250000,500000,1000000):
    est,lo,hi,ni,no=boot(inside,outside,des,name,bp)
    out.append({'design':des,'region':name,'statistic':'D_site_inside_minus_full_chromosome_complement','block_bp':bp,'point_difference':est,'block_bootstrap_low':lo,'block_bootstrap_high':hi,'inside_blocks':ni,'complement_blocks':no,'inside_windows':len(inside),'complement_windows':len(outside)})
 with open(OUT/'corrected_inside_outside_block_uncertainty.tsv','w') as f:
  w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter='\t');w.writeheader();w.writerows(out)
 print('rows',len(out))
if __name__=='__main__':main()
