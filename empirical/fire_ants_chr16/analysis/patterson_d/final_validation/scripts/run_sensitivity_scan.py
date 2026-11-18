#!/usr/bin/env python3
"""Stream the verified VCF once for frozen-Q1/Q2 DP/callability sensitivities."""
import csv,gzip,math,random
from collections import defaultdict
from pathlib import Path
from bisect import bisect_right
ROOT=Path(__file__).resolve().parents[6]; BASE=ROOT/'empirical/fire_ants_chr16/analysis/patterson_d'; OUT=BASE/'final_validation/results'
VCF=ROOT/'data/fire_ants_chr16/raw/genotypes/gt.vcf.gz'; MAP=BASE/'results/window_gene_coordinate_mapping.tsv'; WIN=ROOT/'data/fire_ants_chr16/processed/stage1_window_index.tsv'; MEM=BASE/'config/frozen_sample_membership.tsv'
DP=(8,10,12,15); FRAC=(.5,.6,.7,.8); Q=('Q1_primary_Argentina','Q2_primary_Argentina')
def read(p):
 with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 if rows:
  with open(p,'w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def parse(s,ix):
 z=s.split(':');gt=z[ix['GT']] if 'GT' in ix and ix['GT']<len(z) else '.'
 try:depth=int(z[ix['DP']]) if 'DP' in ix and ix['DP']<len(z) else 0
 except ValueError:depth=0
 if '.' in gt:return None,depth
 if '/' in gt or '|' in gt:
  a=gt.replace('|','/').split('/')
  if len(a)!=2 or a[0]!=a[1]:return None,depth
  gt=a[0]
 return (int(gt) if gt in ('0','1') else None),depth
def acc():return dict(records=0,depth_all_called=0,call_fraction_pass=0,outgroup_fixed=0,polarizable=0,informative_sites=0,ABBA=0.,BABA=0.,derived_REF=0,derived_ALT=0)
def main():
 OUT.mkdir(parents=True,exist_ok=True); wr={int(x['window_index']):x for x in read(WIN)}
 iv=defaultdict(list)
 for x in read(MAP):
  if x['vcf_contig_present']=='1':iv[x['vcf_contig']].append((int(x['vcf_gene_start_1based']),int(x['vcf_gene_end_1based_inclusive']),int(x['window_index'])))
 for c in iv:iv[c].sort()
 starts={c:[x[0] for x in v] for c,v in iv.items()}
 maxends={}
 for c,v in iv.items():
  a=[];m=0
  for _,hi,_ in v:m=max(m,hi);a.append(m)
  maxends[c]=a
 mm=defaultdict(lambda:defaultdict(list))
 for x in read(MEM):
  if x['quartet_id'] in Q:mm[x['quartet_id']][x['role']].append(x['sample_id'])
 samples=sorted({s for g in mm.values() for ids in g.values() for s in ids}); bywin=defaultdict(acc); regions=defaultdict(acc)
 with gzip.open(VCF,'rt') as f:
  h=next(x for x in f if x.startswith('#CHROM')).rstrip().split('\t'); ci={s:i for i,s in enumerate(h[9:],9)}
  if set(samples)-set(ci):raise RuntimeError('frozen sample absent from VCF')
  cols=[ci[s] for s in samples]; si={s:i for i,s in enumerate(samples)}
  roles={q:{r:[si[s] for s in ids] for r,ids in groups.items()} for q,groups in mm.items()}
  for line in f:
   if line.startswith('#'):continue
   i=line.find('\t');j=line.find('\t',i+1);chrom=line[:i]
   if chrom not in iv:continue
   pos=int(line[i+1:j]); hits=[]
   k=bisect_right(starts[chrom],pos)
   j=k-1
   while j>=0 and maxends[chrom][j]>=pos:
    lo,hi,wi=iv[chrom][j]
    if hi>=pos:hits.append(wi)
    j-=1
   if not hits:continue
   z=line.rstrip('\n').split('\t');ref,alt=z[3],z[4]
   if len(ref)!=1 or len(alt)!=1 or alt=='.' or ',' in alt:continue
   ix={x:i for i,x in enumerate(z[8].split(':'))}; raw=[parse(z[k],ix) for k in cols]
   for q,gr in roles.items():
    sizes=[len(gr[r]) for r in ('P1','P2','P3')]; no=len(gr['O'])
    for d in DP:
     gt=[a if a is not None and depth>=d else None for a,depth in raw]
     gs=[[gt[k] for k in gr[r] if gt[k] is not None] for r in ('P1','P2','P3')]; oo=[gt[k] for k in gr['O'] if gt[k] is not None]
     for fcall in FRAC:
      passed=all(len(a)>=math.ceil(n*fcall) for a,n in zip(gs,sizes)); fixed=len(oo)>=min(4,no) and len(set(oo))==1
      for wi in set(hits):
       rg=wr[wi]['region']
       for key,c in (((q,wi,d,fcall),bywin),((q,rg,d,fcall),regions)):
        a=c[key];a['records']+=1;a['depth_all_called']+=int(all(len(g)==n for g,n in zip(gs,sizes)) and len(oo)==no);a['call_fraction_pass']+=int(passed);a['outgroup_fixed']+=int(passed and fixed)
        if passed and fixed:
         anc=oo[0];a['polarizable']+=1;a['derived_REF' if anc==1 else 'derived_ALT']+=1
         p=[sum((v if anc==0 else 1-v) for v in g)/len(g) for g in gs];ab=(1-p[0])*p[1]*p[2];ba=p[0]*(1-p[1])*p[2]
         a['ABBA']+=ab;a['BABA']+=ba;a['informative_sites']+=int(ab+ba>0)
 out=[]
 for (q,reg,d,f),a in sorted(regions.items()):
  out.append(dict(quartet_id=q,region=reg,min_DP=d,min_ingroup_call_fraction=f,**a,D=(a['ABBA']-a['BABA'])/(a['ABBA']+a['BABA']) if a['ABBA']+a['BABA'] else 'NA'))
 write(OUT/'filter_attrition.tsv',out)
 out=[]
 for q in Q:
  for d in DP:
   for f in FRAC:
    for reg in ('chr16_supergene','chr16A','chr16B'):
     cells=[v for (qq,wi,dd,ff),v in bywin.items() if qq==q and dd==d and ff==f and wr[wi]['region']==reg];ab=sum(x['ABBA'] for x in cells);ba=sum(x['BABA'] for x in cells)
     out.append(dict(quartet_id=q,min_DP=d,min_ingroup_call_fraction=f,region=reg,n_windows=sum(wr[i]['region']==reg for i in wr),ABBA=ab,BABA=ba,D=(ab-ba)/(ab+ba) if ab+ba else 'NA',informative_sites=sum(x['informative_sites'] for x in cells),uncertainty='not estimated in this scan'))
 # Save window-level pooled contributions so physical-block uncertainty can be
 # independently reproduced for each threshold combination.
 wout=[]
 for (q,wi,d,f),a in sorted(bywin.items()):
  den=a['ABBA']+a['BABA'];wout.append(dict(quartet_id=q,window_index=wi,region=wr[wi]['region'],mid=wr[wi]['mid'],min_DP=d,min_ingroup_call_fraction=f,**a,D=(a['ABBA']-a['BABA'])/den if den else 'NA'))
 write(OUT/'filter_window_sensitivity.tsv',wout)
 for r in out:
  r['n_informative_physical_blocks']='NA'
  if r['region']=='chr16_supergene':
   other=[x for x in out if x['quartet_id']==r['quartet_id'] and x['min_DP']==r['min_DP'] and x['min_ingroup_call_fraction']==r['min_ingroup_call_fraction'] and x['region'] in ('chr16A','chr16B')]
   ab=sum(x['ABBA'] for x in other);ba=sum(x['BABA'] for x in other);do=(ab-ba)/(ab+ba) if ab+ba else 'NA'
   r['D_chr16A_plus_B']=do;r['Delta_D']=float(r['D'])-do if do!='NA' and r['D']!='NA' else 'NA';r['ABBA_chr16A_plus_B']=ab;r['BABA_chr16A_plus_B']=ba;r['informative_chr16A_plus_B']=sum(x['informative_sites'] for x in other)
   # Block bootstrap for delta D using 250-kb source-window midpoint blocks.
   q,d,f=r['quartet_id'],r['min_DP'],r['min_ingroup_call_fraction']; per={}
   for (qq,wi,dd,ff),a in bywin.items():
    if (qq,dd,ff)!=(q,d,f):continue
    region=wr[wi]['region']
    if region not in ('chr16_supergene','chr16A','chr16B'):continue
    mid=int(wr[wi]['mid']); block=(mid-1)//250000
    key=('inside' if region=='chr16_supergene' else 'outside',block)
    z=per.setdefault(key,[0.,0.,0]);z[0]+=a['ABBA'];z[1]+=a['BABA'];z[2]+=a['informative_sites']
   ins=[v for (reg,_),v in per.items() if reg=='inside' and v[2]>0];outs=[v for (reg,_),v in per.items() if reg=='outside' and v[2]>0]
   r['n_informative_physical_blocks_inside']=len(ins);r['n_informative_physical_blocks_outside']=len(outs)
   r['n_informative_physical_blocks']='%d/%d'%(len(ins),len(outs))
   if len(ins)>=3 and len(outs)>=3:
    rng=random.Random(20261118+DP.index(d)*100+FRAC.index(f)+(0 if q==Q[0] else 1000));vals=[]
    for _ in range(2000):
     aa=bb=cc=dden=0.
     for j in range(len(ins)):
      x=ins[rng.randrange(len(ins))];aa+=x[0];bb+=x[1]
     for j in range(len(outs)):
      x=outs[rng.randrange(len(outs))];cc+=x[0];dden+=x[1]
     if aa+bb and cc+dden:vals.append((aa-bb)/(aa+bb)-(cc-dden)/(cc+dden))
    vals.sort();r['delta_D_block_bootstrap_low']=vals[int(.025*(len(vals)-1))] if vals else 'NA';r['delta_D_block_bootstrap_high']=vals[int(.975*(len(vals)-1))] if vals else 'NA';r['uncertainty']='250-kb physical blocks; 2000 percentile bootstrap replicates'
   else:r['delta_D_block_bootstrap_low']='NA';r['delta_D_block_bootstrap_high']='NA';r['uncertainty']='NA fewer than three informative blocks in one region'
  else:
   for k in ('n_informative_physical_blocks_inside','n_informative_physical_blocks_outside','delta_D_block_bootstrap_low','delta_D_block_bootstrap_high'):r[k]='NA'
   for k in ('D_chr16A_plus_B','Delta_D','ABBA_chr16A_plus_B','BABA_chr16A_plus_B','informative_chr16A_plus_B'):r[k]='NA'
 write(OUT/'genotype_filter_sensitivity.tsv',out)
if __name__=='__main__':main()
