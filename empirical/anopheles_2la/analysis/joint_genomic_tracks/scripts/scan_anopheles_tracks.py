#!/usr/bin/env python3
"""Scan Ag3.10 2L all-site population tracks in bounded 1-Mb chunks.

Adds corrected allele-state-aware Patterson D and Dsuite-style f_dM using the
frozen quartet designs, plus gamb_colu-mask dXY, Hudson F_ST and pi tracks.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,math,os,sys,time
from collections import defaultdict
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[5]
BASE=ROOT/'empirical/anopheles_2la'
DATA=ROOT/'data/anopheles_2la'
OUT=BASE/'analysis/joint_genomic_tracks'
RES=OUT/'results'
GRID=DATA/'processed/stage2r_window_grid.tsv'
INV=(20524058,42165532)
WIN=50000
CHUNK=1000000
sys.path.insert(0,str(BASE/'scripts'))
sys.path.insert(0,str(ROOT/'empirical/scripts'))
import run_patterson_d as pd
from joint_genomic_tracks.fdm import fdm_site
from patterson_d_core import calculate_window

def read_tsv(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(path,rows,fields=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    if fields is None: fields=list(rows[0]) if rows else []
    tmp=path.with_name(path.name+'.partial')
    with tmp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    os.replace(tmp,path)
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def allele_freq(gt, idx, anc):
    z=gt[:,idx,:];ok=z>=0;n=ok.sum(axis=(1,2));cnt=(z==anc[:,None,None]).sum(axis=(1,2))
    return n,np.divide(cnt,n,out=np.zeros(len(n),float),where=n>0)
def corrected_quartet(gt,design,ids,min_mac):
    keys=('P1_sample_ids','P2_sample_ids','P3_sample_ids','outgroup_sample_ids')
    names=('P1','P2','P3','O')
    gi={n:[ids.index(s) for s in design[k].split(';') if s] for n,k in zip(names,keys)}
    calls=gt>=0;complete=calls.all(axis=2);missing=1.0-complete.mean(axis=1)
    nalleles=np.stack([np.any(gt==a,axis=(1,2)) for a in range(4)],axis=1).sum(axis=1).astype(np.int8)
    keep=(missing<=.25)&(nalleles<=2)
    # Population callability and fixed outgroup requirements reproduce the frozen D rules.
    gc={n:gt[:,ix,:] for n,ix in gi.items()}
    n={k:np.sum((v>=0).all(axis=2),axis=1)*2 for k,v in gc.items()}
    for k in names:keep &= n[k]>=2
    out=gc['O'];ofull=(out>=0).all(axis=2).all(axis=1)
    amin=np.min(np.where(out>=0,out,4),axis=(1,2));amax=np.max(np.where(out>=0,out,-1),axis=(1,2))
    fixed=ofull&(amin==amax);keep &= fixed
    ancestral=amin
    total_called=sum(n.values())
    derived_counts={}
    for k in ('P1','P2','P3'):
        vals=gc[k];complete_group=(vals>=0).all(axis=2)
        derived_counts[k]=np.sum(complete_group[:,:,None]&(vals!=ancestral[:,None,None]),axis=(1,2))
    # Outgroup is fixed ancestral by construction; its derived count is zero.
    mac=np.minimum(sum(derived_counts.values()),total_called-sum(derived_counts.values()))
    keep &= mac>=min_mac
    p={k:np.divide(derived_counts[k],n[k],out=np.zeros(len(gt),float),where=n[k]>0) for k in ('P1','P2','P3')}
    a=(1-p['P1'])*p['P2']*p['P3'];b=p['P1']*(1-p['P2'])*p['P3']
    fnum=np.zeros(len(gt));fden=np.zeros(len(gt));included=keep & ~((p['P1']==p['P2'])&(p['P2']==p['P3'])&((p['P1']==0)|(p['P1']==1)))
    # Vectorized Dsuite f_dM denominator algebra, exactly equivalent to fdm_site.
    p1,p2,p3=p['P1'],p['P2'],p['P3'];case1=p1<=p2;case2=p2>p3
    fden=np.where(case1&case2,p2*(p2-p1),fden)
    fden=np.where(case1&~case2,p3*(p3-p1),fden)
    fden=np.where(~case1&(p1>p3),p1*(p1-p2),fden)
    fden=np.where(~case1&(p1<=p3),p3*(p3-p2),fden)
    fnum=a-b
    dkeep=keep
    vals={'abba':a,'baba':b,'fdm_num':fnum,'fdm_den':fden}
    sums={k:float(v[dkeep].sum()) for k,v in vals.items()}
    sums.update(n_pass=int(dkeep.sum()),n_informative=int(((a+b)>0)[dkeep].sum()),n_fdm_included=int(included.sum()))
    denom=sums['abba']+sums['baba'];sums['D']=(sums['abba']-sums['baba'])/denom if denom else None
    sums['fdM']=sums['fdm_num']/sums['fdm_den'] if sums['fdm_den'] else None
    # diagnostics used to compare with the saved legacy estimator
    sums['n_biallelic_or_monomorphic']=int((nalleles<=2).sum());sums['n_legacy_code_le1']=int((~np.any(gt>1,axis=(1,2))).sum())
    return sums,keep,a,b,fnum,fden,included

def pi_site(counts):
    n=counts.sum(axis=1);num=n*n-(counts*counts).sum(axis=1)
    den=n*(n-1)
    return np.divide(num,den,out=np.zeros(len(n),float),where=den>0),n>=2

def scalar_diversity(gt,mask,groups,ia,ib):
    """Independent small-window scalar cross-check for per-site estimators."""
    vals=[]
    for s in np.flatnonzero(mask):
        def count(ix):
            z=gt[s,ix,:].ravel();z=z[z>=0];return np.bincount(z,minlength=4),len(z)
        cx,nx=count(groups[ia]);cy,ny=count(groups[ib])
        if nx<2 or ny<2:continue
        pix=1.0-np.sum(cx*(cx-1))/(nx*(nx-1));piy=1.0-np.sum(cy*(cy-1))/(ny*(ny-1))
        dxy=1.0-np.sum(cx*cy)/(nx*ny)
        vals.append((pix,piy,dxy))
    return vals

def calc_diversity(gt,mask,group_ix):
    counts={}
    pis={}
    for name,idx in group_ix.items():
        z=gt[:,idx,:];c=np.stack([(z==a).sum(axis=(1,2)) for a in range(4)],axis=1).astype(np.int64)
        pi,ok=pi_site(c);counts[name]=c;pis[name]=(pi,ok)
    pi_rows=[]
    for name,(pi,ok) in pis.items():
        use=mask&ok
        pi_rows.append((name,float(pi[use].sum()),int(use.sum()),int(counts[name][mask].sum(axis=1).sum())))
    pair_rows=[]
    pairs=(('G_s','G_i','primary_within_gambiae'),('C_s','C_i','primary_within_coluzzii'),('G_s','C_s','cross_species_standard'),('G_i','C_i','cross_species_inverted'),('G_s','C_i','cross_species_mixed_1'),('G_i','C_s','cross_species_mixed_2'))
    for x,y,label in pairs:
        cx,cy=counts[x],counts[y];nx=cx.sum(axis=1);ny=cy.sum(axis=1)
        pairden=nx*ny;dxy_num=pairden-(cx*cy).sum(axis=1)
        eligible=mask&(nx>0)&(ny>0)
        dxy=np.divide(dxy_num,pairden,out=np.zeros(len(gt),float),where=pairden>0)
        pair_rows.append((x,y,label,float(dxy[eligible].sum()),int(eligible.sum()),int(pairden[eligible].sum()),
            float(dxy_num[eligible].sum()),int(pairden[eligible].sum())))
    return pi_rows,pair_rows,counts,pis

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pilot',action='store_true');ap.add_argument('--chunk-bp',type=int,default=CHUNK);ap.add_argument('--force',action='store_true');args=ap.parse_args()
    import malariagen_data
    RES.mkdir(parents=True,exist_ok=True)
    prefix='anopheles_2la_pilot' if args.pilot else 'anopheles_2la'
    paths={'windows':RES/f'{prefix}_fdm_D_corrected_windows.tsv','pi':RES/f'{prefix}_pi_windows.tsv','dxy':RES/f'{prefix}_dxy_fst_windows.tsv','audit':RES/f'{prefix}_site_mask_audit.tsv','blocks':RES/f'{prefix}_fdm_500kb_blocks.tsv','checkpoint':RES/f'{prefix}_scan_checkpoint.json'}
    grid=read_tsv(GRID);designs=pd.design_rows()
    groups,_=pd.group_ids();sample_ids=sorted({s for d in designs for k in ('P1_sample_ids','P2_sample_ids','P3_sample_ids','outgroup_sample_ids') for s in d[k].split(';') if s})
    group_ix={k:[sample_ids.index(s) for s in v] for k,v in groups.items() if k in ('G_s','G_i','C_s','C_i')}
    ag=malariagen_data.Ag3(url=pd.BUCKET,check_location=False,show_progress=False);meta=ag.sample_metadata(sample_sets=pd.SAMPLE_SET);meta_ix={str(s):i for i,s in enumerate(meta.sample_id)}
    indices=[meta_ix[s] for s in sample_ids]
    if args.pilot:
        grid=pd.read_tsv(DATA/'processed/stage1r_pilot_windows.tsv')
    design_ids={d['test_id']:d for d in designs}
    window_rows=[];pi_rows=[];dxy_rows=[];block_rows=[];mask_audit=[]
    if not args.force:
        for key,var in [('windows','window_rows'),('pi','pi_rows'),('dxy','dxy_rows'),('blocks','block_rows'),('audit','mask_audit')]:
            if paths[key].exists(): locals()[var].extend(read_tsv(paths[key]))
        # Explicit initialization avoids relying on locals() assignment semantics.
        window_rows=read_tsv(paths['windows']) if paths['windows'].exists() else []
        pi_rows=read_tsv(paths['pi']) if paths['pi'].exists() else []
        dxy_rows=read_tsv(paths['dxy']) if paths['dxy'].exists() else []
        block_rows=read_tsv(paths['blocks']) if paths['blocks'].exists() else []
        mask_audit=read_tsv(paths['audit']) if paths['audit'].exists() else []
    done_regions={x['region'] for x in mask_audit}
    chunks=[];cur=[]
    for w in grid:
        if not cur or int(w['end'])-int(cur[0]['start'])+1<=args.chunk_bp:cur.append(w)
        else:chunks.append(cur);cur=[w]
    if cur:chunks.append(cur)
    for ci,ws in enumerate(chunks,1):
        chrom=ws[0]['chrom'];start=int(ws[0]['start']);end=int(ws[-1]['end']);region=f'{chrom}:{start}-{end}'
        if region in done_regions:print(f'[{ci}/{len(chunks)}] skip completed {region}',flush=True);continue
        print(f'[{ci}/{len(chunks)}] query {region}',flush=True);t0=time.time()
        ds=ag.snp_calls(region=region,sample_sets=pd.SAMPLE_SET,sample_indices=indices,site_mask=None,inline_array=True,chunks='native')
        if [str(x) for x in ds.sample_id.values.tolist()]!=sample_ids:raise RuntimeError('sample order mismatch')
        pos=np.asarray(ds.variant_position.values);gt=np.asarray(ds.call_genotype.values,dtype=np.int8)
        mask=np.asarray(ds['variant_filter_pass_gamb_colu'].values,dtype=bool)
        site_mask=ag.site_filters(region=region,mask='gamb_colu',field='filter_pass')
        site_mask=np.asarray(site_mask.compute() if hasattr(site_mask,'compute') else site_mask,dtype=bool)
        if not np.array_equal(mask,site_mask):raise RuntimeError(f'gamb_colu mask misaligned for {region}')
        if len(np.unique(pos))!=len(pos) or (len(pos) and (pos.min()<start or pos.max()>end)):raise RuntimeError(f'coordinate uniqueness/range failed: {region}')
        refseq=np.asarray(ag.genome_sequence(region).compute()).astype(str)
        non_n=np.char.upper(refseq)!='N';expected_pos=np.flatnonzero(non_n)+start
        if not np.array_equal(pos,expected_pos):raise RuntimeError(f'Ag3 sites do not match all non-N reference coordinates: {region}')
        n_observed_alleles=np.stack([np.any(gt==a,axis=(1,2)) for a in range(4)],axis=1).sum(axis=1)
        invariant=int((n_observed_alleles<=1).sum())
        mask_audit.append(dict(region=region,n_positions=len(pos),n_non_n_reference_positions=int(non_n.sum()),n_reference_N_positions=int((~non_n).sum()),n_gamb_colu_pass=int(mask.sum()),n_gamb_colu_fail=int((~mask).sum()),n_invariant_or_unobserved_variable=int(invariant),mean_complete_sample_fraction=float(np.mean(np.all(gt>=0,axis=2))),query_seconds=round(time.time()-t0,3)))
        if args.pilot:
            # independent scalar validation and numerical reference for the corrected quartet estimator
            w=ws[0];d=designs[0];take=[sample_ids.index(s) for k in ('P1_sample_ids','P2_sample_ids','P3_sample_ids','outgroup_sample_ids') for s in d[k].split(';') if s]
            local=gt[:,take,:];off={};q=0
            for k in ('P1_sample_ids','P2_sample_ids','P3_sample_ids','outgroup_sample_ids'):
                n=len(d[k].split(';'));off[k]=list(range(q,q+n));q+=n
            z,keep,a,b,fn,fd,inc=corrected_quartet(local,d,[sample_ids[i] for i in take],0)
            # Check the vectorized denominator against the shared scalar f_dM formula.
            from joint_genomic_tracks.fdm import fdm_site
            cols={k:off[k+'_sample_ids'] for k in ('P1','P2','P3')}
            check_ix=np.flatnonzero(keep)[:1000]
            for si in check_ix:
                pvals=[]
                for k in ('P1','P2','P3'):
                    zz=local[si,cols[k],:];cc=zz[zz>=0];pvals.append(float(np.mean(cc!=int(np.min(local[si,off['outgroup_sample_ids'],:])))))
                ref=fdm_site(*pvals)
                if abs(ref.denominator-fd[si])>1e-12 or abs(ref.numerator-fn[si])>1e-12:raise RuntimeError('vectorized f_dM site terms differ from shared scalar implementation')
            # Scalar f_dM terms independently recomputed via the shared scalar reference function.
            # Compare diversity estimator with direct per-site combinatorics for the primary pair.
            pcheck=mask
            scal=scalar_diversity(gt,pcheck,group_ix,'G_s','G_i')
            cgs=np.stack([(gt[:,group_ix['G_s'],:]==k).sum(axis=(1,2)) for k in range(4)],axis=1)
            cgi=np.stack([(gt[:,group_ix['G_i'],:]==k).sum(axis=(1,2)) for k in range(4)],axis=1)
            pgs,_=pi_site(cgs);pgi,_=pi_site(cgi);ngs=cgs.sum(axis=1);ngi=cgi.sum(axis=1);dx=1-np.divide((cgs*cgi).sum(axis=1),ngs*ngi,out=np.zeros(len(ngs),float),where=(ngs*ngi)>0)
            use=mask&(ngs>=2)&(ngi>=2)
            if len(scal)!=int(use.sum()) or abs(sum(t[0] for t in scal)/len(scal)-pgs[use].mean())>1e-12 or abs(sum(t[2] for t in scal)/len(scal)-dx[use].mean())>1e-12:raise RuntimeError('independent diversity estimator check failed')
            print(f'PILOT_OK corrected D={z["D"]} f_dM={z["fdM"]}; {len(scal)} G_s/G_i sites matched scalar pi/dXY',flush=True)
        # process each 50-kb window in this bounded chunk
        for w in ws:
            lo,hi=int(w['start']),int(w['end']);sel=(pos>=lo)&(pos<=hi);g=gt[sel];p=pos[sel]
            for design in designs:
                take=[sample_ids.index(s) for k in ('P1_sample_ids','P2_sample_ids','P3_sample_ids','outgroup_sample_ids') for s in design[k].split(';') if s]
                sub=g[:,take,:];local_ids=[sample_ids[i] for i in take]
                rowlist=[]
                for mac in (0,2):
                    v,keep,a,b,fn,fd,inc=corrected_quartet(sub,design,local_ids,mac)
                    region_class=w.get('region_class','')
                    outrow=dict(window_id=w['window_id'],chrom=chrom,start=lo,end=hi,region_class=region_class,test_id=design['test_id'],filter_label='primary' if mac==0 else 'MAC2',abba=v['abba'],baba=v['baba'],D_patterson_corrected=v['D'],f_dM_numerator=v['fdm_num'],f_dM_denominator=v['fdm_den'],f_dM=v['fdM'],n_sites_examined=len(p),n_sites_pass=v['n_pass'],n_sites_informative=v['n_informative'],n_fdm_included=v['n_fdm_included'],n_biallelic_or_monomorphic=v['n_biallelic_or_monomorphic'],n_legacy_code_le1=v['n_legacy_code_le1'],allele_state_rule='retain sites with <=2 observed nucleotide states across quartet; derived is non-outgroup state',max_missing_fraction=.25,min_callable_chromosomes=2,mac=mac)
                    rowlist.append(outrow);window_rows.append(outrow)
                    for bstart in range(((lo-1)//500000)*500000+1,hi+1,500000):
                        # 50-kb windows are contained in their corresponding 500-kb blocks.
                        block_rows.append(dict(test_id=design['test_id'],filter_label=outrow['filter_label'],region_class=region_class,block_start=((lo-1)//500000)*500000+1,numerator=v['fdm_num'],denominator=v['fdm_den'],abba=v['abba'],baba=v['baba']))
            # diversity statistics use the published gamb_colu accessibility mask
            m=mask[sel];gpass=g[m]
            pi_s,pair_s,cc,pp=calc_diversity(gpass,np.ones(len(gpass),dtype=bool),group_ix)
            for group,num,nsite,nchrom in pi_s:
                pi_rows.append(dict(window_id=w['window_id'],chrom=chrom,start=lo,end=hi,region_class=w.get('region_class',''),group=group,pi_per_callable_site=num/nsite if nsite else '',pi_sum=num,callable_sites=nsite,called_chromosomes=nchrom,site_mask='gamb_colu'))
            for x,y,label,num,nsite,npairs,raw_num,rawpairs in pair_s:
                dxy_rows.append(dict(window_id=w['window_id'],chrom=chrom,start=lo,end=hi,region_class=w.get('region_class',''),population_X=x,population_Y=y,comparison=label,dxy_per_callable_site=num/nsite if nsite else '',dxy_sum=num,callable_sites=nsite,effective_allele_pairs=npairs,dxy_pair_weighted=raw_num/rawpairs if rawpairs else '',hudson_fst='',hudson_fst_sites=0,pi_X='',pi_Y='',site_mask='gamb_colu'))
            # Hudson FST on pairwise sites callable for >=2 chromosomes in both groups.
            for r in dxy_rows[-len(pair_s):]:
                x,y=r['population_X'],r['population_Y'];cx,cy=cc[x],cc[y];nx=cx.sum(axis=1);ny=cy.sum(axis=1);eligible=np.ones(len(gpass),bool)&(nx>=2)&(ny>=2)
                pix,pokx=pp[x];piy,poky=pp[y];dx=1-np.divide((cx*cy).sum(axis=1),nx*ny,out=np.zeros(len(nx),float),where=nx*ny>0)
                eligible &= pokx&poky
                fstnum=float((dx[eligible]-(pix[eligible]+piy[eligible])/2).sum());fstden=float(dx[eligible].sum())
                r['hudson_fst']=fstnum/fstden if fstden else '';r['hudson_fst_sites']=int(eligible.sum());r['pi_X']=float(pix[eligible].mean()) if eligible.any() else '';r['pi_Y']=float(piy[eligible].mean()) if eligible.any() else ''
        print(f'[{ci}/{len(chunks)}] done positions={len(pos)} mask_pass={int(mask.sum())} elapsed={time.time()-t0:.1f}s; total window rows={len(window_rows)}',flush=True)
        # atomic checkpoint after each completed chunk
        write_tsv(paths['windows'],window_rows)
        write_tsv(paths['pi'],pi_rows)
        write_tsv(paths['dxy'],dxy_rows)
        write_tsv(paths['audit'],mask_audit)
        write_tsv(paths['blocks'],block_rows)
        (paths['checkpoint']).write_text(json.dumps(dict(completed_chunks=ci,last_region=region,completed_windows=len(set(r['window_id'] for r in window_rows)),updated=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())),indent=2)+'\n')
    print('SCAN_COMPLETE',len(window_rows),len(pi_rows),len(dxy_rows),flush=True)
if __name__=='__main__':main()
