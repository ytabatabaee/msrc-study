#!/usr/bin/env python3
"""Additive validation analyses for the Neurospora L10 nucleotide-D package."""
import csv, json, math, random, sys, zipfile
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[6]
BASE=ROOT/'empirical/neurospora/analysis/patterson_d'
OUT=BASE/'final_validation'; RESULTS=OUT/'results'
sys.path.insert(0,str(BASE/'scripts'))
from run_patterson_d import GENO,WINDOWS,LINKAGE_GROUPS,PRIMARY_P1,PRIMARY_P2,P3,read_genotypes,load_windows,vector_site_contributions,safe_div
from validate_outgroup_alignment import alignment_columns,read_fasta_member
SEED=20261118; BASES=set(b'ACGT')

def write_tsv(path,rows):
    if rows:
        with open(path,'w',newline='') as f:
            fields=[]
            for row in rows:
                for key in row:
                    if key not in fields:fields.append(key)
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

def pooled(stats,mask):
    mask=np.asarray(mask,dtype=bool); valid=stats['polarizable']&mask
    a=float(stats['abba'][mask].sum()); b=float(stats['baba'][mask].sum())
    return {'source_sites':int(mask.sum()),'polarizable_sites':int(valid.sum()),'informative_sites':int(((stats['abba']+stats['baba']>0)&valid).sum()),'ABBA':a,'BABA':b,'D':safe_div(a-b,a+b)}

def scalar(p1,p2,p3,out):
    """Independent one-row counter; does not call the production estimator."""
    if out not in BASES:return ('outgroup_gap' if out==45 else 'outgroup_ambiguous',0.,0.)
    gs=[]
    for group in (p1,p2,p3):
        called=[x for x in group if x in BASES]
        if not called or len(called)<math.ceil(len(group)*.5):return ('insufficient_calls',0.,0.)
        gs.append(called)
    alleles=set(x for g in gs for x in g)
    if len(alleles)<2:return ('monomorphic',0.,0.)
    if len(alleles)>2:return ('multiallelic',0.,0.)
    if out not in alleles:return ('outgroup_other',0.,0.)
    der=next(x for x in alleles if x!=out); f=[sum(x==der for x in g)/float(len(g)) for g in gs]
    return ('polarizable',(1-f[0])*f[1]*f[2],f[0]*(1-f[1])*f[2])

def provenance():
    rows=[]
    with open(BASE/'results/sample_mapping_audit.tsv') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['sample_id'].startswith('CJ'):
                rows.append({'sample_id':r['sample_id'],'matrix_present':r['matrix_present'],'mating_type':r['mating_type_assignment'],'source_assignment':r['assignment_source'],'primary_status':r['status_for_primary_group_frequency_test'],'note':'Unresolved and excluded; no ID/tree-based assignment.' if r['mating_type_assignment']=='unknown' else ''})
    write_tsv(RESULTS/'sample_provenance_audit.tsv',rows)

def main():
    OUT.mkdir(exist_ok=True); RESULTS.mkdir(exist_ok=True); provenance()
    samples,pos,calls=read_genotypes(); windows=load_windows(); ix={s:i for i,s in enumerate(samples)}
    if '--samples-only' in sys.argv:
        sample_analysis(samples,pos,calls,ix)
        genealogy_analysis()
        print('sample sensitivity complete')
        return
    p1=[ix[s] for s in PRIMARY_P1]; p2=[ix[s] for s in PRIMARY_P2]; p3=[ix[s] for s in P3]
    winrows=[]; rawrows=[]; unassigned=[]; profile=[]; extent=[]
    with zipfile.ZipFile(ROOT/'data/neurospora/raw/outgroup_alignment.zip') as zf:
      for lg in LINKAGE_GROUPS:
        rec=read_fasta_member(zf,'outgroup_alignment/%s.outgroups.fasta'%lg)
        coords=np.frombuffer(pos[lg],dtype=np.uint32); mat=np.frombuffer(calls[lg],dtype=np.uint8).reshape((-1,len(samples)))
        cols=np.frombuffer(alignment_columns(rec['2509@6'],pos[lg]),dtype=np.uint32); out=np.frombuffer(rec['Ndis@999'],dtype=np.uint8)[cols]
        stats=vector_site_contributions(mat,out,(p1,p2,p3),.5); assigned=np.zeros(len(coords),bool)
        for st,en,mid,nsite in windows[lg]:
            lo=int(np.searchsorted(coords,st,'left')); hi=int(np.searchsorted(coords,en,'right')); assigned[lo:hi]=True; m=np.zeros(len(coords),bool);m[lo:hi]=True; z=pooled(stats,m)
            span=en-st+1; n=int(m.sum()); oc=np.isin(out[lo:hi],list(BASES)); c1=np.isin(mat[lo:hi][:,p1],list(BASES)).sum()/float(max(1,n*len(p1)));c2=np.isin(mat[lo:hi][:,p2],list(BASES)).sum()/float(max(1,n*len(p2)));c3=np.isin(mat[lo:hi][:,p3],list(BASES)).sum()/float(max(1,n*len(p3)))
            winrows.append({'linkage_group':lg,'start':st,'end':en,'midpoint':mid,'source_snps':nsite,'observed_rows':n,'physical_span_bp':span,'snps_per_mb':nsite/(span/1e6),'outgroup_callability':float(oc.mean()) if n else 'NA','P1_callability':c1,'P2_callability':c2,'P3_callability':c3,'polarizable_sites':z['polarizable_sites'],'informative_sites':z['informative_sites'],'informative_per_mb':z['informative_sites']/(span/1e6),'ABBA':z['ABBA'],'BABA':z['BABA'],'D':z['D'],'inside_LGI_1p5_8Mb':int(lg=='LGI' and 1_500_000<=mid<=8_000_000)})
        if (~assigned).any():unassigned.append({'linkage_group':lg,**pooled(stats,~assigned)})
        # Fixed spatial bins, pooled site contributions, no average of window D.
        if lg=='LGI':
            for bw in [50_000,100_000,250_000,500_000]:
                starts=((coords-1)//bw)*bw+1
                for st in np.unique(starts):
                    mask=starts==st; z=pooled(stats,mask)
                    z['Ndis_callability']=float(np.isin(out[mask],list(BASES)).mean()) if mask.sum() else 'NA'
                    z['P1_callability']=float(np.isin(mat[mask][:,p1],list(BASES)).mean()) if mask.sum() else 'NA'
                    z['P2_callability']=float(np.isin(mat[mask][:,p2],list(BASES)).mean()) if mask.sum() else 'NA'
                    profile.append({'linkage_group':lg,'bin_width_bp':bw,'start':int(st),'end':int(st+bw-1),'midpoint':float(st+(bw-1)/2),**z})
            for fw in [100_000,250_000,500_000,1_000_000]:
                zones=[('focal_1.5_8Mb',lambda q:1_500_000<=q<=8_000_000),('left_immediate',lambda q:1_500_000-fw<=q<1_500_000),('left_distant',lambda q:q<1_500_000-fw),('right_immediate',lambda q:8_000_000<q<=8_000_000+fw),('right_distant',lambda q:q>8_000_000+fw),('LGI_complement',lambda q:not (1_500_000<=q<=8_000_000))]
                for name,pred in zones:
                    m=np.fromiter((pred(int(q)) for q in coords),bool,count=len(coords));extent.append({'flank_width_bp':fw,'region':name,**pooled(stats,m)})
        else:
            extent.append({'flank_width_bp':'NA','region':lg,**pooled(stats,np.ones(len(coords),bool))})
        # Stratified deterministic exact comparison of the scalar and vector counters.
        strata={}
        calln=np.isin(mat[:,p1+p2],list(BASES)).sum(axis=1)
        for i,q in enumerate(coords):
            region='LGI_inside' if lg=='LGI' and 1_500_000<=q<=8_000_000 else 'LGI_flank' if lg=='LGI' and (1_250_000<=q<1_500_000 or 8_000_000<q<=8_250_000) else 'LGI_outside' if lg=='LGI' else lg
            state='outgroup_missing' if out[i] not in BASES else 'outgroup_ref' if out[i]==mat[i,ix['2509@6']] else 'outgroup_nonref'
            cstate='low_group_callability' if calln[i] < int(.75*(len(p1)+len(p2))) else 'high_group_callability'
            strata.setdefault((region,state,cstate),[]).append(i)
        rng=random.Random(SEED+sum(map(ord,lg))); selected=set()
        for key,ids in sorted(strata.items()):selected.update(rng.sample(ids,min(5,len(ids))))
        for i in sorted(selected):
            ans=scalar([int(mat[i,j]) for j in p1],[int(mat[i,j]) for j in p2],[int(mat[i,j]) for j in p3],int(out[i]))
            one=vector_site_contributions(mat[i:i+1],out[i:i+1],(p1,p2,p3),.5)
            pstatus=next(k for k in ['outgroup_gap','outgroup_ambiguous','insufficient_calls','monomorphic','multiallelic','outgroup_other','polarizable'] if bool(one[k][0]))
            aa=float(one['abba'][0]);bb=float(one['baba'][0])
            rawrows.append({'linkage_group':lg,'position':int(coords[i]),'stratum':str(next(k for k,v in strata.items() if i in v)),'source_window_assigned':int(assigned[i]),'scalar_status':ans[0],'vector_status':pstatus,'scalar_ABBA':ans[1],'vector_ABBA':aa,'scalar_BABA':ans[2],'vector_BABA':bb,'max_abs_difference':max(abs(ans[1]-aa),abs(ans[2]-bb))})
    write_tsv(RESULTS/'ascertainment_diagnostics.tsv',winrows);write_tsv(RESULTS/'estimator_validation.tsv',rawrows);write_tsv(RESULTS/'unassigned_219_audit.tsv',unassigned);write_tsv(RESULTS/'spatial_extent_analysis.tsv',extent);write_tsv(RESULTS/'spatial_profile_bins.tsv',profile)
    common_support(winrows)
    sample_analysis(samples,pos,calls,ix)
    genealogy_analysis()
    print('SNP rows',sum(len(x) for x in pos.values()),'source windows',sum(len(x) for x in windows.values()),'independent checks',len(rawrows),'unassigned rows',sum(x['source_sites'] for x in unassigned))

def common_support(wins):
    lgi=[x for x in wins if x['linkage_group']=='LGI']; a=[x for x in lgi if str(x['inside_LGI_1p5_8Mb'])=='1']; b=[x for x in lgi if str(x['inside_LGI_1p5_8Mb'])=='0']; cov=['physical_span_bp','outgroup_callability','P1_callability','P2_callability']; out=[]
    for name,vars_ in [('span',['physical_span_bp']),('outgroup_callability',['outgroup_callability']),('P1_callability',['P1_callability']),('P2_callability',['P2_callability']),('informative_density',['informative_per_mb']),('joint',cov)]:
        bounds={}
        for v in vars_:
            av=[float(x[v]) for x in a if x[v]!='NA'];bv=[float(x[v]) for x in b if x[v]!='NA'];bounds[v]=(max(min(av),min(bv)),min(max(av),max(bv)))
        for region,group in [('inside',a),('outside',b)]:
            keep=[x for x in group if all(lo<=float(x[v])<=hi for v,(lo,hi) in bounds.items())]
            aa=sum(float(x['ABBA']) for x in keep);bb=sum(float(x['BABA']) for x in keep)
            out.append({'filter':name,'region':region,'windows_retained':len(keep),'ABBA':aa,'BABA':bb,'D':safe_div(aa-bb,aa+bb),'bounds':json.dumps(bounds)})
    # D-blind nearest-neighbour matching on physical span and raw callability only.
    matchcov=['physical_span_bp','outgroup_callability','P1_callability','P2_callability']
    ia=[x for x in a if int(x['informative_sites'])>0];ob=[x for x in b if int(x['informative_sites'])>0]; scales={v:np.std([float(x[v]) for x in ia+ob]) or 1 for v in matchcov}; used=set();pairs=[]
    for outw in ob:
        choices=[(sum(((float(x[v])-float(outw[v]))/scales[v])**2 for v in matchcov),j,x) for j,x in enumerate(ia) if j not in used]
        if not choices:break
        _,j,inw=min(choices,key=lambda y:y[0]);used.add(j);pairs.append((inw,outw))
    for region,idx in [('inside_matched',0),('outside_matched',1)]:
        group=[p[idx] for p in pairs];aa=sum(float(x['ABBA']) for x in group);bb=sum(float(x['BABA']) for x in group)
        out.append({'filter':'nearest_match_no_D','region':region,'windows_retained':len(group),'ABBA':aa,'BABA':bb,'D':safe_div(aa-bb,aa+bb),'bounds':'nearest standardized distance: physical span, Ndis, P1, P2 callability; no D or informative density used'})
    write_tsv(RESULTS/'ascertainment_common_support.tsv',out)

def sample_analysis(samples,pos,calls,ix):
    p1=[ix[x] for x in PRIMARY_P1];p2=[ix[x] for x in PRIMARY_P2];p3=[ix[x] for x in P3]; rows=[]; draws=[];pairs=[]
    targets=['LGI','LGIV']
    with zipfile.ZipFile(ROOT/'data/neurospora/raw/outgroup_alignment.zip') as zf:
      for lg in targets:
        rec=read_fasta_member(zf,'outgroup_alignment/%s.outgroups.fasta'%lg)
        coords=np.frombuffer(pos[lg],dtype=np.uint32);mat=np.frombuffer(calls[lg],dtype=np.uint8).reshape((-1,len(samples))); cols=np.frombuffer(alignment_columns(rec['2509@6'],pos[lg]),dtype=np.uint32);out=np.frombuffer(rec['Ndis@999'],dtype=np.uint8)[cols]
        regions={'whole':np.ones(len(coords),bool),'inside_LGI':(coords>=1_500_000)&(coords<=8_000_000) if lg=='LGI' else np.zeros(len(coords),bool),'outside_LGI':~((coords>=1_500_000)&(coords<=8_000_000)) if lg=='LGI' else np.zeros(len(coords),bool)}
        full=vector_site_contributions(mat,out,(p1,p2,p3),.5)
        for group,name in [(p1,'mat_a'),(p2,'mat_A')]:
            for j in group:
                reduced=[x for x in group if x!=j]; gs=(reduced,p2,p3) if name=='mat_a' else (p1,reduced,p3);s=vector_site_contributions(mat,out,gs,.5)
                for region,m in regions.items():
                    if region!='whole' and lg!='LGI':continue
                    z=pooled(s,m); base=pooled(full,m)
                    rows.append({'analysis':'leave_one_out','linkage_group':lg,'omitted_sample':samples[j],'subset_n_P1':len(gs[0]),'subset_n_P2':len(gs[1]),'region':region,'D':z['D'],'ABBA':z['ABBA'],'BABA':z['BABA'],'informative_sites':z['informative_sites'],'delta_vs_full_D':z['D']-base['D']})
        # Balanced random subsets, selection independent of D; 20 LGI and 10 LGIV draws per size.
        rng=random.Random(SEED+(1 if lg=='LGI' else 4))
        for n in [2,4,6,8]:
            for rep in range(20 if lg=='LGI' else 10):
                aa=sorted(rng.sample(p1,n));bb=sorted(rng.sample(p2,n));s=vector_site_contributions(mat,out,(aa,bb,p3),.5)
                for region,m in regions.items():
                    if region!='whole' and lg!='LGI':continue
                    z=pooled(s,m); draws.append({'analysis':'balanced_random_subset','linkage_group':lg,'seed':SEED,'replicate':rep+1,'n_per_mating_type':n,'sample_P1':';'.join(samples[j] for j in aa),'sample_P2':';'.join(samples[j] for j in bb),'region':region,**z})
        if lg=='LGIV':
            for a in [ix['CJ04'],ix['CJ16']]:
                for b in [ix['CJ03'],ix['CJ15']]:
                    s=vector_site_contributions(mat,out,([a],[b],p3),1.0);z=pooled(s,regions['whole']);pairs.append({'analysis':'single_matched_pair','linkage_group':lg,'region':'whole','pair_mat_a':samples[a],'pair_mat_A':samples[b],**z})
    write_tsv(RESULTS/'sample_sensitivity.tsv',rows+draws+pairs)

def genealogy_analysis():
    data=[]
    infos={}; calls={}
    # Aggregate source-window diagnostics onto the published 100-kb genealogy bins.
    winbins={}
    with open(RESULTS/'ascertainment_diagnostics.tsv') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['linkage_group']!='LGI':continue
            k=(int(r['midpoint'])//100000)*100000
            x=winbins.setdefault(k,[0,0.0,0])
            x[0]+=int(r['informative_sites']);x[1]+=int(r['observed_rows'])*float(r['outgroup_callability']);x[2]+=int(r['observed_rows'])
    for k,(n,ok,total) in winbins.items():infos[k]=n;calls[k]=ok/total if total else float('nan')
    with open(BASE/'results/matched_nucleotide_genealogy.tsv') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            try:
                st=int(float(r['start']));obj={'start':st,'mid':(float(r['start'])+float(r['end']))/2,'D':float(r['D']),'q3':float(r['q3']),'informative':infos.get(st,0),'Ndis_callability':calls.get(st,float('nan'))}
                if math.isfinite(obj['D']) and math.isfinite(obj['q3']):data.append(obj)
            except Exception:pass
    inside=lambda x:1_500_000<=x['mid']<=8_000_000
    def corr(z):return float(np.corrcoef([x['D'] for x in z],[x['q3'] for x in z])[0,1]) if len(z)>2 else float('nan')
    rows=[];rng=np.random.default_rng(SEED)
    for name,pred in [('whole_LGI',lambda x:True),('inside',inside),('outside',lambda x:not inside(x))]:
        for threshold in [0,10,25,50,100]:
            z=[x for x in data if pred(x) and x['informative']>=threshold]
            # Resample contiguous 500-kb units, retaining all bins within each selected block.
            blockgroups={}
            for x in z:blockgroups.setdefault(x['start']//500_000,[]).append(x)
            ids=sorted(blockgroups); boot=[]
            if len(ids)>=3:
                for _ in range(1000):
                    draw=rng.choice(ids,size=len(ids),replace=True);sample=[]
                    for bid in draw:sample.extend(blockgroups[bid])
                    r=corr(sample)
                    if math.isfinite(r):boot.append(r)
            rows.append({'subset':name,'metric':'D_vs_TWISST_topo3','min_informative_sites':threshold,'n_bins':len(z),'n_500kb_blocks':len(ids),'pearson_r':corr(z),'CI95_low':float(np.quantile(boot,.025)) if boot else 'NA','CI95_high':float(np.quantile(boot,.975)) if boot else 'NA','interpretation':'TWISST q3 is the unrooted split pairing L10 mat-A with N. hispaniola in trees also containing N. crassa; N. discreta is absent.'})
    for threshold in [.8,.9,.95]:
        z=[x for x in data if x['informative']>=1 and x['Ndis_callability']>=threshold]
        rows.append({'subset':'outgroup_callability_ge_%.2f'%threshold,'metric':'D_vs_TWISST_topo3','min_informative_sites':1,'n_bins':len(z),'n_500kb_blocks':'NA','pearson_r':corr(z),'CI95_low':'NA','CI95_high':'NA','interpretation':'Descriptive subset filter based on source-window N. discreta callability averaged into physical bins.'})
    z=[x for x in data if x['informative']>=1]
    dm={b:np.mean([x['D'] for x in z if inside(x)==b]) for b in [True,False]};qm={b:np.mean([x['q3'] for x in z if inside(x)==b]) for b in [True,False]};c=[dict(x,D=x['D']-dm[inside(x)],q3=x['q3']-qm[inside(x)]) for x in z]
    blocks={}
    for x in c:blocks.setdefault(x['start']//500_000,[]).append(x)
    ids=sorted(blocks);boot=[]
    for _ in range(1000):
        draw=rng.choice(ids,size=len(ids),replace=True);sample=[]
        for bid in draw:sample.extend(blocks[bid])
        r=corr(sample)
        if math.isfinite(r):boot.append(r)
    rows.append({'subset':'region_centered','metric':'D_vs_TWISST_topo3','min_informative_sites':1,'n_bins':len(c),'n_500kb_blocks':len(ids),'pearson_r':corr(c),'CI95_low':float(np.quantile(boot,.025)) if boot else 'NA','CI95_high':float(np.quantile(boot,.975)) if boot else 'NA','interpretation':'D and q3 each centered within focal/complement classes; 500-kb moving-block bootstrap. TWISST tree lacks N. discreta.'})
    write_tsv(RESULTS/'genealogy_correspondence.tsv',rows)

if __name__=='__main__':main()
