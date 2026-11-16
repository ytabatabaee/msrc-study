#!/usr/bin/env python3
"""Independent genotype-level Patterson-D, ascertainment, and quartet validation."""
from __future__ import print_function
import csv,gzip,hashlib,json,math,os,statistics,sys,time
from collections import defaultdict,Counter
from pathlib import Path
from datetime import datetime,timezone

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[5]
AN=ROOT/'empirical/atlantic_cod/analysis/patterson_d'
VAL=AN/'validation'; RES=VAL/'results'; CFG=VAL/'config'
RAW=ROOT/'data/atlantic_cod/raw/zenodo/local_downloads'
TREES=ROOT/'data/atlantic_cod/processed/cod_window_trees.tsv'
REGIONS=ROOT/'data/atlantic_cod/metadata/region_manifest.tsv'
LGS=('LG01','LG02','LG07','LG12')
TRANS={('A','G'),('G','A'),('C','T'),('T','C')}

def read_tsv(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write_tsv(p,rows):
    p.parent.mkdir(parents=True,exist_ok=True)
    if not rows: return
    with p.open('w') as f:
        fields=[]
        for row in rows:
            for key in row:
                if key not in fields:fields.append(key)
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n',extrasaction='ignore');w.writeheader()
        for r in rows:w.writerow({k:fmt(v) for k,v in r.items()})
def fmt(x):
    if x is None:return 'NA'
    if isinstance(x,float):return '%.10g'%x if math.isfinite(x) else 'NA'
    return x
def dval(a,b):return (a-b)/(a+b) if a+b>0 else None
def gt_parse(s):
    if not s:return None
    s=s.split(':',1)[0];sep='|' if '|' in s else ('/' if '/' in s else None)
    if not sep:return None
    a=s.split(sep)
    if len(a)!=2 or any(x not in ('0','1') for x in a):return None
    return (int(a[0]),int(a[1]))
def frequency(gts):
    ac=an=nc=0
    for s in gts:
        g=gt_parse(s)
        if g is not None:ac+=sum(g);an+=2;nc+=1
    return (ac/an if an else None,nc)
def contributions(p1,p2,p3,og):
    # Deliberately independent scalar implementation: re-code derived frequency
    # from the fixed homozygous outgroup before applying the equations.
    if og==0: d1,d2,d3=p1,p2,p3
    elif og==1: d1,d2,d3=1-p1,1-p2,1-p3
    else:return None
    return (1-d1)*d2*d3,d1*(1-d2)*d3
def window_id(lg,pos):
    s=((pos-1)//250000)*250000+1
    return '%s_%09d_%09d'%(lg,s,s+249999)
def load_context():
    regions={r['lg']:(int(r['start']),int(r['end'])) for r in read_tsv(REGIONS)}
    windows=defaultdict(dict)
    for r in read_tsv(TREES):
        windows[r['lg']][r['window_id']]={'start':int(r['start']),'end':int(r['end']),'midpoint':int(r['midpoint'])}
    primary=read_tsv(AN/'config/frozen_quartets.tsv')
    extra=read_tsv(CFG/'frozen_validation_quartets.tsv')
    designs=[]
    for r in primary:
        r=dict(r);r['design_class']='primary_arrangement_concordant';designs.append(r)
    designs.extend(extra)
    for r in designs:
        for role in ('P1','P2','P3'):
            r[role+'_ids']=r[role+'_sample_ids'].split(',')
        r['O_ids']=[r.get('O_sample_id','Borsai')]
    classes={}
    for lg,mp in windows.items():
        a,b=regions[lg]
        for wid,w in mp.items():
            if w['start']>=a and w['end']<=b:classes[(lg,wid)]='inside'
            elif w['end']<a or w['start']>b:classes[(lg,wid)]='outside'
            else:classes[(lg,wid)]='boundary_overlap'
    return regions,windows,classes,designs

def get_headers(path):
    with gzip.open(str(path),'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):return line.rstrip('\n').split('\t')[9:]
    return []

def choose_ids(d,role,samples):
    frozen=d[role+'_ids'];present=[s for s in frozen if s in samples]
    if len(present)==len(frozen):return present,False
    pop=d[role]
    avail=sorted(s for s in samples if s.startswith(pop))
    return avail[:int(d['n_'+role])],True

def empty_acc():
    return {'records':0,'pass_records':0,'biallelic_snps':0,'nonpass':0,'multiallelic':0,'non_snp':0,
      'ingroup_missing_exclusions':0,'outgroup_missing':0,'outgroup_heterozygous':0,'outgroup_hom_ref':0,'outgroup_hom_alt':0,
      'polarized_sites':0,'abba_informative_sites':0,'baba_informative_sites':0,'abba_only_sites':0,'baba_only_sites':0,'both_pattern_sites':0,
      'sum_abba':0.0,'sum_baba':0.0,'transition_abba':0.0,'transition_baba':0.0,'transversion_abba':0.0,'transversion_baba':0.0,
      'ref_outgroup_abba':0.0,'ref_outgroup_baba':0.0,'alt_outgroup_abba':0.0,'alt_outgroup_baba':0.0,
      'p1_alt_frequency_sum':0.0,'p2_alt_frequency_sum':0.0,'p3_alt_frequency_sum':0.0,'freq_n':0,
      'polymorphic_p1':0,'polymorphic_p2':0,'polymorphic_p3':0,'ingroup_sample_substitution_sites':0}

def add_site(a,abba,baba,ref,alt,og,p):
    a['polarized_sites']+=1;a['sum_abba']+=abba;a['sum_baba']+=baba
    a['abba_informative_sites']+=int(abba>0);a['baba_informative_sites']+=int(baba>0)
    if abba>0 and baba==0:a['abba_only_sites']+=1
    elif baba>0 and abba==0:a['baba_only_sites']+=1
    elif abba>0 and baba>0:a['both_pattern_sites']+=1
    cls='transition' if (ref,alt) in TRANS else 'transversion'
    a[cls+'_abba']+=abba;a[cls+'_baba']+=baba
    a['ref_outgroup_abba' if og==0 else 'alt_outgroup_abba']+=abba
    a['ref_outgroup_baba' if og==0 else 'alt_outgroup_baba']+=baba
    for i,k in enumerate(('p1_alt_frequency_sum','p2_alt_frequency_sum','p3_alt_frequency_sum')):
        a[k]+=p[i]
    a['freq_n']+=1

def update_d(a):
    a['D']=dval(a['sum_abba'],a['sum_baba'])
    for k in ('p1','p2','p3'):
        a[k+'_mean_alt_frequency']=a[k+'_alt_frequency_sum']/a['freq_n'] if a['freq_n'] else None

def scan_site_vcf(path,source,lg_designs,windows,classes,regions,allow_alt_outgroups=False,
                  site_sets=None,window_acc=None,topheap=None):
    samples=get_headers(path);idx={s:i for i,s in enumerate(samples)}
    selected={}
    for d in lg_designs:
        for role in ('P1','P2','P3'):
            selected[(d['design_id'],role)]=choose_ids(d,role,idx.keys())
    out_ids={d['design_id']:[d.get('O_sample_id','Borsai')] for d in lg_designs}
    if allow_alt_outgroups:
        # Outgroup sensitivity is limited to the original primary design, with
        # the same fixed sample calls and all available external taxa.
        for d in lg_designs:
            if d['design_class']=='primary_arrangement_concordant':out_ids[d['design_id']]=['Borsai','Gadoga','Gadcha','Gadmac']
    acc=defaultdict(empty_acc);wacc=defaultdict(empty_acc);position_sets=set();filter_global=Counter()
    heap=defaultdict(list)
    for line in gzip.open(str(path),'rt'):
        if line.startswith('#'):continue
        v=line.rstrip('\n').split('\t');lg=v[0]
        if lg not in windows:continue
        pos=int(v[1]);wid=window_id(lg,pos)
        if wid not in windows[lg]:continue
        reg=classes[(lg,wid)]
        if reg=='boundary_overlap':continue
        a,b=regions[lg]
        is_in=a<=pos<=b
        if (source in ('SNAPP_inside','DSUITE_inside') and not is_in) or (source in ('SNAPP_outside','DSUITE_outside') and is_in):continue
        position_sets.add((lg,pos))
        fclass='records'
        filter_global['records']+=1
        ref,alt=v[3],v[4]
        if v[6] not in ('PASS','.'):filter_global['nonpass']+=1;continue
        if ',' in alt:filter_global['multiallelic']+=1;continue
        if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT':filter_global['non_snp']+=1;continue
        filter_global['biallelic_snps']+=1
        fields=v[8].split(':')
        try:gt_i=fields.index('GT')
        except ValueError:filter_global['no_gt']+=1;continue
        calls=[x.split(':')[gt_i] if gt_i<len(x.split(':')) else '.' for x in v[9:]]
        for d in lg_designs:
            if d['lg'] != lg:
                continue
            key=d['design_id']
            ids=[];swap=False;good=True
            for role in ('P1','P2','P3'):
                x,changed=selected[(key,role)];ids.append(x);swap|=changed
                if not x:good=False
            if not good:continue
            freqs=[];eligible=True
            for group in ids:
                p,n=frequency(calls[idx[s]] for s in group)
                if p is None or n<len(group):eligible=False;break
                freqs.append(p)
            if not eligible:
                for og in out_ids[key]:acc[(key,reg,og)]['ingroup_missing_exclusions']+=1
                continue
            for og in out_ids[key]:
                if og not in idx:continue
                g=gt_parse(calls[idx[og]])
                z=acc[(key,reg,og)];zw=wacc[(key,reg,og,wid)]
                z['records']+=1;z['pass_records']+=1;z['biallelic_snps']+=1;zw['records']+=1;zw['pass_records']+=1;zw['biallelic_snps']+=1
                if swap:z['ingroup_sample_substitution_sites']+=1
                if g is None:
                    z['outgroup_missing']+=1;zw['outgroup_missing']+=1;z['ingroup_missing_exclusions']+=1;continue
                if g[0]!=g[1]:
                    z['outgroup_heterozygous']+=1;zw['outgroup_heterozygous']+=1;continue
                state=g[0]
                z['outgroup_hom_ref' if state==0 else 'outgroup_hom_alt']+=1
                zw['outgroup_hom_ref' if state==0 else 'outgroup_hom_alt']+=1
                abba,baba=contributions(freqs[0],freqs[1],freqs[2],state)
                add_site(z,abba,baba,ref,alt,state,freqs);add_site(zw,abba,baba,ref,alt,state,freqs)
                if site_sets is not None:site_sets.setdefault((key,reg,og),[]).append((pos,wid,abba,baba,(ref,alt),state,freqs))
                score=abs(abba-baba)
                h=heap[(key,reg,og)]
                if len(h)<15: h.append((score,pos,abba,baba,ref,alt,state,freqs))
                else:
                    j=min(range(len(h)),key=lambda n:h[n][0])
                    if score>h[j][0]:h[j]=(score,pos,abba,baba,ref,alt,state,freqs)
    for x in acc.values():update_d(x)
    for x in wacc.values():update_d(x)
    return acc,wacc,heap,position_sets,filter_global,samples

def main():
    RES.mkdir(parents=True,exist_ok=True)
    regions,windows,classes,designs=load_context()
    bylg={lg:[d for d in designs if d['lg']==lg] for lg in LGS}
    sources=[]
    for lg in LGS:
        sources += [(f'SNAPP_inside_{lg}',RAW/f'gadus_morhua_supergenes_{lg.lower()}_snapp.vcf.gz','SNAPP_inside',lg),
                    (f'DSUITE_inside_{lg}',RAW/f'gadus_morhua_supergene_{lg.lower()}_dsuite.vcf.gz','DSUITE_inside',lg),
                    (f'RELATE_inside_{lg}',RAW/f'gadus_morhua_supergene_{lg.lower()}_relate.vcf.gz','RELATE_inside',lg)]
    sources += [('SNAPP_outside',RAW/'gadus_morhua_outside_supergenes_snapp.vcf.gz','SNAPP_outside',None),
                ('RELATE_outside',RAW/'gadus_morhua_outside_supergenes_relate.vcf.gz','RELATE_outside',None),
                ('DSUITE_outside',RAW/'gadus_morhua_outside_supergenes_dsuite.vcf.gz','DSUITE_outside',None)]
    summaries=[];window_rows=[];toprows=[];filter_rows=[];site_sets=None;source_samples={}
    primary_ids={d['design_id'] for d in designs if d['design_class']=='primary_arrangement_concordant'}
    source_results={};source_windows={};source_heaps={};source_positions={}
    for label,path,kind,lg0 in sources:
        if not path.exists():print('SKIP missing',label);continue
        lgs=[lg0] if lg0 else list(LGS)
        dlist=[d for lg in lgs for d in bylg[lg]]
        # Outside Dsuite is a 3.1-GB whole-genome VCF. Restrict it to the
        # four frozen primary designs; the smaller regional files carry all
        # frozen alternative quartets for robustness.
        if kind=='DSUITE_outside':
            dlist=[d for d in dlist if d['design_class']=='primary_arrangement_concordant']
        alt=(kind=='SNAPP_inside' or kind=='SNAPP_outside')
        acc,wacc,heap,posset,fc,samp=scan_site_vcf(path,kind,dlist,windows,classes,regions,
             allow_alt_outgroups=alt,site_sets=site_sets,topheap=source_heaps)
        source_samples[label]=samp
        source_results[label]=acc;source_windows[label]=wacc;source_heaps[label]=heap;source_positions[label]=posset
        for k,v in fc.items():filter_rows.append({'source':label,'filter_metric':k,'count':v,'input_file':path.name})
        for (did,reg,og),a in acc.items():
            update_d(a)
            summaries.append({'source':label,'dataset':kind,'design_id':did,'lg':next(d['lg'] for d in designs if d['design_id']==did),
              'design_class':next(d['design_class'] for d in designs if d['design_id']==did),'region':reg,'outgroup':og,
              **a,'D':a['D'],'source_file':path.name,'sample_ids_P1':','.join(choose_ids(next(d for d in designs if d['design_id']==did),'P1',set(samp))[0]),
              'sample_ids_P2':','.join(choose_ids(next(d for d in designs if d['design_id']==did),'P2',set(samp))[0]),
              'sample_ids_P3':','.join(choose_ids(next(d for d in designs if d['design_id']==did),'P3',set(samp))[0])})
        for (did,reg,og,wid),a in wacc.items():
            update_d(a);lg_d=next(d['lg'] for d in designs if d['design_id']==did)
            meta=windows.get(lg_d,{}).get(wid)
            if meta is None:
                # A malformed/aliased contig-window key must not be silently
                # assigned to another linkage group.
                continue
            window_rows.append({'source':label,'dataset':kind,'design_id':did,'lg':lg_d,
               'region':reg,'outgroup':og,'window_id':wid,'start':meta['start'],'end':meta['end'],'midpoint':meta['midpoint'],
               'sum_ABBA':a['sum_abba'],'sum_BABA':a['sum_baba'],'D_site':a['D'],'polarized_sites':a['polarized_sites'],
               'informative_ABBA':a['abba_informative_sites'],'informative_BABA':a['baba_informative_sites'],
               'n_SNP':a['biallelic_snps'],'outgroup_heterozygous':a['outgroup_heterozygous'],'outgroup_missing':a['outgroup_missing']})
        if kind.startswith('SNAPP'):
            for (did,reg,og),h in heap.items():
                total=sum(abs(float(a['sum_abba'])-float(a['sum_baba'])) for a in [acc.get((did,reg,og),{})])
                # Attach top numerator sites, ranked by absolute per-site contribution.
                d=next(x for x in designs if x['design_id']==did)
                for rank,x in enumerate(sorted(h,reverse=True),1):
                    toprows.append({'source':label,'design_id':did,'lg':d['lg'],'region':reg,'outgroup':og,'rank':rank,
                      'position':x[1],'ABBA':x[2],'BABA':x[3],'numerator_contribution':x[2]-x[3],
                      'absolute_numerator_contribution':x[0],'REF':x[4],'ALT':x[5],'outgroup_state':x[6],
                      'p1_alt_frequency':x[7][0],'p2_alt_frequency':x[7][1],'p3_alt_frequency':x[7][2],
                      'share_of_abs_numerator':x[0]/total if total else None})
        print('scanned',label,'records',fc.get('records',0),flush=True)
        # per-source checkpoint lets a completed scan survive a later error.
        write_tsv(RES/(label+'_checkpoint.tsv'),[{'design_id':k[0],'region':k[1],'outgroup':k[2],**v,'D':v.get('D')} for k,v in acc.items()])
    write_tsv(RES/'variant_contribution_diagnostics.tsv',summaries)
    write_tsv(RES/'matched_window_statistics.tsv',window_rows)
    write_tsv(RES/'top_contributing_sites.tsv',toprows)
    write_tsv(RES/'source_record_filters.tsv',filter_rows)
    run_overlap(sources,source_positions,source_samples,windows,classes)
    run_independent_comparison(summaries,source_results,source_windows,designs,regions,classes,windows)
    run_quartet_summary(summaries,designs)
    run_outgroup_summary(summaries,designs)
    run_spatial_summaries(summaries,source_windows,designs,regions,windows,classes)
    manifest={'created_utc':datetime.now(timezone.utc).isoformat(),'designs':'12 additional quartets frozen from Stage-3 baseline splits plus four primary designs','inputs':[],'coordinate_system':'1-based inclusive gadMor2 linkage-group positions; D windows matched to published 250-kb windows','primary_filters':'PASS/. biallelic SNP; complete genotype calls for chosen individuals; fixed homozygous outgroup; no REF-as-ancestral assumption','no_remote_index':True}
    for _,p,kind,_ in sources:
        if p.exists():
            h=hashlib.md5()
            with p.open('rb') as f:
                for b in iter(lambda:f.read(1<<20),b''):h.update(b)
            manifest['inputs'].append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'md5':h.hexdigest(),'kind':kind})
    (RES/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Validation scans complete.',flush=True)

def run_overlap(sources,position_sets,samples,windows,classes):
    # Site-position overlap of deposited datasets; genotype-complete subsets are
    # separately inferred from VCF call states during the source scan.
    rows=[]
    by={x[0]:x for x in sources}
    for lg in LGS:
        names=['SNAPP_inside_'+lg,'DSUITE_inside_'+lg,'RELATE_inside_'+lg]
        for i,n1 in enumerate(names):
            for n2 in names[i+1:]:
                a=position_sets.get(n1,set());b=position_sets.get(n2,set());aa={p for c,p in a if c==lg};bb={p for c,p in b if c==lg}
                rows.append({'lg':lg,'source_a':n1,'source_b':n2,'n_a':len(aa),'n_b':len(bb),'n_intersection':len(aa&bb),
                    'n_union':len(aa|bb),'jaccard':len(aa&bb)/float(len(aa|bb)) if aa|bb else None,
                    'a_fraction_shared':len(aa&bb)/float(len(aa)) if aa else None,'b_fraction_shared':len(aa&bb)/float(len(bb)) if bb else None,
                    'same_sample_roster':set(samples.get(n1,[]))==set(samples.get(n2,[]))})
        n1='SNAPP_outside';n2='RELATE_outside'
        a=position_sets.get(n1,set());b=position_sets.get(n2,set());aa={p for c,p in a if c==lg};bb={p for c,p in b if c==lg}
        rows.append({'lg':lg,'source_a':n1,'source_b':n2,'n_a':len(aa),'n_b':len(bb),'n_intersection':len(aa&bb),'n_union':len(aa|bb),
           'jaccard':len(aa&bb)/float(len(aa|bb)) if aa|bb else None,'a_fraction_shared':len(aa&bb)/float(len(aa)) if aa else None,
           'b_fraction_shared':len(aa&bb)/float(len(bb)) if bb else None,'same_sample_roster':set(samples.get(n1,[]))==set(samples.get(n2,[]))})
    write_tsv(RES/'variant_overlap_audit.tsv',rows)

def run_independent_comparison(summaries,source_results,source_windows,designs,regions,classes,windows):
    old=read_tsv(AN/'results/inside_outside_summary.tsv');oldmap={r['lg']:r for r in old}
    rows=[]
    for d in designs:
        if d['design_class']!='primary_arrangement_concordant':continue
        lg=d['lg'];did=d['design_id']
        for source in ('SNAPP_inside_'+lg,'SNAPP_outside'):
            for og in ('Borsai','Gadoga','Gadcha','Gadmac'):
                x=source_results.get(source,{}).get((did,'inside' if source.startswith('SNAPP_inside') else 'outside',og))
                if not x:continue
                reg='inside' if source.startswith('SNAPP_inside') else 'outside'
                rows.append({'lg':lg,'design_id':did,'source':source,'outgroup':og,'D_independent_scalar':dval(x['sum_abba'],x['sum_baba']),
                  'sum_ABBA':x['sum_abba'],'sum_BABA':x['sum_baba'],'polarized_sites':x['polarized_sites'],
                  'previous_D':oldmap[lg]['inside_D' if reg=='inside' else 'outside_D'] if og=='Borsai' else None,
                  'absolute_difference_from_previous':abs(dval(x['sum_abba'],x['sum_baba'])-float(oldmap[lg]['inside_D' if reg=='inside' else 'outside_D'])) if og=='Borsai' and dval(x['sum_abba'],x['sum_baba']) is not None else None,
                  'comparison':'Independent direct genotype parsing and scalar frequency estimator; exact old region pooling rule.'})
    write_tsv(RES/'independent_estimator_validation.tsv',rows)

def run_quartet_summary(summaries,designs):
    primary={d['lg']:d for d in designs if d['design_class']=='primary_arrangement_concordant'}
    smap={(r['source'],r['design_id'],r['region'],r['outgroup']):r for r in summaries}
    out=[];bylg=defaultdict(list)
    for d in designs:
        lg=d['lg'];source='SNAPP_inside_'+lg;og='Borsai';di=smap.get((source,d['design_id'],'inside',og),{});do=smap.get(('SNAPP_outside',d['design_id'],'outside',og),{})
        if not di or not do:continue
        expected='ABBA' if d['design_class'] in ('primary_arrangement_concordant','same_arrangement_q3') else ('BABA' if d['design_class']=='alternative_arrangement_q2' else 'none_baseline_pair')
        row={'lg':lg,'design_id':d['design_id'],'design_class':d['design_class'],'P1':d['P1'],'P2':d['P2'],'P3':d['P3'],
          'arrangement_state_pattern':d.get('arrangement_state_pattern',''.join('A' if d[f'P{i}_arrangement']=='ancestral' else 'D' for i in (1,2,3))),
          'expected_arrangement_topology':expected,'inside_D':di['D'],'outside_D':do['D'],
          'inside_minus_outside_D':di['D']-do['D'] if di['D'] is not None and do['D'] is not None else None,
          'inside_ABBA':di['sum_abba'],'inside_BABA':di['sum_baba'],'outside_ABBA':do['sum_abba'],'outside_BABA':do['sum_baba'],
          'inside_informative_sites':di['abba_informative_sites']+di['baba_informative_sites'],'outside_informative_sites':do['abba_informative_sites']+do['baba_informative_sites'],
          'inside_samples':'2/2/2 + Borsai','note':'Different quartet orientations are not averaged as signed D.'}
        out.append(row);bylg[lg].append(row)
    write_tsv(RES/'additional_quartet_designs.tsv',out)
    qsum=[]
    for lg,rr in bylg.items():
        for cls in ('primary_arrangement_concordant','same_arrangement_q3','alternative_arrangement_q2','baseline_pair_control'):
            xs=[r for r in rr if r['design_class']==cls]
            if not xs:continue
            qsum.append({'lg':lg,'design_class':cls,'n_designs':len(xs),'fraction_inside_absD_gt_outside_absD':sum(abs(r['inside_D'])>abs(r['outside_D']) for r in xs)/float(len(xs)),
             'mean_abs_inside_D':statistics.mean(abs(r['inside_D']) for r in xs if r['inside_D'] is not None),
             'expected_side_sign_fraction':sum((r['inside_D'] is not None and ((r['inside_D']>0) if r['expected_arrangement_topology']=='ABBA' else ((r['inside_D']<0) if r['expected_arrangement_topology']=='BABA' else True))) for r in xs)/float(len(xs)),
             'note':'One design per class/LG; descriptive only, quartets share populations and are not independent replicates.'})
    write_tsv(RES/'quartet_robustness.tsv',qsum)

def run_outgroup_summary(summaries,designs):
    rows=[];primary=[d for d in designs if d['design_class']=='primary_arrangement_concordant']
    smap={(r['source'],r['design_id'],r['region'],r['outgroup']):r for r in summaries}
    for d in primary:
        for reg,source in (('inside','SNAPP_inside_'+d['lg']),('outside','SNAPP_outside')):
            for og in ('Borsai','Gadoga','Gadcha','Gadmac'):
                x=smap.get((source,d['design_id'],reg,og))
                if not x:continue
                rows.append({'lg':d['lg'],'design_id':d['design_id'],'region':reg,'outgroup':og,'outgroup_sample_count':1,
                  'n_fixed_ref':x['outgroup_hom_ref'],'n_fixed_alt':x['outgroup_hom_alt'],'n_heterozygous':x['outgroup_heterozygous'],
                  'n_missing':x['outgroup_missing'],'n_polarized':x['polarized_sites'],'D':x['D'],
                  'callability':'VCF-filtered globally complete sites; sitewise homozygous polarization only'})
    write_tsv(RES/'outgroup_sensitivity.tsv',rows)

def run_spatial_summaries(summaries,source_windows,designs,regions,windows,classes):
    # Spatial effect summaries from both filtered complete-call and missing-call
    # Dsuite sources. Delete-block uncertainty uses ABBA/BABA pooled contributions.
    out=[];sw={}
    for src, per_source in source_windows.items():
        for (did,reg,og,wid),v in per_source.items():
            sw.setdefault((src,did,og),{})[wid]=v
    for d in designs:
        for label,src_in,src_out in (
                ('SNAPP','SNAPP_inside_'+d['lg'],'SNAPP_outside'),
                ('DSUITE','DSUITE_inside_'+d['lg'],'DSUITE_outside')):
            a=dict(sw.get((src_in,d['design_id'],'Borsai'),{}))
            a.update(sw.get((src_out,d['design_id'],'Borsai'),{}))
            if not a:continue
            totals={r:[0.,0.,0] for r in ('inside','outside')}
            blocks=defaultdict(lambda:{'inside':[0.,0.,0],'outside':[0.,0.,0]})
            for wid,v in a.items():
                lg=d['lg'];reg=classes[(lg,wid)]
                if reg not in totals:continue
                start=windows[lg][wid]['start'];block=(start-1)//250000
                x=totals[reg];x[0]+=v['sum_abba'];x[1]+=v['sum_baba'];x[2]+=v['polarized_sites']
                y=blocks[block][reg];y[0]+=v['sum_abba'];y[1]+=v['sum_baba'];y[2]+=v['polarized_sites']
            di=dval(totals['inside'][0],totals['inside'][1]);do=dval(totals['outside'][0],totals['outside'][1])
            if di is None or do is None:continue
            out.append({'lg':d['lg'],'design_id':d['design_id'],'source':label,'inside_D':di,'outside_D':do,'inside_minus_outside_D':di-do,
              'inside_polarized_sites':totals['inside'][2],'outside_polarized_sites':totals['outside'][2],
              'inside_windows':sum(1 for w in windows[d['lg']] if classes[(d['lg'],w)]=='inside'),
              'outside_windows':sum(1 for w in windows[d['lg']] if classes[(d['lg'],w)]=='outside'),
              'boundary_overlap_excluded':sum(1 for w in windows[d['lg']] if classes[(d['lg'],w)]=='boundary_overlap')})
            for size in (250000,500000,1000000):
                blockvals=defaultdict(lambda:{'inside':[0.,0.],'outside':[0.,0.]})
                for wid,v in a.items():
                    reg=classes[(d['lg'],wid)]
                    if reg in ('inside','outside'):
                        bl=(windows[d['lg']][wid]['start']-1)//size
                        blockvals[bl][reg][0]+=v['sum_abba'];blockvals[bl][reg][1]+=v['sum_baba']
                bs=sorted(blockvals);jack=[]
                for b0 in bs:
                    x={r:(totals[r][0],totals[r][1]) for r in ('inside','outside')}
                    for r in ('inside','outside'):
                        x[r]=(x[r][0]-blockvals[b0][r][0],x[r][1]-blockvals[b0][r][1])
                    d1=dval(*x['inside']);d2=dval(*x['outside'])
                    if d1 is not None and d2 is not None:jack.append(d1-d2)
                mean=sum(jack)/len(jack) if jack else 0
                se=math.sqrt((len(jack)-1)/len(jack)*sum((z-mean)**2 for z in jack)) if len(jack)>1 else None
                out.append({'lg':d['lg'],'design_id':d['design_id'],'source':label,'block_bp':size,'n_blocks':len(bs),
                    'inside_minus_outside_D_jackknife_SE':se,'CI_low':di-do-1.96*se if se is not None else None,
                    'CI_high':di-do+1.96*se if se is not None else None})
    write_tsv(RES/'inside_outside_validation.tsv',out)
    write_tsv(RES/'spatial_robustness.tsv',out)

def finalize_existing():
    """Rebuild post-scan summaries after an interrupted postprocessing step."""
    regions,windows,classes,designs=load_context()
    summaries=read_tsv(RES/'variant_contribution_diagnostics.tsv')
    for r in summaries:
        if 'outgroup_sample_swaps' in r:
            r['ingroup_sample_substitution_sites']=r.pop('outgroup_sample_swaps')
    write_tsv(RES/'variant_contribution_diagnostics.tsv',summaries)
    if summaries:
        # Record substitutions explicitly: the Dsuite VCF has Gadmor_lfo3
        # where the primary SNAPP quartet uses Gadmor_lfo2.
        sample_rows=[]
        for r in summaries:
            if r['source'].startswith('DSUITE') and r['design_id'].endswith('_PRIMARY'):
                sample_rows.append({'source':r['source'],'lg':r['lg'],'design_id':r['design_id'],
                    'P1':r['sample_ids_P1'],'P2':r['sample_ids_P2'],'P3':r['sample_ids_P3'],
                    'ingroup_sample_substitution_sites':r.get('ingroup_sample_substitution_sites',r.get('outgroup_sample_swaps','NA')),
                    'note':'Source roster substitution is explicit in sample IDs; only Gadmor_lfo2 -> Gadmor_lfo3 occurs in these primary designs.'})
        write_tsv(RES/'source_sample_mapping.tsv',sample_rows)
    source_windows=defaultdict(dict)
    for r in read_tsv(RES/'matched_window_statistics.tsv'):
        source_windows[r['source']][(r['design_id'],r['region'],r['outgroup'],r['window_id'])]={
            'sum_abba':float(r['sum_ABBA']),'sum_baba':float(r['sum_BABA']),
            'polarized_sites':int(r['polarized_sites'])}
    run_spatial_summaries(summaries,source_windows,designs,regions,windows,classes)
    manifest={'created_utc':datetime.now(timezone.utc).isoformat(),
      'designs':'12 additional quartets frozen from Stage-3 baseline splits plus four primary designs',
      'inputs':[],'coordinate_system':'1-based inclusive gadMor2 linkage-group positions; D windows matched to published 250-kb windows',
      'primary_filters':'PASS/. biallelic SNP; complete genotype calls for chosen individuals; fixed homozygous outgroup; no REF-as-ancestral assumption',
      'no_remote_index':True,'postprocessing':'Rebuilt from completed per-source scan summaries; see per-source checkpoints'}
    for p in sorted(RAW.glob('*.vcf.gz')):
        h=hashlib.md5()
        with p.open('rb') as f:
            for b in iter(lambda:f.read(1<<20),b''):h.update(b)
        manifest['inputs'].append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'md5':h.hexdigest()})
    (RES/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('Post-scan summaries rebuilt from completed local scans.')

if __name__=='__main__':
    if '--finalize-existing' in sys.argv:finalize_existing()
    else:main()
