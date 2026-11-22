"""Exploratory cod Tajima D from deposited Relate SNP VCFs, population by LG."""
from pathlib import Path
from collections import defaultdict
import gzip,csv,sys,json
import numpy as np,pandas as pd,matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[5]; BASE=ROOT/'empirical/atlantic_cod/analysis/tajima_d'; DATA=ROOT/'data/atlantic_cod'
sys.path.insert(0,str(ROOT/'empirical/analysis/tajimas_d')); from tajimas import tajima_from_counts
D=DATA/'raw/zenodo/local_downloads'; TREE=pd.read_csv(DATA/'processed/cod_window_trees.tsv',sep='\t'); ARR=pd.read_csv(DATA/'processed/population_arrangements.tsv',sep='\t').set_index('canonical_population'); REG=pd.read_csv(DATA/'metadata/region_manifest.tsv',sep='\t').set_index('lg')
POP={'avc':'More_stationary','avo':'More_migratory','bat':'Labrador','bor':'Bornholm_Basin','icc':'Iceland_coastal','ico':'Iceland_oceanic','kie':'Kiel_Bight','lfc':'Lofoten_coastal','lfo':'Lofoten_oceanic','low':'Lowestoft','twc':'Trondheim_coastal','two':'Trondheim_oceanic'}

def main():
    windows={lg:df.sort_values('start').to_dict('records') for lg,df in TREE[TREE.lg.isin(['LG01','LG02','LG07','LG12'])].groupby('lg')}
    acc=defaultdict(list); candidate_counts=defaultdict(int); audit=[]; records=0; sampleinfo={}
    files=[(lg,'inside',D/f'gadus_morhua_supergene_{lg.lower()}_relate.vcf.gz') for lg in ['LG01','LG02','LG07','LG12']]+[('outside','outside',D/'gadus_morhua_outside_supergenes_relate.vcf.gz')]
    for lgpart,region,path in files:
        with gzip.open(path,'rt') as f:
            ids=[]
            for line in f:
                if line.startswith('#CHROM'):
                    ids=line.rstrip('\n').split('\t')[9:];break
            if not ids: raise RuntimeError('no VCF sample header: '+str(path))
            sampleinfo[str(path.name)]=ids
            groups={code:[i for i,s in enumerate(ids) if s.startswith('Gadmor_'+code)] for code in POP}
            fmt_seen=set(); rec_ct=0; miss_ct=0; mult_ct=0; nonpass=0
            for line in f:
                if line.startswith('#'):continue
                z=line.rstrip('\n').split('\t'); chrom=z[0]; pos=int(z[1]); llg=chrom if chrom.startswith('LG') else ''
                if lgpart!='outside' and llg!=lgpart: continue
                if lgpart=='outside' and llg not in ['LG01','LG02','LG07','LG12']: continue
                # Keep only the source's outside-supergene records; require the marker's midpoint to be outside frozen interval.
                b=int((pos-1)//250000); wlist=[w for w in windows.get(llg,[]) if int(w['start'])<=pos<=int(w['end'])]
                if not wlist: continue
                w=wlist[0]
                inside=(str(w['inside_inversion']).lower()=='true')
                if (region=='inside') != inside: continue
                if z[6] not in ('PASS','.'): nonpass+=1; continue
                als=[z[3]]+z[4].split(',')
                if len(als)!=2 or any(x not in ('A','C','G','T') for x in als): mult_ct+=1;continue
                fm=z[8].split(':'); gt_ix=fm.index('GT') if 'GT' in fm else -1
                if gt_ix<0:continue
                gts=[]; failed=False
                for j in range(9,len(z)):
                    gt=z[j].split(':')[gt_ix].replace('|','/').split('/')
                    if len(gt)!=2 or any(x=='.' for x in gt): failed=True;break
                    if any(x not in ('0','1') for x in gt):failed=True;break
                    gts.append([int(x) for x in gt])
                if failed:miss_ct+=1;continue
                rec_ct+=1; records+=1; fmt_seen.add(z[8])
                for code,ix in groups.items():
                    if len(ix)<2:continue
                    candidate_counts[(llg,w['window_id'],code)]+=1
                    alleles=[a for j in ix for a in gts[j]]; c0=alleles.count(0);c1=alleles.count(1)
                    if c0 and c1: acc[(llg,w['window_id'],code)].append([c0,c1])
            audit.append({'source_file':path.name,'subset':region,'LG':lgpart,'n_sample_columns':len(ids),'n_population_groups_n_ge_2':sum(len(v)>=2 for v in groups.values()),'n_records_used':rec_ct,'n_missing_or_malformed_records':miss_ct,'n_multiallelic_or_non_ACGT':mult_ct,'n_nonpass':nonpass,'FORMAT_seen':';'.join(sorted(fmt_seen)),'samples':';'.join(ids)})
    rows=[]
    for lg,ww in windows.items():
      for w in ww:
       inside=str(w['inside_inversion']).lower()=='true'
       for code in POP:
        ids=sampleinfo.get('gadus_morhua_supergene_'+lg.lower()+'_relate.vcf.gz',sampleinfo.get('gadus_morhua_outside_supergenes_relate.vcf.gz',[]))
        nind=sum(s.startswith('Gadmor_'+code) for s in ids); n=2*nind
        result=tajima_from_counts(acc[(lg,w['window_id'],code)],n) if n>=4 else {'S':0,'pi_total':0,'theta_W':None,'tajimas_D':None,'D_undefined_reason':'n_haplotypes<=3','n_sites_examined':0}
        canon=POP[code]; state=ARR.loc[canon,lg+'_state'] if canon in ARR.index else 'unmapped'
        state_conf=ARR.loc[canon,lg+'_confidence'] if canon in ARR.index else 'unmapped'
        rows.append({'LG':lg,'window_id':w['window_id'],'start':w['start'],'end':w['end'],'midpoint':w['midpoint'],'window_bp':w['window_size_bp'],'region_class':'inside_supergene' if inside else 'outside','inversion_id':w['inversion_id'],'population_code':code,'population_label':POP[code],'population_arrangement_state':state,'population_arrangement_confidence':state_conf,'n_individuals':nind,'n_haplotypes':n,'ploidy':'diploid; phase ignored for site allele counts','S':result['S'],'pi_total':result['pi_total'],'theta_W':result['theta_W'],'tajimas_D':result['tajimas_D'],'D_undefined_reason':result['D_undefined_reason'],'n_sites_examined':candidate_counts[(lg,w['window_id'],code)],'n_segregating_records_in_population':result['S'],'n_multiallelic_excluded_per_window':'','multiallelic_note':'source VCF prefiltered biallelic; per-window discarded multiallelic count unavailable','n_missing_or_malformed_retained_records':0,'callable_bases':'','per_base_unavailable_reason':'variant-only VCF; no invariant callable mask','ascertainment':'Beagle GT VCF produced from masked GATK calls; source split filters biallelic sites and removes any missing genotype; D is conditional on retained SNP callset, not established complete for segregating sites'})
    out=pd.DataFrame(rows); BASE.joinpath('results').mkdir(parents=True,exist_ok=True);out.to_csv(BASE/'results/atlantic_cod_tajima_d_250kb.tsv',sep='\t',index=False,na_rep='NA')
    pd.DataFrame(audit).drop(columns=['samples']).to_csv(BASE/'results/atlantic_cod_tajima_d_source_audit.tsv',sep='\t',index=False,na_rep='NA')
    summaries=[];rng=np.random.default_rng(20261121)
    for (lg,code),df in out.groupby(['LG','population_code']):
        for m in ['tajimas_D']:
            a=df[df.region_class=='inside_supergene'].dropna(subset=[m]); b=df[df.region_class=='outside'].dropna(subset=[m]);a=a.copy();b=b.copy();a['block']=a.midpoint//1_000_000;b['block']=b.midpoint//1_000_000
            av=a.groupby('block')[m].mean().values;bv=b.groupby('block')[m].mean().values;boot=[]
            if len(av)>1 and len(bv)>1:
                for _ in range(3000):boot.append(rng.choice(av,len(av),replace=True).mean()-rng.choice(bv,len(bv),replace=True).mean())
            summaries.append({'LG':lg,'population_code':code,'population':POP[code],'inside_mean':a[m].mean(),'outside_mean':b[m].mean(),'difference':a[m].mean()-b[m].mean(),'n_inside_windows':len(a),'n_outside_windows':len(b),'n_inside_1Mb_blocks':len(av),'n_outside_1Mb_blocks':len(bv),'block_95_low':np.quantile(boot,.025) if boot else np.nan,'block_95_high':np.quantile(boot,.975) if boot else np.nan,'status':'exploratory; no cross-population pooling'})
    pd.DataFrame(summaries).to_csv(BASE/'results/atlantic_cod_tajima_d_inside_outside.tsv',sep='\t',index=False,na_rep='NA')
    fig,axs=plt.subplots(4,1,figsize=(12,9),sharex=False,constrained_layout=True)
    for ax,lg in zip(axs,['LG01','LG02','LG07','LG12']):
        dd=out[(out.LG==lg)&(out.population_code=='avc')].sort_values('midpoint'); ax.axvspan(int(REG.loc[lg,'start'])/1e6,int(REG.loc[lg,'end'])/1e6,color='#e9a743',alpha=.15);ax.scatter(dd.midpoint/1e6,dd.tajimas_D,s=12,color='#29627d');ax.axhline(0,color='.5',lw=.6);ax.set_ylabel(lg+"\nTajima's D\nMore stationary");ax.grid(axis='y',alpha=.2)
    axs[-1].set_xlabel('Linkage-group coordinate (Mb; 250-kb source window axis)');fig.suptitle("Atlantic cod population-specific Tajima's D\nInside intervals only: local outside Relate VCF has no LG01/02/07/12 records; exploratory masked-SNP VCF",fontsize=11)
    for ext in ('pdf','svg','png'):fig.savefig(BASE/'figures'/f'atlantic_cod_tajima_d.{ext}',dpi=300)
    print('records_used',records,'rows',len(rows),'D_estimable',int(out.tajimas_D.notna().sum()));print(pd.DataFrame(summaries).to_string(index=False))
if __name__=='__main__':main()
