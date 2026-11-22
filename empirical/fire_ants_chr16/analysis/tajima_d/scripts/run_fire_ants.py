from pathlib import Path
import sys, csv, math
from collections import defaultdict
import gzip
import numpy as np, pandas as pd, matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[5]
sys.path.insert(0,str(ROOT/'empirical/analysis/tajimas_d'))
from tajimas import tajima_from_counts, weighted_delta
BASE=ROOT/'empirical/fire_ants_chr16/analysis/tajima_d'; DATA=ROOT/'data/fire_ants_chr16'
VCF=DATA/'raw/genotypes/gt.vcf.gz'; MAP=ROOT/'empirical/fire_ants_chr16/analysis/patterson_d/results/window_gene_coordinate_mapping.tsv'; MEM=ROOT/'empirical/fire_ants_chr16/analysis/patterson_d/config/frozen_sample_membership.tsv'; INV=ROOT/'empirical/fire_ants_chr16/analysis/patterson_d/results/vcf_sample_inventory.tsv'

def main():
    mp=pd.read_csv(MAP,sep='\t'); mp=mp[(mp.tree_chromosome=='chr16') & mp.vcf_contig_present.eq(1)].copy()
    win=mp.drop_duplicates('window_index')[['window_index','tree_window_start','tree_window_end','tree_window_mid','supplement_genomic_region']].sort_values('window_index')
    memb=pd.read_csv(MEM,sep='\t'); inv=pd.read_csv(INV,sep='\t')
    q=memb[memb.quartet_id=='Q1_primary_Argentina'].merge(inv,left_on='sample_id',right_on='vcf_sample_id')
    q=q[(q.species=='S.invicta')&(q.caste=='Male')]
    sb=q[q.supergene_variant=='SB'].vcf_sample_id.tolist(); sb2=q[q.supergene_variant=='Sb'].vcf_sample_id.tolist()
    # Parse the existing BGZF VCF directly so the pipeline has no pysam dependency.
    ids=[]
    with gzip.open(VCF,'rt') as f:
        for line in f:
            if line.startswith('#CHROM'):
                ids=line.rstrip('\n').split('\t')[9:]; break
    ix1=[ids.index(s) for s in sb]; ix2=[ids.index(s) for s in sb2]
    if (len(sb),len(sb2))!=(78,30): raise ValueError(f'expected 78/30, got {len(sb)}/{len(sb2)}')
    bycontig=defaultdict(list)
    for r in mp.itertuples(): bycontig[r.vcf_contig].append((int(r.vcf_gene_start_1based),int(r.vcf_gene_end_1based_inclusive),int(r.window_index)))
    site_counts=defaultdict(lambda:defaultdict(list)); matched_counts=defaultdict(lambda:defaultdict(list)); examined_counts=defaultdict(lambda:defaultdict(int)); diag=defaultdict(lambda:defaultdict(int)); nrec=0; n_target=0
    with gzip.open(VCF,'rt') as f:
      for line in f:
        if line.startswith('#'): continue
        fields=line.rstrip('\n').split('\t')
        if len(fields)<10: continue
        contig,pos=fields[0],int(fields[1]); spans=bycontig.get(contig)
        if not spans: continue
        ws={w for lo,hi,w in spans if lo<=pos<=hi}
        if not ws: continue
        n_target+=1
        alleles=[fields[3]]+fields[4].split(',')
        if len(alleles)!=2 or any(a not in ('A','C','G','T') for a in alleles):
            for w in ws: diag[w]['excluded_non_biallelic_or_non_ACGT']+=1
            continue
        fmt=fields[8].split(':'); gi=fmt.index('GT') if 'GT' in fmt else -1; di=fmt.index('DP') if 'DP' in fmt else -1
        if gi<0 or di<0:
            for w in ws: diag[w]['excluded_missing_GT_or_DP_format']+=1
            continue
        vals={}
        valid=True
        for label,indices in [('SB',ix1),('Sb',ix2)]:
            a=[]; failed=False
            for w in ws: diag[w][label+'_genotype_opportunities']+=len(indices)
            for j in indices:
                parts=fields[9+j].split(':')
                gttext=parts[gi] if gi<len(parts) else '.'; dptext=parts[di] if di<len(parts) else '.'
                gt=gttext.replace('|','/').split('/')
                if len(gt)!=2 or any(v=='.' for v in gt):
                    for w in ws: diag[w][label+'_missing_or_partial_GT']+=1
                    failed=True; continue
                if gt[0]!=gt[1]:
                    for w in ws: diag[w][label+'_heterozygous_GT']+=1
                    failed=True; continue
                try: dp=int(dptext)
                except ValueError: dp=0
                if dp<3:
                    for w in ws: diag[w][label+'_DP_lt_3']+=1
                    failed=True; continue
                try: a.append(int(gt[0]))
                except ValueError: failed=True; break
            if not failed: vals[label]=a
        groups=[]
        if 'SB' in vals: groups.append(('SB',vals['SB']))
        if 'Sb' in vals: groups.append(('Sb',vals['Sb']))
        if len(vals)==2: groups.append(('pooled',vals['SB']+vals['Sb']))
        for label,a in groups:
            for w in ws: examined_counts[w][label]+=1
            if len(set(a))<2: continue
            row=[a.count(0),a.count(1)]
            for w in ws:
                site_counts[w][label].append(row); diag[w][label+'_complete_variant_sites']+=1
        # Use the same sites, complete across all 108 samples, for the pooled
        # versus SB/Sb contrast so it is not driven by group-specific masks.
        if 'SB' in vals and 'Sb' in vals:
            for label,a in [('pooled',vals['SB']+vals['Sb']),('SB',vals['SB']),('Sb',vals['Sb'])]:
                row=[a.count(0),a.count(1)]
                for w in ws: matched_counts[w][label].append(row)
        for w in ws: diag[w]['target_biallelic_records']+=1
    rows=[]
    for r in win.itertuples():
        w=int(r.window_index); row={'window_index':w,'start':r.tree_window_start,'end':r.tree_window_end,'mid':r.tree_window_mid,'interval_class':r.supplement_genomic_region,'coordinate_system':'local chr16 tree-window coordinate; BUSCO VCF contig intervals; no chromosome-scale transform','n_SB':78,'n_Sb':30,'ploidy':'male haploid; homozygous diploid GT recoded as one allele','ascertainment':'variant-only FreeBayes VCF, four selected BUSCO gene bodies/window; observed-panel conditional; not validated as complete segregating-site record'}
        for label in ['pooled','SB','Sb']:
            result=tajima_from_counts(site_counts[w][label],{'pooled':108,'SB':78,'Sb':30}[label])
            for k,v in result.items():
                if k=='n_sites_examined': row[label+'_n_segregating_records_used_for_D']=v
                else: row[label+'_'+k]=v
            row[label+'_callable_bases']=''; row[label+'_per_base_unavailable_reason']='no invariant/callable base record; gene-body SNP panel only'
            row[label+'_n_biallelic_records_examined_complete_group']=examined_counts[w][label]
        d,wd=weighted_delta(row['pooled_tajimas_D'],row['SB_tajimas_D'],row['Sb_tajimas_D'],78,30)
        row['Delta_D_pooled_minus_background_mean']=d; row['Delta_D_sample_weighted']=wd
        matched={g:tajima_from_counts(matched_counts[w][g],{'pooled':108,'SB':78,'Sb':30}[g]) for g in ['pooled','SB','Sb']}
        row['n_sites_matched_complete_108']=matched['pooled']['n_sites_examined']
        for g in ['pooled','SB','Sb']:
            row['matched_'+g+'_S']=matched[g]['S']; row['matched_'+g+'_pi_total']=matched[g]['pi_total']; row['matched_'+g+'_theta_W']=matched[g]['theta_W']; row['matched_'+g+'_tajimas_D']=matched[g]['tajimas_D']; row['matched_'+g+'_undefined_reason']=matched[g]['D_undefined_reason']
        md,mwd=weighted_delta(row['matched_pooled_tajimas_D'],row['matched_SB_tajimas_D'],row['matched_Sb_tajimas_D'],78,30)
        row['matched_Delta_D_unweighted']=md; row['matched_Delta_D_weighted']=mwd
        for k,v in diag[w].items(): row['audit_'+k]=v
        for label in ['SB','Sb']:
            denom=diag[w][label+'_genotype_opportunities']
            row[label+'_missing_or_partial_GT_fraction']=diag[w][label+'_missing_or_partial_GT']/denom if denom else np.nan
            row[label+'_fraction_DP_lt_3']=diag[w][label+'_DP_lt_3']/denom if denom else np.nan
            row[label+'_fraction_heterozygous_GT']=diag[w][label+'_heterozygous_GT']/denom if denom else np.nan
        pooled_denom=diag[w]['SB_genotype_opportunities']+diag[w]['Sb_genotype_opportunities']
        row['pooled_missing_or_partial_GT_fraction']=(diag[w]['SB_missing_or_partial_GT']+diag[w]['Sb_missing_or_partial_GT'])/pooled_denom if pooled_denom else np.nan
        rows.append(row)
    out=pd.DataFrame(rows); BASE.joinpath('results').mkdir(parents=True,exist_ok=True); out.to_csv(BASE/'results/fire_ant_chr16_tajima_d_windows.tsv',sep='\t',index=False,na_rep='NA')
    # Regional summaries: window means, with 1-Mb block bootstrap over local tree coordinate.
    out['inside']=out.interval_class.eq('chr16nr'); out['block']=(pd.to_numeric(out.mid)//1_000_000).astype(int)
    metrics=['pooled_tajimas_D','SB_tajimas_D','Sb_tajimas_D','Delta_D_pooled_minus_background_mean','Delta_D_sample_weighted','matched_Delta_D_unweighted','matched_Delta_D_weighted']
    summaries=[]; rng=np.random.default_rng(20261121)
    for metric in metrics:
        inside=out[out.inside & out[metric].notna()]; outside=out[(~out.inside)&out[metric].notna()]
        a=inside.groupby('block')[metric].mean().values; b=outside.groupby('block')[metric].mean().values
        delta=inside[metric].mean()-outside[metric].mean(); boots=[]
        if len(a)>1 and len(b)>1:
            for _ in range(10000): boots.append(rng.choice(a,len(a),replace=True).mean()-rng.choice(b,len(b),replace=True).mean())
        lo,hi=np.quantile(boots,[.025,.975]) if boots else (np.nan,np.nan)
        summaries.append({'statistic':metric,'n_inside_windows':len(inside),'n_outside_windows':len(outside),'inside_mean':inside[metric].mean(),'outside_mean':outside[metric].mean(),'difference_inside_minus_outside':delta,'block_bp':1_000_000,'n_inside_blocks':len(a),'n_outside_blocks':len(b),'block_bootstrap_95_low':lo,'block_bootstrap_95_high':hi,'interpretation':'exploratory observed-panel contrast; not neutrality-test p-value'})
    pd.DataFrame(summaries).to_csv(BASE/'results/fire_ant_chr16_tajima_d_inside_outside.tsv',sep='\t',index=False,na_rep='NA')
    # Existing diversity and tree tracks are joined on the same frozen window_index.
    div=pd.read_csv(ROOT/'empirical/fire_ants_chr16/window_diversity/results/sinvicta_chr16_tree_diversity_join.tsv',sep='\t')
    joined=out.merge(div,on='window_index',how='left',suffixes=('','_existing'))
    joined.to_csv(BASE/'results/fire_ant_chr16_tajima_d_joined_tracks.tsv',sep='\t',index=False,na_rep='NA')
    # Figure: preserve gaps; show D tracks and background contrast.
    fig,axs=plt.subplots(4,1,figsize=(12,9),sharex=True,constrained_layout=True)
    for ax,col,label in zip(axs,['pooled_tajimas_D','SB_tajimas_D','Sb_tajimas_D','Delta_D_pooled_minus_background_mean'],['Pooled Tajima\'s D (n=108)','SB Tajima\'s D (n=78)','Sb Tajima\'s D (n=30)','Pooled minus mean background D']):
        ax.axvspan(11.680438,27.917498,color='#e7a34a',alpha=.16,lw=0); ax.scatter(out.mid/1e6,out[col],s=14,color='#28627b',alpha=.82); ax.axhline(0,color='.5',lw=.6); ax.set_ylabel(label); ax.grid(axis='y',alpha=.2)
    axs[-1].set_xlabel('Local chr16 tree-window coordinate (Mb; BUSCO mapped VCF intervals)'); fig.suptitle("S. invicta Argentine males: Tajima's D on observed BUSCO SNP panel\nExploratory; variant-only and gene-ascertained; no callable-base denominator",fontsize=12)
    for ext in ('pdf','png'): fig.savefig(BASE/'figures'/f'fire_ant_chr16_tajima_d.{ext}',dpi=300)
    audit={'vcf_records_in_target_intervals':n_target,'groups':{'SB':78,'Sb':30,'pooled':108},'windows':len(rows),'windows_with_D':{g:int(out[g+'_tajimas_D'].notna().sum()) for g in ['pooled','SB','Sb']},'mapping':'exact BUSCO gene-body intervals to existing local chr16 tree windows; overlapping intervals deduplicated by window per site','site_filter':'biallelic A/C/G/T; DP>=3; male diploid GT must be homozygous and is recoded haploid; require complete calls in the entire calculation group at the site','ascertainment':'variant-only FreeBayes callset; records include target gene bodies only and upstream source filter completeness is not demonstrated; outputs are observed-panel exploratory D, not asserted unbiased classical SFS','callability':'no invariant/callable positions'}
    import json; (BASE/'results/fire_ant_chr16_tajima_d_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print(json.dumps(audit,indent=2)); print(pd.DataFrame(summaries).to_string(index=False))
if __name__=='__main__': main()
