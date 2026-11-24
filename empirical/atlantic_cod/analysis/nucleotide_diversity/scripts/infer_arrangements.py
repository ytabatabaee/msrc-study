#!/usr/bin/env python3
"""Inspect individual VCF genotypes and exploratory inversion-region PCA.

PCA/GMM are diagnostic only; arrangement labels are never assigned from cluster
number. Published homokaryotype selection and documented population/ecotype
arrangement assignments are kept separate from PCA/GMM diagnostics.
"""
import csv, gzip, re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture

ROOT=Path(__file__).resolve().parents[5]
BASE=ROOT/'empirical/atlantic_cod/analysis/nucleotide_diversity'
VCFDIR=ROOT/'data/atlantic_cod/raw/zenodo/local_downloads'
POPFILE=ROOT/'data/atlantic_cod/processed/population_arrangements.tsv'
OUT=BASE/'results'; FIG=BASE/'figures'
LGFILES={lg:VCFDIR/f'gadus_morhua_supergenes_{lg.lower()}_snapp.vcf.gz' for lg in ['LG01','LG02','LG07','LG12']}


def sample_pop(s):
    if s.startswith('Gadmor_'):
        return s[:10]
    return 'non_cod:'+s


def read_pop_states():
    with POPFILE.open() as f:
        return {r['tree_label'].replace('_spc',''):r for r in csv.DictReader(f,delimiter='\t')}


def parse_gt(field):
    gt=field.split(':',1)[0]
    if gt in ('.','./.','.|.') or '.' in gt: return None
    sep='|' if '|' in gt else '/'
    alleles=gt.split(sep)
    if len(alleles)!=2 or any(a not in ('0','1') for a in alleles): return None
    return int(alleles[0])+int(alleles[1])


def analyze(lg,path,popstates):
    samples=[]; rows=[]; raw=qual=0; raw_missing=qual_missing=0; gq_vals=[]; dp_vals=[]; siteq=[]
    with gzip.open(path,'rt') as fh:
        for line in fh:
            if line.startswith('#CHROM'):
                h=line.rstrip().split('\t'); samples=h[9:]; continue
            if line.startswith('#'): continue
            a=line.rstrip().split('\t')
            if a[0]!=lg or len(a[4].split(','))!=1 or len(a[3])!=1 or len(a[4])!=1: continue
            raw+=1; fmt=a[8].split(':'); vals=[]; gqs=[]; dps=[]
            for sample in a[9:]:
                parts=sample.split(':'); d=dict(zip(fmt,parts)); gt=parse_gt(sample)
                if gt is None: raw_missing+=1
                if d.get('GQ','.') not in ('.',''):
                    try:gqs.append(float(d['GQ']))
                    except ValueError:pass
                if d.get('DP','.') not in ('.',''):
                    try:dps.append(float(d['DP']))
                    except ValueError:pass
                vals.append((gt,float(d['GQ']) if d.get('GQ','.') not in ('.','') else np.nan,float(d['DP']) if d.get('DP','.') not in ('.','') else np.nan))
            if a[6] not in ('PASS','.') or (a[5] not in ('.','') and float(a[5])<30): continue
            qual+=1; siteq.append(float(a[5]) if a[5] not in ('.','') else np.nan)
            gq_vals.extend(gqs);dp_vals.extend(dps)
            call=[]
            for gt,gq,dp in vals:
                if gt is None or not np.isfinite(gq) or gq<20 or not np.isfinite(dp) or dp<8:
                    call.append(np.nan); qual_missing+=1
                else: call.append(gt)
            rows.append(call)
    X=np.asarray(rows,dtype=float) if rows else np.empty((0,len(samples)))
    # Keep Atlantic cod only for PCA; outgroups are not arrangement candidates.
    cod_idx=[i for i,s in enumerate(samples) if s.startswith('Gadmor_')]
    Xc=X[:,cod_idx]
    if len(Xc):
        rate=np.mean(np.isfinite(Xc),axis=1); af=np.nanmean(Xc,axis=1)/2
        maf=np.minimum(af,1-af); keep=(rate>=.9)&(maf>=.02)
        M=Xc[keep]
        if M.shape[0]>=3:
            means=np.nanmean(M,axis=1); M=np.where(np.isfinite(M),M,means[:,None])
            sd=M.std(axis=1); use=sd>0; Z=((M[use]-M[use].mean(axis=1)[:,None])/M[use].std(axis=1)[:,None]).T
            pca=PCA(n_components=min(5,len(cod_idx)-1,Z.shape[0]),svd_solver='full')
            scores=pca.fit_transform(Z); exp=pca.explained_variance_ratio_
        else: scores=np.empty((len(cod_idx),0)); exp=[]
    else: keep=[];M=np.empty((0,0));scores=np.empty((0,0));exp=[]
    # Candidate mixtures are descriptive exploratory partitions, never karyotype labels.
    bic={}; best_k=None
    if scores.shape[1]>=2:
        pcs=scores[:,:min(5,scores.shape[1])]
        for k in (1,2,3,4):
            if len(pcs)>=k+1:
                gm=GaussianMixture(n_components=k,random_state=2026,n_init=10,reg_covar=1e-4).fit(pcs)
                bic[k]=float(gm.bic(pcs))
        if bic: best_k=min(bic,key=bic.get)
    samples_out=[]
    for j,i in enumerate(cod_idx):
        s=samples[i]; called=np.isfinite(Xc[:,j]); calls=Xc[called,j]
        state=popstates.get(sample_pop(s),{}).get(f'{lg}_state','unknown')
        conf=popstates.get(sample_pop(s),{}).get(f'{lg}_confidence','unknown')
        samples_out.append({'lg':lg,'sample_id':s,'population':sample_pop(s),'population_arrangement_state':state,
            'population_state_confidence':conf,'raw_variant_records':raw,
            'post_sitefilter_called_genotypes':int(called.sum()),'post_sitefilter_missing_genotypes':int((~called).sum()),
            'post_sitefilter_call_rate':float(called.mean()) if len(called) else 'NA',
            'heterozygous_fraction_among_calls':float(np.mean(calls==1)) if len(calls) else 'NA',
            'putative_arrangement_genotype':('XX' if state=='ancestral' else 'YY' if state=='derived' else 'unassigned'),'assignment_confidence':('high' if conf=='direct_outgroup_alignment' else 'moderate' if state in ('ancestral','derived') else 'NA'),
            'assignment_reason':('published homokaryotype-selected cohort; sample mapped by VCF population/ecotype label; arrangement orientation from paper and frozen metadata' if state in ('ancestral','derived') else 'no population arrangement state')})
    # Population-colored PCA. Unknown/unmapped populations are listed in legend too.
    if scores.shape[1]>=2:
        fig,ax=plt.subplots(figsize=(10,7)); pops=[sample_pop(samples[i]) for i in cod_idx]
        unique=sorted(set(pops)); cmap=plt.get_cmap('tab20',max(20,len(unique)))
        for z,pop in enumerate(unique):
            ids=[i for i,p in enumerate(pops) if p==pop]
            ax.scatter(scores[ids,0],scores[ids,1],label=pop.replace('Gadmor_',''),s=50,color=cmap(z),edgecolor='black',linewidth=.35)
            for i in ids: ax.annotate(samples[cod_idx[i]].replace('Gadmor_',''),(scores[i,0],scores[i,1]),xytext=(3,2),textcoords='offset points',fontsize=6)
        ax.set_xlabel(f'PC1 ({exp[0]*100:.1f}%)');ax.set_ylabel(f'PC2 ({exp[1]*100:.1f}%)')
        ax.set_title(f'{lg} inversion-region PCA; filtered biallelic SNPs\nPopulation colors are descriptive, not arrangement assignments')
        ax.legend(title='Population',bbox_to_anchor=(1.02,1),loc='upper left',fontsize=8);fig.tight_layout();fig.savefig(FIG/f'{lg}_pca.png',dpi=180);plt.close(fig)
    gts={}
    for row in samples_out:
        gts[row['sample_id']]=row
    return {'lg':lg,'samples':samples_out,'raw_snp_records':raw,'site_filtered_records':qual,
      'quality_filtered_missing_fraction':qual_missing/(qual*len(samples)) if qual and samples else 'NA',
      'raw_missing_fraction':raw_missing/(raw*len(samples)) if raw and samples else 'NA',
      'mean_site_QUAL':float(np.nanmean(siteq)) if siteq and np.isfinite(siteq).any() else 'NA',
      'mean_genotype_GQ':float(np.mean(gq_vals)) if gq_vals else 'NA','mean_genotype_DP':float(np.mean(dp_vals)) if dp_vals else 'NA',
      'samples_in_vcf':len(samples),'Atlantic_cod_samples':len(cod_idx),'pca_snps':int(sum(keep)) if len(keep) else 0,
      'pc1_variance':float(exp[0]) if len(exp)>0 else 'NA','pc2_variance':float(exp[1]) if len(exp)>1 else 'NA',
      'gmm_bic_by_k':bic,'best_exploratory_k':best_k,'interpretation':'exploratory genotype structure only; arrangement classes assigned from published homokaryotype cohort and population/ecotype metadata'}


def main():
    FIG.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True); states=read_pop_states(); summaries=[]; sample_rows=[]
    for lg,p in LGFILES.items():
        r=analyze(lg,p,states); samples=r.pop('samples'); sample_rows+=samples; summaries.append(r)
    def write(path,rows,fields):
        with path.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
    write(OUT/'genotype_pca_lg_summary.tsv',summaries,['lg','samples_in_vcf','Atlantic_cod_samples','raw_snp_records','site_filtered_records','raw_missing_fraction','quality_filtered_missing_fraction','mean_site_QUAL','mean_genotype_GQ','mean_genotype_DP','pca_snps','pc1_variance','pc2_variance','gmm_bic_by_k','best_exploratory_k','interpretation'])
    write(OUT/'individual_genotype_classification.tsv',sample_rows,['lg','sample_id','population','population_arrangement_state','population_state_confidence','raw_variant_records','post_sitefilter_called_genotypes','post_sitefilter_missing_genotypes','post_sitefilter_call_rate','heterozygous_fraction_among_calls','putative_arrangement_genotype','assignment_confidence','assignment_reason'])
    with (OUT/'genotype_pca_audit.md').open('w') as f:
        f.write('# Individual genotype and arrangement identifiability audit\n\n')
        f.write('The inversion VCFs contain diploid individual genotype calls (GT), genotype quality (GQ), read depth (DP), and site-level QUAL; they are not population summaries. This corrects the earlier local audit, which had not inspected the available headers. Arrangement classes are assigned from the published homokaryotype-selected cohort and documented population/ecotype mapping; PCA/GMM do not assign classes.\n\n')
        for r in summaries:
            f.write(f"## {r['lg']}\n\n- Samples: {r['samples_in_vcf']} total VCF columns; {r['Atlantic_cod_samples']} Gadmor cod individuals used in PCA.\n- Biallelic SNP records in inversion VCF: {r['raw_snp_records']}; PASS/QUAL≥30: {r['site_filtered_records']}; PCA SNPs after ≥90% call rate and MAF≥0.02 with GQ≥20, DP≥8: {r['pca_snps']}.\n- Raw missing genotype fraction: {r['raw_missing_fraction']}; after genotype filters: {r['quality_filtered_missing_fraction']}; mean QUAL: {r['mean_site_QUAL']}; mean GQ: {r['mean_genotype_GQ']}; mean DP: {r['mean_genotype_DP']}.\n- PC variance: {r['pc1_variance']}, {r['pc2_variance']}; exploratory GMM BIC (PC1–PC5): {r['gmm_bic_by_k']}; lowest BIC k={r['best_exploratory_k']}.\n- Interpretation: PCA/clustering describes SNP structure. All 24 cod individuals are part of the published cohort: 22 additional individuals were selected as homokaryotypic for all four supergenes, and the two assembly specimens were selected/identified as complementary homozygotes. We therefore map their VCF population/ecotype IDs to published arrangement classes; PCA is not used to define those classes. No XY class is expected because the cohort was ascertained to exclude arrangement heterokaryotypes.\n\n")
        f.write('All 24 cod samples map to the published cohort: 22 selected homokaryotypes plus the two complementary assembly specimens. The cohort has no expected XY class because heterokaryotypes were excluded by study design. The PCA does not recover three clusters; its lowest GMM BIC among k=1–4 is k=4 for each LG, likely reflecting population structure, so it is not interpreted as arrangement count. SNP-conditional tracks are estimated separately in estimate_tracks.py; no callable invariant mask exists, and windows outside or crossing suppressed-region boundaries remain unassigned.\n')

if __name__=='__main__':main()
