#!/usr/bin/env python3
"""Scan the full mouse chr17 for sitewise pseudo-t SNP-panel estimates.
Tree coverage is indicated in the output; all statistics remain SNP-only.
"""
import os, certifi, json
os.environ['SSL_CERT_FILE']=certifi.where(); os.environ['REQUESTS_CA_BUNDLE']=certifi.where()
from pathlib import Path
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
ROOT=Path(__file__).resolve().parents[4]
OUT=ROOT/'empirical/house_mouse_t_complex/window_diversity'
URL='https://wwwuser.gwdguser.de/~evolbio/evolgen/wildmouse/vcf/AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz'
TBI=ROOT/'AllMouse.vcf_90_recalibrated_snps_raw_indels_reheader_PopSorted.PASS.vcf.gz.tbi'
MEM=ROOT/'empirical/house_mouse_t_complex/analysis/patterson_d/config/frozen_quartet_members.tsv'
CHECK=OUT/'results/mouse_full_chr17_checkpoint_5Mb.json'
START=0; WIN=5_000; CHUNK=500_000
INVERSIONS=[(5_747_331,6_957_919),(7_435_308,13_016_510),(14_068_229,18_051_507),(19_401_000,40_049_964)]

def gt_pair(gt):
    if gt is None or len(gt)!=2 or any(a is None for a in gt): return None
    return int(gt[0]),int(gt[1])

def main():
    mem=pd.read_csv(MEM,sep='\t')
    designs={}
    for des,g in mem[mem.design.str.startswith(('D1_Mmd','D2_Mmm','D3_Mmc'))].groupby('design'):
        designs[des]=(g[g.role=='P1'].sample_id.tolist(),g[g.role=='P2'].sample_id.tolist())
    vf=pysam.VariantFile(URL,index_filename=str(TBI)); samples=set(vf.header.samples)
    END=int(vf.header.contigs['chr17'].length)
    for des,(t,std) in designs.items():
        if not set(t+std)<=samples: raise ValueError(f'{des}: sample mismatch')
    union=sorted(set(x for groups in designs.values() for g in groups for x in g))
    starts=list(range(START,END,WIN)); tree=pd.read_csv(ROOT/'empirical/house_mouse_t_complex/results/stage5_unified_tracks/window_unified_tracks.tsv',sep='\t')
    tree_starts=set(tree.start_bp.astype(int));
    state=json.loads(CHECK.read_text()) if CHECK.exists() else {'processed':[],'acc':{},'variant_records':0,'ambiguous':{d:0 for d in designs},'opportunities':{d:0 for d in designs}}
    if 'completed_ranges' not in state:
        state['completed_ranges']=[[int(c),min(int(c)+5_000_000,END)] for c in state['processed']]
        CHECK.write_text(json.dumps(state,separators=(',',':')))
    completed=[(int(a),int(b)) for a,b in state['completed_ranges']]
    for chunk in range(START,END,CHUNK):
        upper=min(chunk+CHUNK,END)
        if any(a<=chunk and upper<=b for a,b in completed): continue
        print(f'Fetching chr17:{chunk:,}-{upper:,}',flush=True)
        for rec in vf.fetch('chr17',chunk,upper):
            if len(rec.alleles)!=2 or not rec.ref or len(rec.ref)!=1 or rec.ref not in 'ACGT' or not rec.alts or len(rec.alts[0])!=1 or rec.alts[0] not in 'ACGT': continue
            st=((rec.pos-1)//WIN)*WIN
            if st<START or st>=END: continue
            state['variant_records']+=1
            gtmap={sample:gt_pair(rec.samples[sample].get('GT')) for sample in union}
            for des,(ts,ss) in designs.items():
                std_calls=[gtmap[s] for s in ss]
                if any(x is None for x in std_calls): continue
                std_alleles=[a for gt in std_calls for a in gt]; present=set(std_alleles); t_alleles=[]; amb=0
                for sample in ts:
                    gt=gtmap[sample]
                    if gt is None: amb+=1; continue
                    if gt[0]==gt[1]: t_alleles.append(gt[0]); continue
                    unique=[a for a in set(gt) if a not in present]
                    if len(unique)==1: t_alleles.append(unique[0])
                    else: amb+=1
                state['ambiguous'][des]+=amb; state['opportunities'][des]+=len(ts)
                if len(t_alleles)<2: continue
                nt=len(t_alleles); pt=np.bincount(np.asarray(t_alleles),minlength=2)/nt
                sa=np.asarray(std_alleles); ns=len(sa); ps=np.bincount(sa,minlength=2)/ns
                pi_t=nt/(nt-1)*(1-float((pt*pt).sum())); pi_s=ns/(ns-1)*(1-float((ps*ps).sum()))
                dxy=1-float((pt*ps).sum()); num=dxy-(pi_t+pi_s)/2
                key=f'{des}|{st}'; a=state['acc'].setdefault(key,[0.,0.,0.,0.,0,0,0])
                a[0]+=pi_t; a[1]+=pi_s; a[2]+=dxy; a[3]+=num; a[4]+=1; a[5]+=nt; a[6]+=amb
        state['processed'].append(chunk); state['completed_ranges'].append([chunk,upper]); completed.append((chunk,upper)); CHECK.write_text(json.dumps(state,separators=(',',':')))
    vf.close()
    rows=[]
    for des,(ts,ss) in designs.items():
        for st in starts:
            a=state['acc'].get(f'{des}|{st}',[0.,0.,0.,0.,0,0,0]); n=int(a[4])
            mid=(st+min(st+WIN,END)-1)/2
            in_inversion=any(lo<=mid<=hi for lo,hi in INVERSIONS)
            rows.append(dict(design=des,start_bp=st,end_bp=min(st+WIN,END)-1,midpoint_bp=mid,
                region='inversion_union' if in_inversion else 'outside_inversions',has_local_tree=st in tree_starts,n_variant_sites_reconstructed=n,
                Hudson_FST_SNP_panel=(a[3]/a[2] if n and a[2]>0 else np.nan),dXY_per_observed_SNP=(a[2]/n if n else np.nan),
                pi_t_pseudohap_per_observed_SNP=(a[0]/n if n else np.nan),pi_standard_per_observed_SNP=(a[1]/n if n else np.nan),
                mean_reconstructed_t_carriers_per_site=(a[5]/n if n else np.nan),mean_ambiguous_carrier_calls_per_site=(a[6]/n if n else np.nan),
                n_t_carriers=len(ts),n_standard_diploids=len(ss)))
    d=pd.DataFrame(rows); d.to_csv(OUT/'results/mouse_t_pseudohap_distal_flank_window_stats.tsv',sep='\t',index=False,na_rep='NA')
    metrics=['Hudson_FST_SNP_panel','dXY_per_observed_SNP','pi_t_pseudohap_per_observed_SNP','pi_standard_per_observed_SNP']
    rng=np.random.default_rng(20261120); summ=[]
    for des,g in d.groupby('design'):
        g=g.copy(); g['block']=g.midpoint_bp//500_000
        for col in metrics:
            x=g[g.region=='inversion_union'][col].dropna(); y=g[g.region=='outside_inversions'][col].dropna(); delta=x.mean()-y.mean()
            xb=g[g.region=='inversion_union'].groupby('block')[col].mean().dropna().to_numpy(); yb=g[g.region=='outside_inversions'].groupby('block')[col].mean().dropna().to_numpy(); boot=[]
            if len(xb)>1 and len(yb)>1:
                for _ in range(10000): boot.append(rng.choice(xb,len(xb),replace=True).mean()-rng.choice(yb,len(yb),replace=True).mean())
            ci=np.quantile(boot,[.025,.975]) if boot else [np.nan,np.nan]
            summ.append(dict(design=des,statistic=col,n_inversion_union_windows=len(x),n_outside_windows=len(y),inversion_union_mean=x.mean(),outside_mean=y.mean(),difference_inversion_minus_outside=delta,block_width_bp=500000,n_inversion_blocks=len(xb),n_outside_blocks=len(yb),block_bootstrap_95_low=ci[0],block_bootstrap_95_high=ci[1]))
    pd.DataFrame(summ).to_csv(OUT/'results/mouse_t_pseudohap_chr17_full_summary.tsv',sep='\t',index=False,na_rep='NA')
    fig,axs=plt.subplots(4,1,figsize=(12,9),sharex=True)
    fig.subplots_adjust(left=.14,right=.99,top=.91,bottom=.11,hspace=.12)
    for ax,col,label in zip(axs,metrics,['Hudson FST (unitless)','dXY / PASS SNP','π t allele / PASS SNP','π standard / PASS SNP']):
        for ilo,ihi in INVERSIONS:
            ax.axvspan(ilo/1e6,ihi/1e6,color='#e7a34a',alpha=.14,lw=0)
            ax.axvline(ilo/1e6,color='#a75c0b',lw=.65,alpha=.8); ax.axvline(ihi/1e6,color='#a75c0b',lw=.65,alpha=.8)
        ax.axvline(5.0,color='#444444',lw=.8,ls='--'); ax.axvline(37.055,color='#444444',lw=.8,ls='--')
        for des,g in d.groupby('design'): ax.scatter(g.midpoint_bp/1e6,g[col],s=2,alpha=.38,label=des.replace('_t_vs_standard',''),rasterized=True)
        ax.set_ylabel(label); ax.grid(axis='y',alpha=.2)
    hh,ll=axs[0].get_legend_handles_labels(); hh.append(Patch(facecolor='#e7a34a',alpha=.25,edgecolor='none')); ll.append('Inversion intervals')
    axs[0].legend(hh,ll,fontsize=7,ncol=4); axs[-1].set_xlabel('mm10/GRCm38 chr17 (Mb; local-tree windows cover 5.0–37.055 Mb; dashed lines mark limits)')
    fig.suptitle(f'EXPLORATORY — mouse t pseudohaplotype vs standard: full chr17 (0–{END/1e6:.3f} Mb), PASS-SNP panel',fontsize=13)
    fig.text(.5,.025,'dXY and π are per observed PASS SNP, not per callable nucleotide.',ha='center',fontsize=9)
    d.to_csv(OUT/'results/mouse_t_pseudohap_chr17_full_window_stats.tsv',sep='\t',index=False,na_rep='NA')
    fig.savefig(OUT/'figures/mouse_t_pseudohap_chr17_full_tracks.pdf');fig.savefig(OUT/'figures/mouse_t_pseudohap_chr17_full_tracks.png',dpi=300)
    print('variant_records',state['variant_records'],'windows',len(starts),'windows_with_any_estimate',int((d.n_variant_sites_reconstructed>0).sum()))
    print(pd.DataFrame(summ).to_string(index=False))
if __name__=='__main__': main()
