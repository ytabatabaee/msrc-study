#!/usr/bin/env python3
"""Dependency-light SVG figures and final interpretation report."""
import csv, html, shutil, subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[5]; BASE=ROOT/"empirical/neoaves_chr4/analysis/patterson_d_finite_sites"; RES=BASE/"results"; FIG=BASE/"figures"
INT=(('PNAS_1',25030000,32670000),('PNAS_2',33510000,34470000),('PNAS_3',44130000,56810000))
def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def n(x):
    try:return float(x)
    except:return None
def reg(s):
    for q,a,b in INT:
        if a<=s<=b:return q
    return 'outside_all_PNAS'
def svg(name,body,w=1500,h=800):
    s=f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial,sans-serif;font-size:14px}} .title{{font-size:22px;font-weight:bold}} .axis{{stroke:#555;stroke-width:1}} .line{{fill:none;stroke-width:2}}</style>{body}</svg>'
    p=FIG/(name+'.svg');p.write_text(s)
    c=shutil.which('rsvg-convert')
    if c:
        subprocess.run([c,'-f','png','-o',str(FIG/(name+'.png')),str(p)],check=False)
        subprocess.run([c,'-f','pdf','-o',str(FIG/(name+'.pdf')),str(p)],check=False)
def xmap(x,x0,x1,left=100,right=1450):return left+(x-x0)/(x1-x0)*(right-left)
def plot_tracks():
    w=[r for r in read(RES/'window_model_comparison.tsv') if r['resolution']=='250kb']; x0=0;x1=max(n(r['window_end']) for r in w); body='<text x="100" y="30" class="title">Neoaves chr4 N61 finite-sites predictions</text>'
    specs=(('Observed D','observed_D','#d62728'),('JC69 D','JC69_D','#1f77b4'),('HKY85 D','HKY85_D','#2ca02c'),('GTR+Gamma D','GTR_Gamma_D','#9467bd'))
    for j,(lab,key,col) in enumerate(specs):
        top=55+j*135; body+=f'<text x="10" y="{top+55}">{lab}</text><line class="axis" x1="100" y1="{top+105}" x2="1450" y2="{top+105}"/>'
        for _,a,b in INT:body+=f'<rect x="{xmap(a,x0,x1):.1f}" y="{top}" width="{xmap(b,x0,x1)-xmap(a,x0,x1):.1f}" height="105" fill="#f4b183" opacity=".25"/>'
        pts=[]
        for r in w:
            v=n(r[key]);
            if v is not None:pts.append(f'{xmap(n(r["midpoint"]),x0,x1):.1f},{top+105*(1-(v+1)/2):.1f}')
        body+=f'<polyline class="line" stroke="{col}" points="{" ".join(pts)}"/>'
    body+='<text x="100" y="610">Coordinate on Stiller/GalGal6 axis (Mb; published PNAS intervals shown as approximate mask)</text>'
    svg('neoaves_N61_finite_sites_diagnostics',body,1500,640)
def barfig():
    reg=read(RES/'regional_model_comparison.tsv'); regions=['PNAS_1','PNAS_2','PNAS_3','all_PNAS','outside_all_PNAS']; r=[q for q in reg if q['design']=='N61_primary' and q['filter']=='all']; body='<text x="100" y="30" class="title">N61 observed versus finite-sites regional D</text>'; cols=(('#d62728','observed_D','Observed'),('#1f77b4','JC69_D','JC69'),('#2ca02c','HKY85_D','HKY85'),('#9467bd','GTR_Gamma_D','GTR+Gamma'))
    for i,q in enumerate(regions):
        x=140+i*260; body+=f'<text x="{x}" y="780">{q}</text>'
        for j,(c,k,l) in enumerate(cols):
            v=n(next(z[k] for z in r if z['region']==q)); y=600-v*380; body+=f'<rect x="{x+j*42}" y="{min(y,600):.1f}" width="30" height="{abs(y-600):.1f}" fill="{c}"/>'
    body+='<line class="axis" x1="100" y1="600" x2="1450" y2="600"/><text x="100" y="810">D=(ABBA-BABA)/(ABBA+BABA)</text><text x="100" y="840">Legend: Observed (red) | JC69 (blue) | HKY85 (green) | GTR+Gamma (purple)</text>';svg('neoaves_N61_finite_sites_regional_comparison',body,1500,880)
def residualfig():
    z=[r for r in read(RES/'observed_predicted_site_counts_with_residuals.tsv') if r['design']=='N61_primary']; body='<text x="100" y="30" class="title">N61 finite-sites residuals versus branch lengths</text>'
    for j,(key,lab) in enumerate((('internal_branch','internal branch'),('terminal_mean','mean terminal branch'))):
        top=70+j*300; body+=f'<text x="40" y="{top}">{lab}</text><line class="axis" x1="100" y1="{top+250}" x2="1450" y2="{top+250}"/>'
        vals=[n(r[key]) for r in z if n(r[key]) is not None]; mx=max(vals) if vals else 1
        for r in z:
            a=n(r[key]); b=n(r['GTR_Gamma_D_residual']);
            if a is not None and b is not None:body+=f'<circle cx="{100+1350*a/mx:.1f}" cy="{top+125-80*b:.1f}" r="1.8" fill="#555" opacity=".35"/>'
    svg('neoaves_N61_finite_sites_residuals',body,1500,680)
def topologyfig():
    z=[r for r in read(RES/'topology_site_pattern_audit.tsv') if r['design']=='N61_primary' and r['topology'] in ('q1','q2','q3')]; qs=['PNAS_1','PNAS_2','PNAS_3','all_PNAS','outside_all_PNAS']; body='<text x="100" y="30" class="title">N61 inferred topology composition</text>'; colors={'q1':'#7f7f7f','q2':'#ff7f0e','q3':'#2ca02c'}
    for i,q in enumerate(qs):
        x=150+i*260; total=sum(n(next(r['n_loci'] for r in z if r['region']==q and r['topology']==t)) or 0 for t in ('q1','q2','q3')); bottom=600
        for t in ('q1','q2','q3'):
            v=n(next(r['n_loci'] for r in z if r['region']==q and r['topology']==t)) or 0; h=400*v/total if total else 0; bottom-=h; body+=f'<rect x="{x}" y="{bottom:.1f}" width="100" height="{h:.1f}" fill="{colors[t]}"/><text x="{x+105}" y="{bottom+h/2:.1f}">{t}</text>'
        body+=f'<text x="{x}" y="650">{q}</text>'
    svg('neoaves_N61_topology_composition',body,1500,720)
def controls():
    reg=read(RES/'regional_model_comparison.tsv'); ds=(('Columbea_primary','Columbea'),('N61_primary','N61'),('N62_primary','N62')); body='<text x="100" y="30" class="title">Observed and GTR+Gamma D across clades</text>'
    for i,(d,l) in enumerate(ds):
        for j,q in enumerate(('all_PNAS','outside_all_PNAS')):
            r=next(x for x in reg if x['design']==d and x['filter']=='all' and x['region']==q);v=n(r['observed_D']);x=180+i*400+j*120;y=600-v*400;body+=f'<rect x="{x}" y="{min(y,600):.1f}" width="70" height="{abs(y-600):.1f}" fill="{"#d62728" if j==0 else "#999"}"/><text x="{x}" y="650">{l} {"in" if j==0 else "out"}</text>'
    body+='<line class="axis" x1="100" y1="600" x2="1400" y2="600"/>';svg('neoaves_all_clades_finite_sites',body,1500,720)
def report():
    r=read(RES/'regional_model_comparison.tsv');g=lambda d,q:next(x for x in r if x['design']==d and x['filter']=='all' and x['region']==q);a=g('N61_primary','all_PNAS');o=g('N61_primary','outside_all_PNAS'); lines=['# Finite-sites interpretation','','## Numerical result','', '| Region | observed D | JC69 | HKY85 | GTR+Gamma | observed ABBA/BABA | GTR expected ABBA/BABA |','|---|---:|---:|---:|---:|---:|---:|']
    for label,x in [('all PNAS',a),('outside',o)]:lines.append(f"| {label} | {x['observed_D']} | {x['JC69_D']} | {x['HKY85_D']} | {x['GTR_Gamma_D']} | {x['observed_ABBA']}/{x['observed_BABA']} | {float(x['GTR_Gamma_expected_ABBA']):.1f}/{float(x['GTR_Gamma_expected_BABA']):.1f} |")
    lines.append(f"| inside minus outside | {float(a['observed_D'])-float(o['observed_D']):.6f} | {float(a['JC69_D'])-float(o['JC69_D']):.6f} | {float(a['HKY85_D'])-float(o['HKY85_D']):.6f} | {float(a['GTR_Gamma_D'])-float(o['GTR_Gamma_D']):.6f} | | |")
    lines += ['', 'The near−1 infinite-sites predictions arise from the extreme inferred discordant-topology composition inside the PNAS mask. Complete finite-sites pruning, including terminal branches, recurrent substitutions, and the empirical base composition, substantially attenuates this ratio. GTR+Gamma predicts the N61 shift at approximately −0.1988, versus observed −0.1776; JC69 and HKY85 predict −0.3022 and −0.2558. Thus terminal-branch/homoplasy effects explain much of the mismatch, while model choice changes the residual modestly.', '', 'Expected absolute counts remain lower than observed counts, so agreement in D is not a full likelihood fit. The trees and alignments are not independent: the published gene trees were estimated from these alignments. The original inference model and a faithful parametric re-estimation experiment were not recovered, so gene-tree error cannot be quantified here.', '', 'Columbea remains positive and N62 near zero, unlike N61. These results are compatible with N61-specific genealogy asymmetry, but do not distinguish MSRC from introgression, ancestral structure, or other causes. PNAS intervals are the published GalGal4-derived mask applied to the Stiller GalGal6 axis; exact boundary lifting remains unresolved.', '', 'See `statistical_validation.tsv` for 500-kb, 1-Mb, and 2-Mb block bootstrap intervals.']
    (RES/'finite_sites_interpretation.md').write_text('\n'.join(lines)+'\n')
if __name__=='__main__':
    FIG.mkdir(parents=True,exist_ok=True);plot_tracks();barfig();residualfig();topologyfig();controls();report()
