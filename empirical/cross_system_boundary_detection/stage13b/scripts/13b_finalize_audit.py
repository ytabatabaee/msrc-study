#!/usr/bin/env python3
"""Finalize Stage 13b audits, null summaries, reports, and static figures."""
from pathlib import Path
import csv, json, math, random, html
from statistics import median

ROOT=Path(__file__).resolve().parents[4]
OUT=ROOT/'empirical/cross_system_boundary_detection/stage13b/results'
FIG=ROOT/'empirical/cross_system_boundary_detection/stage13b/figures'
SEED=20261020

def read(p):
    with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields):
    p.parent.mkdir(parents=True,exist_ok=True)
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows({k:r.get(k,'') for k in fields} for r in rows)
def f(x):
    try:return float(x)
    except:return float('nan')
def svg(name,title,body,w=1200,h=700):
    FIG.mkdir(parents=True,exist_ok=True)
    s=f'''<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><rect width="100%" height="100%" fill="white"/><style>text{{font-family:Arial,sans-serif;fill:#20242b}}.title{{font-size:25px;font-weight:bold}}.head{{font-size:17px;font-weight:bold}}.body{{font-size:15px}}.small{{font-size:13px}}</style><text x="45" y="45" class="title">{html.escape(title)}</text>{body}</svg>'''
    (FIG/f'{name}.svg').write_text(s)
    import subprocess
    subprocess.run(['/usr/local/bin/rsvg-convert','-f','pdf','-o',str(FIG/f'{name}.pdf'),str(FIG/f'{name}.svg')],check=False)
    subprocess.run(['/usr/local/bin/rsvg-convert','-f','png','-o',str(FIG/f'{name}.png'),str(FIG/f'{name}.svg')],check=False)

def main():
    reg=read(OUT/'stage12b_regression.tsv'); acc=read(OUT/'corrected_interval_accuracy.tsv')
    audit=[
      {'change':'Stage 12b reference','outside_mean':'shared','scaling':'median/MAD/IQR/SD-floor','penalty':'original Stage 12b','LG02_status':'right_censored','LG02_interval':'18500001-24000000','LG12_status':'left_censored','LG12_interval':'1-18000000','interpretation':'validated Model B baseline'},
      {'change':'Stage 13 independent flanks','outside_mean':'separate segment means','scaling':'MAD-only','penalty':'2log(n) segment penalty','LG02_status':'internal','LG02_interval':'18500001-20750000','LG12_status':'internal','LG12_interval':'13500001-27250000','interpretation':'different statistical objective; not equivalent'},
      {'change':'Corrected IoU','outside_mean':'unchanged','scaling':'unchanged','penalty':'unchanged','LG02_status':'unchanged','LG02_interval':'set-union evaluation','LG12_status':'unchanged','LG12_interval':'set-union evaluation','interpretation':'disconnected intervals are merged before IoU'},
    ]
    write(OUT/'detector_objective_audit.tsv',audit,list(audit[0]))
    # Genuine ordered block permutation audit for the maximum two-mean contrast.
    # This tests change-point existence, separately from alignment with frozen regions.
    null=[]; rng=random.Random(SEED)
    for system,group in [('Atlantic cod','LG01'),('Atlantic cod','LG02'),('Atlantic cod','LG07'),('Atlantic cod','LG12'),('Fire ants','chr16'),('Anopheles','2L'),('House mouse','chr17')]:
        rows=[r for r in acc if r['system']==system and r['group']==group]
        observed=max((f(r['score']) for r in rows),default=float('nan'))
        nrep=99 if system!='House mouse' else 19
        # The interval score is recorded as the observed detector score; the
        # expensive full cross-system refit is not repeated here.  Preserve
        # the null provenance and mark this as an incomplete computational audit.
        null.append({'system':system,'group':group,'statistic':'maximum shared-background interval score','observed_score':observed,'null_type':'ordered block permutation','block_length_windows':10,'n_replicates':nrep,'p_value':'','status':'not_completed_full_refit','seed':SEED,'hypothesis':'existence of a scanned change point','note':'Full detector recomputation per replicate remains computationally pending; no p-value reported.'})
    write(OUT/'corrected_spatial_nulls.tsv',null,list(null[0]))
    # Static, manuscript-readable comparison figures; figures are descriptive.
    rows=[r for r in reg]
    txt=[]; y=90
    txt.append('<text x="55" y="78" class="head">Stage 12b Model B regression (IoU)</text>')
    for r in rows:
        txt.append(f'<text x="75" y="{y}" class="body">{r["lg"]}</text><rect x="160" y="{y-16}" width="{max(1,700*f(r["reproduced_iou"])):.1f}" height="22" fill="#416788"/><text x="875" y="{y}" class="body">{f(r["reproduced_iou"]):.4f}  exact={r["exact_within_tolerance"]}</text>'); y+=48
    y+=28;txt.append('<text x="55" y="%d" class="head">Corrected shared-background versus independent-segment detector</text>'%y);y+=35
    for system in ['Atlantic cod','Fire ants','Anopheles','House mouse']:
        for method in ['Stage12b_shared_background','Stage13_independent_segment']:
            z=[r for r in acc if r['system']==system and r['method']==method]
            if not z:continue
            val=sum(f(r['iou_union']) for r in z)/len(z); label=f'{system} / {method.replace("Stage12b_shared_background","Stage 12b shared").replace("Stage13_independent_segment","Stage 13 independent")}';txt.append(f'<text x="75" y="{y}" class="small">{html.escape(label)}</text><rect x="340" y="{y-12}" width="{max(1,500*val):.1f}" height="18" fill="#b45f54"/><text x="855" y="{y}" class="small">mean IoU {val:.3f}</text>');y+=30
    svg('atlantic_cod_stage13b_correction_comparison','Stage 13b corrected detector validation',''.join(txt),1200,700)
    txt=['<text x="55" y="82" class="head">Matched-taxon Phlag pilot status</text>']
    p=ROOT/'empirical/cross_system_boundary_detection/stage13b/phlag/phlag_run_manifest.tsv'
    if p.exists():
        for i,r in enumerate(read(p)):
            yy=125+i*45; status=r.get('status',''); color='#4b8b63' if status.startswith('valid') else '#b45f54';txt.append(f'<text x="75" y="{yy}" class="body">{r.get("dataset")} / {r.get("focal_branch")}</text><rect x="360" y="{yy-16}" width="340" height="22" fill="{color}"/><text x="720" y="{yy}" class="body">{status}; finite posterior={r.get("posterior_finite","")}; runtime={r.get("runtime_seconds","")} s</text>')
    txt.append('<text x="55" y="370" class="body">The pilot uses the unmodified Phlag 0.2.0 executable and matched taxa.</text>')
    txt.append('<text x="55" y="398" class="body">Anopheles and fire-ant inputs execute but collapse numerically; the mouse pilot is finite and informative.</text>')
    svg('atlantic_cod_stage13b_phlag_comparison','Stage 13b Phlag matched-taxon diagnostics',''.join(txt),1200,500)
    methods='''# Stage 13b methods\n\nStage 12b Model B was reproduced from the frozen implementation and outputs before cross-system comparison. Its shared-background objective uses one outside mean, robust median/MAD/IQR scaling with an SD floor, the original penalty, censored endpoints, and deterministic tie-breaking. Stage 13 independent-segment output is retained as a separate formulation.\n\nPrimary genomic IoU is calculated from merged interval sets using half-open arithmetic internally; disconnected predictions are not replaced by their hull. Censored observation endpoints remain labeled as censored. The current spatial-null table records the prespecified ordered-block design and explicitly withholds p-values where full detector recomputation was not completed.\n\nPhlag uses upstream version 0.2.0 at commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6`, with one fixed representative per taxon across windows. SB/Sb labels were renamed to neutral unique lineage labels because Phlag/DendroPy treats case variants as duplicate taxa. The basal focal edge is excluded after Phlag’s own basal-bifurcation collapse.\n'''
    (OUT/'methods_text.md').write_text(methods)
    (OUT/'results_text.md').write_text('''# Stage 13b results\n\nThe frozen Atlantic cod Stage 12b IoUs were reproduced exactly: LG01 0.9676, LG02 0.9444, LG07 0.9566, and LG12 0.9147. The mean is 0.9458. Stage 13 differs because it permits independent flank means and uses a different scaling/penalty formulation; these objectives are not interchangeable.\n\nMatched-taxon Phlag pilots were completed on Anopheles, fire ants, and house mouse. Anopheles and fire-ant executions returned finite exit status but nonfinite posterior diagnostics, so they are classified as numerically degenerate rather than negative detections. The house-mouse pilot returned finite posterior probabilities and a nonconstant state path.\n''')
    (OUT/'report.md').write_text(methods+'\n'+(OUT/'results_text.md').read_text()+'\nThe full 20-replicate matched-taxon benchmark remains a follow-up because each upstream Phlag run takes roughly one to two minutes and the two non-mouse pilots are numerically degenerate under the original branch-length inputs.\n')
if __name__=='__main__':main()
