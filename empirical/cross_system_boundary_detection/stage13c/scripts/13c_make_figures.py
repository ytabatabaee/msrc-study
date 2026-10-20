#!/usr/bin/env python3
from pathlib import Path
import csv, html, subprocess
ROOT=Path(__file__).resolve().parents[4]; RES=ROOT/'empirical/cross_system_boundary_detection/stage13c/results'; FIG=ROOT/'empirical/cross_system_boundary_detection/stage13c/figures'
def read(p):
 with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def make(name,title,rows,subtitle=''):
 FIG.mkdir(parents=True,exist_ok=True); y=90; body=[f'<text x="50" y="45" class="title">{html.escape(title)}</text>',f'<text x="50" y="70" class="sub">{html.escape(subtitle)}</text>']
 for r in rows:
  label=' / '.join(str(r.get(k,'')) for k in ('system','group','method') if r.get(k,'')); val=r.get('value',r.get('score',r.get('iou',''))); body.append(f'<text x="60" y="{y}" class="body">{html.escape(label[:80])}</text><text x="900" y="{y}" class="body">{html.escape(str(val))}</text>');y+=28
 s='''<svg xmlns="http://www.w3.org/2000/svg" width="1200" height="%d"><rect width="100%%" height="100%%" fill="white"/><style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:25px;font-weight:bold}.sub{font-size:14px;fill:#555}.body{font-size:15px}</style>%s</svg>'''%(max(500,y+40),''.join(body));p=FIG/(name+'.svg');p.write_text(s);subprocess.run(['/usr/local/bin/rsvg-convert','-f','pdf','-o',str(FIG/(name+'.pdf')),str(p)],check=False);subprocess.run(['/usr/local/bin/rsvg-convert','-f','png','-o',str(FIG/(name+'.png')),str(p)],check=False)
def main():
 reg=read(RES/'stage12b_regression.tsv'); prof=read(RES/'lg02_detector_discrepancy.tsv')
 make('lg02_detector_score_profiles','Atlantic cod LG02: detector score profiles',prof[:40],'Top candidates from the exact Stage 12b path and the distinct cross-system D+A implementation; coordinates are Mb.')
 make('segmentation_model_selection','Single versus multi-interval model selection',read(RES/'segmentation_model_selection.tsv'),'BIC-style selection allows K=0..4; reference coordinates are evaluation-only.')
 fire=read(RES/'robust_loss_sensitivity.tsv');make('fire_ant_robust_distance','Fire-ant distance outlier sensitivity',fire,'Huber loss is applied without deleting windows; source units are substitution distances.')
 ph=read(RES/'phlag_replicate_accuracy.tsv');make('phlag_house_mouse_replicates','House-mouse Phlag valid pilot replicates',ph,'Only finite, informative runs are scored; adjacent windows are not treated as independent replicates.')
 diag=read(RES/'phlag_numerical_diagnostics.tsv');make('phlag_numerical_diagnostics','Phlag numerical validity by dataset',diag,'Nonfinite posterior runs are classified as degenerate and excluded from accuracy metrics.')
 acc=read(RES/'cross_system_summary.tsv');make('cross_system_detector_accuracy','Cross-system detector comparison',acc,'Continuous tracks retain their source units; this is a descriptive comparison of frozen and exploratory analyses.')
if __name__=='__main__':main()
