#!/usr/bin/env python3
"""Collect matched-taxon Phlag pilot outputs and score only valid runs."""
from pathlib import Path
import csv, math, re
ROOT=Path(__file__).resolve().parents[4]; P=ROOT/'empirical/cross_system_boundary_detection/stage13b/phlag'; RAW=P/'raw_results'
def read(p):
    with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows,fields):
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows({k:r.get(k,'') for k in fields} for r in rows)
def parse(path):
    lines=[x.strip() for x in path.read_text(errors='replace').splitlines() if x.strip() and not x.startswith('#') and ',' in x]
    state=next((x for x in lines if all(y in {'0','1'} for y in x.split(','))), '')
    probs=lines[lines.index(state)+1].split(',') if state and lines.index(state)+1<len(lines) else []
    finite=bool(probs) and all(math.isfinite(float(x)) for x in probs)
    cls='numerically_degenerate' if not finite else ('valid_no_detection' if len(set(state.split(',')))<=1 else 'valid_informative')
    return cls,state.split(',') if state else [],probs
def main():
    runs=[]; windows=[]; acc=[]
    coord_sources={
      'anopheles_2la':ROOT/'empirical/anopheles_2la/results/stage12_distance/stage12_window_distance.tsv',
      'fire_ants_chr16':ROOT/'empirical/fire_ants_chr16/results/stage9_distance/stage9_window_patristic_distance.tsv',
      'house_mouse_t_complex':ROOT/'empirical/house_mouse_t_complex/results/stage4_distance/stage4_window_distance.tsv'}
    refs={'anopheles_2la':[(20524058,42165532)],'fire_ants_chr16':[(6614000,10400000)]}
    mi=ROOT/'empirical/house_mouse_t_complex/results/t_complex_inversion_intervals.tsv'
    if mi.exists(): refs['house_mouse_t_complex']=[(int(x['start_bp']),int(x['end_bp'])) for x in read(mi)]
    for dataset in ('anopheles_2la','fire_ants_chr16','house_mouse_t_complex'):
      for path in sorted((RAW/'runs'/dataset).glob('replicate_*/*.txt')):
        rep=int(path.parent.name.split('_')[-1]); focal=path.stem; cls,states,probs=parse(path)
        mapping=read(RAW/'matched_inputs'/dataset/f'replicate_{rep:03d}'/'window_mapping.tsv')
        r={'dataset':dataset,'replicate':rep,'focal_branch':focal,'status':cls,'n_windows':len(mapping),'n_state1':sum(x=='1' for x in states),'posterior_finite':str(all(math.isfinite(float(x)) for x in probs)).lower() if probs else 'false','posterior_max':max((float(x) for x in probs),default=float('nan')),'output':str(path)};runs.append(r)
        if cls.startswith('valid'):
          coords=read(coord_sources[dataset]); coords=coords[:len(mapping)]
          pred=[]; tp=0; refbp=sum(b-a+1 for a,b in refs.get(dataset,[]))
          for i,m in enumerate(mapping):
            c=coords[i] if i<len(coords) else {}; start=int(float(c.get('start_bp',c.get('start',0)) or 0)); end=int(float(c.get('end_bp',c.get('end',0)) or 0)); state=states[i] if i<len(states) else '0';
            if start and end and state=='1': pred.append((start,end))
            windows.append({'dataset':dataset,'replicate':rep,'focal_branch':focal,**m,'start':start,'end':end,'predicted_state':state,'posterior_anomaly':probs[i] if i<len(probs) else ''})
          pred=sorted(pred); merged=[]
          for a,b in pred:
            if merged and a<=merged[-1][1]+1: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
            else: merged.append((a,b))
          p_bp=sum(b-a+1 for a,b in merged); inter=sum(max(0,min(b,d)-max(a,c)+1) for a,b in merged for c,d in refs.get(dataset,[])); union=p_bp+refbp-inter; precision=inter/p_bp if p_bp else 0; recall=inter/refbp if refbp else 0; f1=2*precision*recall/(precision+recall) if precision+recall else 0
          acc.append({'dataset':dataset,'replicate':rep,'focal_branch':focal,'status':cls,'accuracy_status':'scored','predicted_intervals':';'.join(f'{a}-{b}' for a,b in merged),'interval_iou':inter/union if union else 0,'f1':f1,'note':'Default state-1 call; pilot representative and focal edge only.'})
    fields=list(runs[0]) if runs else ['dataset','replicate','focal_branch','status'];write(P/'phlag_run_manifest.tsv',runs,fields);write(P/'phlag_numerical_diagnostics.tsv',runs,fields);write(P/'window_level_predictions.tsv',windows,list(windows[0]) if windows else ['dataset','replicate','focal_branch','window_id','predicted_state','posterior_anomaly']);write(P/'phlag_accuracy.tsv',acc,list(acc[0]) if acc else ['dataset','replicate','focal_branch','status','accuracy_status','note'])
if __name__=='__main__':main()
