#!/usr/bin/env python3
"""Stage 13 unified cross-system spatial boundary analysis."""
from __future__ import annotations
import csv, hashlib, json, math, random, sys
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[3]; OUT=ROOT/'empirical/cross_system_boundary_detection/results'; FIG=ROOT/'empirical/cross_system_boundary_detection/figures'; SHARED=ROOT/'empirical/shared'; OUT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
sys.path.insert(0,str(SHARED)); from spatial_boundary_detector import segment_track, intervals_to_coordinates
SEED=20261020

def read(p):
 with open(p,newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
def write(p,rows,fields):
 with open(p,'w',newline='') as h:
  w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def num(x):
 try:return float(x)
 except:return float('nan')
def good(x):return math.isfinite(num(x))
def md5(p):
 h=hashlib.sha256();
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
def overlap(a,b,c,d):return max(0,min(b,d)-max(a,c)+1)
def eval_intervals(pred,refs,rows):
 if not pred:return {'status':'no_detection','n_predicted_intervals':0,'iou_union':0,'f1':0,'precision':0,'recall':0,'false_positive_bp':0,'mean_boundary_error_kb':float('nan')}
 p_union=(min(a for a,b in pred),max(b for a,b in pred)); r_union=(min(a for a,b in refs),max(b for a,b in refs)); inter=overlap(*p_union,*r_union); union=p_union[1]-p_union[0]+1+r_union[1]-r_union[0]+1-inter
 pred_pos=sum((b-a+1) for a,b in pred); ref_pos=sum((b-a+1) for a,b in refs); tp=sum(overlap(a,b,c,d) for a,b in pred for c,d in refs); precision=tp/pred_pos if pred_pos else 0; recall=tp/ref_pos if ref_pos else 0; f1=2*precision*recall/(precision+recall) if precision+recall else 0
 errors=[]
 for a,b in pred:
  c,d=min(refs,key=lambda x:abs(x[0]-a)+abs(x[1]-b));errors += [abs(a-c)/1000,abs(b-d)/1000]
 return {'status':'detected','n_predicted_intervals':len(pred),'iou_union':inter/union,'f1':f1,'precision':precision,'recall':recall,'false_positive_bp':max(0,pred_pos-tp),'mean_boundary_error_kb':float(np.mean(errors)) if errors else float('nan')}
def shifts(rows,feature,refs,nperm=999):
 z=[r for r in rows if good(r.get(feature,''))]; z=sorted(z,key=lambda r:float(r['midpoint'])); vals=np.array([num(r[feature]) for r in z]); mask=np.array([any(a<=float(r['midpoint'])<=b for a,b in refs) for r in z]); obs=float(vals[mask].mean()-vals[~mask].mean()) if mask.any() and (~mask).any() else float('nan'); rng=random.Random(SEED); null=[];n=len(vals);k=int(mask.sum())
 for _ in range(nperm):
  start=rng.randrange(n);m=np.zeros(n,bool);m[(start+np.arange(k))%n]=True;null.append(float(vals[m].mean()-vals[~m].mean()))
 p=(1+sum(abs(x)>=abs(obs) for x in null))/(nperm+1) if good(obs) else float('nan')
 return obs,p,nperm
def quick_change_score(rows,feature):
 vals=np.array([num(r[feature]) for r in rows if good(r.get(feature,''))],float)
 if len(vals)<6:return float('nan')
 med=np.median(vals);mad=np.median(np.abs(vals-med));scale=1.4826*mad or np.std(vals) or 1.0;y=(vals-med)/scale;cs=np.cumsum(y-y.mean());return float(np.max(cs*cs)/len(y))
def load_systems():
 systems={}
 # Atlantic cod: SNAPP population-tree time contrast, joined to frozen quartet D.
 t=read(ROOT/'empirical/atlantic_cod/results/stage7_window_divergence_time_signal.tsv'); q=read(ROOT/'data/atlantic_cod/processed/stage4a_window_topology_scores.tsv'); qm={(r['lg'],r['window_id']):r for r in q}; rr=read(ROOT/'data/atlantic_cod/metadata/region_manifest.tsv')
 refs=defaultdict(list)
 for r in rr:
  if r['lg'] in {'LG01','LG02','LG07','LG12'}:refs[r['lg']].append((int(r['start']),int(r['end'])))
 cod=defaultdict(list)
 for r in t:
  if r['lg'] not in refs:continue
  x=dict(r);q0=qm.get((r['lg'],r['window_id']),{});x.update({'topology':q0.get('D_arrangement_minus_baseline',''),'continuous':r['A_opposite_minus_same'],'blind':r['mean_time_all']});cod[r['lg']].append(x)
 systems['Atlantic cod']={'groups':cod,'refs':refs,'units':'SNAPP population-tree time units'}
 # Fire ants: observed chromosome-16 BUSCO windows only.
 f=read(ROOT/'empirical/fire_ants_chr16/results/stage9_distance/stage9_window_patristic_distance.tsv'); fm=read(ROOT/'data/fire_ants_chr16/metadata/region_manifest.tsv'); fire_refs=[(int(r['coordinate_start']),int(r['coordinate_end'])) for r in fm if r['region_id']=='chr16_supergene']; fire=[dict(r,midpoint=r['mid'],topology=r['D'],continuous=r['A_distance'],blind=r['mean_all_focal_pair_distances']) for r in f if r['chrom']=='chr16'];systems['Fire ants']={'groups':{'chr16':fire},'refs':{'chr16':fire_refs},'units':'RAxML patristic substitution units'}
 # Anopheles: 2La local NJ trees.
 a=read(ROOT/'empirical/anopheles_2la/results/stage12_distance/stage12_window_distance.tsv'); an=[dict(r,topology=r['D'],continuous=r['A_distance'],blind=r['mean_pair_distance']) for r in a];systems['Anopheles']={'groups':{'2L':an},'refs':{'2L':[(20524058,42165532)]},'units':'NJ patristic substitution units'}
 # Mouse: retain all local trees as observed windows; reference has four intervals.
 m=read(ROOT/'empirical/house_mouse_t_complex/results/stage4_distance/stage4_window_distance.tsv'); mi=read(ROOT/'empirical/house_mouse_t_complex/results/t_complex_inversion_intervals.tsv'); mouse_refs=[(int(r['start_bp']),int(r['end_bp'])) for r in mi]; mouse=[dict(r,start=r['start_bp'],end=r['end_bp'],midpoint=r['midpoint_bp'],topology=r['q_t_alt'],continuous=r['A_distance'],blind=r['mean_pair_distance_no_spretus']) for r in m];systems['House mouse']={'groups':{'chr17':mouse},'refs':{'chr17':mouse_refs},'units':'IQ-TREE ML patristic substitution units'}
 return systems
def main():
 systems=load_systems(); audits=[]; inventories=[];fixed=[];multi=[];accuracy=[];spnull=[];cpnull=[];robust=[];valid=[]
 for system,s in systems.items():
  for group,rows in s['groups'].items():
   rows=sorted(rows,key=lambda r:float(r['midpoint'])); refs=s['refs'][group]; audits.append({'system':system,'group':group,'n_windows':len(rows),'first_start':rows[0]['start'],'last_end':rows[-1]['end'],'n_reference_intervals':len(refs),'units':s['units'],'missing_coordinate_policy':'retain observed rows; no interpolation'})
   inventories.append({'system':system,'group':group,'feature':'topology','source':'frozen local quartet/topology track','units':'dimensionless','n_nonmissing':sum(good(r.get('topology','')) for r in rows)})
   inventories.append({'system':system,'group':group,'feature':'continuous','source':'arrangement-informed A_distance','units':s['units'],'n_nonmissing':sum(good(r.get('continuous','')) for r in rows)})
   inventories.append({'system':system,'group':group,'feature':'arrangement_blind','source':'mean pairwise distance/time summary','units':s['units'],'n_nonmissing':sum(good(r.get('blind','')) for r in rows)})
   configs=[('topology',['topology']),('continuous',['continuous']),('joint',['topology','continuous']),('arrangement_blind',['blind'])]
   for method,fs in configs:
    z=[r for r in rows if all(good(r.get(f,'')) for f in fs)]; vals=np.array([[num(r[f]) for f in fs] for r in z]); seg=segment_track(vals,fs,min_size=3,max_intervals=1); pred=[]
    for x in intervals_to_coordinates(seg,z):pred.append((int(float(x['start'])),int(float(x['end']))))
    ev=eval_intervals(pred,refs,z); fixed.append({'system':system,'group':group,'method':method,'status':seg.status,'n_predicted_intervals':len(pred),'predicted_intervals':';'.join(f'{a}-{b}' for a,b in pred),'score':seg.score,'raw_improvement':seg.raw_improvement,'penalty':seg.penalty,'n_segments':seg.n_segments,'feature_contributions':json.dumps(seg.feature_contributions),'iou_union':ev['iou_union'],'f1':ev['f1'],'precision':ev['precision'],'recall':ev['recall'],'false_positive_bp':ev['false_positive_bp'],'mean_boundary_error_kb':ev['mean_boundary_error_kb']})
   for f in ('topology','continuous','blind'):
    obs,p,n=shifts(rows,f,refs,999);spnull.append({'system':system,'group':group,'feature':f,'observed_inside_minus_outside':obs,'p_circular_shift':p,'n_replicates':n,'hypothesis':'alignment with frozen reference interval(s)'})
    rng=random.Random(SEED); reps=19 if len(rows)>1000 else 99; scores=[]
    for _ in range(reps):
     perm=rows[:];rng.shuffle(perm);scores.append(quick_change_score(perm,f))
    observed=quick_change_score(rows,f);cpnull.append({'system':system,'group':group,'feature':f,'observed_max_score':observed,'n_replicates':reps,'p_block_permutation':(1+sum(x>=observed for x in scores))/(len(scores)+1),'hypothesis':'existence of a scanned one-interval change point; heuristic permutation null'})
   # fixed robustness summaries: penalty multiplier and minimum segment size, no coordinates used.
   for min_size in (3,5,10):
    seg=segment_track(np.array([[num(r['continuous'])] for r in rows if good(r.get('continuous',''))]),['continuous'],min_size=min_size,max_intervals=1);robust.append({'system':system,'group':group,'method':'continuous','min_size':min_size,'status':seg.status,'n_intervals':len(seg.intervals),'score':seg.score})
   z=[r for r in rows if good(r.get('continuous',''))];mseg=segment_track(np.array([[num(r['continuous'])] for r in z]),['continuous'],min_size=3,max_intervals=4);mpred=[(int(float(x['start'])),int(float(x['end']))) for x in intervals_to_coordinates(mseg,z)];mev=eval_intervals(mpred,refs,z)
   multi.append({'system':system,'group':group,'method':'continuous','model':'exact_dynamic_programming_over_0_to_4_intervals','status':mseg.status,'n_predicted_intervals':len(mpred),'predicted_intervals':';'.join(f'{a}-{b}' for a,b in mpred),'score':mseg.score,'raw_improvement':mseg.raw_improvement,'penalty':mseg.penalty,'n_segments':mseg.n_segments,'feature_contributions':json.dumps(mseg.feature_contributions),'iou_union':mev['iou_union'],'f1':mev['f1'],'precision':mev['precision'],'recall':mev['recall'],'false_positive_bp':mev['false_positive_bp'],'mean_boundary_error_kb':mev['mean_boundary_error_kb']})
   valid.append({'system':system,'group':group,'validation':'fixed common detector; no reference coordinates used in fitting','status':'complete'})
 fields=list(audits[0]);write(OUT/'input_audit.tsv',audits,fields);write(OUT/'feature_inventory.tsv',inventories,list(inventories[0]));write(OUT/'fixed_model_predictions.tsv',fixed,list(fixed[0]));write(OUT/'multi_interval_predictions.tsv',multi,list(multi[0]));write(OUT/'boundary_accuracy.tsv',fixed,list(fixed[0]));write(OUT/'spatial_null_results.tsv',spnull,list(spnull[0]));write(OUT/'change_point_null_results.tsv',cpnull,list(cpnull[0]));write(OUT/'robustness_summary.tsv',robust,list(robust[0]));write(OUT/'cross_system_validation.tsv',valid,list(valid[0]));
 # Main summary by system and method.
 summary=[]
 for system in systems:
  for method in ('topology','continuous','joint','arrangement_blind'):
   x=[r for r in fixed if r['system']==system and r['method']==method];summary.append({'system':system,'method':method,'mean_iou':float(np.mean([num(r['iou_union']) for r in x])),'mean_f1':float(np.mean([num(r['f1']) for r in x])),'detected_fraction':sum(r['status']!='no_detection' for r in x)/len(x),'groups':len(x)})
 write(OUT/'main_summary_table.tsv',summary,list(summary[0]))
 # Figure 1: one row per system, tracks on observed coordinates; gaps remain gaps.
 fig,axes=plt.subplots(4,1,figsize=(13,11),constrained_layout=True); colors={'topology':'#405a8d','continuous':'#c44e52'}
 for ax,(system,s) in zip(axes,systems.items()):
  group=list(s['groups'])[0];rows=s['groups'][group];x=np.array([float(r['midpoint'])/1e6 for r in rows]);ax.plot(x,[num(r['topology']) for r in rows],'.',ms=1.5,color=colors['topology'],label='topology');ax2=ax.twinx();ax2.plot(x,[num(r['continuous']) for r in rows],'.',ms=1.5,color=colors['continuous'],label='continuous');ax.set_title(f'{system} ({s["units"]})',loc='left');ax.set_ylabel('topology');ax2.set_ylabel('continuous');
  for a,b in s['refs'][group]:ax.axvspan(a/1e6,b/1e6,color='#999999',alpha=.12)
 axes[-1].set_xlabel('genomic position (Mb)');fig.savefig(FIG/'stage13_spatial_signals.pdf');fig.savefig(FIG/'stage13_spatial_signals.png',dpi=220);plt.close(fig)
 # Figure 2 accuracy.
 fig,ax=plt.subplots(1,2,figsize=(12,5),constrained_layout=True);systems_order=list(systems);methods=['topology','continuous','joint','arrangement_blind'];xx=np.arange(len(systems_order));w=.2
 for j,method in enumerate(methods):
  vals=[np.mean([num(r['iou_union']) for r in fixed if r['system']==sy and r['method']==method]) for sy in systems_order];ax[0].bar(xx+(j-1.5)*w,vals,w,label=method)
  vals=[np.mean([num(r['f1']) for r in fixed if r['system']==sy and r['method']==method]) for sy in systems_order];ax[1].bar(xx+(j-1.5)*w,vals,w,label=method)
 for a,title in zip(ax,['union IoU','window F1']):a.set_xticks(xx,systems_order,rotation=25,ha='right');a.set_ylim(0,1);a.set_title(title);a.legend(fontsize=7)
 fig.savefig(FIG/'stage13_boundary_accuracy.pdf');fig.savefig(FIG/'stage13_boundary_accuracy.png',dpi=220);plt.close(fig)
 # Figure 3 case studies: cod LG07, Anopheles, mouse, fire ants.
 fig,axes=plt.subplots(4,1,figsize=(13,10),constrained_layout=True)
 cases=[('Atlantic cod','LG07'),('Anopheles','2L'),('House mouse','chr17'),('Fire ants','chr16')]
 for ax,(sy,gr) in zip(axes,cases):
  rs=systems[sy]['groups'][gr];x=np.array([float(r['midpoint'])/1e6 for r in rs]);ax.plot(x,[num(r['topology']) for r in rs],'.',ms=1.2,label='topology');ax.plot(x,[num(r['continuous']) for r in rs],'.',ms=1.2,label='continuous');
  for a,b in systems[sy]['refs'][gr]:ax.axvspan(a/1e6,b/1e6,color='#999999',alpha=.12)
  ax.set_title(f'{sy} {gr}',loc='left');ax.set_ylabel('signal')
 axes[-1].set_xlabel('genomic position (Mb)');axes[0].legend();fig.savefig(FIG/'stage13_case_studies.pdf');fig.savefig(FIG/'stage13_case_studies.png',dpi=220);plt.close(fig)
 inputs=[]
 for s in systems.values():
  for rs in s['groups'].values():
   pass
 manifest={'seed':SEED,'detector':'empirical/shared/spatial_boundary_detector.py','penalty':'2 log(n) per added free-mean segment after robust MAD scaling','min_size':3,'max_intervals':4,'input_hashes':{str(p.relative_to(ROOT)):md5(p) for p in [ROOT/'empirical/atlantic_cod/results/stage7_window_divergence_time_signal.tsv',ROOT/'data/atlantic_cod/processed/stage4a_window_topology_scores.tsv',ROOT/'empirical/fire_ants_chr16/results/stage9_distance/stage9_window_patristic_distance.tsv',ROOT/'empirical/anopheles_2la/results/stage12_distance/stage12_window_distance.tsv',ROOT/'empirical/house_mouse_t_complex/results/stage4_distance/stage4_window_distance.tsv']}}
 (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 txt='''# Stage 13 unified spatial boundary detection\n\nThis analysis applies one robustly standardized, penalized piecewise-constant detector to frozen topology and continuous tracks from Atlantic cod, fire ants, Anopheles, and house mouse. The fixed penalty is `2 log(n)` per added free-mean segment after MAD scaling. Exact dynamic programming optimizes each fixed segment count; the selected number of intervals is chosen from zero through four. Observed rows retain their coordinates and are not interpolated across gaps.\n\nThe circular-shift tests assess alignment with the frozen reference intervals. The block-permutation change-point tests assess existence of a scanned change point and are reported separately; they are heuristic under strong spatial autocorrelation. Branch-length quantities are kept in their source units. Atlantic cod uses SNAPP population-tree time summaries; the other distance tracks use substitution-distance units.\n''';(OUT/'methods_text.md').write_text(txt);(OUT/'report.md').write_text(txt);(OUT/'results_text.md').write_text(txt)
if __name__=='__main__':main()
