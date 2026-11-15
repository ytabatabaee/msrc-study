#!/usr/bin/env python3
"""Audit and analyze Short & Streisfeld monkeyflower TWISST weights."""
from __future__ import annotations
import csv, hashlib, math
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
RAW = ROOT / 'empirical/monkeyflowers/data/raw'
OUT = ROOT / 'empirical/monkeyflowers/results'
FIG = ROOT / 'empirical/monkeyflowers/figures'
DATA_NAME = 'ari.par.isl_lon.lon.cle.output.phyml_w_100000.data.tsv'
WEIGHT_NAME = 'ari.par.isl_lon.lon.cle.output.weights.w_100000.csv'
COMBINED_NAME = 'ari.par.isl_lon.lon.cle.output.weights.w_100000.window_data_100kb_weights.csv'
COORD_NAME = '100kb_win_data_nonoverlap_plos_bio_coords.csv'

# All genomic intervals in these inputs are 1-based inclusive on the M. aurantiacus
# reference assembly used by Stankowski et al. (2019).
# q1 = species history (topo3); q2 = recent introgression set (topo9+12+13);
# q3 = ancient introgression topology (topo10). Other 11 classes are retained.
SPECIES = ('topo3',)
RECENT = ('topo9', 'topo12', 'topo13')
ANCIENT = ('topo10',)
ALL_TOPO = tuple(f'topo{i}' for i in range(1, 16))
C5_PRIMARY = (24_700_001, 24_900_000)  # two 100-kb intervals explicitly noted in the paper
C5_BROAD = (24_700_001, 25_000_000)   # ~300 kb region; approximate boundary sensitivity


def read_delimited(path, delimiter='\t'):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def key(row):
    return (row['scaffold'], int(row['start']), int(row['end']))


def topology_class(topology):
    if topology in SPECIES: return 'species'
    if topology in RECENT: return 'recent'
    if topology in ANCIENT: return 'ancient'
    return 'other'


def validate_coordinate_keys(sequences):
    """Require unique, identical (scaffold,start,end) keys across coordinate files."""
    keysets=[]
    for rows in sequences:
        values=list(map(key,rows))
        if len(values)!=len(set(values)):
            raise ValueError('Duplicate genomic coordinate keys')
        keysets.append(set(values))
    if any(k!=keysets[0] for k in keysets[1:]):
        raise ValueError('Coordinate-key mismatch across coordinate-bearing tables')
    return keysets[0]


def read_weight_file(path):
    """Read TWISST count CSV with four topology-definition comment rows."""
    with open(path, newline='') as f:
        lines = f.readlines()
    defs = {}
    for line in lines:
        if not line.startswith('#topo'):
            break
        label, tree = line[1:].strip().split(' ', 1)
        defs[label] = tree.rstrip(';')
    rows = list(csv.DictReader(lines[len(defs):], delimiter='\t'))
    parsed = []
    for r in rows:
        vals = {c: float(r[c]) for c in ALL_TOPO}
        total = sum(vals.values())
        if total <= 0 or not math.isfinite(total):
            raise ValueError(f'Invalid topology count total: {total}')
        parsed.append({c: v / total for c, v in vals.items()})
    return defs, parsed


def validate_and_load():
    data = read_delimited(RAW / DATA_NAME)
    # Despite its .csv suffix, this deposited coordinate file is tab-delimited.
    coords = read_delimited(RAW / COORD_NAME, '\t')
    combined = read_delimited(RAW / COMBINED_NAME, ',')
    defs, raw_weights = read_weight_file(RAW / WEIGHT_NAME)
    coordinate_keys = validate_coordinate_keys((data, coords, combined))
    if len(raw_weights) != len(combined):
        raise ValueError('Count-only weight file and joined export differ in row count')
    bykey = {key(r): r for r in data}
    cbykey = {key(r): r for r in combined}
    # The deposited count-only file has no genomic coordinates. The combined
    # table is the coordinate-bearing TWISST export; validate the count file
    # against every combined row before assigning any coordinate to it.
    if len(raw_weights) != len(combined):
        raise ValueError('Count-only topology file and coordinate-bearing table differ in row count')
    for i, c in enumerate(combined):
        raw = raw_weights[i]
        if max(abs(float(c[name]) - raw[name]) for name in ALL_TOPO) > 2e-6:
            raise ValueError(f'Count-only weights disagree with coordinate-bearing export at row {i+1}')
    rows = []
    for k in sorted(coordinate_keys, key=lambda x: (x[0], x[1], x[2])):
        counts = {name: float(cbykey[k][name]) for name in ALL_TOPO}
        c = cbykey[k]
        vals_combined = {name: float(c[name]) for name in ALL_TOPO}
        if abs(sum(vals_combined.values()) - 1) > 1e-7:
            raise ValueError(f'Combined weights do not sum to one at {k}')
        if max(abs(vals_combined[name] - counts[name]) for name in ALL_TOPO) > 2e-6:
            raise ValueError(f'Combined weights disagree with normalized raw counts at {k}')
        r = dict(bykey[k])
        r.update(counts)
        r['q1_species'] = sum(counts[x] for x in SPECIES)
        r['q2_recent'] = sum(counts[x] for x in RECENT)
        r['q3_ancient'] = sum(counts[x] for x in ANCIENT)
        r['other_weight'] = sum(counts[x] for x in ALL_TOPO if x not in set(SPECIES+RECENT+ANCIENT))
        r['D_recent_minus_species'] = r['q2_recent'] - r['q1_species']
        r['entropy15'] = -sum(v * math.log(v) for v in counts.values() if v > 0) / math.log(15)
        r['dominant_topology'] = max(ALL_TOPO, key=counts.get)
        r['dominant_class'] = topology_class(r['dominant_topology'])
        r['is_informative'] = int(float(r['sites']) > 0 and all(math.isfinite(v) for v in counts.values()))
        rows.append(r)
    return data, coords, combined, defs, rows


def physical_runs(rows, field, threshold):
    """Runs of windows with weight > threshold, requiring exact physical adjacency."""
    output = []
    bychr = defaultdict(list)
    for r in rows:
        bychr[r['scaffold']].append(r)
    for chrom, cr in bychr.items():
        cr.sort(key=lambda x: int(x['start']))
        run = []
        def finish():
            if not run: return
            output.append({'chromosome': chrom, 'start': int(run[0]['start']),
                           'end': int(run[-1]['end']),
                           'run_length_bp': int(run[-1]['end']) - int(run[0]['start']) + 1,
                           'n_windows': len(run), 'mean_weight': float(np.mean([x[field] for x in run])),
                           'threshold': threshold, 'support_field': field})
            run.clear()
        for r in cr:
            good = float(r[field]) > threshold and int(r['is_informative']) == 1
            adjacent = not run or int(r['start']) == int(run[-1]['end']) + 1
            if good and adjacent:
                run.append(r)
            else:
                finish()
                if good: run.append(r)
        finish()
    return output


def transition_stats(rows, field, threshold):
    ntrans = 0; span = 0
    bychr = defaultdict(list)
    for r in rows: bychr[r['scaffold']].append(r)
    for cr in bychr.values():
        cr.sort(key=lambda r: int(r['start']))
        for a,b in zip(cr,cr[1:]):
            if int(a['end']) + 1 != int(b['start']): continue
            if not (a['is_informative'] and b['is_informative']): continue
            span += int(b['start']) - int(a['start'])
            ntrans += int((float(a[field]) > threshold) != (float(b[field]) > threshold))
    return ntrans, span / 1e6, ntrans / (span / 1e6) if span else float('nan')


def autocorrelation(rows, field, max_dist=2_000_000, step=100_000):
    """Descriptive chromosome-stratified Pearson correlations by physical lag bins."""
    bychr = defaultdict(list)
    for r in rows:
        try: finite_value=math.isfinite(float(r[field]))
        except (TypeError,ValueError): finite_value=False
        if r['is_informative'] and finite_value: bychr[r['scaffold']].append(r)
    bins = defaultdict(list)
    for chrom, cr in bychr.items():
        cr.sort(key=lambda r:int(r['start']))
        segments=[]; segment=[]
        for r in cr:
            if segment and int(r['start']) != int(segment[-1]['end'])+1:
                segments.append(segment); segment=[]
            segment.append(r)
        if segment: segments.append(segment)
        for seg in segments:
          for i,a in enumerate(seg):
            xa=float(a[field]); ma=(int(a['start'])+int(a['end']))/2
            for b in seg[i+1:]:
                mb=(int(b['start'])+int(b['end']))/2; d=mb-ma
                if d > max_dist: break
                bi=max(1,int(round(d/step)))
                bins[bi].append((xa,float(b[field])))
    result=[]
    for bi,pairs in sorted(bins.items()):
        x=np.array([p[0] for p in pairs]); y=np.array([p[1] for p in pairs])
        corr=float(np.corrcoef(x,y)[0,1]) if len(pairs)>=3 and np.std(x)>0 and np.std(y)>0 else float('nan')
        result.append({'distance_start_bp':max(0,(bi-.5)*step), 'distance_end_bp':(bi+.5)*step,
                       'mid_distance_bp':bi*step,'n_pairs':len(pairs),'correlation':corr})
    return result


def write_tsv(path, rows, fields=None):
    if fields is None: fields=list(rows[0]) if rows else []
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore')
        w.writeheader()
        for row in rows:
            w.writerow({k: ('NA' if isinstance(row.get(k),float) and not math.isfinite(row[k]) else row.get(k,'')) for k in fields})


def save_figure(fig, stem, dpi=300, bbox=None):
    for ext in ('png','pdf','svg'):
        path=FIG/f'{stem}.{ext}'
        kwargs={'dpi':dpi}
        if bbox: kwargs['bbox_inches']=bbox
        fig.savefig(path,**kwargs)
        if ext=='svg':
            # Matplotlib emits harmless trailing spaces in multi-line SVG path data.
            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines())+'\n')


def summarize(rows, label):
    vals={f:np.array([float(r[f]) for r in rows],dtype=float) for f in ('q1_species','q2_recent','q3_ancient','other_weight','D_recent_minus_species','entropy15')}
    n=len(rows)
    dom={x:sum(r['dominant_class']==x for r in rows) for x in ('species','recent','ancient','other')}
    return {'subset':label,'n_windows':n,'mean_q1_species':float(vals['q1_species'].mean()) if n else float('nan'),
            'mean_q2_recent_set':float(vals['q2_recent'].mean()) if n else float('nan'),
            'mean_q3_ancient':float(vals['q3_ancient'].mean()) if n else float('nan'),
            'mean_other_11_topologies':float(vals['other_weight'].mean()) if n else float('nan'),
            'mean_D_recent_minus_species':float(vals['D_recent_minus_species'].mean()) if n else float('nan'),
            'mean_entropy15':float(vals['entropy15'].mean()) if n else float('nan'),
            'dominant_species_windows':dom['species'],'dominant_recent_windows':dom['recent'],
            'dominant_ancient_windows':dom['ancient'],'dominant_other_windows':dom['other']}


def plot_figures(rows, runs):
    plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':9,'xtick.labelsize':7,'ytick.labelsize':7,'svg.fonttype':'none'})
    colors={'q1_species':'#b2182b','q2_recent':'#ef8a62','q3_ancient':'#1b7837'}
    bychr=defaultdict(list)
    for r in rows: bychr[r['scaffold']].append(r)
    chroms=sorted(bychr,key=lambda x:int(''.join(c for c in x if c.isdigit()) or 0))
    fig,axes=plt.subplots(len(chroms),1,figsize=(8,1.05*len(chroms)),sharey=True,constrained_layout=True)
    for ax,ch in zip(np.atleast_1d(axes),chroms):
        rr=sorted(bychr[ch],key=lambda r:int(r['start']))
        segments=[]; segment=[]
        for r in rr:
            if segment and int(r['start']) != int(segment[-1]['end'])+1:
                segments.append(segment); segment=[]
            segment.append(r)
        if segment: segments.append(segment)
        for si,seg in enumerate(segments):
            x=np.array([(int(r['start'])+int(r['end']))/2/1e6 for r in seg])
            for f in colors:
                y=np.array([float(r[f]) for r in seg]); ax.plot(x,y,lw=.8,color=colors[f],label=f.replace('_',' ') if si==0 else None)
        ax.set_ylabel(ch); ax.set_ylim(0,1); ax.grid(axis='y',lw=.3,alpha=.3)
    axes[-1].set_xlabel('Position on M. aurantiacus reference (Mb)')
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,ncol=3,loc='upper center',bbox_to_anchor=(.5,1.035),frameon=False)
    save_figure(fig,'figure1_chromosome_topology_tracks',bbox='tight')
    plt.close(fig)

    fig,axes=plt.subplots(2,1,figsize=(7.2,4.2),sharex=True,constrained_layout=True)
    rr=sorted([r for r in rows if r['scaffold']=='LG5' and int(r['end'])>=24_300_000 and int(r['start'])<=25_300_000],key=lambda r:int(r['start']))
    x=np.array([(int(r['start'])+int(r['end']))/2/1e6 for r in rr])
    for f in colors: axes[0].plot(x,[r[f] for r in rr],marker='.',ms=2,lw=.9,color=colors[f],label=f.replace('_',' '))
    axes[0].axvspan(C5_BROAD[0]/1e6,C5_BROAD[1]/1e6,color='#fdae61',alpha=.17,label='approx. 300-kb candidate')
    axes[0].axvspan(C5_PRIMARY[0]/1e6,C5_PRIMARY[1]/1e6,color='#fdae61',alpha=.35,label='two reported 100-kb windows')
    axes[0].set_ylabel('Topology weight'); axes[0].set_ylim(0,1); axes[0].legend(ncol=2,frameon=False)
    axes[1].plot(x,[r['D_recent_minus_species'] for r in rr],color='#762a83',lw=1)
    axes[1].axhline(0,color='#555',lw=.6); axes[1].axvspan(C5_BROAD[0]/1e6,C5_BROAD[1]/1e6,color='#fdae61',alpha=.17); axes[1].axvspan(C5_PRIMARY[0]/1e6,C5_PRIMARY[1]/1e6,color='#fdae61',alpha=.35)
    axes[1].set_ylabel('D = recent − species'); axes[1].set_xlabel('LG5 position (Mb)')
    save_figure(fig,'figure2_chromosome5_topology')
    plt.close(fig)

    fig,ax=plt.subplots(figsize=(6.6,4.2),constrained_layout=True)
    for field,label,color in [('q2_recent','Recent topology set','#ef8a62'),('q3_ancient','Ancient topology','#1b7837')]:
        rruns=[r for r in runs if r['support_field']==field and r['threshold']==.5]
        if rruns:
            xx=np.sort(np.array([r['run_length_bp']/1e6 for r in rruns])); yy=np.arange(1,len(xx)+1)/len(xx)
            ax.step(xx,yy,where='post',label=f'{label} (n={len(rruns)})',color=color,lw=1.4)
    ax.set_xlabel('Contiguous run length (Mb)');ax.set_ylabel('Cumulative fraction of runs');ax.legend(frameon=False)
    save_figure(fig,'figure3_spatial_run_lengths')
    plt.close(fig)

    # Cross-system 100-kb block comparison, direct feature definitions noted in caption/report.
    neuro=[]
    p=ROOT/'empirical/neurospora/results/lgi_100kb_blocks.tsv'
    if p.exists(): neuro=read_delimited(p)
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.7),constrained_layout=True)
    for field,label,color in [('q2_recent','Monkeyflower recent set','#ef8a62'),('q3_ancient','Monkeyflower ancient','#1b7837')]:
        rruns=[r for r in runs if r['support_field']==field and r['threshold']==.5]
        if rruns:
            xx=np.sort([r['run_length_bp']/1e6 for r in rruns]); axes[0].step(xx,np.arange(1,len(xx)+1)/len(xx),where='post',label=label,color=color)
    nruns=[]
    if neuro:
        # The established Neurospora block file has 100-kb means; use q3 as its documented introgression-associated topology.
        vals=[]
        for r in neuro:
            vals.append({'scaffold':'LGI','start':int(r['start'])+1,'end':int(r['end'])+1,
                         'is_informative':int(float(r.get('n_windows',0))>0),'q3n':float(r['q3'])})
        nruns=physical_runs(vals,'q3n',.5)
        if nruns:
            xx=np.sort([r['run_length_bp']/1e6 for r in nruns]); axes[0].step(xx,np.arange(1,len(xx)+1)/len(xx),where='post',label='Neurospora q3',color='#2166ac')
    axes[0].set_xlabel('Run length (Mb), 100-kb blocks');axes[0].set_ylabel('Cumulative fraction');axes[0].legend(frameon=False,fontsize=7)
    acm=autocorrelation(rows,'q3_ancient',max_dist=2_000_000)
    if acm: axes[1].plot([r['mid_distance_bp']/1e6 for r in acm],[r['correlation'] for r in acm],label='Monkeyflower ancient',color='#1b7837')
    if neuro:
        # Calculate correlogram for the existing 100-kb block means.
        na=[]
        for r in neuro:
            na.append({'scaffold':'LGI','start':int(r['start'])+1,'end':int(r['end'])+1,'is_informative':int(float(r.get('n_windows',0))>0),'q3':float(r['q3'])})
        acn=autocorrelation(na,'q3',max_dist=2_000_000)
        if acn: axes[1].plot([r['mid_distance_bp']/1e6 for r in acn],[r['correlation'] for r in acn],label='Neurospora q3',color='#2166ac')
    axes[1].set_xlabel('Physical separation (Mb)');axes[1].set_ylabel('Pearson correlation');axes[1].legend(frameon=False,fontsize=7)
    save_figure(fig,'figure4_monkeyflower_neurospora')
    plt.close(fig)


def run():
    OUT.mkdir(parents=True,exist_ok=True); FIG.mkdir(parents=True,exist_ok=True)
    data,coords,combined,defs,rows=validate_and_load()
    write_tsv(OUT/'topology_weight_tracks.tsv',rows)
    run_records=[]
    for f in ('q2_recent','q3_ancient'):
        for t in (.5,.7,.9): run_records.extend(physical_runs(rows,f,t))
    write_tsv(OUT/'spatial_run_lengths.tsv',run_records)
    ac=[]
    for f in ('q2_recent','q3_ancient','D_recent_minus_species'):
        for r in autocorrelation(rows,f): ac.append({'metric':f,**r})
    write_tsv(OUT/'autocorrelation_summary.tsv',ac)
    trans=[]
    for f in ('q2_recent','q3_ancient'):
        for t in (.5,.7,.9):
            n,mb,rate=transition_stats(rows,f,t);trans.append({'support_field':f,'threshold':t,'transitions':n,'adjacent_covered_Mb':mb,'transitions_per_Mb':rate})
    write_tsv(OUT/'topology_transition_summary.tsv',trans)
    inside=[]; flanks=[]; chr5=[]; outside=[]
    for r in rows:
        if r['scaffold']=='LG5':
            chr5.append(r)
            if int(r['start'])>=C5_PRIMARY[0] and int(r['end'])<=C5_PRIMARY[1]: inside.append(r)
            elif int(r['end'])==C5_PRIMARY[0]-1 or int(r['start'])==C5_PRIMARY[1]+1: flanks.append(r)
        if not (r['scaffold']=='LG5' and int(r['start'])>=C5_PRIMARY[0] and int(r['end'])<=C5_PRIMARY[1]): outside.append(r)
    c5=[]
    for label,subset in [('chr5_candidate_24.7_24.9Mb',inside),('immediate_flanks',flanks),('other_LG5', [r for r in chr5 if r not in inside and r not in flanks]),('genome_outside_candidate',outside)]:
        c5.append({**summarize(subset,label),'mean_entropy15':float(np.mean([r['entropy15'] for r in subset])) if subset else float('nan')})
    broad=[r for r in rows if r['scaffold']=='LG5' and int(r['start'])>=C5_BROAD[0] and int(r['end'])<=C5_BROAD[1]]
    outside_broad=[r for r in rows if r not in broad]
    c5.append({**summarize(broad,'chr5_approx_300kb_candidate_24.7_25.0Mb'),'mean_entropy15':float(np.mean([r['entropy15'] for r in broad])) if broad else float('nan')})
    write_tsv(OUT/'chromosome5_region_comparison.tsv',c5)
    # Genomewide removal sensitivity: independently reported 200-kb window pair and 300-kb sensitivity.
    # Genomewide removal sensitivity: 200-kb focal interval and broad 300-kb sensitivity.
    sens=[summarize(rows,'all_windows'),summarize(outside,'excluding_published_200kb_candidate'),summarize(outside_broad,'excluding_approx_300kb_candidate'),summarize([r for r in rows if r['dominant_class']=='species'],'species_dominant'),summarize([r for r in rows if r['dominant_class']=='recent'],'recent_dominant'),summarize([r for r in rows if r['dominant_class']=='ancient'],'ancient_dominant'),summarize(broad,'approx_300kb_candidate_sensitivity')]
    write_tsv(OUT/'genomewide_topology_sensitivity.tsv',sens)
    # Compare fragmentation at common 100-kb block resolution with Neurospora's committed blocks.
    neuro_path=ROOT/'empirical/neurospora/results/lgi_100kb_blocks.tsv'; neuro=read_delimited(neuro_path)
    nr=[]
    for r in neuro:
        nr.append({'scaffold':'LGI','start':int(r['start'])+1,'end':int(r['end'])+1,'is_informative':int(float(r.get('n_windows',0))>0),'q1':float(r['q1']),'q2':float(r['q2']),'q3':float(r['q3'])})
    nrun=physical_runs(nr,'q3',.5)
    nrun70=physical_runs(nr,'q3',.7); nrun90=physical_runs(nr,'q3',.9)
    comparison=[]
    for system,field,rlist in [('monkeyflower_recent','q2_recent',run_records),('monkeyflower_ancient','q3_ancient',run_records),('neurospora_q3','q3',nrun+nrun70+nrun90)]:
        selected=[r for r in rlist if r['support_field']==field and r['threshold']==.5]
        if system=='neurospora_q3':
            ntrans,nmb,nrate=transition_stats(nr,'q3',.5)
            acsys=autocorrelation(nr,'q3',max_dist=2_000_000)
        else:
            ntrans,nmb,nrate=transition_stats(rows,field,.5)
            acsys=autocorrelation(rows,field,max_dist=2_000_000)
        lag1=next((z['correlation'] for z in acsys if z['mid_distance_bp']==100000),float('nan'))
        scale=next((z['mid_distance_bp'] for z in acsys if math.isfinite(z['correlation']) and math.isfinite(lag1) and z['correlation'] <= lag1/math.e),float('nan'))
        comparison.append({'system_metric':system,'window_resolution_bp':100000,'threshold':.5,'n_runs':len(selected),'total_run_bp':sum(r['run_length_bp'] for r in selected),'median_run_bp':float(np.median([r['run_length_bp'] for r in selected])) if selected else float('nan'),'longest_run_bp':max([r['run_length_bp'] for r in selected],default=0),'max_run_windows':max([r['n_windows'] for r in selected],default=0),'transitions':ntrans,'adjacent_covered_Mb':nmb,'transitions_per_Mb':nrate,'lag_100kb_autocorrelation':lag1,'approx_decay_distance_bp':scale})
    write_tsv(OUT/'neurospora_spatial_comparison.tsv',comparison)
    plot_figures(rows,run_records)
    manifest=[]
    for name in (DATA_NAME,WEIGHT_NAME,COMBINED_NAME,COORD_NAME):
        p=RAW/name;b=p.read_bytes();manifest.append({'file':str(p.relative_to(ROOT)),'size_bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'role':{'data':'window likelihood/coordinates','weights':'15-topology count weights','combined':'joined coordinate and normalized topology weights','coords':'100kb nonoverlapping/recombination coordinate grid'}[{'ari.par.isl_lon.lon.cle.output.phyml_w_100000.data.tsv':'data','ari.par.isl_lon.lon.cle.output.weights.w_100000.csv':'weights','ari.par.isl_lon.lon.cle.output.weights.w_100000.window_data_100kb_weights.csv':'combined','100kb_win_data_nonoverlap_plos_bio_coords.csv':'coords'}[name]]})
    write_tsv(OUT/'input_manifest.tsv',manifest)
    # Return summary to create a human report from exact computed values.
    print('windows',len(rows),'scaffolds',len(set(r['scaffold'] for r in rows)),'informative',sum(r['is_informative'] for r in rows))
    for s in sens: print(s)
    print('run summary')
    for f in ('q2_recent','q3_ancient'):
        for t in (.5,.7,.9):
            rr=[r for r in run_records if r['support_field']==f and r['threshold']==t]
            print(f,t,'n',len(rr),'median',np.median([r['run_length_bp'] for r in rr]) if rr else None,'max',max([r['run_length_bp'] for r in rr],default=0),'maxwindows',max([r['n_windows'] for r in rr],default=0))
    print('c5',c5)
    print('neuro',comparison)
    print('autocorr sample',ac[:5])

if __name__=='__main__': run()
