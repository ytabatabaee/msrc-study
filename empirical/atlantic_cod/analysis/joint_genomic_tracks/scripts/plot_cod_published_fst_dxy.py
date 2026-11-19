#!/usr/bin/env python3
"""Plot Matschiner et al. published gadMor2 F_ST and d_XY source tracks."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[5]
AN = ROOT / 'empirical/atlantic_cod/analysis/joint_genomic_tracks'
RESULTS = AN / 'results'
FIG = AN / 'figures'

def read_tsv(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def plot(metric, column, ylabel, filename, title, color):
    rows = [r for r in read_tsv(RESULTS/'atlantic_cod_published_edfig3_stats.tsv') if r['reference']=='gadMor2']
    bounds = {r['lg']:(int(r['start']), int(r['end'])) for r in read_tsv(ROOT/'data/atlantic_cod/metadata/region_manifest.tsv')}
    fig, axes = plt.subplots(4, 1, figsize=(11, 8.2), sharex=False)
    fig.subplots_adjust(left=.09, right=.985, top=.92, bottom=.085, hspace=.32)
    for ax, lg in zip(axes, ('LG01','LG02','LG07','LG12')):
        data = sorted((r for r in rows if r['lg']==lg), key=lambda r:int(r['window_center']))
        x = [int(r['window_center'])/1e6 for r in data]
        y = [float(r[column]) for r in data]
        lo, hi = bounds[lg]
        ax.axvspan(lo/1e6, hi/1e6, color='#8c79b8', alpha=.16, lw=0, zorder=0)
        # Points only: source windows can be missing, so connecting lines could imply interpolation.
        ax.scatter(x, y, s=11, color=color, linewidths=0, zorder=2)
        ax.set_ylabel(ylabel)
        ax.set_title(lg, loc='left', fontsize=10, fontweight='bold')
        ax.grid(axis='y', color='0.88', linewidth=.6)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        ax.set_xlim(0, max(x)+.15 if x else 1)
    axes[-1].set_xlabel('gadMor2 linkage group coordinate (Mb)')
    fig.suptitle(title, fontsize=13, fontweight='bold')
    FIG.mkdir(parents=True, exist_ok=True)
    for ext in ('pdf','svg','png'):
        fig.savefig(FIG/f'{filename}.{ext}', dpi=300 if ext=='png' else None)
    plt.close(fig)

plot('f_st','f_st',r"$F_{ST}$",'atlantic_cod_published_fst_tracks',r'Atlantic cod: published $F_{ST}$ by linkage group','#246a8d')
plot('d_xy','d_xy_physical_window_snp_ascertained',r"Published $d_{XY}^{*}$",'atlantic_cod_published_dxy_tracks',r'Atlantic cod: published SNP-ascertained $d_{XY}^{*}$ by linkage group','#b45538')
