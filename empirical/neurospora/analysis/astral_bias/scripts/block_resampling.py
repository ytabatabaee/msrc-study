#!/usr/bin/env python3
"""Linkage-group-stratified block bootstrap and matched block-removal controls."""
import argparse, csv, random, re, subprocess
from collections import defaultdict
from pathlib import Path
from Bio import Phylo
ROOT=Path(__file__).resolve().parents[5]
HERE=Path(__file__).resolve().parents[1]
COORD=ROOT/'data/neurospora/raw/neurospora/trees/cras_sit_hist_tet_L10.min20MAC1.phyml_GTR_bionj_w50sites.data.tsv.gz'
ASTRAL=Path('/Users/ytabatabaee/Desktop/ASTER/bin/astral4')

def fields(cl): return dict(re.findall(r'([A-Za-z0-9_]+)=([^;\]]+)',cl.name or ''))
def split_key(side,taxa):
 a=frozenset(side); b=frozenset(taxa)-a
 if len(a)>len(b) or (len(a)==len(b) and tuple(sorted(a))>tuple(sorted(b))): a=b
 return '|'.join(sorted(a))
def read_splits(path):
 t=Phylo.read(str(path),'newick'); taxa={x.name for x in t.get_terminals()}; out={}
 for cl in t.get_nonterminals():
  side={x.name for x in cl.get_terminals()}
  if len(side)>=2 and len(taxa-side)>=2:
   f=fields(cl); out[split_key(side,taxa)]={'cu':float(f.get('CULength',cl.branch_length or 0)), 'su':float(f.get('SULength',0))}
 return out
def astral_run(design, trees, label, seed):
 d=HERE/'logs'/'block_resampling'/design; d.mkdir(parents=True,exist_ok=True)
 inp=d/f'{label}.trees'; out=d/f'{label}.nwk'; log=d/f'{label}.log'
 inp.write_text('\n'.join(trees)+'\n')
 cmd=[str(ASTRAL),'-R','-u','2','-t','8','--seed',str(seed),'--length','CULength','--genelength','50','--sulengthtype','1','-a',str(HERE/'inputs'/design/'taxon_mapping.tsv'),'-i',str(inp),'-o',str(out),'-v','1']
 cp=subprocess.run(cmd,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True); log.write_text(cp.stdout)
 return read_splits(out)
def block_sets(coords, inside_only=None, scaffold_filter=None, block=250000):
 groups=defaultdict(list)
 for i,r in enumerate(coords):
  inn=r['scaffold']=='LGI' and 1500000<=int(r['mid'])<=8000000
  if inside_only is not None and inn!=inside_only: continue
  if scaffold_filter is True and r['scaffold']=='LGI': continue
  if scaffold_filter is False and r['scaffold']!='LGI': continue
  k=(r['scaffold'],int(r['mid'])//block)
  groups[k].append(i)
 return groups
def draw_blocks(groups,rng):
 bysc=defaultdict(list)
 for k in groups:bysc[k[0]].append(k)
 selected=[]
 for scaffold,keys in sorted(bysc.items()):
  selected.extend(rng.choice(keys) for _ in keys)
 return [i for k in selected for i in groups[k]]
def main():
 ap=argparse.ArgumentParser(); ap.add_argument('--reps',type=int,default=20); ap.add_argument('--block-bp',type=int,default=250000); ap.add_argument('--seed',type=int,default=88421); a=ap.parse_args()
 import gzip
 with gzip.open(COORD,'rt') as f: coords=list(csv.DictReader(f,delimiter='\t'))
 valid=[i for i,r in enumerate(coords) if (HERE/'inputs/species_level/all.trees').exists() and True]
 # Each prepared input contains only valid trees; retain source row indices by synchronously filtering NA Newicks.
 treesraw=[]
 treearc=ROOT/'data/neurospora/raw/neurospora/trees/cras_sit_hist_tet_L10.min20MAC1.phyml_GTR_bionj_w50sites.trees.gz'
 with gzip.open(treearc,'rt') as f: nwks=[x.rstrip('\r\n') for x in f]
 valid=[i for i,n in enumerate(nwks) if n not in ('NA','')]
 coords=[coords[i] for i in valid]
 rng=random.Random(a.seed); rows=[]
 for design in ('species_level','background_lineage'):
  with (HERE/'inputs'/design/'all.trees').open() as f: alltrees=[x.rstrip() for x in f]
  with (HERE/'inputs'/design/'region_excluded.trees').open() as f: outtrees=[x.rstrip() for x in f]
  if len(alltrees)!=len(coords): raise ValueError('coordinate and valid tree counts disagree')
  inside=[j for j,r in enumerate(coords) if r['scaffold']=='LGI' and 1500000<=int(r['mid'])<=8000000]
  outside=[j for j in range(len(coords)) if j not in set(inside)]
  outtrees_by_idx=dict(zip(outside,outtrees))
  # Baseline all-window reference branches
  base=read_splits(HERE/design/'all'/'astral_culength_type1.nwk')
  g_all=block_sets(coords,block=a.block_bp); g_out=block_sets(coords,inside_only=False,block=a.block_bp)
  g_in=block_sets(coords,inside_only=True,block=a.block_bp)
  for rep in range(a.reps):
   # Region contrast: independently resample genomic blocks within scaffolds for all and region-excluded sets.
   for condition,groups,source in [('boot_all',g_all,alltrees),('boot_excluded',g_out,outtrees_by_idx)]:
    inds=draw_blocks(groups,rng); result=astral_run(design,[source[j] for j in inds],f'{condition}_{rep:03d}',a.seed+rep*17+(0 if design=='species_level' else 100000))
    for split,b in result.items():
     bb=base.get(split,{})
     rows.append([design,'bootstrap_branch',rep,condition,a.block_bp,len(inds),None,None,split,b.get('cu'),b.get('su'),b.get('cu',0)-bb.get('cu',0) if bb else None,b.get('su',0)-bb.get('su',0) if bb else None,split in base])
   # Matched control: remove whole outside blocks until removed count approximates focal interval's valid-window count.
   keys=list(g_out); rng.shuffle(keys); removed=[]
   for key in keys:
    if len(removed)>=len(inside): break
    removed.extend(g_out[key])
   removed=set(removed); retained=[j for j in range(len(alltrees)) if j not in removed]
   result=astral_run(design,[alltrees[j] for j in retained],f'matched_remove_{rep:03d}',a.seed+rep*31+(1 if design=='species_level' else 200000))
   rows.append([design,'matched_removal',rep,'outside_blocks_removed',a.block_bp,len(retained),len(removed),len(inside),None,None,None,None])
   for split,b in result.items():
    bb=base.get(split,{}); rows.append([design,'matched_removal_branch',rep,'outside_blocks_removed',a.block_bp,len(retained),len(removed),len(inside),split,b.get('cu'),b.get('su'),b.get('cu',0)-bb.get('cu',0) if bb else None,b.get('su',0)-bb.get('su',0) if bb else None,split in base])
  # Note: bootstrap deltas are relative to full-data baseline; conditions remain separate so contrast can be paired by replicate.
 out=HERE/'results'; out.mkdir(exist_ok=True)
 with (out/'block_resampling.tsv').open('w') as f:
  w=csv.writer(f,delimiter='\t');w.writerow(['design','analysis','replicate','condition','block_bp','n_trees','n_removed','target_region_n','split','CU','SU_per_site','delta_CU_vs_all','delta_SU_vs_all','split_in_full_reference']);w.writerows(rows)
 (out/'block_resampling_parameters.json').write_text(__import__('json').dumps({'reps':a.reps,'block_bp':a.block_bp,'seed':a.seed,'focal_valid_windows':5767,'method':'within-scaffold moving-block bootstrap; blocks assigned by midpoint; random outside-block removal until at least focal count; ASTRAL-IV with same mappings and settings'},indent=2)+'\n')
 print(f'completed {a.reps} bootstrap replicates and matched-removal replicates per design')
if __name__=='__main__':main()
