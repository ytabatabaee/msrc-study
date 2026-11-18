"""Summarize group-mean patristic distances and four-point edge contrasts."""
import csv,io,math
from collections import defaultdict
from pathlib import Path
from Bio import Phylo
ROOT=Path(__file__).resolve().parents[6];B=ROOT/'empirical/fire_ants_chr16/analysis/patterson_d';R=B/'results';O=B/'final_validation/results'
def rd(p):
 with open(p,newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def main():
 groups=defaultdict(lambda:defaultdict(set))
 for r in rd(B/'config/frozen_sample_membership.tsv'):
  if r['quartet_id'] in ('Q1_primary_Argentina','Q2_primary_Argentina'):groups[r['quartet_id']][r['role']].add(r['sample_id'])
 out=[]
 for r in rd(ROOT/'data/fire_ants_chr16/processed/stage1_window_trees.tsv'):
  if r['chrom']!='chr16':continue
  tree=Phylo.read(io.StringIO(r['tree_newick']),'newick');tips={x.name for x in tree.get_terminals()};
  parent={};depth={tree.root:0.0};stack=[(tree.root,None)]
  while stack:
   node,par=stack.pop();parent[node]=par
   for child in node.clades:depth[child]=depth[node]+(child.branch_length or 0.0);stack.append((child,node))
  chains={};sets={}
  for tip in tree.get_terminals():
   chain=[];node=tip
   while node is not None:chain.append(node);node=parent[node]
   chains[tip.name]=chain;sets[tip.name]=set(chain)
  def dist(a,b):
   for node in chains[a]:
    if node in sets[b]:return depth[a_tip[a]]+depth[b_tip[b]]-2*depth[node]
   return float('nan')
  a_tip={x.name:x for x in tree.get_terminals()};b_tip=a_tip
  for q,g in groups.items():
   ids={k:sorted(v&tips) for k,v in g.items()};means={}
   for a,b in (('P1','P2'),('P1','P3'),('P1','O'),('P2','P3'),('P2','O'),('P3','O')):
    vals=[dist(x,y) for x in ids[a] for y in ids[b]]
    means[a+b]=sum(vals)/len(vals) if vals else float('nan')
   d12,d13,d14,d23,d24,d34=(means[x] for x in ('P1P2','P1P3','P1O','P2P3','P2O','P3O'))
   s1=d12+d34;s2=d13+d24;s3=d14+d23
   out.append({'quartet_id':q,'window_index':r['window_index'],'region':r['region'],'D12_P1P2':d12,'D13_P1P3':d13,'D14_P1O':d14,'D23_P2P3':d23,'D24_P2O':d24,'D34_P3O':d34,'four_point_edge_q1':(s2+s3-2*s1)/4,'four_point_edge_q2':(s1+s3-2*s2)/4,'four_point_edge_q3':(s1+s2-2*s3)/4,'edge_note':'contrasts of group-mean patristic distances; can be negative under population mixture and are not mean positive quartet branch lengths'})
 with open(O/'branch_length_sensitivity.tsv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=out[0].keys(),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(out)
if __name__=='__main__':main()
