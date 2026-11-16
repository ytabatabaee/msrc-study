#!/usr/bin/env python3
"""Render six figure families in PNG, PDF and SVG from frozen Stage 2 outputs."""
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[5];BASE=ROOT/"empirical/fire_ants_chr16/analysis/patterson_d";R=BASE/"results";F=BASE/"figures"
F.mkdir(exist_ok=True)
def read(p):
    with open(p,newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def save(fig,name):
    fig.tight_layout()
    for ext in ("png","pdf","svg"):fig.savefig(F/f"{name}.{ext}",dpi=300,bbox_inches="tight")
    plt.close(fig)
def num(x):
    try:return float(x)
    except:return float("nan")
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":120})
windows=read(R/"window_patterson_d.tsv");tops=read(R/"induced_quartet_topologies.tsv");regions=read(R/"inside_outside_summary.tsv");corr=read(R/"spatial_correlations.tsv")
old={r["window_index"]:r for r in read(ROOT/"data/fire_ants_chr16/processed/stage4a_window_quartet_support.tsv")}
primary=("Q1_primary_Argentina","Q2_primary_Argentina")
colors={"Q1_primary_Argentina":"#2878B5","Q2_primary_Argentina":"#D95319"}

# A: pooled nucleotide D along original linkage-group-16 coordinates; exact tree window intervals.
fig,axs=plt.subplots(2,1,figsize=(11,6),sharex=False)
for ax,region in zip(axs,("chr16","chr1")):
  for q in primary:
    rr=sorted([r for r in windows if r["quartet_id"]==q and r["chrom"]==region],key=lambda r:num(r["mid"]))
    x=[num(r["mid"])/1e6 for r in rr];y=[num(r["D"]) for r in rr]
    ax.scatter(x,y,s=14,color=colors[q],alpha=.85,label=q.replace("_primary_Argentina",""))
  if region=="chr16":
    s=[r for r in windows if r["region"]=="chr16_supergene"]
    ax.axvspan(min(num(r["start"]) for r in s)/1e6,max(num(r["end"]) for r in s)/1e6,color="#f2c14e",alpha=.22,label="published supergene windows")
  ax.axhline(0,color="black",lw=.6);ax.set_ylabel("Patterson's D");ax.set_title(region);ax.legend(frameon=False,ncol=2)
axs[-1].set_xlabel("Author-window midpoint (Mb)")
save(fig,"figure_A_nucleotide_D_tracks")

# B: original TWISST supports and induced topology tracks alongside site D.
fig,axs=plt.subplots(3,1,figsize=(11,8),sharex=True)
for ax,field,label in zip(axs,("D","q_S","q_H","q_3"),("D_site","q_S (species split)","q_H (haplotype split)","q_3 (third split)")):
  pass
# Use four rows because D_tree_old is distinct and preserve original definitions.
plt.close(fig)
fig,axs=plt.subplots(4,1,figsize=(11,9),sharex=True)
for ax,key,label in zip(axs,("D","q_S","q_H","q_3"),("D_site","q_S","q_H","q_3")):
  if key=="D":
    for q in primary:
      rr=sorted([r for r in windows if r["quartet_id"]==q and r["chrom"]=="chr16"],key=lambda r:num(r["mid"]))
      ax.scatter([num(r["mid"])/1e6 for r in rr],[num(r["D"]) for r in rr],s=13,color=colors[q],label=q.split("_")[0])
  else:
    rr=sorted([r for r in windows if r["quartet_id"]==primary[0] and r["chrom"]=="chr16"],key=lambda r:num(r["mid"]))
    ax.scatter([num(r["mid"])/1e6 for r in rr],[num(old[r["window_index"]][key]) for r in rr],s=13,color={"q_S":"#1b9e77","q_H":"#d95f02","q_3":"#7570b3"}[key])
  ax.axhline(0,color="black",lw=.5);ax.set_ylabel(label)
  if key=="D":ax.legend(frameon=False,ncol=2)
s=[r for r in windows if r["region"]=="chr16_supergene"]
for ax in axs:ax.axvspan(min(num(r["start"]) for r in s)/1e6,max(num(r["end"]) for r in s)/1e6,color="#f2c14e",alpha=.16)
axs[-1].set_xlabel("Author-window midpoint (Mb)")
save(fig,"figure_B_D_and_TWISST_tracks")

# C: pooled inside/outside contrast across each frozen quartet, paired with flank components.
fig,ax=plt.subplots(figsize=(10,5));qids=sorted({r["quartet_id"] for r in regions});x=np.arange(len(qids));width=.36
for off,field,label,color in ((-.18,"D_inside","Supergene","#D95319"),(.18,"D_outside","chr16A + chr16B","#2878B5")):
 vals=[]
 for q in qids:
  rr=next((r for r in regions if r["quartet_id"]==q and r["outside_region"]=="chr16A+chr16B"),None);vals.append(num(rr[field]) if rr else np.nan)
 ax.bar(x+off,vals,width,label=label,color=color)
ax.set_xticks(x);ax.set_xticklabels([q.replace("_","\n") for q in qids],fontsize=7);ax.set_ylabel("Pooled Patterson's D");ax.axhline(0,color="black",lw=.6);ax.legend(frameon=False)
save(fig,"figure_C_inside_outside_quartets")

# D: exact induced D_topology against matched nucleotide D by local window.
top={(r["quartet_id"],r["window_index"]):r for r in tops}
fig,axs=plt.subplots(1,2,figsize=(10,4),sharex=True,sharey=True)
for ax,q in zip(axs,primary):
 for r in windows:
  if r["quartet_id"]!=q or r["chrom"]!="chr16" or r["D"]=="NA":continue
  t=top.get((q,r["window_index"]));
  if t:ax.scatter(num(r["D"]),num(t["D_topology"]),s=16,c="#D95319" if r["region"]=="chr16_supergene" else "#2878B5",alpha=.65)
 ax.axhline(0,color="black",lw=.5);ax.axvline(0,color="black",lw=.5);ax.set_title(q.split("_")[0]);ax.set_xlabel("D_site")
axs[0].set_ylabel("D_topology = q3 − q2")
save(fig,"figure_D_site_vs_exact_topology")

# E: contributions, informative SNP density and callability.
fig,axs=plt.subplots(3,1,figsize=(11,8),sharex=True)
q="Q1_primary_Argentina";rr=sorted([r for r in windows if r["quartet_id"]==q and r["chrom"]=="chr16"],key=lambda r:num(r["mid"]))
x=[num(r["mid"])/1e6 for r in rr]
axs[0].scatter(x,[num(r["abba"]) for r in rr],label="ABBA",color="#D95319",s=14);axs[0].scatter(x,[num(r["baba"]) for r in rr],label="BABA",color="#2878B5",s=14);axs[0].legend(frameon=False);axs[0].set_ylabel("Pooled contribution")
axs[1].scatter(x,[int(r["informative"]) for r in rr],color="#7A9E7E",s=14);axs[1].set_ylabel("Informative SNPs")
axs[2].scatter(x,[num(r.get("genotype_call_rate_at_variant_sites","nan")) for r in rr],color="#674EA7",s=14);axs[2].set_ylabel("Call rate");axs[2].set_xlabel("Author-window midpoint (Mb)")
s=[r for r in rr if r["region"]=="chr16_supergene"]
for ax in axs:ax.axvspan(min(num(r["start"]) for r in s)/1e6,max(num(r["end"]) for r in s)/1e6,color="#f2c14e",alpha=.15)
save(fig,"figure_E_contributions_density_callability")

# F: sensitivity of pooled inside-outside delta; colors distinguish designs.
fig,ax=plt.subplots(figsize=(10,6));rows=[r for r in regions if r["outside_region"]=="chr16A+chr16B"];rows.sort(key=lambda r:r["quartet_id"])
items=[(r["quartet_id"].replace("Q1_","Q1 ").replace("Q2_","Q2 ").replace("_"," "),num(r["Delta_D"]),"#2878B5") for r in rows]
for r in read(R/"transversions_inside_outside.tsv"):
  if r["quartet_id"] in primary:items.append((r["quartet_id"].split("_")[0]+" transversions",num(r["Delta_D_supergene_minus_chr16outside"]),"#7A9E7E"))
for q in primary:items.append((q.split("_")[0]+" strict DP15 / 80%",float("nan"),"#9E9E9E"))
labels=[x[0] for x in items];vals=[x[1] for x in items];cols=[x[2] for x in items];y=np.arange(len(items))
for yi,(v,c) in enumerate(zip(vals,cols)):
  if np.isfinite(v):ax.barh(yi,v,color=c)
  else:ax.text(.01,yi,"insufficient informative sites",va="center",ha="left",fontsize=7,color="#555555")
ax.set_yticks(y);ax.set_yticklabels(labels,fontsize=7);ax.axvline(0,color="black",lw=.6);ax.set_xlabel("D_supergene − D_chr16A+chr16B")
save(fig,"figure_F_quartet_and_outgroup_sensitivity")
print("wrote 18 figure files")
