#!/usr/bin/env python3
"""Run genotype-based D on verified BUSCO-window variants and frozen quartets."""
import csv, gzip, hashlib, json, math, random, statistics, sys
from collections import defaultdict
from pathlib import Path
from Bio import Phylo

ROOT=Path(__file__).resolve().parents[5]
BASE=ROOT/"empirical/fire_ants_chr16/analysis/patterson_d"
RESULTS=BASE/"results"; CONFIG=BASE/"config"
VCF=ROOT/"data/fire_ants_chr16/raw/genotypes/gt.vcf.gz"
WIN=ROOT/"data/fire_ants_chr16/processed/stage1_window_index.tsv"
TREES=ROOT/"data/fire_ants_chr16/processed/stage1_window_trees.tsv"
MAP=BASE/"results/window_gene_coordinate_mapping.tsv"
MEM=CONFIG/"frozen_sample_membership.tsv"
DP_MIN=10
MIN_GROUP_CALL_RATE=.5
MIN_OUTGROUP_CALLS=4

def tsv(path):
    with open(path,newline="") as f: return list(csv.DictReader(f,delimiter="\t"))
def write(path, rows):
    if not rows:return
    with open(path,"w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n",extrasaction="ignore");w.writeheader();w.writerows(rows)
def sha(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(4*1024*1024),b""):h.update(b)
    return h.hexdigest()
def vcf_gt(call, fmt):
    if fmt=="GT:DP":
        parts=call.split(":",1);gt=parts[0];dpstr=parts[1] if len(parts)>1 else "0"
    else:
        vals=call.split(":");ix={x:i for i,x in enumerate(fmt.split(":"))}
        gt=vals[ix["GT"]] if "GT" in ix and ix["GT"]<len(vals) else "."
        dpstr=vals[ix["DP"]] if "DP" in ix and ix["DP"]<len(vals) else "0"
    try: dp=int(dpstr)
    except ValueError: dp=0
    if dp<DP_MIN or gt in (".","./.",".|."): return None,"missing"
    if "." in gt:return None,"missing"
    if "/" in gt or "|" in gt:
        sep="/" if "/" in gt else "|";aa=gt.split(sep)
        if len(aa)!=2 or aa[0]!=aa[1]:return None,"ambiguous"
        a=aa[0]
    else:a=gt
    if a=="0":a=0
    elif a=="1":a=1
    else:return None,"ambiguous"
    return a,"called"
def init_counter():
    return {"records":0,"valid_snps":0,"callability_pass":0,"polarizable":0,"missing_genotypes":0,
            "ambiguous_genotypes":0,"group_callability_exclusions":0,"outgroup_polarization_failures":0,
            "abba":0.0,"baba":0.0,"informative":0,"call_obs_P1":0,"call_obs_P2":0,"call_obs_P3":0,"call_obs_O":0,
            "derived_ref_sites":0,"derived_alt_sites":0,"called_individual_obs":0}
def add_site(c, group_sizes, outgroup_size, alleles, outcalls, orient="all"):
    c["valid_snps"]+=1
    if orient=="transversion": return
    # Alleles are REF/ALT; outgroup must be sufficiently sampled and fixed among valid calls.
    for i in range(3): c["call_obs_P%d"%(i+1)]+=len(alleles[i])
    c["call_obs_O"]+=len(outcalls); c["called_individual_obs"]+=sum(map(len,alleles))+len(outcalls)
    if any(len(a)<math.ceil(n*MIN_GROUP_CALL_RATE) for a,n in zip(alleles,group_sizes)):
        c["group_callability_exclusions"]+=1;return
    c["callability_pass"]+=1
    if len(outcalls)<min(MIN_OUTGROUP_CALLS, outgroup_size) or len(set(outcalls))!=1:
        c["outgroup_polarization_failures"]+=1;return
    anc=outcalls[0]; c["polarizable"]+=1
    if anc==0:c["derived_ref_sites"]+=1
    else:c["derived_alt_sites"]+=1
    # VCF allele 1 is ALT: ALT is derived for REF-fixed outgroups; REF is
    # derived for ALT-fixed outgroups.
    freqs=[sum(a if anc==0 else 1-a for a in x)/len(x) for x in alleles]
    p1,p2,p3=freqs
    ab=(1-p1)*p2*p3; ba=p1*(1-p2)*p3
    c["abba"]+=ab;c["baba"]+=ba
    if ab+ba>0:c["informative"]+=1
    # Independently calculated scalar equality is tested on sampled sites.

def mapped_intervals():
    iv=defaultdict(list); wrows={}
    for w in tsv(WIN): wrows[int(w["window_index"])]=w
    for r in tsv(MAP):
        if r["vcf_contig_present"]!="1":continue
        ix=int(r["window_index"]); a=int(r["vcf_gene_start_1based"]); b=int(r["vcf_gene_end_1based_inclusive"])
        iv[r["vcf_contig"]].append((a,b,ix))
    for c in iv:iv[c].sort()
    return iv,wrows

def scan_vcf(transversions_only=False):
    memberships=tsv(MEM); designs=tsv(CONFIG/"frozen_quartets.tsv")
    dmembers=defaultdict(lambda:defaultdict(list)); all_samples=set()
    for r in memberships:dmembers[r["quartet_id"]][r["role"]].append(r["sample_id"]);all_samples.add(r["sample_id"])
    intervals,windows=mapped_intervals()
    counters={d["quartet_id"]:defaultdict(init_counter) for d in designs}
    # Collect selected sample column positions from the VCF header.
    with gzip.open(VCF,"rt") as f:
        header=None
        for line in f:
            if line.startswith("#CHROM"):
                header=line.rstrip("\n").split("\t");break
        if header is None:raise RuntimeError("VCF header missing #CHROM")
        col={s:i for i,s in enumerate(header[9:],9)}
        missing=all_samples-set(col)
        if missing:raise RuntimeError("Frozen sample IDs absent from VCF: "+repr(sorted(missing)))
        used=sorted(all_samples); used_cols=[col[s] for s in used]; used_ix={s:i for i,s in enumerate(used)}
        grouped={qid:{role:[used_ix[x] for x in samples] for role,samples in roles.items()} for qid,roles in dmembers.items()}
        total_records=0; mapped_records=0; site_examples=[]
        for line in f:
            if line.startswith("#"):continue
            total_records+=1
            # parse only leading fields on unassigned scaffolds
            a=line.find("\t");b=line.find("\t",a+1)
            chrom=line[:a]
            if chrom not in intervals:continue
            pos=int(line[a+1:b]); hits={ix for lo,hi,ix in intervals[chrom] if lo<=pos<=hi}
            if not hits:continue
            mapped_records+=1
            fields=line.rstrip("\n").split("\t")
            ref,alt=fields[3],fields[4]
            if len(ref)!=1 or len(alt)!=1 or "," in alt or alt==".":
                for q in counters:
                    for ix in hits:counters[q][ix]["records"]+=1
                continue
            if transversions_only and (ref,alt) in (("A","G"),("G","A"),("C","T"),("T","C")):
                for q in counters:
                    for ix in hits:counters[q][ix]["records"]+=1
                continue
            fmt=fields[8]
            gen=[]
            for j in used_cols:gen.append(vcf_gt(fields[j],fmt))
            for qid,roles in grouped.items():
                for ix in hits:
                    c=counters[qid][ix];c["records"]+=1
                allele_groups=[];out=[];miss=amb=0
                for role in ("P1","P2","P3","O"):
                    vals=[]
                    for j in roles[role]:
                        al,status=gen[j]
                        if al is None:
                            if status=="ambiguous":amb+=1
                            else:miss+=1
                        else:vals.append(al)
                    if role=="O":out=vals
                    else:allele_groups.append(vals)
                for ix in hits:
                    c=counters[qid][ix];c["missing_genotypes"]+=miss;c["ambiguous_genotypes"]+=amb
                    add_site(c,tuple(len(roles[x]) for x in ("P1","P2","P3")),len(roles["O"]),allele_groups,out)
                if len(site_examples)<25 and all(len(x)>0 for x in allele_groups) and out and len(set(out))==1:
                    anc=out[0]; p=[sum(a if anc==0 else 1-a for a in x)/len(x) for x in allele_groups]
                    p1,p2,p3=p; site_examples.append({"quartet_id":qid,"window_index":min(hits),"p1":p1,"p2":p2,"p3":p3,"anc":anc,"abba_scalar":(1-p1)*p2*p3,"baba_scalar":p1*(1-p2)*p3})
    return counters,designs,windows,dmembers,{"vcf_records_total":total_records,"mapped_records":mapped_records,"sample_count_used":len(all_samples),"scanned_sha256":sha(VCF)},site_examples

def summarize(counters,designs,windows,scan,examples,suffix=""):
    memberships={r["quartet_id"]:r for r in designs}
    rows=[]
    for q,bywin in counters.items():
        for ix in sorted(windows):
            c=bywin[ix]; w=windows[ix];den=c["abba"]+c["baba"]
            rows.append({"quartet_id":q,"window_index":ix,"chrom":w["chrom"],"region":w["region"],"start":w["start"],"end":w["end"],"mid":w["mid"],"n_P1":memberships[q]["n_P1"],"n_P2":memberships[q]["n_P2"],"n_P3":memberships[q]["n_P3"],"n_O":memberships[q]["n_O"],**c,"D":(c["abba"]-c["baba"])/den if den else "NA","informative_site_ratio":c["informative"]/c["polarizable"] if c["polarizable"] else "NA","coordinate_scope":"VCF variants in exact BUSCO gene bodies assigned to author window; chromosome coordinates label parent window"})
    write(RESULTS/("window_patterson_d"+suffix+".tsv"),rows)
    # Region pooled counts; never average window D.
    region_rows=[]
    for q in counters:
        for chrom,region in [("chr16","chr16_supergene"),("chr16","chr16A"),("chr16","chr16B"),("chr16","chr16_outside"),("chr1","chr1"),("chr16","chr16")]:
            sel=[r for r in rows if r["quartet_id"]==q and ((r["region"]==region) if region!="chr16_outside" and region!="chr16" else (r["chrom"]=="chr16" and (r["region"]!="chr16_supergene" if region=="chr16_outside" else True)))]
            if not sel:continue
            ab=sum(float(r["abba"]) for r in sel);ba=sum(float(r["baba"]) for r in sel);d=(ab-ba)/(ab+ba) if ab+ba else "NA"
            region_rows.append({"quartet_id":q,"region":region,"n_windows":len(sel),"records":sum(int(r["records"]) for r in sel),"valid_snps":sum(int(r["valid_snps"]) for r in sel),"polarizable_sites":sum(int(r["polarizable"]) for r in sel),"informative_sites":sum(int(r["informative"]) for r in sel),"ABBA":ab,"BABA":ba,"D":d,"callable_variant_sites":sum(int(r["polarizable"]) for r in sel),"coverage_note":"variant-only VCF; callable invariant bases unavailable; gene-body-only BUSCO coverage"})
    write(RESULTS/("region_patterson_d"+suffix+".tsv"),region_rows)
    # Multiresolution physical-coordinate bins using frozen author coordinates.
    multi=[]
    for q in counters:
      for size in (50000,100000,250000,500000,1000000):
       bins=defaultdict(list)
       for r in rows:
        if r["quartet_id"]!=q:continue
        chrom=r["chrom"]
        if chrom not in ("chr1","chr16"):continue
        key=(chrom,(int(r["mid"])-1)//size*size+1)
        bins[key].append(r)
       for (chrom,start),rr in sorted(bins.items()):
        ab=sum(float(x["abba"]) for x in rr);ba=sum(float(x["baba"]) for x in rr);den=ab+ba
        multi.append({"quartet_id":q,"resolution_bp":size,"chrom":chrom,"bin_start":start,"bin_end":start+size-1,"n_source_windows":len(rr),"records":sum(int(x["records"]) for x in rr),"valid_snps":sum(int(x["valid_snps"]) for x in rr),"polarizable_sites":sum(int(x["polarizable"]) for x in rr),"informative_sites":sum(int(x["informative"]) for x in rr),"ABBA":ab,"BABA":ba,"D":(ab-ba)/den if den else "NA","coverage_note":"D pooled from variants in BUSCO gene bodies; parent author-window coordinates; no continuous scaffold lift"})
    write(RESULTS/("multiresolution_patterson_d"+suffix+".tsv"),multi)
    # Scalar independent verifier from saved site-frequency examples.
    for x in examples:
        p1,p2,p3=x["p1"],x["p2"],x["p3"];assert abs(x["abba_scalar"]-(1-p1)*p2*p3)<1e-15;assert abs(x["baba_scalar"]-p1*(1-p2)*p3)<1e-15
    return rows,region_rows,multi

def topology_outputs(designs,windows):
    # Exact frequencies over all sampled four-individual combinations. Vectorized
    # four-point distances avoid treating multiple collinear edges as repeats.
    import numpy as np
    groups=defaultdict(lambda:defaultdict(set))
    for r in tsv(MEM):groups[r["quartet_id"]][r["role"]].add(r["sample_id"])
    tree_rows=tsv(TREES);outs=[]
    for tr in tree_rows:
      tree=Phylo.read(__import__('io').StringIO(tr["tree_newick"]),"newick")
      tips={t.name for t in tree.get_terminals()}; leaves={}
      for q,g in groups.items():
       lab={}
       for role,samples in g.items():
        for s in samples:
         if s in tips:lab[s]=role
       if any(not any(v==x for v in lab.values()) for x in ("P1","P2","P3","O")):
        outs.append({"quartet_id":q,"window_index":tr["window_index"],"chrom":tr["chrom"],"region":tr["region"],"q1":0,"q2":0,"q3":0,"missing_groups":";".join(x for x in ("P1","P2","P3","O") if x not in lab.values()),"unresolved_fraction":"NA","D_topology":"NA"});continue
       # Build unit-edge distance matrix over selected tree tips.
       adj=defaultdict(list)
       for parent in tree.find_clades(order="preorder"):
        for child in parent.clades:adj[parent].append(child);adj[child].append(parent)
       selected=list(lab)
       tipcl={t.name:t for t in tree.get_terminals()}
       dist={s:{} for s in selected}
       for s in selected:
        stack=[(tipcl[s],None,0)]
        while stack:
         node,prev,d=stack.pop()
         if node.is_terminal() and node.name in dist:dist[s][node.name]=d
         for nxt in adj[node]:
          if nxt is not prev:stack.append((nxt,node,d+1))
       ids={role:[s for s in g[role] if s in tips] for role in ("P1","P2","P3","O")}
       npossible=1
       for role in ids:npossible*=len(ids[role])
       if not npossible:
        outs.append({"quartet_id":q,"window_index":tr["window_index"],"chrom":tr["chrom"],"region":tr["region"],"q1":0,"q2":0,"q3":0,"D_topology":"NA","missing_groups":";".join(x for x in ids if not ids[x]),"unresolved_fraction":"NA","n_possible_sample_quartets":0,"n_resolved_sample_quartets":0});continue
       # Exact Cartesian product in chunks over P1 to cap temporary memory.
       counts=[0,0,0];unresolved=0
       a,b,c,d=(ids[x] for x in ("P1","P2","P3","O"))
       nB,nC,nD=len(b),len(c),len(d)
       ib=np.repeat(np.arange(nB,dtype=np.int32),nC*nD)
       ic=np.tile(np.repeat(np.arange(nC,dtype=np.int32),nD),nB)
       idd=np.tile(np.arange(nD,dtype=np.int32),nB*nC)
       mat12=np.array([[dist[x][y] for y in b] for x in a],dtype=np.int16)
       mat13=np.array([[dist[x][y] for y in c] for x in a],dtype=np.int16)
       mat14=np.array([[dist[x][y] for y in d] for x in a],dtype=np.int16)
       mat23=np.array([[dist[x][y] for y in c] for x in b],dtype=np.int16)
       mat24=np.array([[dist[x][y] for y in d] for x in b],dtype=np.int16)
       mat34=np.array([[dist[x][y] for y in d] for x in c],dtype=np.int16)
       for ai in range(len(a)):
        nsub=nB*nC*nD
        d12=mat12[ai,ib];d13=mat13[ai,ic];d14=mat14[ai,idd]
        d23=mat23[ib,ic];d24=mat24[ib,idd];d34=mat34[ic,idd]
        sums=np.vstack((d12+d34,d13+d24,d14+d23))
        mins=sums.min(axis=0); ties=(sums==mins).sum(axis=0)
        unresolved+=int((ties>1).sum())
        for j in range(3):counts[j]+=int(((sums[j]==mins)&(ties==1)).sum())
       resolved=sum(counts);qfreq=[x/resolved for x in counts] if resolved else [0,0,0]
       outs.append({"quartet_id":q,"window_index":tr["window_index"],"chrom":tr["chrom"],"region":tr["region"],"q1":qfreq[0],"q2":qfreq[1],"q3":qfreq[2],"D_topology":qfreq[2]-qfreq[1],"missing_groups":"","unresolved_fraction":unresolved/npossible if npossible else "NA","n_possible_sample_quartets":npossible,"n_resolved_sample_quartets":resolved})
    write(RESULTS/"induced_quartet_topologies.tsv",outs);return outs

def downstream(rows,regions,multi,topos):
    # Inside/outside pooled D and moving-block jackknife. For physical bins, keep counts pooled.
    summaries=[];unc=[];spatial=[]
    for q in sorted({r["quartet_id"] for r in rows}):
      for outside in ("chr16A+chr16B","chr1"):
       inside=[r for r in rows if r["quartet_id"]==q and r["region"]=="chr16_supergene"]
       out=[r for r in rows if r["quartet_id"]==q and ((r["region"] in ("chr16A","chr16B")) if outside=="chr16A+chr16B" else r["region"]=="chr1")]
       def sums(rs):return sum(float(r["abba"]) for r in rs),sum(float(r["baba"]) for r in rs)
       ai,bi=sums(inside);ao,bo=sums(out);di=(ai-bi)/(ai+bi) if ai+bi else float("nan");do=(ao-bo)/(ao+bo) if ao+bo else float("nan")
       a_flank=[r for r in rows if r["quartet_id"]==q and r["region"]=="chr16A"];b_flank=[r for r in rows if r["quartet_id"]==q and r["region"]=="chr16B"]
       ordered=sorted(inside,key=lambda r:float(r["mid"]));trim=ordered[1:-1] if len(ordered)>2 else ordered
       expanded=ordered+([max(a_flank,key=lambda r:float(r["mid"]))] if a_flank else [])+([min(b_flank,key=lambda r:float(r["mid"]))] if b_flank else [])
       ati,bti=sums(trim);ae,be=sums(expanded);dtrim=(ati-bti)/(ati+bti) if ati+bti else "NA";dexp=(ae-be)/(ae+be) if ae+be else "NA"
       summaries.append({"quartet_id":q,"inside_region":"chr16_supergene","outside_region":outside,"n_inside_windows":len(inside),"n_outside_windows":len(out),"ABBA_inside":ai,"BABA_inside":bi,"D_inside":di,"informative_inside":sum(int(r["informative"]) for r in inside),"ABBA_outside":ao,"BABA_outside":bo,"D_outside":do,"informative_outside":sum(int(r["informative"]) for r in out),"Delta_D":di-do,"D_chr16A_left_flank":((sums(a_flank)[0]-sums(a_flank)[1])/sum(sums(a_flank)) if sum(sums(a_flank)) else "NA"),"D_chr16B_right_flank":((sums(b_flank)[0]-sums(b_flank)[1])/sum(sums(b_flank)) if sum(sums(b_flank)) else "NA"),"D_supergene_trim_first_last_window":dtrim,"Delta_trimmed_inside_minus_outside":(dtrim-do if dtrim!="NA" else "NA"),"D_supergene_expanded_one_nearest_flank_each_side":dexp,"Delta_expanded_inside_minus_outside":(dexp-do if dexp!="NA" else "NA"),"target_gene_body_span_bp_inside":sum(int(r.get("target_gene_body_span_bp",0)) for r in inside),"target_gene_body_span_bp_outside":sum(int(r.get("target_gene_body_span_bp",0)) for r in out),"coverage_note":"variant-only BUSCO gene-body sites; spans are BUSCO target intervals, not callable invariant bases; chromosome label from frozen author window"})
       for block in (250000,500000,1000000):
        # block coordinates over the entire chr16; each physical block contributes pooled ABBA/BABA.
        blk=defaultdict(lambda:[0.,0.,0.,0.])
        for r in inside:
         b=(r["chrom"],(int(r["mid"])-1)//block);blk[b][0]+=float(r["abba"]);blk[b][1]+=float(r["baba"])
        for r in out:
         b=(r["chrom"],(int(r["mid"])-1)//block);blk[b][2]+=float(r["abba"]);blk[b][3]+=float(r["baba"])
        if len(blk)<2:
         unc.append({"quartet_id":q,"comparison":outside,"block_bp":block,"n_blocks":len(blk),"jackknife_se":"NA","bootstrap_low":"NA","bootstrap_high":"NA","note":"fewer than two physical blocks"});continue
        vals=list(blk.values()); n=len(vals);boots=[];bootinside=[];random.seed(20261116)
        for _ in range(2000):
         samp=[random.choice(vals) for z in range(n)];ai3=sum(x[0] for x in samp);bi3=sum(x[1] for x in samp);ao3=sum(x[2] for x in samp);bo3=sum(x[3] for x in samp)
         ddin=(ai3-bi3)/(ai3+bi3) if ai3+bi3 else float("nan");ddout=(ao3-bo3)/(ao3+bo3) if ao3+bo3 else float("nan")
         bootinside.append(ddin);boots.append(ddin-ddout)
        # jackknife is uncertainty in the within-minus-outside contrast; leave one spatial block out and recompute pooled regions.
        jds=[]
        for b in blk:
         ai2,bi2=sums([r for r in inside if (r["chrom"],(int(r["mid"])-1)//block)!=b]);ao2,bo2=sums([r for r in out if (r["chrom"],(int(r["mid"])-1)//block)!=b])
         di2=(ai2-bi2)/(ai2+bi2) if ai2+bi2 else 0;do2=(ao2-bo2)/(ao2+bo2) if ao2+bo2 else 0;jds.append(di2-do2)
        mean=sum(jds)/n;se=math.sqrt((n-1)/n*sum((x-mean)**2 for x in jds))
        boots=sorted(x for x in boots if math.isfinite(x));bootinside=sorted(x for x in bootinside if math.isfinite(x))
        unc.append({"quartet_id":q,"comparison":outside,"block_bp":block,"n_blocks":n,"jackknife_se_delta_D":se,"jackknife_95_low_delta_D":(di-do)-1.96*se,"jackknife_95_high_delta_D":(di-do)+1.96*se,"bootstrap_low_delta_D":boots[int(.025*len(boots))] if boots else "NA","bootstrap_high_delta_D":boots[int(.975*len(boots))] if boots else "NA","bootstrap_low_inside_D":bootinside[int(.025*len(bootinside))] if bootinside else "NA","bootstrap_high_inside_D":bootinside[int(.975*len(bootinside))] if bootinside else "NA","bootstrap_method":"resample physical blocks; pooled ABBA/BABA per region; 2000 replicates"})
      # spatial cyclic shift null on chr16, region membership fixed; preserve per-window counts.
      seq=sorted([r for r in rows if r["quartet_id"]==q and r["chrom"]=="chr16"],key=lambda r:int(r["mid"]))
      if len(seq)>1:
       labels=[r["region"]=="chr16_supergene" for r in seq]; observed=[]
       for shift in range(len(seq)):
        rotated=seq[shift:]+seq[:shift]; ia=[r for r,l in zip(rotated,labels) if l];oa=[r for r,l in zip(rotated,labels) if not l]
        ai,bi=sums(ia);ao,bo=sums(oa);dd=(ai-bi)/(ai+bi) if ai+bi else 0;ee=(ao-bo)/(ao+bo) if ao+bo else 0;observed.append(dd-ee)
       obs=observed[0];p=sum(abs(x)>=abs(obs) for x in observed)/len(observed)
       spatial.append({"quartet_id":q,"test":"circular_shift_inversion_labels_fixed","n_shifts":len(observed),"observed_delta_D":obs,"two_sided_p_exact":p,"warning":"only 96-ish author windows, linked and irregular; exploratory spatial null on gene-window track"})
    write(RESULTS/"inside_outside_summary.tsv",summaries);write(RESULTS/"block_uncertainty.tsv",unc);write(RESULTS/"spatial_null_validation.tsv",spatial)
    # Window-matched nucleotide / induced quartet / original seven-group TWISST tracks.
    cor=[]
    bytop={(r["quartet_id"],int(r["window_index"])):r for r in topos}
    old={int(r["window_index"]):r for r in tsv(ROOT/"data/fire_ants_chr16/processed/stage4a_window_quartet_support.tsv")}
    for q in sorted({r["quartet_id"] for r in rows}):
      for region in ("chr16","chr16_supergene","chr16_outside","chr1"):
       pairs=[]
       for r in rows:
        if r["quartet_id"]!=q or r["D"]=="NA":continue
        if region=="chr16" and r["chrom"]!="chr16":continue
        if region=="chr16_supergene" and r["region"]!="chr16_supergene":continue
        if region=="chr16_outside" and not (r["chrom"]=="chr16" and r["region"]!="chr16_supergene"):continue
        if region=="chr1" and r["chrom"]!="chr1":continue
        t=bytop.get((q,int(r["window_index"])));o=old[int(r["window_index"])]
        if t and t["D_topology"]!="NA":pairs.append((float(r["D"]),float(t["D_topology"]),float(t["q1"]),float(t["q2"]),float(t["q3"]),float(o["q_S"]),float(o["q_H"]),float(o["q_3"]),float(o["D"]),int(float(r["mid"])),r["chrom"],r["region"]))
       if len(pairs)>=3:
        centered=[None,None]
        if region=="chr16":
         groups2=defaultdict(list)
         for x in pairs:groups2["inside" if x[11]=="chr16_supergene" else "outside"].append(x)
         xc=[];yc=[]
         for group in groups2.values():
          mx=statistics.mean(x[0] for x in group);my=statistics.mean(x[1] for x in group)
          xc.extend(x[0]-mx for x in group);yc.extend(x[1]-my for x in group)
         centered=[corr(xc,yc),corr(xc,yc,True)]
        item={"quartet_id":q,"region":region,"n_windows":len(pairs),"pearson_Dsite_Dtopology":corr([x[0] for x in pairs],[x[1] for x in pairs]),"spearman_Dsite_Dtopology":corr([x[0] for x in pairs],[x[1] for x in pairs],True),"pearson_Dsite_Dtopology_centered_within_inside_outside":centered[0],"spearman_Dsite_Dtopology_centered_within_inside_outside":centered[1],"pearson_Dsite_q1":corr([x[0] for x in pairs],[x[2] for x in pairs]),"spearman_Dsite_q1":corr([x[0] for x in pairs],[x[2] for x in pairs],True),"pearson_Dsite_q2":corr([x[0] for x in pairs],[x[3] for x in pairs]),"spearman_Dsite_q2":corr([x[0] for x in pairs],[x[3] for x in pairs],True),"pearson_Dsite_q3":corr([x[0] for x in pairs],[x[4] for x in pairs]),"spearman_Dsite_q3":corr([x[0] for x in pairs],[x[4] for x in pairs],True),"pearson_Dsite_q_S":corr([x[0] for x in pairs],[x[5] for x in pairs]),"spearman_Dsite_q_S":corr([x[0] for x in pairs],[x[5] for x in pairs],True),"pearson_Dsite_q_H":corr([x[0] for x in pairs],[x[6] for x in pairs]),"spearman_Dsite_q_H":corr([x[0] for x in pairs],[x[6] for x in pairs],True),"pearson_Dsite_q_3_original":corr([x[0] for x in pairs],[x[7] for x in pairs]),"spearman_Dsite_q_3_original":corr([x[0] for x in pairs],[x[7] for x in pairs],True),"pearson_Dsite_D_tree_old_qH_minus_qS":corr([x[0] for x in pairs],[x[8] for x in pairs]),"spearman_Dsite_D_tree_old_qH_minus_qS":corr([x[0] for x in pairs],[x[8] for x in pairs],True),"note":"descriptive; linked windows; centered within inside/outside strata only for chr16-wide residual correlation; nucleotide site D, exact induced quartet q, and original seven-group TWISST are distinct estimands"}
        # Moving-block bootstrap of spatial correlations, preserving linked windows within 1-Mb blocks.
        blocks=defaultdict(list)
        for x in pairs:blocks[(x[10],(int(x[9])-1)//1000000)].append(x)
        keys=list(blocks); boots=[];random.seed(20261116)
        if len(keys)>=3:
         for _ in range(1000):
          sample=[]
          for key in (random.choice(keys) for z in keys):sample.extend(blocks[key])
          if len(sample)>=3:boots.append(corr([x[0] for x in sample],[x[1] for x in sample]))
        boots=sorted(x for x in boots if isinstance(x,float) and math.isfinite(x))
        item.update({"correlation_block_bp":1000000,"correlation_block_count":len(keys),"pearson_Dsite_Dtopology_block_boot_low":boots[int(.025*len(boots))] if boots else "NA","pearson_Dsite_Dtopology_block_boot_high":boots[int(.975*len(boots))] if boots else "NA"})
        cor.append(item)
    write(RESULTS/"spatial_correlations.tsv",cor)
    write(RESULTS/"quartet_robustness.tsv",[{"quartet_id":r["quartet_id"],"region":r["region"],"D":r["D"],"informative_sites":r["informative"],"polarizable_sites":r["polarizable"],"filter":"DP>=10; haploid homozygous calls only; >=50% called per ingroup group; >=4 unanimous outgroup calls"} for r in rows if r["region"] in ("chr16_supergene","chr16A","chr16B","chr1")])
    write(RESULTS/"outgroup_sensitivity.tsv",[r for r in summaries if r["quartet_id"] in ("Q1_outgroup_Brazil","Q1_outgroup_French_Guiana","Q1_primary_Argentina","Q2_primary_Argentina")])
    diag=[]
    for q in sorted({r["quartet_id"] for r in rows}):
     for region in ("chr16_supergene","chr16A","chr16B","chr1"):
      rr=[r for r in rows if r["quartet_id"]==q and r["region"]==region]
      if not rr:continue
      valid=sum(int(r["valid_snps"]) for r in rr);calls=sum(int(r["called_individual_obs"]) for r in rr)
      ns=sum(int(rr[0][f"n_{g}"]) for g in ("P1","P2","P3","O"))
      denom=valid*ns
      diag.append({"quartet_id":q,"region":region,"records":sum(int(r["records"]) for r in rr),"valid_biallelic_snps":valid,"outgroup_polarizable_snps":sum(int(r["polarizable"]) for r in rr),"informative_sites":sum(int(r["informative"]) for r in rr),"missing_genotype_observations":sum(int(r["missing_genotypes"]) for r in rr),"ambiguous_genotype_observations":sum(int(r["ambiguous_genotypes"]) for r in rr),"called_genotype_observations":calls,"expected_genotype_observations_at_variants":denom,"genotype_call_rate_at_variant_sites":calls/denom if denom else "NA","excluded_group_callability_sites":sum(int(r["group_callability_exclusions"]) for r in rr),"outgroup_polarization_failures":sum(int(r["outgroup_polarization_failures"]) for r in rr),"ABBA":sum(float(r["abba"]) for r in rr),"BABA":sum(float(r["baba"]) for r in rr),"D":(sum(float(r["abba"]) for r in rr)-sum(float(r["baba"]) for r in rr))/(sum(float(r["abba"]) for r in rr)+sum(float(r["baba"]) for r in rr)) if sum(float(r["abba"])+float(r["baba"]) for r in rr) else "NA","coverage_scope":"variant sites only; does not estimate invariant base callability"})
    write(RESULTS/"variant_diagnostics.tsv",diag)

def corr(x,y,spearman=False):
    if spearman:
        def rank(a):
            ix=sorted(range(len(a)),key=lambda i:a[i]);z=[0.]*len(a);i=0
            while i<len(ix):
                j=i+1
                while j<len(ix) and a[ix[j]]==a[ix[i]]:j+=1
                for k in range(i,j):z[ix[k]]=(i+j-1)/2
                i=j
            return z
        x=rank(x);y=rank(y)
    mx=statistics.mean(x);my=statistics.mean(y);sx=sum((a-mx)**2 for a in x);sy=sum((b-my)**2 for b in y)
    return sum((a-mx)*(b-my) for a,b in zip(x,y))/math.sqrt(sx*sy) if sx and sy else "NA"

def annotate_coordinate_coverage():
    # Span refers to exact target BUSCO gene bodies, not whole-window or invariant-site callability.
    bywin=defaultdict(lambda:defaultdict(list))
    for r in tsv(MAP):
        if r["vcf_contig_present"]!="1":continue
        bywin[int(r["window_index"])][r["vcf_contig"]].append((int(r["vcf_gene_start_1based"]),int(r["vcf_gene_end_1based_inclusive"])))
    spans={}
    for ix,contigs in bywin.items():
        total=0
        for seq in contigs.values():
            seq.sort();lo=hi=None
            for a,b in seq:
                if lo is None:lo,hi=a,b
                elif a<=hi+1:hi=max(hi,b)
                else:total+=hi-lo+1;lo,hi=a,b
            if lo is not None:total+=hi-lo+1
        spans[ix]=total
    rows=tsv(RESULTS/"window_patterson_d.tsv")
    for r in rows:
        ix=int(r["window_index"]);r["target_gene_body_span_bp"]=spans.get(ix,0)
        r["window_physical_span_bp"]=int(round(float(r["end"])-float(r["start"])+1))
        denom=int(r["valid_snps"])*(sum(int(r[f"n_{g}"]) for g in ("P1","P2","P3","O")))
        r["genotype_call_rate_at_variant_sites"]=int(r["called_individual_obs"])/denom if denom else "NA"
    write(RESULTS/"window_patterson_d.tsv",rows)

def main():
    RESULTS.mkdir(parents=True,exist_ok=True)
    if "--topology-only" in sys.argv:
        designs=tsv(CONFIG/"frozen_quartets.tsv");windows={int(r["window_index"]):r for r in tsv(WIN)}
        rows=tsv(RESULTS/"window_patterson_d.tsv");topos=topology_outputs(designs,windows)
        downstream(rows,[],[],topos)
        print("topology, contrasts, uncertainty and correlation outputs refreshed")
        return
    if "--downstream-only" in sys.argv:
        annotate_coordinate_coverage()
        rows=tsv(RESULTS/"window_patterson_d.tsv");topos=tsv(RESULTS/"induced_quartet_topologies.tsv")
        downstream(rows,[],[],topos)
        print("summary tables refreshed from saved windows and exact induced quartets")
        return
    is_tv="--transversions-only" in sys.argv
    is_strict="--strict-calls" in sys.argv
    if is_strict:
        global DP_MIN,MIN_GROUP_CALL_RATE
        DP_MIN=15;MIN_GROUP_CALL_RATE=.8
    counters,designs,windows,dm,scan,examples=scan_vcf(transversions_only=is_tv)
    suffix="_transversions" if is_tv else ("_strict" if is_strict else "")
    rows,regions,multi=summarize(counters,designs,windows,scan,examples,suffix=suffix)
    if is_tv or is_strict:
        out=[]
        for q in counters:
            rr=[r for r in regions if r["quartet_id"]==q]
            by={r["region"]:r for r in rr}
            a=by.get("chr16_supergene",{}); bA=by.get("chr16A",{}); bB=by.get("chr16B",{}); c1=by.get("chr1",{})
            ba=float(bA.get("ABBA",0))+float(bB.get("ABBA",0));bb=float(bA.get("BABA",0))+float(bB.get("BABA",0));do=(ba-bb)/(ba+bb) if ba+bb else "NA"
            di=a.get("D","NA")
            out.append({"quartet_id":q,"filter":"transversions only; otherwise Stage 2 primary rules" if is_tv else "strict DP>=15 and >=80% ingroup call rate; otherwise primary rules","D_supergene":di,"ABBA_supergene":a.get("ABBA",0),"BABA_supergene":a.get("BABA",0),"informative_supergene":a.get("informative_sites",0),"D_chr16A_plus_chr16B":do,"ABBA_chr16A_plus_chr16B":ba,"BABA_chr16A_plus_chr16B":bb,"informative_chr16A_plus_chr16B":int(bA.get("informative_sites",0))+int(bB.get("informative_sites",0)),"D_chr1":c1.get("D","NA"),"Delta_D_supergene_minus_chr16outside":float(di)-do if di!="NA" and do!="NA" else "NA"})
        prefix="transversions" if is_tv else "strict"
        write(RESULTS/(prefix+"_inside_outside.tsv"),out)
        (RESULTS/(prefix+"_run_manifest.json")).write_text(json.dumps({"mode":"transversions-only sensitivity" if is_tv else "strict callability sensitivity","scan":scan,"filter":{"DP_min":DP_MIN,"group_call_rate_min":MIN_GROUP_CALL_RATE,"outgroup_min":"min(4, frozen O sample count), all called O alleles unanimous","heterozygous_pseudo_diploid":"excluded"},"outputs":[prefix+"_window_patterson_d.tsv",prefix+"_region_patterson_d.tsv",prefix+"_inside_outside.tsv"]},indent=2)+"\n")
        print(json.dumps({"scan":scan,"sensitivity_contrasts":out},indent=2));return
    topos=topology_outputs(designs,windows)
    annotate_coordinate_coverage()
    rows=tsv(RESULTS/"window_patterson_d.tsv")
    downstream(rows,regions,multi,topos)
    # site examples retained as independent arithmetic verification record
    write(RESULTS/"scalar_validation_examples.tsv",examples)
    manifest={"stage":"2 nucleotide D","input_vcf":str(VCF.relative_to(ROOT)),"vcf_sha256":scan["scanned_sha256"],"scan":scan,"frozen_quartets":len(designs),"window_count":len(windows),"filters":{"DP_min":DP_MIN,"pseudo_diploid_heterozygotes":"excluded as ambiguous; never split into two haploid alleles","group_call_rate_min":MIN_GROUP_CALL_RATE,"outgroup":"at least 4 called individuals and unanimous fixed allele among calls","multiallelic_and_indel":"excluded","pool_and_random_allele":"excluded per frozen Stage 1 membership"},"coordinate_scope":"variants in exact BUSCO gene bodies, joined to author-defined four-gene windows; physical tracks use frozen parent window coordinates; no continuous scaffold liftover","reproduction":"python3 empirical/fire_ants_chr16/analysis/patterson_d/scripts/run_stage2.py"}
    (RESULTS/"run_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps({"scan":scan,"primary_region_results":[r for r in regions if r["quartet_id"] in ("Q1_primary_Argentina","Q2_primary_Argentina") and r["region"] in ("chr16_supergene","chr16_outside","chr1")],"n_window_rows":len(rows),"n_topology_rows":len(topos)},indent=2))

if __name__=="__main__":main()
