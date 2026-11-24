"""Rebuild dataset-specific sample-group audit TSVs from frozen metadata."""
import csv,gzip
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
def read(p):return list(csv.DictReader(open(p),delimiter='\t'))
def write(p,rows):
 p.parent.mkdir(parents=True,exist_ok=True)
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
# Anopheles: population groups and homozygous background subsets.
rows=[r for r in read(ROOT/'data/anopheles_2la/metadata/sample_manifest.tsv') if r['species'] in ('gambiae','coluzzii')]; groups={}
for sp in ('gambiae','coluzzii'):
 sr=[r for r in rows if r['species']==sp];groups['species_pooled_'+sp]=sr
 for st,lab in [('A0_homozygous','standard'),('A1_homozygous','inverted')]:groups['species_'+lab+'_'+sp]=[r for r in sr if r['2La_state']==st]
 for loc in sorted(set(r['locality'] for r in sr)):
  lr=[r for r in sr if r['locality']==loc];groups['population_pooled_'+sp+'_'+loc]=lr
  for st,lab in [('A0_homozygous','standard'),('A1_homozygous','inverted')]:
   sub=[r for r in lr if r['2La_state']==st]
   if len(sub)>=2:groups['population_'+lab+'_'+sp+'_'+loc]=sub
out=[]
for g,rs in groups.items():
 if len(rs)>=2:
  for r in rs:out.append({'group':g,'n_individuals':len(rs),'sample_id':r['sample_id'],'species':r['species'],'population':r['population'],'locality':r['locality'],'2La_state':r['2La_state'],'sample_set':r['sample_set'],'reference_build':r['reference_build']})
write(ROOT/'empirical/anopheles_2la/analysis/tajima_d/results/anopheles_tajima_d_sample_groups.tsv',out)
# Fire ants: frozen Q1 membership intersected with VCF sample inventory.
mem=read(ROOT/'empirical/fire_ants_chr16/analysis/patterson_d/config/frozen_sample_membership.tsv');inv=read(ROOT/'empirical/fire_ants_chr16/analysis/patterson_d/results/vcf_sample_inventory.tsv');info={x['vcf_sample_id']:x for x in inv};out=[]
for m in mem:
 if m['quartet_id']=='Q1_primary_Argentina' and m['sample_id'] in info:
  x=info[m['sample_id']]
  if x['species']=='S.invicta' and x['caste']=='Male' and x['supergene_variant'] in ('SB','Sb'):
   out.append({'group':x['supergene_variant'],'sample_id':m['sample_id'],'species':x['species'],'country':x['country'],'caste':x['caste'],'fraction_missing':x['fraction_missing'],'sample_quality_flag':x.get('sample_quality_flag',''),'assignment_source':'frozen Q1 Argentina membership + VCF inventory'})
write(ROOT/'empirical/fire_ants_chr16/analysis/tajima_d/results/fire_ant_tajima_d_sample_groups.tsv',out)
# Cod populations and source arrangement labels; the VCF contains 24 diploid fish.
pop={'avc':'More_stationary','avo':'More_migratory','bat':'Labrador','bor':'Bornholm_Basin','icc':'Iceland_coastal','ico':'Iceland_oceanic','kie':'Kiel_Bight','lfc':'Lofoten_coastal','lfo':'Lofoten_oceanic','low':'Lowestoft','twc':'Trondheim_coastal','two':'Trondheim_oceanic'};arr={x['canonical_population']:x for x in read(ROOT/'data/atlantic_cod/processed/population_arrangements.tsv')}
vcf=ROOT/'data/atlantic_cod/raw/zenodo/local_downloads/gadus_morhua_supergene_lg01_relate.vcf.gz'
with gzip.open(vcf,'rt') as f:
 ids=[]
 for line in f:
  if line.startswith('#CHROM'):ids=line.rstrip().split('\t')[9:];break
out=[]
for code,label in pop.items():
 ss=[s for s in ids if s.startswith('Gadmor_'+code)]
 for sample in ss:
  a=arr.get(label,{});r={'population_code':code,'population':label,'sample_id':sample,'n_individuals':len(ss),'ploidy':'diploid','phase_present_in_source':True}
  for lg in ('LG01','LG02','LG07','LG12'):
   r[lg+'_arrangement_state']=a.get(lg+'_state','');r[lg+'_confidence']=a.get(lg+'_confidence','')
  out.append(r)
write(ROOT/'empirical/atlantic_cod/analysis/tajima_d/results/atlantic_cod_tajima_d_sample_groups.tsv',out)
print('Anopheles rows',len(rows),'groups',len(groups),'Fire ants',len([x for x in mem if x['quartet_id']=='Q1_primary_Argentina']),'Cod fish',len(ids))
