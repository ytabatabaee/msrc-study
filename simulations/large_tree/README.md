# MV root simulated data
---

Date belonging to the following paper:

Mai, Uyen, Erfan Sayyari, and Siavash Mirarab. “Minimum Variance Rooting of Phylogenetic Trees and Implications for Species Tree Reconstruction.” Edited by Gabriel Moreno-Hagelsieb. PLOS ONE 12, no. 8 (August 11, 2017): e0182238. https://doi.org/10.1371/journal.pone.0182238.

## Description of the data and file structure


### D1 dataset

D1—30-taxon heterogeneous dataset: Here, the number of ingroup species was fixed to 30. We simulated 100 replicates, each with a different species tree and 500 gene trees. This dataset is used for extensive analysis of all methods.

Each of D1 files includes the results for all replciates of all of the following models, given as subdirectories:
* outgroup.0.25.species.1.5.genes.1.5
* outgroup.0.5.species.1.5.genes.1.5
* outgroup.0.species.0.15.genes.0.15
* outgroup.0.species.1.5.genes.1.5
* outgroup.0.species.5.genes.5
* outgroup.1.species.0.15.genes.0.15
* outgroup.1.species.1.5.genes.1.5
* outgroup.1.species.5.genes.5
* outgroup.2.species.1.5.genes.1.5
* outgroup.4.species.1.5.genes.1.5

* D1-true.tar.gz: True gene trees and species trees. 
  * For each model and replicate, we have:
	* `control.txt`: Indelible control file
	* `s_tree.trees`: True species tree file
	* `truegenetrees`: True gene trees
	* `s_tree.trees.score`: The quartet score of the true species tree with respect to the true gene trees
  * In addition, we include the simphy scripts to do the simulations: `post_stidsim.pl`, `run-analysis.sh`, `simulate-alignments.sh`, `generateData.sh`
  * `LOG.txt`, `trees.txt`: the log files of simphy 
  * `mean-seq-length.txt`: mean sequence length per model condition and replicate
* D1-estimated.tar.gz: The input gene tree files (true and estimated)
  For each model and replicate, we have:
	* truegenetrees: true gene trees
	* estimatedgenetre.gtr: FastTree gene tree using GTR
	* estimatedgenetre.jc: FastTree gene tree using JC
	* estimatedgenetrees.tar.gz: estimated gene trees (fasttree)
	* estimatedgenetre.jc.info: FastTree log file for JC gene tree estimation
	* estimatedgenetre.gtr.info: FastTree lof file for GTR gene tree estimation
	* score-estimategenetre-gtr.sc: The accuracy of GTR gene trees
	* score-estimategenetre-jc.sc: The accuracy of JC gene trees
	* all-genes.phylip: The alignment file for all genes
* D1-estimated-rooted.tar.gz
	* results of methods MP, MV, Outgroup rooting (OG), and best posibile rooting by triplet criterion (`estimatedgenetre.gtr_MP`, `estimatedgenetre.gtr_MV`, `estimatedgenetre.gtr_OG`, `estimatedgenetre.gtr_ideal`). 
* D1-true-rooted.tar.gz: Similar to `D1-estimated-rooted.tar.gz` but with true gene trees used as input

### D2 dataset

D2—Large heterogeneous dataset: This dataset includes two subsets, one with 2000 and another with 5000 taxa, and is used for testing performance on large datasets. 
For both datasets we created 20 replicates with different species trees and 50 gene trees.

Each of D2 files includes the results for all replciates of all of the following models, given as subdirectories:
* sp-2000/outgroup.0.25.species.1.5.genes.1.5
* sp-2000/outgroup.0.5.species.1.5.genes.1.5
* sp-2000/outgroup.0.species.0.15.genes.0.15
* sp-2000/outgroup.0.species.1.5.genes.1.5
* sp-2000/outgroup.0.species.5.genes.5
* sp-2000/outgroup.1.species.0.15.genes.0.15
* sp-2000/outgroup.1.species.1.5.genes.1.5
* sp-2000/outgroup.1.species.5.genes.5
* sp-2000/outgroup.2.species.1.5.genes.1.5
* sp-2000/outgroup.4.species.1.5.genes.1.5
* sp-5000/outgroup.0.25.species.1.5.genes.1.5
* sp-5000/outgroup.0.5.species.1.5.genes.1.5
* sp-5000/outgroup.0.species.0.15.genes.0.15
* sp-5000/outgroup.0.species.1.5.genes.1.5
* sp-5000/outgroup.0.species.5.genes.5
* sp-5000/outgroup.1.species.0.15.genes.0.15
* sp-5000/outgroup.1.species.1.5.genes.1.5
* sp-5000/outgroup.1.species.5.genes.5
* sp-5000/outgroup.2.species.1.5.genes.1.5
* sp-5000/outgroup.4.species.1.5.genes.1.5

The data are divided into multiple files:

* D2-true.tar.gz: The true species trees and gene trees for 2000 (sp-2000.tar.gz) and 5000 (sp-5000.tar.gz) taxon datasets 
  * Each directory (as noted above) includes:
	* `s_tree.trees`: True species tree file
	* `truegenetrees`: True gene trees
	* `s_tree.trees.score`: The quartet score of the true species tree with respect to the true gene trees
* D2-estimated.tar.gz: Estimate gene tree file (`estimatedgenetre.gtr`), one per model condition / replicate combination 
* D2-estimated-rooted.tar.gz
	* results of methods MP, MV, Outgroup rooting (OG) (`estimatedgenetre.gtr_MP`, `estimatedgenetre.gtr_MV`, `estimatedgenetre.gtr_OG`). 
* D2-true-rooted.tar.gz: Similar to `D2-estimated-rooted.tar.gz`  but with true gene trees used as input

### D3 dataset

D3—ASTRAL-II dataset: We reused a previously published dataset (http://www.ncbi.nlm.nih.gov/pubmed/26072508) to investigate performance for intermediate number of species (i.e, 10, 50, 100, 200, 500, and 1000).

* D3-rooted.tar.gz: We simply provide the rooted trees (produced in this paper) here and refer the reader to the original publication for the original data
	* results of methods MP, MV, Outgroup rooting (OG) (`estimatedgenetre.gtr_MP`, `estimatedgenetre.gtr_MV`, `estimatedgenetre.gtr_OG`). 


## Sharing/Access information

See more on 
	* [more info](https://uym2.github.io/MinVar-Rooting/)
	* [Github page](https://github.com/esayyari/simphy-uyen) for scripts and pipelines.

## Code/Software

These are generated using Simphy
