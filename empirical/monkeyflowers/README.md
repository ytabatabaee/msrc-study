# Channel Island monkeyflowers

Topology-only spatial analysis of Short & Streisfeld (2025), *Molecular Ecology* 34:e17778, using locally supplied Dryad TWISST outputs. The input files are under `data/raw/` and are excluded from Git. See `results/stage0_data_audit.md` for file validation and `results/stage1_topology_analysis_report.md` for results, caveats, and the conditional-go recommendation.

Run the analysis from the repository root:

```sh
python empirical/monkeyflowers/scripts/analyze_topologies.py
python -m unittest discover -s empirical/monkeyflowers/tests -v
```

The analysis preserves the 15 rooted TWISST classes. It defines the species-history topology as topo3, recent-introgression support as the sum of topo9/topo12/topo13, ancient-introgression support as topo10, and retains the other 11 classes. The count-only weights file has no coordinates; its rowwise normalized values are checked against the combined coordinate-bearing export, then the export is reconciled to metadata and coordinate grid by genomic keys.

**Scope limitation:** no local Newick trees are present in the available inputs. Results support topology-weight and spatial analyses only. No A(w), B(w), internal branch, CU, or gene-tree species-tree inference is reported. The topology-support runs are not introgressed ancestry-tract lengths.
