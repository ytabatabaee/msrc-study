# House mouse chr17 Patterson D

The analysis fetches indexed regions across the full chr17 contig (mm10/GRCm38 POS 1–94,987,271) directly from the Harr et al. AllMouse PASS VCF; it does not download the 55.7 GB compressed VCF. Inversion/collinear/flank contrasts focus on chr17:5.0–42.05 Mb; exact tree overlap is 5.0–37.055 Mb. `config/frozen_quartets.tsv` records the contrasts, and `config/frozen_quartet_members.tsv` records each sample’s pre-scan P1/P2/P3/outgroup assignment. The source Kelemen & Vicoso archive is inspected separately and contains gene trees/coordinate annotations rather than genotype data.

## Reproduce

From the repository root, using `.venv-malaria310` with `pysam` and `certifi` installed:

```sh
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/scripts/run_patterson_d.py
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/scripts/finalize_patterson_d.py
.venv-malaria310/bin/python empirical/house_mouse_t_complex/analysis/patterson_d/tests/test_synthetic_counts.py
```

The scan checkpoints counts in `results/checkpoint_counts.json` every 500 kb and resumes completed chunks. Delete that file only to force a clean rescan. Remote indexed access requires network connectivity to the published VCF host.

For each biallelic one-base A/C/G/T SNP with all eight *M. spretus* samples called and fixed for REF, population ALT frequencies p1–p3 yield ABBA=(1-p1)p2p3 and BABA=p1(1-p2)p3. Missing individuals are excluded from that population's frequency mean; sites with an entirely missing population are tallied. Sites with any outgroup ALT frequency are unpolarized and excluded. D is formed from pooled contributions, not mean window D. Invariant sites are not in the SNP-only source. Indels, MNPs, multiallelic and ambiguous records are excluded.

The designs estimate exploratory arrangement-associated population allele-frequency asymmetry. They do not isolate t-haplotype chromosomes because calls are unphased and carrier samples are diploid. The inversion/collinear contrasts are descriptive; geography, subspecies ancestry, SNP ascertainment, admixture, ancestral structure and incomplete lineage sorting remain alternatives. Robustness scans for transversions, alternative outgroups and SNP filters have not been run.


Final validation bundle: see `final_validation/README.md` for exact reproduction commands, completed robustness scan, independent scalar tests, new figures, and explicit distinction between the matched physical span and exact tree-matched windows. The sample-level scan in `final_validation/scripts/sample_sensitivity.py` was attempted but the endpoint truncated BGZF at byte offset 48,567,409,620 before outputs completed; no sample-level estimates are reported.
