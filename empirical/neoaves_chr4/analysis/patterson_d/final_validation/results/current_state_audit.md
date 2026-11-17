# Current-state audit

## Existing verified work

- The original Stiller 63,430 intergenic alignment archive and associated local tree archives are present under `empirical/neoaves_chr4/external/stiller2024/`; the prior input audit recorded 6,247 anchored chr4 loci and exact locus-ID matches to the alignment archive.
- The saved nucleotide locus table contains ABBA/BABA counts for six frozen combinations (three primary quartets and the same ingroup quartet with *Meleagris gallopavo* outgroup) and two filters (all informative patterns and transversions only). Pooled regional results are reproducible from those rows without re-reading 3.6 GB of alignments.
- Existing 100-kb, 250-kb, 500-kb and 1-Mb display tables and 500-kb, 1-Mb and 2-Mb physical block-jackknife summaries exist. The 250-kb output is reused for Figure A.
- The source tree archive includes the exact quartet taxa for 5,727 Columbea, 5,648 N61 and 5,442 N62 chr4 trees. Exact induced quartet frequencies, matched 1-Mb D tracks and branch-length summaries were newly computed from the named source trees.
- Existing tests cover the original estimator and track generation. The new final-validation tests independently check scalar site-pattern cases and result integration.

## Superseded or incomplete evidence

- The top-level `interpretation_report.md` predates the independent rooted-quartet audit and should not be read as validating D=0 for all three frozen designs. The frozen source tables are retained unchanged.
- The `recheck_20261108/` directory has four pre-existing deleted files in the working tree (input audit, locus mapping, locus results and summary). They were already deleted at task start and were not restored or modified. Top-level locus-level results were used where present.
- Stage 4D is explicitly marked as incomplete reproduction; its chr4 species-tree/topology output is not used to set the D null.
- The published 48-taxon tree is available as supplement/figure evidence, but an exact local Newick and *Ciconia maguari* tip are unavailable. The 48-taxon exact quartet cannot be reconstructed without inventing a taxon mapping.

## Remaining scientific weaknesses

- Columbea's frozen P1/P2 relationship conflicts with the induced quartet in the published Stiller 363-taxon reference; Jarvis 2014 supports broad Columbea but lacks the exact *Ciconia* tip. Its conventional D=0 interpretation is disputed.
- N61 and N62 are consistent with the Stiller reference quartet, but the exact Jarvis quartet cannot be checked for *Ciconia*. Nonzero genome-wide D means the simple D=0 baseline does not hold empirically.
- The alternative outgroup sensitivity is a second Galliform, not a distant independent clade. Transversion-only results reduce transition-driven finite-sites concern but are much sparser; they do not constitute a fitted finite-sites demographic null.
- Ingroup representative sensitivity beyond the three frozen quartets was not run. The saved nucleotide locus outputs encode no additional ingroup representatives; doing so requires a new alignment scan. Quartet occupancy varies from about 86% to 92% of chr4 loci across frozen designs.
- The independent scalar estimator tests use synthetic alignments and compare against the production routine; no sampled raw archived alignment was independently re-decoded. Alignment-length equality, duplicate FASTA identifiers, and per-taxon callability were not audited from raw alignments in this finalization.
- Nucleotide site columns and nearby 10-kb loci are linked. Uncertainty is therefore described with genomic block jackknife/bootstrap, not site-level binomial tests. No spatial randomization p-value was generated.
- Exact nucleotide-to-genealogy matching is possible for the exact four species, but only on loci with a full four-tip local tree. Those complete-tree subsets differ by design and are not every nucleotide locus.
