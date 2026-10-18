# Stage 10 legacy validity audit

The existing Stage 10 analysis treated per-window SNAPP population-tree
internal node ages from `stage7_pairwise_mrca_times.tsv` as realized lineage
coalescence times drawn from an MSRC density. That observation model is not
valid: SNAPP integrates over gene trees, so a population-tree divergence
estimate is not a realized gene-lineage coalescence time.

The legacy construction found 11,673 of 27,786 audited non-boundary
observations with `T < d` (42.010%), which is a symptom of this
mismatch. Its lambda and m values are exploratory only and should not be
reported as empirical MSRC parameter estimates.

Stage 10b does not use Stage-7 MRCA times, the legacy Stage-10 scaffold, or
the legacy parameter estimates. It fits the validated conditional 2:2
quartet model directly to quartet outcomes and evaluates held-out
contiguous genomic blocks. The frozen Stage-10 files were hash-audited and
remain unchanged.

Legacy file hashes are recorded in the Stage-10b manifest. The primary
Stage-10b estimands are effective reduced-model parameters and rho=m/lambda;
they are not meiotic recombination rates or full biological MSRC estimates.
