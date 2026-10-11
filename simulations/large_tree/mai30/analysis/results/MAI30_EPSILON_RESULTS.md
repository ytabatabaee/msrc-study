# Mai30 epsilon-mixture results

## Design

The epsilon mixtures use exactly 500 gene trees per dataset and three
source-pool conditional mixture replicates at each epsilon. The source pools
are matched MSC controls and complete MSRC collections. Each MSRC source pool
contains 25 independent WF histories × 20 loci; the mixture manifest retains
that clustering. The incomplete 074 and 005 MSRC conditions are excluded.

Datasets and manifests are under
`datasets/processed/mai30_epsilon_mixtures/`. The ASTRAL summary is
`analysis/results/mai30_epsilon_astral.tsv`.

## Interpretation

The epsilon experiment is conditional on six realized source pools. Resampling
source trees creates mixture sampling variation, while the manifest allows
history-level uncertainty to be analyzed separately. It does not create new
biological histories.

The topology results should be summarized by replicate, suppression treatment,
and epsilon with replicate-level intervals. Any apparent trend must be checked
against the three mixture resamples and the 25-history source clustering.
The ASTRAL CULength results use the audited common-split method and include
signed, absolute, and median errors. Long true branches in replicate 042 are
subject to quartet-estimator saturation, as documented in
`MAI30_BRANCH_LENGTH_AUDIT.md`.

No result from this source-pool experiment is labeled statistical
inconsistency. A valid population failure claim requires independent history
expansion and analysis of the expected ASTRAL quartet objective.

## Reproduction

```bash
/Users/ytabatabaee/Desktop/msrc-sim/.venv/bin/python \
  simulations/large_tree/mai30/analysis/scripts/build_epsilon_mixtures.py \
  --replicates 3
/Users/ytabatabaee/Desktop/msrc-sim/.venv/bin/python \
  simulations/large_tree/mai30/analysis/scripts/analyze_epsilon_astral.py
/Users/ytabatabaee/Desktop/msrc-sim/.venv/bin/python \
  simulations/large_tree/mai30/analysis/scripts/analyze_epsilon_quartets.py
```
