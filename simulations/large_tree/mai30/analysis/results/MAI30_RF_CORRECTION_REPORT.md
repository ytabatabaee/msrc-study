# Mai30 RF-distance correction report

## What changed

The earlier Mai analysis used an explicit split symmetric difference for
binary trees, but it was limited to species-tree discordance summaries and the
small set of pilot collections. It did not provide strict duplicate-taxon
checks, an unresolved-tree normalization rule, a 300-replicate matching
manifest, or within-pool and cross-pool RF distributions. The prior
MSRC/original comparison also treated distributional summaries as the main
result without reporting the requested within-original and within-MSRC
pairwise baselines.

The corrected implementation is `analysis/scripts/rf_metrics.py` plus
`analysis/scripts/analyze_mai_rf.py`. It:

- parses unrooted topology through nontrivial bipartitions and ignores root
  placement and branch lengths;
- requires identical taxon sets and rejects duplicate or malformed labels;
- uses `2(n-3)` for two fully resolved `n`-taxon trees;
- uses `|S_1|+|S_2|` for an unresolved comparison, so the normalized value is
  bounded by one and the resolution status is retained;
- preserves replicate IDs as three-digit strings;
- matches condition, replicate, and MSRC treatment before reading trees;
- computes reproducible random cross-pool RF pairs, within-original pairs,
  within-MSRC pairs, species-tree RF, and descriptive index-paired RF;
- records missing, partial, and invalid collections explicitly.

The RF unit tests cover identical trees, rerooting, branch-length changes,
quartet NNI differences, unresolved trees, duplicate taxa, mismatched taxa,
and an independent DendroPy calculation. The file-matching tests verify
leading-zero replicate IDs and prevent substitution across conditions.

## Validation and completion status

All 300 original no-outgroup Mai collections were found and validated as
500-tree, 30-taxon collections across the three main conditions:

| Condition | Original collections | Original trees | Complete MSRC collections | Partial MSRC collections | Missing MSRC collections |
|---|---:|---:|---:|---:|---:|
| species 0.15 | 100 | 50,000 | 0 | 0 | 100 per treatment |
| species 1.5 | 100 | 50,000 | 2 per treatment | 1 per treatment | 97 per treatment |
| species 5 | 100 | 50,000 | 0 | 1 per treatment | 99 per treatment |

The six complete MSRC collections are replicates 007 and 042 under the three
available treatments. Replicates 005 and 074 are partial 20-tree smoke
collections. Therefore the requested 300-replicate MSRC analysis is **not
complete** and no 100-replicate MSRC aggregate is claimed.

## Corrected RF results

The per-replicate files are:

- `mai30_rf_replicate_results.csv` and `.tsv`
- `mai30_rf_validation_manifest.csv` and `.tsv`
- `mai30_rf_aggregate_results.csv` and `.tsv`

The original gene-tree-to-species-tree means, calculated over all 100 valid
original replicates in each condition, are:

| Condition | Mean normalized RF |
|---|---:|
| species 0.15 | 0.4511 |
| species 1.5 | 0.4506 |
| species 5 | 0.4520 |

The currently available complete MSRC gene-tree-to-species-tree means are only
two-replicate summaries for the species-1.5 condition:

| Treatment | Complete-replicate mean nRF | Complete replicates |
|---|---:|---:|
| no suppression | 0.4910 | 007, 042 |
| moderate suppression | 0.4914 | 007, 042 |
| strong suppression | 0.4910 | 007, 042 |

The original Mai means and the two-replicate MSRC values are descriptive. The
large confidence intervals in `mai30_rf_aggregate_results.csv` reflect the
fact that only two complete MSRC replicates are available for that condition.

## Effect on previously reported results

The corrected unrooted binary RF values for the six complete pilot
collections agree with the previously reported species-tree RF values. The
correction therefore does not change those six gene-tree-to-species-tree
means. It adds the missing distributional baselines and makes the incomplete
300-replicate coverage explicit. The earlier 12-row pilot tables remain
preserved.

Index-paired RF is included only as a descriptive diagnostic. Original Mai
and MSRC gene trees were independently simulated, so equal gene indices do
not identify the same genealogy. The primary cross-pool metric uses
deterministic random pairing, and the within-pool metrics quantify baseline
ILS and finite-sample variability.

## MSC validation and remaining issues

The existing ordinary MSC validation remains a distributional comparison for
the four selected Mai replicates; it is not a zero-RF paired-tree test. The
new RF implementation confirms that the existing species-tree RF summaries
are topologically well-defined for the fully resolved inputs. It does not
resolve the separate scientific limitations of the scaled-WF approximation or
the lack of 300 complete MSRC collections.

No changes were made to `msrc-sim`. Suspected simulator behavior remains
documented in the existing compatibility and validation reports.

Reproduce the corrected analysis with:

```bash
MPLCONFIGDIR=/tmp/mai30-mpl \
/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria/bin/python \
  simulations/large_tree/mai30/analysis/scripts/analyze_mai_rf.py
```

Tests:

```bash
MPLCONFIGDIR=/tmp/mai30-mpl \
/Users/ytabatabaee/Desktop/msrc-paper/.venv-malaria/bin/python -m pytest -q \
  simulations/large_tree/mai30/analysis/tests
```
