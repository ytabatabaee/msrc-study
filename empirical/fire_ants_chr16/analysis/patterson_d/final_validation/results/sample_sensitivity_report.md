# Sample-level sensitivity

The original frozen alternative sampling analyses include all-geography and
*S. invicta*-only focal sets, plus Brazil-only and French-Guiana-only
*S. saevissima* outgroups. The original report finds that these preserve the
positive supergene and near-zero outside direction. They address geography and
outgroup composition, not individual leverage.

For Q1, existing alternative geographic/outgroup supergene estimates fall near
+0.802 to +0.832, with chr16A+B near zero; Q2 all-geography and
*S. invicta*-only estimates are +0.844 and +0.850. Exact per-design estimates
and counts are retained in `sample_sensitivity.tsv`.

Individual leave-one-out and balanced subsampling are not present in saved
source outputs. The per-individual genotype scan is required because removing an
individual changes group call thresholds and site-specific outgroup unanimity.
No claim that no single individual dominates is made here. `sample_sensitivity.tsv`
marks these analyses explicitly as not computed; they remain a paper-readiness
limitation.
