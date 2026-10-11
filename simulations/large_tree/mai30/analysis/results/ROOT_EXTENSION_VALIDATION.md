# Root-extension validation

## What the current code does

`msrc-sim` represents the root as a population branch with
`older_age = root_age + root_extension`. However, both
`simulate_msc_genealogy` and the final coalescent loop in
`simulate_genealogy` continue ordinary pairwise coalescence whenever all
remaining lineages are in the root population. There is no finite-endpoint
truncation, forced coalescence, or unresolved-genealogy return at
`older_age`. The `max_time` variable in the MSRC loop is currently not used as
a stopping condition.

Thus, for ordinary MSC controls, the ancestral population is continued
indefinitely. This is consistent with the standard MSC treatment of the
population above the species-tree root, and it is preferable to imposing an
arbitrary finite truncation. The `root_extension` value does not affect MSC
genealogy generation in the current implementation.

The forward rearrangement history is a separate issue. Its records stop at the
configured root extension and `FrequencyHistory.frequency_at` holds the oldest
state beyond the final record. That is a finite forward-history approximation;
it should not be described as an infinite rearrangement history merely because
the coalescent continuation is infinite.

## Empirical same-seed check

For Mai condition `outgroup.0.species.1.5.genes.1.5`, replicate `042`, I
generated 100 ordinary MSC trees at `R = 5, 10, 20, 40` coalescent units using
the same seed. All four Newick collections have the same SHA-256. The result
is recorded in:

`analysis/results/root_extension_same_seed.tsv`

This confirms that the current ordinary-MSC implementation has no finite root
endpoint whose multiple-lineage frequency can be measured. The endpoint field
is therefore recorded as “not applicable,” rather than assigning a misleading
completion rate.

## Generation conversion

For a replicate with SimPhy haploid size `N_H`, the exact mapped diploid
coalescent size is `N_H/2`, so an extension of `R` coalescent units is
`R*N_H` generations. Finite Wright–Fisher histories use the documented integer
mapping in `population_mapping.py`, giving `R*2*N_msrc_integer` generations.

No finite root extension is required for the ordinary MSC validation. A
separate MSRC sensitivity analysis should report the forward-history extension
explicitly and test whether it changes rearrangement-frequency outcomes.
