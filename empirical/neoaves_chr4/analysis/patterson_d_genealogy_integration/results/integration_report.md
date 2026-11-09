# Neoaves Patterson's D / genealogy integration

## Data correspondence
The Stage-4D manifest contains 6247 chr4 loci. Patterson, QQS, and distance inputs were joined through exact manifest `locus_id` mappings from the integer `tree_index` fields; no positional join was used. The full audit is in `locus_mapping_audit.tsv`. The QQS table contributed 6132 exact N61 locus/clade matches (115 N61 clade records are absent from the 6247-locus manifest join); the distance table has 6140 exact chr4 locus matches. Patterson N61 all-site records cover 5648 unique loci, with 1221 inside and 4427 outside the PNAS mask; the remaining manifest loci have no N61 all-site record. The complete per-locus missingness audit is retained in `locus_mapping_audit.tsv`.

The published PNAS intervals are GalGal4 window-start intervals, while Stiller local-tree coordinates are GalGal6. Existing repository provenance supplies approximate GalGal6 counterparts but no exact lift; this analysis uses the previously validated published-predicate mask and does not claim nucleotide-resolution boundary equivalence.

## Main N61 result

At 250 kb, nucleotide D is compared with q2-q3, because negative Patterson D (BABA > ABBA) corresponds to q2 > q3 under the frozen N61 taxon mapping. The existing q2-q1 track is retained separately. Correlations use nonoverlapping 250-kb windows; block bootstrap intervals and circular-shift spatial nulls are reported in `spatial_correspondence.tsv`.

The original locus-level pooled Patterson analysis gives N61 ΔD = -0.177615 (inside D = -0.258194; outside D = -0.080579). The 250-kb window classification gives ΔD = -0.174987 because windows overlapping a published boundary are assigned to the interval; this small difference is a classification consequence, not a change to the frozen result.

N61 all-PNAS inside-minus-outside effects (D is pooled from ABBA/BABA counts before division):
- Patterson_D: Δ inside−outside = -0.175, block-jackknife SE = 0.0113, inside/outside windows = 87/274
- q2_minus_q3: Δ inside−outside = 0.8596, block-jackknife SE = 0.02174, inside/outside windows = 86/274
- topology_shift_q2_minus_q1: Δ inside−outside = 0.8739, block-jackknife SE = 0.02573, inside/outside windows = 86/274
- B_observed: Δ inside−outside = 0.009679, block-jackknife SE = 0.006749, inside/outside windows = 86/274
- quartet_internal_branch: Δ inside−outside = 0.006875, block-jackknife SE = 0.0006853, inside/outside windows = 86/274

N61 250-kb spatial correspondence:
- D vs q2_minus_q3: r = -0.579, 95% 1-Mb-block bootstrap CI [-0.7035, -0.4237], circular-shift p = 0.005525
- D vs topology_shift_q2_minus_q1: r = -0.5354, 95% 1-Mb-block bootstrap CI [-0.6766, -0.3795], circular-shift p = 0.01381
- D vs B_observed: r = -0.0395, 95% 1-Mb-block bootstrap CI [-0.2005, 0.1126], circular-shift p = 0.605
- D vs quartet_internal_branch: r = -0.4706, 95% 1-Mb-block bootstrap CI [-0.5951, -0.3308], circular-shift p = 0.005525

The q2-q3 and D relationship should be interpreted as a sign-reversed discordance comparison, not as the same estimand. B_observed is a separate arrangement-blind substitution-distance summary.

## Branch-length feasibility

The local named trees yielded 6247 chr4 N61 quartet records. Status counts are {'resolved': 5648, 'missing_taxon': 599}. Quartet-level topology, internal branch, total quartet length, and pairwise distances are in `quartet_branch_lengths.tsv`. Missing taxa, unresolved induced quartets, and parse failures are retained as status fields.

## Interpretation

A spatially coincident nucleotide asymmetry and discordant-topology contrast would be compatible with a shared local genealogy regime, but neither statistic uniquely distinguishes MSRC from introgression, incomplete lineage sorting, mutation-rate variation, alignment error, or other causes. Branch lengths remain substitution distances. The joint analysis is therefore additive evidence of spatial association, not a mechanistic identification of MSRC.

Candidate Results wording: The N61 nucleotide site-pattern asymmetry was spatially aligned with the independently defined chromosome-4 outlier intervals and was compared with the corresponding q2-q3 discordant-topology contrast. Because ABBA and BABA map to q3 and q2 respectively under the frozen quartet, the sign-reversed contrasts were directionally comparable. This association is compatible with a shared recombination-partitioned genealogy regime, but does not by itself discriminate MSRC from introgression or incomplete lineage sorting.
