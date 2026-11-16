# Frozen quartet design audit

Quartets were defined before any nucleotide D values were computed. Population/state membership comes from the explicit Supplementary Data 1 `Supergene Variant` field and its independently reproduced `stage6c_individual_to_group.tsv` mapping, not from sample-name substrings or local topologies. Candidate individuals are exact local-tree tips, male, have source missing fraction below 0.05, and are not marked as pooled workers or random-allele treated. Every chosen ID is in the released VCF.

The independent background tree is the published chromosome 1–15 ASTRAL history in `data/fire_ants_chr16/processed/stage3_background_tree.nwk`, rooted to *S. geminata* in the source study. For each retained quartet, the script checked that the P1+P2 most recent common ancestor excludes P3 and O and that the P1+P2+P3 ancestor excludes O. All eight designs pass. Cross-species P1/P2 arrangement-control proposals were rejected because the baseline tree does not support their P1/P2 sister assumption.

## Primary population-history / arrangement contrasts

- `Q1_primary_Argentina`: invicta/macdonaghi SB (94), invicta/macdonaghi Sb (40), richteri Sb (32), *S. saevissima* O (6).
- `Q2_primary_Argentina`: richteri SB (28), richteri Sb (32), invicta/macdonaghi Sb (40), *S. saevissima* O (6).

Focal samples are restricted to Argentina to limit geographic imbalance. Q1 compares arrangements within the invicta/macdonaghi group against the richteri Sb population; Q2 compares arrangements within richteri against invicta/macdonaghi Sb. The independent background tree supports the grouping of the two invicta/macdonaghi arrangement strata and of the two richteri strata, while supporting the focal ingroup relative to the selected outgroup.

## Sensitivity designs

Two S. invicta-only focal sensitivities address the invicta/macdonaghi mixture; two all-geography versions evaluate geographic sampling; and Brazil (4 individuals) and French Guiana (2 individuals) subsets evaluate the only clean outgroup species' geographic composition. The latter are outgroup sensitivities, not independent outgroup-species replications. Exact member IDs and counts are in `config/frozen_sample_membership.tsv`; design rationale and topology checks are in `config/frozen_quartets.tsv`.

## Outgroup interpretation

The primary outgroup is *S. saevissima*, represented by six individual males: four Brazil and two French Guiana. They are outside the focal ingroup clade in the independent background tree and are not in the source pipeline's random-allele or pooled-worker sets. *S. geminata* has only a worker-pool genotype and is not used as an individual-level primary outgroup. *S. pusillignis* has one worker pool and one male from the random-allele list; neither provides a clean individual-level alternative. Outgroup genotypes still require per-site fixation/missingness filtering and cannot be assumed to polarize every site.

No design was chosen based on observed nucleotide patterns. This file and both frozen TSVs must remain unchanged during the subsequent D run unless a documented data-integrity issue requires a versioned amendment.
