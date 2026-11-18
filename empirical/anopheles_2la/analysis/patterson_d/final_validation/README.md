# Anopheles 2La Patterson's D final validation

This additive package audits the completed Ag3.10 analysis without altering its frozen designs, source tables, or primary results. The derivation script reads the retained 15,808-row window table and Stage 2R topology table. It writes new result tables and four figure sets.

Run from the repository root:

```bash
python empirical/anopheles_2la/analysis/patterson_d/final_validation/scripts/derive_validation.py
python -m pytest -q empirical/anopheles_2la/tests/test_patterson_d.py empirical/anopheles_2la/analysis/patterson_d/final_validation/tests
```

The script does not access the network or refetch genotypes. Individual-level leave-one-out, balanced subsampling, direct geography-adjusted estimates, and an expanded raw-site audit require Ag3.10 genotype access. The ordinary interpreter in the validation session did not have `malariagen_data`; the existing `.venv-malaria310` environment is documented by the original workflow and should be checked before a future raw-data run. Do not use a different release.

The topology track is the Stage 2R species/arrangement-haplotype quartet, not an exact match to the nucleotide population quartets. Figure 2 is explicitly descriptive. Correlation P values are not inferential because genomic positions are dependent; use physical-block uncertainty for future formal association tests.
