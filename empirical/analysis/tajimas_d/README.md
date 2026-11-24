# Empirical Tajima's D

The shared estimator in `tajimas.py` implements classical fixed-n Tajima's D, pairwise diversity totals, Watterson theta, and per-callable-base summaries only when the caller supplies a validated denominator. Variant-only panels are never treated as callable sequence. The estimator rejects unsupported ploidy and labels undefined n/S/variance cases. Eight unit tests pass; the implementation was independently checked against R `pegas::tajima.test` on a four-sequence hand calculation.

## Reproduce

From repository root:

```bash
PYTHONPATH=empirical/analysis/tajimas_d python3 -m unittest discover -s empirical/analysis/tajimas_d/tests -v
python3 empirical/analysis/tajimas_d/scripts/write_sample_group_manifests.py
python3 empirical/fire_ants_chr16/analysis/tajima_d/scripts/run_fire_ants.py
python3 empirical/fire_ants_chr16/analysis/tajima_d/scripts/summarize_fire_ants.py
python3 empirical/atlantic_cod/analysis/tajima_d/scripts/run_cod.py
python3 empirical/atlantic_cod/analysis/tajima_d/scripts/join_cod_tracks.py
SSL_CERT_FILE="$PWD/.venv-malaria310/lib/python3.10/site-packages/certifi/cacert.pem" REQUESTS_CA_BUNDLE="$PWD/.venv-malaria310/lib/python3.10/site-packages/certifi/cacert.pem" .venv-malaria310/bin/python empirical/anopheles_2la/analysis/tajima_d/scripts/run_anopheles.py
.venv-malaria310/bin/python empirical/anopheles_2la/analysis/tajima_d/scripts/calculate_background_contrasts.py
.venv-malaria310/bin/python empirical/anopheles_2la/analysis/tajima_d/scripts/analyze_track_relationships.py
python3 empirical/analysis/tajimas_d/scripts/build_cross_dataset_summary.py
python3 empirical/analysis/tajimas_d/scripts/make_summary.py
```

Anopheles performs bounded 2-Mb Ag3.10 API queries, checks each returned position against non-N AgamP4 reference positions and verifies the `gamb_colu` site mask. It serializes sufficient statistics before output generation. If only output formatting needs to be rerun, use `--summary-only` on `run_anopheles.py`.

The Neurospora analysis was completed previously and is preserved; this follow-up does not rerun or alter it.

## Interpretation

Read `results/analysis_report.md` and `results/dataset_feasibility.tsv` before interpreting estimates. `VALIDATED` means the group SFS and callability rules were validated for that data source; it does not mean the biological comparison is causal. `EXPLORATORY` marks SNP-panel estimates with unverified completeness or gene ascertainment. `BLOCKED` means no suitable population-level nucleotide SFS was available. No classical Tajima's D neutrality p-values are interpreted under linkage/structure. Spatial uncertainty uses contiguous blocks.
