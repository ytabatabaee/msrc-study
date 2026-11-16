# Atlantic cod Patterson's D validation

This validation supplements the frozen primary analysis at `../`. It freezes 12 additional population triples before scanning their D values, then independently reparses deposited SNAPP, Dsuite, and Relate VCFs. Only SNAPP and Dsuite contain the external Borsai and Gadmor outgroup genotypes required for the rooted fixed-outgroup estimator. Relate is used for source/sample/position diagnostics only; it lacks an external outgroup and must not be used to claim rooted D.

## Reproduce

From the repository root:

```sh
python3 empirical/atlantic_cod/analysis/patterson_d/validation/scripts/freeze_validation_quartets.py
python3 empirical/atlantic_cod/analysis/patterson_d/validation/scripts/run_validation.py
python3 -m unittest discover -s empirical/atlantic_cod/analysis/patterson_d/validation/tests -v
```

The VCFs are large and intentionally excluded from Git. The run manifest records source file hashes. The scanner reports pooled ABBA/BABA counts, site contributions and per-window summaries. Dsuite `.full` inputs allow missing calls; this runner requires complete calls among each frozen ingroup sample set and a fixed homozygous outgroup at each included site. Genotypes are frequency averaged. No REF allele is assumed ancestral. Windows overlapping frozen inversion boundaries are omitted from inside/outside summaries.

Quartet effects are descriptive because designs share populations and are not independent replicates. The validation does not equate allele-sharing asymmetry with introgression or an MSRC-specific mechanism.

The deposited Relate colinear VCF is specifically restricted to LG03–LG06, LG08–LG11, and LG13–LG23; it contains no records from LG01, LG02, LG07, or LG12. Its 24 cod samples also omit external taxa. It therefore cannot validate regional nucleotide D, despite being a less aggressively filtered genotype source. Source overlap is reported explicitly rather than imputing missing regions.

Regenerate the validation figures after a complete scan with:

```sh
MPLCONFIGDIR=/tmp/cod-matplotlib python3 empirical/atlantic_cod/analysis/patterson_d/validation/scripts/make_validation_figures.py
```

`run_validation.py --finalize-existing` rebuilds block summaries and the run manifest from the completed per-source checkpoint and per-window tables without rescanning the large VCFs.
