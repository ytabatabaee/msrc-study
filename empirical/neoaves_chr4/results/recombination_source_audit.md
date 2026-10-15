# Neoaves chromosome 4 recombination-source audit

## Outcome

The audit found published quantitative avian recombination maps, including two
chicken pedigree linkage-map sources that contain marker-defined local rates:

* Groenen et al. (2009), DOI [10.1101/gr.086538.108](https://doi.org/10.1101/gr.086538.108), reports a 9,268-marker chicken linkage map and approximately 1-Mb marker-defined rate bins. It includes sex-average and sex-specific maps and explicitly includes GGA4.
* Elferink et al. (2010), DOI [10.1186/1471-2156-11-11](https://doi.org/10.1186/1471-2156-11-11), reports a 13,340-SNP map from two pedigree populations and 1,819 approximately 500-kb bins with sex-average and sex-specific cM/Mb rates. Its supplementary S1, S2, and S5 files are the most promising raw numeric route for the focal interval.

The Elferink source is suitable in principle for a genuine local cM/Mb track,
but its public supplementary endpoint returned HTML protection pages rather
than the XLS payload in this environment. The source therefore was not treated
as downloaded numeric data. Its coordinates are also reported on WASHUC2, so a
GalGal6 coordinate conversion must be validated before any interval can be
aligned to the frozen Neoaves outputs.

No validated marker-level source file was recovered locally. Consequently this
work does not create a coordinate-resolved recombination table, a cM/Mb figure,
an inside/outside rate comparison, or a significance test. The requested
`neoaves_chr4_recombination_rates.png` and `.pdf` files are intentionally absent.

## Source assessment

The full candidate inventory, coordinate systems, marker counts, measurement
types, sex/population scope, coverage, download status, and decision for each
source is in [recombination_source_audit.tsv](recombination_source_audit.tsv).

Rabie et al. (2004) is a chromosome 4 radiation-hybrid/BAC physical map with
genetic anchors. It is useful for correspondence context but is not a
recombination-rate map. Stapley et al. (2010) provides a pedigree-based zebra
finch map and 1-Mb cM/Mb estimates, but synteny alone does not establish that a
zebra finch interval is the same physical interval as GalGal6 chr4. The original
Neoaves PNAS study provides the approximately 21-Mb GalGal6-referenced focal
region and an inference of historical recombination suppression, but no direct
quantitative cM/Mb measurements.

The published GGA4 chromosome averages illustrate why a local track cannot be
substituted by a scalar: Elferink et al. report 195.2 cM / 94.2 Mb = 2.07 cM/Mb
for their combined map, and population-specific averages of 237.6/94.2 = 2.52
and 183.7/94.2 = 1.95 cM/Mb. These are whole-chromosome summaries, not rates in
the approximately 21-Mb focal region, and are retained only as audit metadata.

## Required recovery path

1. Download the original binary `1471-2156-11-11-S1.XLS`, `S2.XLS`, and `S5.XLS` files from the publisher URLs recorded in the TSV. Confirm with `file` and an XLS parser that each payload is an actual workbook, then record SHA-256 checksums.
2. Extract the GGA4 marker and 500-kb-bin records. Retain the combined, female, and male population columns separately.
3. Validate WASHUC2 marker positions against a sequence-based marker mapping or a documented WASHUC2-to-GalGal6 chain. Reject markers with unresolved, duplicated, chromosome-inconsistent, or order-inconsistent mappings.
4. For consecutive validated markers or published bins, calculate `r = delta_cM / delta_Mb`. Preserve source bin boundaries and missing intervals; do not interpolate over large gaps.
5. Join the resulting GalGal6 intervals to the frozen structural and quartet tables. The frozen structural coordinates are in the GalGal6 chr4 axis and include the approximately 21-Mb Neoaves focal context; no rearrangement intervals should be redefined from quartet support.

If the conversion succeeds, the follow-up should report marker counts, interval
spacing, covered fraction, sex/population variation, and boundary-proximal
coverage before any comparison. Present-day chicken recombination cannot be
interpreted as ancestral recombination suppression between alternative Neoaves
arrangements, especially where an inversion is fixed within a species. The
MSRC parameter `m` must not be estimated from cM/Mb without a population-genetic
model connecting the measurements.

## Validation performed

The repository-side audit checks were run with
`tests/test_recombination_audit.sh`. They verify the audit schema, the published
whole-chromosome arithmetic, the GalGal6 coordinate bounds of the frozen
structural/quartet inputs, and that no recombination figure is emitted while a
validated local rate table is absent. The test does not claim that the external
XLS files were recovered.
