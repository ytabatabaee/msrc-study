# t-complex inversion interval provenance

## Result

The four historical t-complex inversions are shown in the updated spatial
figure as approximate intervals: `In(17)1` (5–12 Mb), `In(17)2` (12–25 Mb),
`In(17)3` (25–31 Mb), and `In(17)4` (31–40 Mb). These are schematic
placements rounded to approximately 0.5-Mb display precision, recorded in
`empirical/house_mouse_t_complex/results/t_complex_inversion_intervals.tsv`.
The rectangles are hatched and labeled approximate in the figure. This is a
rounding/display convention, not a confidence interval: the actual boundary
uncertainty is not reported by the schematic and may be larger than 0.5 Mb.

## Evidence

Kelemen & Vicoso (2018, *Genetics* 208:365–375, DOI
10.1534/genetics.117.300513) describe the t-haplotype as a roughly 40-Mb
chr17 region and reproduce the four putative nonoverlapping inversion
placements in Figure 1. Their text explicitly calls the locations putative;
the paper and its supplementary material do not provide a machine-readable
base-pair breakpoint table. The primary structural literature establishes the
four nonoverlapping inversions genetically and cytologically, but predates a
reference-assembly coordinate system (Hammer et al. 1989; Herrmann et al.
1986; Hammer et al. 1991; Willison & Lyon 2000).

The newer long-read study by Swanepoel et al. (2026, *Nature Communications*
17:21, DOI 10.1038/s41467-025-66616-9) compares a haplotype-resolved
`t^w5` assembly with GRCm38/mm10 and experimentally validates five inversion
junctions, including a newly recognized overlapping inversion. Its methods
and article text establish that numerical breakpoint coordinates exist in the
supplementary table, but that linked table could not be retrieved in this
environment. The newer study therefore does not supply numerical coordinates
to this repository, and its five-inversion nomenclature is not silently
collapsed onto the four historical intervals.

## Assembly compatibility and analysis gate

The 2017 Kelemen & Vicoso archive stores 5-kb tree filenames as offsets into
the published chr17 interval; the local audit applies the documented +5-Mb
offset and recovers non-overlapping windows from chr17:5,000,000–37,054,999.
The archive does not record an assembly accession in the recovered metadata,
and the local figure has no independent assembly manifest. The intervals in
the TSV are therefore labeled `GRCm38/mm10` only as the reference frame of
the published 2018 schematic, not as a validated lift-over into the archive.

Because a chain-file or sequence-based coordinate conversion was not
established, `comparison_eligible=false` for every interval. No inside/outside
quartet-support or recombination-history comparison is reported. The panel is
descriptive and should not be read as exact breakpoint evidence.

## Sources

- Kelemen RK & Vicoso B. 2018. *Complex history and differentiation patterns
  of the t-haplotype, a mouse meiotic driver*. Genetics 208:365–375.
  https://doi.org/10.1534/genetics.117.300513
- Hammer MF, Schimenti J, Silver LM. 1989. *Evolution of mouse chromosome 17
  and the origin of inversions associated with t haplotypes*. PNAS 86:3261–3265.
  https://doi.org/10.1073/pnas.86.9.3261
- Herrmann BG, Barlow DP, Lehrach H. 1987. *A large inverted duplication
  allows homologous recombination between chromosomes heterozygous for the
  proximal t complex inversion*. Cell 48:813–825.
- Hammer MF, Bliss S, Silver LM. 1991. *Genetic exchange across a paracentric
  inversion of the mouse t complex*. Genetics 128:799–812.
  https://doi.org/10.1093/genetics/128.4.799
- Swanepoel CM et al. 2026. *Acquisition of ampliconic sequences marks a
  selfish mouse t-haplotype*. Nature Communications 17:21.
  https://doi.org/10.1038/s41467-025-66616-9
