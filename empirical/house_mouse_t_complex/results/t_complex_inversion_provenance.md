# t-complex inversion interval provenance

## Result

The four historical t-complex inversions are shown in the updated spatial
figure using the validated mm10 wild-type coordinates from Swanepoel et al.
(2026): `In(17)1` (5,747,331–6,957,919), `In(17)2` (7,435,308–13,016,510),
`In(17)3` (14,068,229–18,051,507), and `In(17)4` (19,401,000–40,049,964).
These coordinates are recorded in
`empirical/house_mouse_t_complex/results/t_complex_inversion_intervals.tsv`.
The boundary evidence is sequence-based synteny with PCR/Sanger validation;
the supplementary table also reports the PCR junction-region widths.

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
`t^w5` assembly with GRCm38/mm10 and experimentally validates five inversions.
Supplementary Table 1 supplies the wild-type mm10 coordinates and Table 2
supplies the PCR junction regions. Its Inv1, Inv3, Inv4, and Inv5 correspond
to the four historical inversions used here; its novel overlapping Inv2 is
recorded below but is not silently added to the historical four-track panel.

The novel fifth inversion is Inv2, with mm10 wild-type coordinates
6,957,919–7,435,308. It overlaps the rearranged t-haplotype structure in the
authors' comparison and is therefore treated as a separate structural finding,
not forced into the four historical labels.

## Assembly compatibility and analysis gate

The 2017 Kelemen & Vicoso archive stores 5-kb tree filenames as offsets into
the published chr17 interval; the local audit applies the documented +5-Mb
offset and recovers non-overlapping windows from chr17:5,000,000–37,054,999.
The upstream Harr et al. genomic resources state that reads and variants were
mapped to the mm10 reference, and Kelemen & Vicoso explicitly used the
GRCm38.p4 reference for the corresponding analysis. GRCm38 and mm10 are the
same major coordinate assembly name. No chain-file conversion is therefore
needed; the exact mm10 intervals are directly compatible with the archived
window coordinates. Inversion 4 extends beyond the available tree coverage.

The supplementary descriptive comparison is in
`t_complex_inversion_interval_summary.tsv`. It uses only complete 5-kb windows,
excludes windows intersecting each boundary, and reports recombination-class
counts/fractions separately for domesticus, musculus, and castaneus. The
windows remain linked and the table is descriptive rather than an independent-
replicate test.

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
