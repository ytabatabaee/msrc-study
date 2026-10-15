# Recombination-suppression source audit

## Result

Direct quantitative evidence was recovered, but it is coarse and arrangement-specific rather than coordinate-resolved. Silver & Artzt (1981; [10.1038/290068a0](https://doi.org/10.1038/290068a0)) report a conventional T–tf map distance of 7–12 cM and 0.1–0.5 cM in naturally occurring t/+ heterozygotes. Hammer, Bliss & Silver (1991; [10.1093/genetics/128.4.799](https://doi.org/10.1093/genetics/128.4.799)) describe the four inversions as spanning about 20 cM and summarize 100–200-fold suppression. These are genetic-map or region-wide statements, not rates assigned to the 5-kb Kelemen & Vicoso windows.

The Silver & Artzt comparison is arrangement-specific: the low value is measured in naturally occurring t/+ heterozygotes and the 7–12 cM comparison is the conventional map. It is not a matched, modern, coordinate-resolved +/+ versus +/t map with confidence intervals. The accessible record does not provide progeny denominators or interval-specific confidence limits for the quoted range. The Hammer et al. 100–200-fold statement is a published region-wide estimate, not a continuous map.

The full source inventory, design, units, uncertainty, arrangement status, and coordinate compatibility are in recombination_suppression_source_audit.tsv.

## Coordinate validation

The frozen 2017 archive uses GRCm38.p4/mm10. Window names are offsets into the published chromosome-17 interval, and the local audit recovers chr17 = archive coordinate + 5,000,000, yielding retained windows from approximately chr17:5,000,000–37,054,999. GRCm38 and mm10 refer to the same major assembly. No lift-over is required for the archived windows or sequence-validated wild-type inversion intervals.

The four historical inversion intervals are In(17)1 5,747,331–6,957,919; In(17)2 7,435,308–13,016,510; In(17)3 14,068,229–18,051,507; and In(17)4 19,401,000–40,049,964. These are the Swanepoel et al. (2026) GRCm38/mm10 wild-type coordinates already recorded in t_complex_inversion_intervals.tsv. Junction regions were PCR/Sanger validated, but the older literature provides only historical/cytogenetic placements. The novel overlapping Inv2 is documented separately.

The T and tf marker identities can be anchored approximately to modern mm10 orthologs (T at chr17:8,434,423–8,442,496 and tf/Itpr3 at chr17:27,275,976–27,341,197). This does not convert the old genetic map into a coordinate-resolved recombination profile: the original marker assays predate mm10 and the published cM values are for the marker pair.

## Figure and interpretation

figures/house_mouse_t_complex_suppression.{png,pdf} is a schematic four-panel figure. Panel A shows the four published inversion intervals on the shared chr17 axis. Panel B shows the T–tf marker-pair observations as interval-level cM ranges, modern marker anchors, and the region-wide 20-cM / 100–200-fold summary without turning it into a window track. Panel C reuses the frozen all-tip quartet support profiles. Panel D shows the frozen D(w)=q_t_alt(w)-q_species(w) contrast. Quartet panels retain raw 5-kb points and existing 250-kb summaries; no direct recombination rate is imputed to any window.

The independent suppression observations have one marker-pair resolution and region-wide summaries, lack coordinate-resolved uncertainty, and are not measurements on the Kelemen & Vicoso sample backgrounds. Therefore no formal spatial association test was run. The existing inversion inside/outside summaries and phylogeny-based classification remain descriptive and are not direct suppression rates. A coordinate-resolved analysis would require matched +/+ and +/t pedigrees or linkage maps with marker positions and uncertainty on the same reference assembly.

## Preservation

The frozen quartet scan, local trees, ASTRAL results, biological sample classifications, and source archives were not modified. The historical classification figure remains a separate supplementary result.
