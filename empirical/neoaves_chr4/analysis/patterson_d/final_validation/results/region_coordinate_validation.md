# Region and coordinate validation

The frozen region rule is the published inclusive GalGal4 `chr4` **10-kb window-start** predicate: PNAS_1 25,030,000–32,670,000; PNAS_2 33,510,000–34,470,000; PNAS_3 44,130,000–56,810,000. The saved locus result classification exactly matches that predicate on the stored source window-start coordinate for all included designs. The independent pooling audit found no locus-membership or ABBA/BABA sum discrepancies.

The Stage-4D manifest retains source 10-kb window coordinates and 1-kb alignment offsets separately. The structural/rearrangement track scripts describe their physical axis as GalGal6, while the published PNAS membership is defined on GalGal4. The local package documents approximate GalGal6 counterparts from the paper, but the final validation did not recover and apply an independent base-level GalGal4-to-GalGal6 chain to every alignment/locus. Therefore:

- New figures use the Stiller source window-start coordinate and label the PNAS bands as the frozen published membership predicate.
- They do not copy GalGal4 breakpoints onto GalGal6 structural tracks or imply that interval endpoints are structural breakpoints.
- Direct locus-ID matching of nucleotide counts to local trees is valid and does not require coordinate liftover; spatial plotting against external GalGal6 annotations does.
- Boundary sensitivity is computed on the same Stiller source window-start axis used for the frozen predicate and is explicitly not a GalGal6 lift.

The original 10-kb windows and selected 1-kb alignments are locus units, not continuous nucleotide-resolution breakpoint observations. A fully assembly-normalized overlay requires an explicit, verified coordinate crosswalk and remains unresolved.
