# Stage 11B repair report

The previous committed manifest was empty (header only), while its raw output contained 6,096 states/posteriors. That positional mismatch made the zero Phlag interval invalid. Stage 11B archives those files under `provenance_invalid_stage11/`, rewrites the input and manifest together, and requires exact count equality before parsing Phlag.

The regenerated input contains **5619** trees and the manifest contains **5619** rows; 477 loci were excluded after exact four-taxon pruning. Official Phlag v0.2.0 was run without modifying its algorithm. The output parser requires exactly one state sequence and one posterior sequence, each equal in length to the manifest.

The four-taxon detector is a separately labelled categorical adaptation because one local tree yields one resolved topology indicator, whereas the original multi-taxon QQS values are fractional summaries. The old 6,096-locus detector remains an unmatched exploratory comparison.
