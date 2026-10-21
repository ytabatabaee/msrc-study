# Stage 11 original multi-taxon detector analysis

This directory preserves the detector artifacts committed by `3b71262`
(`phlag on neoaves`). They use the original Columbea locus-level D2/D3
representation with 6,096 source loci, and the topology-plus-B configuration
uses the additional `B_observed` feature. These analyses are scientifically
distinct from the official Phlag benchmark and from the matched four-taxon
detector, which use the validated 5,619-entry four-taxon manifest.

The copied `inconsistent_5619_entry_manifest.tsv`,
`inconsistent_phlag_run_manifest.tsv`, and
`inconsistent_stage11_manifest.json` are retained as provenance only: their
5,619-entry/Phlag metadata do not establish a coordinate manifest for the
6,096-locus detector outputs. They must not be used to evaluate the
6,096-locus predictions.
