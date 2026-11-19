# Anopheles 2La f_dM input status

**Status: unavailable; Phase 1 remains open.** The local checkout contains the
existing Patterson's D summaries and prior validation products, but no
individual-level Ag3.10 genotype array/VCF that can be used to reconstruct
sitewise allele frequencies. The available D results are aggregated; f_dM
cannot be recovered from ABBA/BABA totals. The raw-data directory
`data/anopheles_2la/raw/malariagen_ag3_10/` has no local genotype-call files.

A previous metadata access attempt through `malariagen_data` failed because
DNS resolution for the Google Cloud Storage and OAuth endpoints was
unavailable. Therefore no new genotype download, f_dM scan, or f_dM result was
produced. The supported next step is to rerun the frozen Stage 1A design against
the same Ag3.10 genotype release once that source is reachable, then validate
sitewise D against the existing scan before publishing f_dM. No zero, missing,
or proxy estimate should be interpreted as a biological result.
