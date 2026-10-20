# Phlag source audit

The primary software is the upstream Phlag v0.2.0 checkout at `/private/tmp/phlag-upstream`, commit `76597d0e4ebd0e507fdafebf7f5163ef634a36d6`. It accepts a Newick species tree with labeled internal focal nodes and an ordered one-Newick-per-line local gene-tree file. It computes or reads quadripartition quartet scores and fits its HMM. The bundled example was retained in the earlier cross-system audit and is the software validity check.

The local Neoaves pilot uses the original 6,140 chromosome-4 local trees, pruned to one fixed representative from each Columbea role. It is a matched-taxon pilot, not a reproduction of the published all-taxon avian analysis.
