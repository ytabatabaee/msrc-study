# Audit of the 317 unresolved collapsed records

The audit compares each collapsed record with its header-derived named locus and indexes all named loci by exact taxon set. Candidate taxon-set matches are mathematical alternatives only; they are not locus assignments.

- Published collapsed records: **63,430**.
- Header-derived taxon-set matches: **63,113**.
- Unresolved records audited: **317**.
- Rows with at least one exact taxon-set candidate: **117**.
- Rows with no exact taxon-set candidate: **200**.
- Maximum one-to-one matching among unresolved rows: **117**.

## Interpretation

No unresolved record is promoted to a verified locus. A record is treatment-equivalent only when every exact taxon-set candidate has the same membership for all eight treatments; otherwise the identity remains unresolved. Candidate matches do not establish ordering, genomic proximity, or shared upstream provenance. Because every candidate is already a proposed header locus for another collapsed row, candidate reassignment would require a complete global permutation, which this audit does not claim.

## Distribution

- Proposed chromosome counts: `{'chr10': 9, 'chr15': 1, 'chr17': 2, 'chr18': 4, 'chr1': 59, 'chr20': 4, 'chr21': 2, 'chr22': 4, 'chr23': 3, 'chr24': 2, 'chr25': 2, 'chr26': 5, 'chr27_JH375221_random': 1, 'chr2': 55, 'chr3': 38, 'chr4': 47, 'chr5': 21, 'chr6': 12, 'chr7': 10, 'chr8': 6, 'chr9': 6, 'chrLGE22C19W28_E50C23': 1, 'chrUn_AADN03019083': 1, 'chrZ': 22}`.
- Symmetric taxon-set difference counts: `{12: 10, 17: 13, 5: 16, 13: 13, 36: 4, 18: 11, 14: 11, 10: 7, 6: 9, 52: 2, 45: 3, 19: 6, 47: 2, 23: 12, 43: 3, 32: 7, 7: 12, 37: 5, 20: 12, 15: 13, 3: 5, 27: 8, 24: 7, 11: 8, 2: 9, 46: 3, 16: 9, 42: 4, 22: 5, 9: 6, 35: 3, 29: 9, 26: 5, 41: 2, 4: 6, 21: 6, 31: 5, 55: 3, 28: 5, 40: 2, 30: 5, 48: 1, 1: 5, 8: 9, 44: 3, 34: 3, 49: 1, 33: 2, 54: 1, 39: 2, 25: 3, 51: 1}`.
- Single-taxon differences: **5**.
- Collapsed records containing a polytomy: **317**.
- Naming-normalization-only differences: **0**.
- Classification counts: `{'source discrepancy': 200, 'unresolved': 91, 'treatment-equivalent mapping': 26}`.

Under the proposed header-derived identities, 311 records fall inside the
PNAS outlier mask and 6 fall outside it. The corresponding proposed exclusion
counts are 315 for structural primary, 316 for structural stringent, and 313
for structural inclusive. These are descriptive header-locus memberships;
they are not verified memberships for the collapsed records.

The source-side producer of the collapsed archive remains absent from the cached repository snapshot. Ordinary branch collapsing cannot explain changed leaf sets, so any taxon-set discrepancy remains a source or preprocessing discrepancy until an upstream file or author metadata resolves it.
