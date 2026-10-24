# Patterson's D methods

Patterson's D is calculated directly from unphased diploid nucleotide genotypes using population allele frequencies and a fixed-outgroup ancestral allele. ABBA and BABA are summed across sites before taking the ratio. Sites with multiallelic calls, unsuitable missingness, insufficient group callability, missing or polymorphic outgroups, and sensitivity MAC failures are retained in counts and excluded from the estimator. The primary filter uses no MAC threshold; MAC >= 2 is a sensitivity analysis matching Stage 2R.
