# Nucleotide D and exact quartet genealogy

The exact topology convention is q1=(P1,P2)|(P3,O), q2=(P1,P3)|(P2,O),
q3=(P2,P3)|(P1,O), with D_topology=q3−q2. Existing results use all
available frozen group combinations per tree and report frequencies within
each window; combinations within one tree are not independent genomic
replicates. Window IDs match the frozen nucleotide window grid.

Existing whole-chr16 Pearson correlations are 0.838 (Q1) and 0.762 (Q2); within
the 52 supergene windows they are 0.544 and 0.466. At 1-Mb physical-block
bootstrap, the within-supergene intervals are 0.158–0.752 and 0.249–0.668.
After centering both tracks within inside/outside strata, whole-chr16
correlations remain 0.441 and 0.356, showing that the whole-axis association is
not solely the shared regional step. Within-supergene values are descriptive
and based on linked windows.

In(16)1 midpoint-class windows have high pooled nucleotide D for both designs;
In(16)2 has only four windows and In(16)3 none, so the latter two cannot support
equivalent spatial correlation inference. Topology frequencies describe
topological shifts, not branch-length magnitude. `branch_length_sensitivity.tsv`
separately reports group-mean patristic distances and four-point edge contrasts
derived from source Newick branch lengths. These contrasts can be negative under
population mixtures and are not the mean positive internal branch length of
resolved individual quartets. Within the supergene, D_site correlates only
weakly with the group-mean q3−q2 edge contrast (r=0.178 for Q1 and 0.132 for
Q2; on chr16A, r=0.010 and 0.045). A strong topology-frequency shift therefore
does not imply a proportional shift in the group-average branch contrast. This
is descriptive, not a fitted coalescent branch-length model.

The site and tree data are based on related selected BUSCO loci. Concordance is
therefore not independent replication and cannot distinguish introgression,
ancestral structure, persistent arrangements, or MSRC by itself.
