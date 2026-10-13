# Stage 1B reproduction report: shared_CM027530.1

Stage 1B did not reach raw-read processing, ANGSD SAF generation, winsfs pairwise SFS estimation, dXY calculation, or NJ plotting in this environment.

The frozen Stage-1A design validated: 22 representatives, 11 species, AA = 11, BB = 11, one AA and one BB representative per species, and the frozen interval `CM027530.1:574805-874807` on `GCA_001746935.2`.

Execution is blocked because the required command-line tools were absent from PATH: prefetch, fasterq-dump, bwa, samtools, angsd, winsfs. `Rscript` is present, but the upstream raw-read/SFS stack is not.

No samples were substituted, no focal interval was changed, and no background/control regions were analyzed.

Because no complete dXY matrix or NJ tree exists, the published AA/BB-across-species signal has not yet been reproduced. The exact topology cannot be compared, and no anomalous species/sample can be identified from Stage 1B output.

The frozen sample accession field contains SAMN BioSample identifiers rather than run-level SRR accessions. Before downloading reads, the next executable run must resolve each SAMN accession to its public run accession(s) and record that mapping.
