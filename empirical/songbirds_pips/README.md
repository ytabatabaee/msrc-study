# Songbirds PIPs empirical analysis

Stage-0 feasibility/provenance audit for Pegan and Winger (2025) North American songbird putative inversion polymorphisms. Stage 0 and Stage 1A are frozen. Stage 1B is a minimal raw-read pilot for the shared `CM027530.1` PIP, but it is currently blocked before read download/alignment because the local environment does not provide the required SRA/alignment/ANGSD/winsfs toolchain.

## Stage 1B status

Focal candidate: `shared_CM027530.1`, interval `CM027530.1:574805-874807` on `GCA_001746935.2`, using the 22 frozen author representatives from `data/songbirds_pips/processed/stage1a_shared_pip_samples.tsv`.

Validated frozen design:

- 22 samples
- 11 species
- AA = 11, BB = 11
- one AA and one BB representative per species
- 231 expected sample pairs
- no candidate-set changes
- no background/control analysis run

Execution blocker:

- missing required tools on PATH: `prefetch`, `fasterq-dump`, `bwa`, `samtools`, `angsd`, `winsfs`
- the frozen accession column contains SAMN BioSample accessions; the executable run must resolve these to run-level SRR accessions before download

No FASTQ, BAM, SAF, SFS, dXY, NJ tree, or NJ figure was produced. Placeholder/status tables under `results/` and `data/songbirds_pips/processed/` are explicitly marked `not_run_preflight_blocked_missing_tools`; they are not biological results.

Key Stage 1B records:

- `empirical/songbirds_pips/results/stage1b_design_validation.json`
- `empirical/songbirds_pips/results/stage1b_manifest.json`
- `empirical/songbirds_pips/results/stage1b_reproduction_report.md`
- `empirical/songbirds_pips/results/stage1b_stage1c_preanalysis_plan.md`
- `data/songbirds_pips/intermediate/stage1b_shared_pip/logs/stage1b_preflight_tool_versions.tsv`
- `data/songbirds_pips/intermediate/stage1b_shared_pip/logs/stage1b_planned_commands.sh`

Stage 1C is not justified until Stage 1B is rerun with the missing toolchain, all 22 representatives are downloaded or otherwise supplied, focal-region QC is usable, ANGSD/winsfs completes for all or nearly all pairs, and a complete or nearly complete dXY matrix reproduces the published AA/BB-across-species signal.
