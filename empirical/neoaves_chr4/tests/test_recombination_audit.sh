#!/bin/sh
set -eu

ROOT=$(CDPATH= cd -- "$(dirname "$0")/../../.." && pwd)
AUDIT="$ROOT/empirical/neoaves_chr4/results/recombination_source_audit.tsv"
REPORT="$ROOT/empirical/neoaves_chr4/results/recombination_source_audit.md"
STRUCT="$ROOT/empirical/neoaves_chr4/results/structural_intervals.tsv"
QUARTET="$ROOT/empirical/neoaves_chr4/results/segment_quartet_support.tsv"

test "$(awk -F '\t' 'NR==1 {print NF}' "$AUDIT")" -eq 14
test "$(awk -F '\t' 'NR>1 {n++} END {print n+0}' "$AUDIT")" -ge 5
awk -F '\t' 'NR>1 {if ($1=="" || $2=="" || $8=="" || $14=="") exit 1}' "$AUDIT"
awk -F '\t' 'NR>1 && $1=="neoaves_pnas2024" {found=1; if ($8 !~ /QQS|not cM\/Mb/) exit 1} END {if (!found) exit 1}' "$AUDIT"
awk 'BEGIN {a=195.2/94.2; b=237.6/94.2; c=183.7/94.2; if (a<2.07 || a>2.08 || b<2.52 || b>2.53 || c<1.95 || c>1.96) exit 1}'
awk -F '\t' 'NR>1 {if ($3!="chr4" || $4<0 || $5<=$4) exit 1}' "$STRUCT"
awk -F '\t' 'NR>1 {if ($5<0 || $6<=$5) exit 1}' "$QUARTET"
test -s "$REPORT"
test ! -e "$ROOT/empirical/neoaves_chr4/figures/neoaves_chr4_recombination_rates.png"
test ! -e "$ROOT/empirical/neoaves_chr4/figures/neoaves_chr4_recombination_rates.pdf"
printf '%s\n' 'recombination audit tests passed'
