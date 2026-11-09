#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../../../../../" && pwd)"
out="$repo_root/empirical/neoaves_chr4/external/stiller2024/63430.alns.tar.gz"
url='https://sid.erda.dk/share_redirect/ENhZODU9YE/01_alignments_and_gene_trees/intergenic_regions/63430.alns.tar.gz'
expected='87a2a643083bb567022683d99b4f5446cf9d52b1da0b0d694bd52aca55f37e43'
if [[ ! -s "$out" ]]; then curl -L --fail --retry 2 --retry-delay 3 "$url" -o "$out"; fi
actual="$(shasum -a 256 "$out" | awk '{print $1}')"
[[ "$actual" == "$expected" ]] || { echo "checksum mismatch: $actual" >&2; exit 1; }
echo "$out $actual"
