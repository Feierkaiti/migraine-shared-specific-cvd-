#!/bin/bash
# MAGMA v1.10 gene, tissue (GTEx v8) and drug-target gene-set analyses for one GWAS.
# Usage: bash m01_magma.sh <root> <name> <pval_file with SNP P N columns>
set -euo pipefail
root="$1"; name="$2"; pv="$3"
MAGMA="$root/data/ref/magma/bin/magma"
REF="$root/data/ref/g1000/g1000_eur"
A="$root/data/ref/annot"
OUT="$root/results/magma"; mkdir -p "$OUT"
if [ ! -f "$OUT/g1000.genes.annot" ]; then
  # gene windows: 35 kb upstream, 10 kb downstream (as in FUMA defaults)
  "$MAGMA" --annotate window=35,10 --snp-loc "$REF.bim" --gene-loc "$root/data/ref/genes/NCBI37.3.gene.loc" --out "$OUT/g1000"
fi
"$MAGMA" --bfile "$REF" --pval "$pv" use=SNP,P ncol=N --gene-annot "$OUT/g1000.genes.annot" --out "$OUT/$name"
"$MAGMA" --gene-results "$OUT/$name.genes.raw" --gene-covar "$A/gtex_v8_ts_avg_log2TPM_54.txt" \
  --model condition-hide=Average direction=greater --out "$OUT/${name}_gtex"
"$MAGMA" --gene-results "$OUT/$name.genes.raw" --set-annot "$A/genesets.txt" --out "$OUT/${name}_drugsets"
# migraine drug targets split into those shared with depression drugs and the rest (set file from m02_split_drugsets.py)
if [ -f "$A/genesets_split.txt" ]; then
  "$MAGMA" --gene-results "$OUT/$name.genes.raw" --set-annot "$A/genesets_split.txt" --out "$OUT/${name}_drugsplit"
fi
