#!/bin/bash
# Reverse-direction Mendelian randomization: cardiovascular outcomes as exposures, migraine as outcome.
# Usage: bash r13_reverse_mr.sh <root>
root="$1"; cd "$root"; PY="python"
mkdir -p results/mr/reverse
for t in CAD IS HF AF; do
  $PY src/p07_clump.py "$root" data/std/$t.tsv.gz results/mr/reverse/inst_$t.tsv 5e-8 0.001 10000 > logs/p07_rev_$t.log 2>&1
  $PY src/p08_mr.py "$root" $t results/mr/reverse/inst_$t.tsv results/mr/reverse/mr_${t}_to_MIG.csv MIG,FG_MIG,MVP_MIG > logs/p08_rev_$t.log 2>&1
  cat logs/p08_rev_$t.log
done
