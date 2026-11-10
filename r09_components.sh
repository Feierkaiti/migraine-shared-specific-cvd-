#!/bin/bash
# After the SNP-level decomposition: component summary statistics, LDSC with outcomes, loci, instrument decomposition, MR.
# Usage: bash r09_components.sh <root> <stages: stage1 | stage1,stage2>
set -e
root="$1"; st="$2"; cd "$root"; PY="python"
$PY src/c01_components.py "$root" "$st" > logs/c01_components.log 2>&1
$PY src/p04_hm3.py "$root" COMP_S COMP_M > logs/p04_comp.log 2>&1
Rscript src/r01_munge.R "$root" COMP_S COMP_M > logs/r01_munge_comp.log 2>&1
Rscript src/r02_ldsc.R "$root" comp COMP_S,COMP_M,MIG,MDD25,BACK,FIB,CAD,IS,HF,AF,SBP,AS,UKB_T2D,TG,HDL,BMI > logs/r02_ldsc_comp.log 2>&1
Rscript src/r07_rg_table.R "$root" results/ldsc_comp.rds results/ldsc_comp_rg.csv > logs/r07_comp.log 2>&1
mv results/rg_mig_se.csv results/rg_mig_se_comp_tmp.csv && Rscript src/r07_rg_table.R "$root" results/ldsc_full.rds results/rg_full.csv >> logs/r07_comp.log 2>&1 && rm results/rg_mig_se_comp_tmp.csv
mkdir -p results/loci
for c in COMP_S COMP_M; do $PY src/p07_clump.py "$root" results/gwas/primary/${c}_pos.tsv.gz results/loci/${c}_clump.tsv 5e-8 0.1 250 > logs/p07_${c}.log 2>&1; done
[ -f results/loci/MIG_clump.tsv ] || $PY src/p07_clump.py "$root" data/std/MIG.tsv.gz results/loci/MIG_clump.tsv 5e-8 0.1 250 > logs/p07_MIG.log 2>&1
$PY src/c02_loci.py "$root" > logs/c02_loci.log 2>&1
$PY src/p10_mr_components.py "$root" "$st" 0.75 > logs/p10_mr_components.log 2>&1
$PY src/p08_mr.py "$root" MIG results/mr/inst_M.tsv results/mr/mr_MIG_inst_M.csv CAD,IS,HF,AF > logs/p08_inst_M.log 2>&1
$PY src/p08_mr.py "$root" MIG results/mr/inst_mixed.tsv results/mr/mr_MIG_inst_mixed.csv CAD,IS,HF,AF > logs/p08_inst_mixed.log 2>&1
# M-dominant instruments without PHACTR1 (rs9349379), and without PHACTR1 and LRP1 (rs11172113)
grep -v rs9349379 results/mr/inst_M.tsv > results/mr/inst_M_noPHACTR1.tsv
grep -v -E "rs9349379|rs11172113" results/mr/inst_M.tsv > results/mr/inst_M_noPHACTR1_LRP1.tsv
for f in inst_M_noPHACTR1 inst_M_noPHACTR1_LRP1; do $PY src/p08_mr.py "$root" MIG results/mr/$f.tsv results/mr/mr_MIG_$f.csv CAD,IS,HF,AF > logs/p08_$f.log 2>&1; done
echo done
