#!/bin/bash
# All model-level analyses from the full LDSC run (results/ldsc_full.rds). Usage: bash r04_run_models.sh <root>
root="$1"; R="Rscript $root/src/r03_models.R $root $root/results/ldsc_full.rds"
PRIM=CAD,IS,HF,AF,SBP
SEC=AS,LAS,CES,SVS,DBP,UKB_T2D,HDL,TG,BMI,UKB_IHD,FG_CHD,FG_STR,FG_T2D,FG_HTN
SH=MDD25,BACK,FIB
$R primary MIG $SH $PRIM,$SEC
$R adj_edu_smk MIG $SH $PRIM EDU,SMK
for m in FG_MIG MVP_MIG FG_MIGTRIP FG_MA FG_MO MVP_MA CHOQ; do
  $R def_$m $m $SH $PRIM,LAS,CES,SVS
done
$R sh_dep MIG MDD25 $PRIM
$R sh_dep_back MIG MDD25,BACK $PRIM
$R sh_wide MIG MDD25,BACK,FIB,FG_ANX,FG_IBS $PRIM
$R sh_depclin MIG MDD25_CLIN,BACK,FIB $PRIM
$R sh_nofib MIG MDD25,BACK,FG_ANX,FG_IBS $PRIM
