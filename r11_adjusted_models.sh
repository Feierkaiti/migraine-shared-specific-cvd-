#!/bin/bash
# Further model-level analyses: cardiometabolic adjustment of the outcomes, a shared component defined by chronic pain only,
# and a three-component model. Usage: bash r11_adjusted_models.sh <root>
root="$1"; R="Rscript $root/src/r03_models.R $root $root/results/ldsc_full.rds"
PRIM=CAD,IS,HF,AF,SBP
SH=MDD25,BACK,FIB
$R adj_bmi MIG $SH $PRIM BMI
$R adj_bmi_t2d MIG $SH $PRIM BMI,UKB_T2D
$R adj_all MIG $SH $PRIM EDU,SMK,BMI,UKB_T2D
$R sh_pain MIG BACK,FIB $PRIM
Rscript "$root/src/r12_dep_pain_model.R" "$root" CAD,IS,HF,AF,SBP,UKB_T2D,BMI
