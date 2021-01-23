#!/bin/bash
read -p "Enter your working directory (must not end with /): " workdir
cd -- "$workdir/MasterProject/conquer_framework"

module load r/3.5.0 nixpkgs/16.09  gcc/5.4.0 intel/2016.4

Rscript merge_sub_result.R "GSE45719sim123" "NsLogit" "FALSE" $workdir
Rscript merge_sub_result.R "GSE60749-GPL13112sim123" "NsLogit" "FALSE" $workdir
Rscript merge_sub_result.R "GSE74596sim123" "NsLogit" "FALSE" $workdir

Rscript merge_sub_result.R "GSE45719sim123" "NsLogit" "TRUE" $workdir
Rscript merge_sub_result.R "GSE60749-GPL13112sim123" "NsLogit" "TRUE" $workdir
Rscript merge_sub_result.R "GSE74596sim123" "NsLogit" "TRUE" $workdir

Rscript merge_sub_result.R "GSE45719sim123" "LR" "FALSE" $workdir
Rscript merge_sub_result.R "GSE60749-GPL13112sim123" "LR" "FALSE" $workdir
Rscript merge_sub_result.R "GSE74596sim123" "LR" "FALSE" $workdir

Rscript merge_sub_result.R "GSE45719sim123" "LR" "TRUE" $workdir
Rscript merge_sub_result.R "GSE60749-GPL13112sim123" "LR" "TRUE" $workdir
Rscript merge_sub_result.R "GSE74596sim123" "LR" "TRUE" $workdir

