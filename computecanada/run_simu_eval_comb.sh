#!/bin/bash
read -p "Enter your working directory (must not end with /): " workdir
cd -- "$workdir/MasterProject/scDD_framework"

module load r/3.5.0 nixpkgs/16.09  gcc/5.4.0 intel/2016.4

for mt in "DEsingle" "edgeRLRT" "edgeRLRTcensus" "edgeRLRTdeconv" "edgeRLRTrobust" "edgeRQLF" "edgeRQLFDetRate" "limmatrend" "MASTcpm" "MASTcpmDetRate" "MASTtpmDetRate" "MASTtpm" "metagenomeSeq" "monocle" "monoclecensus" "monoclecount" "ROTScpm" "ROTStpm" "ROTSvoom" "SAMseq" "scDD" "SeuratBimod" "SeuratBimodIsExpr2" "SeuratBimodnofilt" "SeuratTobit" "SeuratTobitnofilt" "voomlimma" "Wilcoxon" "ttest" "BPSC" "LR" "NsLogit"
do
  Rscript simu_eval_comb.R $mt $workdir
done
