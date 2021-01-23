#!/bin/bash
read -p "Enter your working directory (must not end with /): " workdir
cd -- "$workdir/MasterProject/"

module load r/3.5.0 nixpkgs/16.09  gcc/5.4.0 intel/2016.4

for mt in "DEsingle" "edgeRLRT" "edgeRLRTcensus" "edgeRLRTdeconv" "edgeRLRTrobust" "edgeRQLF" "edgeRQLFDetRate" "limmatrend" "MASTcpm" "MASTcpmDetRate" "MASTtpmDetRate" "MASTtpm" "metagenomeSeq" "monocle" "monoclecensus" "monoclecount" "ROTScpm" "ROTStpm" "ROTSvoom" "SAMseq" "scDD" "SeuratBimod" "SeuratBimodIsExpr2" "SeuratBimodnofilt" "SeuratTobit" "SeuratTobitnofilt" "voomlimma" "Wilcoxon" "ttest" "BPSC" "LR" "NsLogit"

do
  for nm in "GSE74596" "GSE60749-GPL13112" "GSE45719"
  do
    for se in 123 456
    do 
      for sz in 50 100 500 1000 1500 2000 2500 3000 3500
      do 
        sbatch --job-name=$sz.$se.$mt ./computecanada/simu_eval.sh $se $sz $mt $nm $workdir
      done
    done
  done
done
