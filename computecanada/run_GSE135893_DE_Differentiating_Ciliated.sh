#!/bin/bash
read -p "Enter your working directory (must not end with /): " workdir
cd -- "$workdir/MasterProject/"

module load r/3.5.0 nixpkgs/16.09  gcc/5.4.0 intel/2016.4

for mt in "BPSC" "DEsingle" "edgeRLRT" "edgeRLRTdeconv" "edgeRLRTrobust" "edgeRQLF" "edgeRQLFDetRate" "limmatrend" "LR" "MASTcpm" "MASTcpmDetRate" "metagenomeSeq" "monoclecount" "ROTScpm" "ROTSvoom" "SAMseq" "scDD" "SeuratBimod" "SeuratBimodIsExpr2" "SeuratBimodnofilt" "voomlimma"
do
	sbatch --job-name="real_data".$mt ./computecanada/GSE135893_DE.sh "Differentiating_Ciliated" $mt $workdir
done

