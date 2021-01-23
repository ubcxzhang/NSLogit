#!/bin/bash
read -p "Enter your working directory (must not end with /): " workdir
cd -- "$workdir/MasterProject/"

module load r/3.5.0 nixpkgs/16.09  gcc/5.4.0 intel/2016.4

for file in "AT1" "AT2" "B_Cells" "Basal" "cDCs" "Ciliated" "Differentiating_Ciliated" "Endothelial_Cells" "Fibroblasts" "HAS1_High_Fibroblasts" "KRT5-_KRT17+" "Lymphatic_Endothelial_Cells" "Macrophages" "Mast_Cells" "Mesothelial_Cells" "Monocytes" "MUC5AC+_High" "MUC5B+" "Myofibroblasts" "NK_Cells" "pDCs" "Plasma_Cells" "PLIN2+_Fibroblasts" "Proliferating_Epithelial_Cells" "Proliferating_Macrophages" "Proliferating_T_Cells" "SCGB3A2+" "SCGB3A2+_SCGB1A1+" "Smooth_Muscle_Cells" "T_Cells" "Transitional_AT2"
  do
  for mt in "DEsingle" "edgeRLRT" "edgeRLRTcensus" "edgeRLRTdeconv" "edgeRLRTrobust" "edgeRQLF" "edgeRQLFDetRate" "limmatrend" "MASTcpm" "MASTcpmDetRate" "MASTtpmDetRate" "MASTtpm" "metagenomeSeq" "monocle" "monoclecensus" "monoclecount" "ROTScpm" "ROTStpm" "ROTSvoom" "SAMseq" "scDD" "SeuratBimod" "SeuratBimodIsExpr2" "SeuratBimodnofilt" "SeuratTobit" "SeuratTobitnofilt" "voomlimma" "Wilcoxon" "ttest" "BPSC" "LR" "NsLogit"
    do
    sbatch --job-name="real_data".$mt ./computecanada/GSE135893_DE.sh $file $mt $workdir
  done
done

