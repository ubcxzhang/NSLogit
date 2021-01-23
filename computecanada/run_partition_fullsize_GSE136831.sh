#!/bin/bash
read -p "Enter your working directory (must not end with /): " workdir
cd -- "$workdir/MasterProject/"

module load r/3.5.0 nixpkgs/16.09  gcc/5.4.0 intel/2016.4

for file in "ATI" "ATII" "B" "Basal" "Ciliated" "Fibroblast" "Macrophage" "Mast" "Mesothelial" "Myofibroblast" "NK" "pDC" "SMC" "T" 
do
	sbatch --job-name=$file ./computecanada/partition_fullsize_GSE136831.sh $file $workdir
done

