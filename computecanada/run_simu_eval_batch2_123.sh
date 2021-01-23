#!/bin/bash
read -p "Enter your working directory (must not end with /): " workdir
cd -- "$workdir/MasterProject/"

module load r/3.5.0 nixpkgs/16.09  gcc/5.4.0 intel/2016.4

for mt in "edgeRLRTrobust" "metagenomeSeq" "SeuratBimodnofilt" "SeuratTobit" "ROTScpm" "ROTStpm"
do
  for nm in "GSE74596" "GSE60749-GPL13112" "GSE45719"
  do
    for se in 123
    do 
      for sz in 50 100 500
      do 
        sbatch --job-name=$sz.$se.$mt ./computecanada/simu_eval.sh $se $sz $mt $nm $workdir
      done
    done
  done
done
