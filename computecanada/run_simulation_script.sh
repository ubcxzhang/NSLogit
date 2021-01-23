#!/bin/bash
read -p "Enter your working directory (must not end with /): " workdir
cd -- "$workdir/MasterProject/"

for data in "GSE74596" "GSE45719" "GSE60749-GPL13112"
do
  for se in 123 456
  do 
    for ns in 50 100 500 1000 2500 1500 2000 3000 3500
    do 
      sbatch --job-name=$ns.$se.$data ./computecanada/simulation_script.sh $data $ns $se $workdir
    done
  done
done

