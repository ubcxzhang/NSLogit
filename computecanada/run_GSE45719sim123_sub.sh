#!/bin/bash
read -p "Enter your working directory (must not end with /): " workdir
cd -- "$workdir/MasterProject/"

nm="GSE45719sim123"
jn="45719"

for mt in "NsLogit" "LR"
  do
  for filt in "FALSE" "TRUE"
    do
    sbatch --job-name=$mt.$jn.50.1.$filt ./computecanada/experiment.sh $nm $mt "50" "1" $filt $workdir
    for sz in 24 12 6
      do 
      for i in 1 2 3 4 5
        do 
	sbatch --job-name=$mt.$jn.$sz.$i.$filt ./computecanada/experiment.sh $nm $mt $sz $i $filt $workdir
      done
    done
  done
done

