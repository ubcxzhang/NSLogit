#!/bin/bash
read -p "Enter your working directory (must not end with /): " workdir
cd -- "$workdir/MasterProject/"

nm="GSE60749-GPL13112sim123"
jn="60749"

for mt in "NsLogit" "LR"
  do
  for filt in "FALSE" "TRUE"
    do
    sbatch --job-name=$mt.$jn.90.1.$filt ./computecanada/experiment.sh $nm $mt "90" "1" $filt $workdir

    for sz in 48 24 12
      do 
      for i in 1 2 3 4 5
        do 
        sbatch --job-name=$mt.$jn.$sz.$i.$filt ./computecanada/experiment.sh $nm $mt $sz $i $filt $workdir
      done
    done

  done
done

