#!/bin/bash
read -p "Enter your working directory (must not end with /): " workdir
cd -- "$workdir/MasterProject/"

nm="GSE74596sim123"
jn="74596"

for mt in "NsLogit" "LR"
  do
  for filt in "FALSE" "TRUE"
    do
    sbatch --job-name=$mt.$jn.40.1.$filt ./computecanada/experiment.sh $nm $mt "44" "1" $filt $workdir

    for sz in 22 12 6
      do 
      for i in 1 2 3 4 5
        do 
        sbatch --job-name=$mt.$jn.$sz.$i.$filt ./computecanada/experiment.sh $nm $mt $sz $i $filt $workdir
      done
    done

  done
done

