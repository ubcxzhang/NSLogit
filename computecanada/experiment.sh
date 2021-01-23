#!/bin/bash
#SBATCH -J BV1
#SBATCH --account=rrg-ubcxzh
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -t 0-00:30

module load r/3.5.0 nixpkgs/16.09  gcc/5.4.0 intel/2016.4
Rscript ./conquer_framework/experiment.R $1 $2 $3 $4 $5 $6

