#!/bin/bash
#SBATCH -J BV1
#SBATCH --account=def-ubcxzh
#SBATCH --cpus-per-task=18
#SBATCH --mem=150G
#SBATCH -t 0-50:00

module load r/3.5.0 nixpkgs/16.09  gcc/5.4.0 intel/2016.4
Rscript --max-ppsize=500000 ./scDD_framework/simu_eval.R $1 $2 $3 $4 $5

