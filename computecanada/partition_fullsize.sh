#!/bin/bash
#SBATCH -J BV1
#SBATCH --account=rrg-ubcxzh
#SBATCH --cpus-per-task=1
#SBATCH --mem=125G
#SBATCH -t 0-03:00

module load r/3.5.0 nixpkgs/16.09  gcc/5.4.0 intel/2016.4
Rscript --max-ppsize=500000  ./real_data/partition_fullsize.R $1 $2

