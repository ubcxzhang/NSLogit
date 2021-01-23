#!/bin/bash
#SBATCH -J BV1
#SBATCH --account=rrg-ubcxzh
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH -t 0-00:20

module load r/3.5.0 nixpkgs/16.09  gcc/5.4.0 intel/2016.4
Rscript --max-ppsize=500000 ./real_data/real_data_analysis_new_filter.R $1


