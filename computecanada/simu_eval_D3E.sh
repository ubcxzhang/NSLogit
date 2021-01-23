#!/bin/bash
#SBATCH -J BV1
#SBATCH --account=def-ubcxzh
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH -t 0-15:00

module load r/3.5.0 nixpkgs/16.09  gcc/5.4.0 intel/2016.4

module load python/2.7
virtualenv --no-download ../ENV
source ../ENV/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r ../py_requirements.txt

Rscript --max-ppsize=500000 ./scDD_framework/simu_eval.R $1 $2 $3 $4 $5

