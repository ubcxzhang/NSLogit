#!/bin/bash
read -p "Enter your working directory (must not end with /): " workdir
cd -- "$workdir/MasterProject/"

sbatch --job-name="preplotsim123" ./computecanada/pre_plot_sim123.sh $workdir

