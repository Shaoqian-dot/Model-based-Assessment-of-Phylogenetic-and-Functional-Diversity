#!/bin/bash

#SBATCH --job-name=phyloSim
#SBATCH --array=1-16
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%A_%a.out

PARAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" config/params.csv)

Rscript scripts/main_Test_Katana.R ...