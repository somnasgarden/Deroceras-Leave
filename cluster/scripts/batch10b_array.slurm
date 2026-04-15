#!/bin/bash
#SBATCH --job-name=nme_array
#SBATCH --partition=defq
#SBATCH --array=1-31
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=cluster/logs/nme_%A_%a.out
#SBATCH --error=cluster/logs/nme_%A_%a.err

# Each array task processes ONE chromosome across all 4 samples
# After all 31 finish, run batch10b_merge.R to combine + produce figures

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

CHR_NUM=$SLURM_ARRAY_TASK_ID
echo "Processing chromosome: chr${CHR_NUM}"

Rscript cluster/scripts/batch10b_per_chr.R $CHR_NUM
