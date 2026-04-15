#!/bin/bash
#SBATCH --job-name=annotate_wg
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=cluster/logs/annotate_wg_%j.out
#SBATCH --error=cluster/logs/annotate_wg_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave
Rscript cluster/scripts/annotate_wg_motifs.R
