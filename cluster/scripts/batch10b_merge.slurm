#!/bin/bash
#SBATCH --job-name=nme_merge
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=cluster/logs/nme_merge_%j.out
#SBATCH --error=cluster/logs/nme_merge_%j.err

PROJECT_DIR="${PROJECT_DIR:-/mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave}"
cd "$PROJECT_DIR"
Rscript cluster/scripts/batch10b_merge.R
