#!/bin/bash
#SBATCH --job-name=batch09
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=cluster/logs/batch09_%j.out
#SBATCH --error=cluster/logs/batch09_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

echo "=== Batch 09 rerun: TF + Motif + GENIE3 Integration ==="
echo "Start: $(date)"
echo ""

Rscript methylation_pipeline/batch09/code/09_tf_motif_integration.R

echo ""
echo "=== Batch 09 complete: $(date) ==="
