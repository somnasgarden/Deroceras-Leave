#!/bin/bash
#SBATCH --job-name=batch11
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=cluster/logs/batch11_%j.out
#SBATCH --error=cluster/logs/batch11_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

echo "=== Batch 11: Motif x NME x Expression ==="
echo "Start: $(date)"
echo ""

Rscript methylation_pipeline/batch11/code/11_motif_entropy_expression.R

echo ""
echo "=== Batch 11 complete: $(date) ==="
