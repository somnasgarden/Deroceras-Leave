#!/bin/bash
#SBATCH --job-name=batch09
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=cluster/logs/batch09_%j.out
#SBATCH --error=cluster/logs/batch09_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

echo "=== Rerunning batch09 (TF + motif + concordant) ==="
Rscript methylation_pipeline/batch09/code/09_tf_motif_integration.R
echo "=== batch09 exit code: $? ==="

echo "=== Regenerating HTML report ==="
Rscript generate_report.R
echo "=== Done ==="
