#!/bin/bash
#SBATCH --job-name=fix_04_09
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=cluster/logs/fix_04_09_%j.out
#SBATCH --error=cluster/logs/fix_04_09_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

echo "=== Rerunning batch04 ==="
Rscript methylation_pipeline/batch04/code/04_baseline_methylation_landscape.R 2>&1
echo "=== batch04 exit code: $? ==="

echo "=== Rerunning batch09 ==="
Rscript methylation_pipeline/batch09/code/09_tf_motif_integration.R 2>&1
echo "=== batch09 exit code: $? ==="

echo "=== Generating HTML report ==="
Rscript generate_report.R 2>&1
echo "=== Done ==="
