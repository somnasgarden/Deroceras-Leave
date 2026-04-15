#!/bin/bash
#SBATCH --job-name=fix_05_09
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=cluster/logs/fix_05_09_%j.out
#SBATCH --error=cluster/logs/fix_05_09_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

echo "=== Rerunning batch05 (TE plots fix) ==="
Rscript methylation_pipeline/batch05/code/05_te_methylation_age.R 2>&1
echo "=== batch05 exit code: $? ==="

echo "=== Rerunning batch09 (volcano fix) ==="
Rscript methylation_pipeline/batch09/code/09_tf_motif_integration.R 2>&1
echo "=== batch09 exit code: $? ==="

echo "=== Regenerating HTML report ==="
Rscript generate_report.R 2>&1
echo "=== Done ==="
