#!/bin/bash
#SBATCH --job-name=fix_06_09
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=cluster/logs/fix_06_09_%j.out
#SBATCH --error=cluster/logs/fix_06_09_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

echo "=== Rerunning batch06 (DMR volcano fix) ==="
Rscript methylation_pipeline/batch06/code/06_differential_methylation.R
echo "=== batch06 exit code: $? ==="

echo "=== Rerunning batch09 ==="
Rscript methylation_pipeline/batch09/code/09_tf_motif_integration.R
echo "=== batch09 exit code: $? ==="

echo "=== Regenerating HTML report ==="
Rscript generate_report.R
echo "=== Done ==="
