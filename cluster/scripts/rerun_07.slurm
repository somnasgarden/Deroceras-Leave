#!/bin/bash
#SBATCH --job-name=fix_07
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=cluster/logs/fix_07_%j.out
#SBATCH --error=cluster/logs/fix_07_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

echo "=== Rerunning batch07 ==="
Rscript methylation_pipeline/batch07/code/07_methylation_expression_decoupling.R
echo "=== batch07 exit code: $? ==="

echo "=== Regenerating HTML report ==="
Rscript generate_report.R
echo "=== Done ==="
