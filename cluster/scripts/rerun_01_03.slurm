#!/bin/bash
#SBATCH --job-name=fix_01_03
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=cluster/logs/fix_01_03_%j.out
#SBATCH --error=cluster/logs/fix_01_03_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

echo "=== Rerunning batch01 ==="
Rscript methylation_pipeline/batch01/code/01_genomic_cpg_annotation.R
echo "=== batch01 exit code: $? ==="

echo "=== Rerunning batch03 ==="
Rscript methylation_pipeline/batch03/code/03_promoter_cpg_classification.R
echo "=== batch03 exit code: $? ==="

echo "=== Regenerating HTML report ==="
Rscript generate_report.R
echo "=== Done ==="
