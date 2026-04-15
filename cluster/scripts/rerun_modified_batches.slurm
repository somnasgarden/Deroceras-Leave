#!/usr/bin/env bash
#SBATCH --job-name=rerun_mods
#SBATCH --output=cluster/logs/rerun_mods_%j.out
#SBATCH --error=cluster/logs/rerun_mods_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --partition=defq

set -euo pipefail
cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

echo "========================================"
echo "Rerun modified batches: 02, 03, 04, 05, 08, 13"
echo "Job ID: $SLURM_JOB_ID | Node: $(hostname)"
echo "Started: $(date)"
echo "========================================"

for batch_script in \
  methylation_pipeline/batch02/code/02_methylation_machinery.R \
  methylation_pipeline/batch03/code/03_promoter_cpg_classification.R \
  methylation_pipeline/batch04/code/04_baseline_methylation_landscape.R \
  methylation_pipeline/batch05/code/05_te_methylation_age.R \
  methylation_pipeline/batch08/code/08_wgcna_methylation.R \
  methylation_pipeline/batch13/code/13_entropy_metagene.R; do

  echo ""
  echo "========================================"
  echo "=== $(basename $batch_script) ==="
  echo "=== Started: $(date) ==="
  echo "========================================"
  Rscript "$batch_script" || {
    echo "ERROR: $(basename $batch_script) failed. Continuing."
    continue
  }
  echo "=== Completed: $(date) ==="
done

echo ""
echo "========================================"
echo "Regenerating HTML report..."
echo "========================================"
Rscript generate_report.R || echo "WARNING: Report generation failed."

echo ""
echo "========================================"
echo "All done: $(date)"
echo "========================================"
