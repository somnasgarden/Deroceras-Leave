#!/bin/bash
#SBATCH --job-name=rerun4567
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=cluster/logs/rerun4567_%j.out
#SBATCH --error=cluster/logs/rerun4567_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

echo "=== Rerun batches 04, 05, 06, 07 + regenerate HTML report ==="
echo "Start: $(date)"
echo ""

for b in batch04 batch05 batch06 batch07; do
  script=$(ls methylation_pipeline/${b}/code/*.R 2>/dev/null | head -1)
  echo ""
  echo "==================== ${b}: $(basename $script) ===================="
  echo "Start: $(date)"
  Rscript "$script" || {
    echo "ERROR: ${b} failed. Continuing to next batch."
    continue
  }
  echo "${b} done: $(date)"
done

echo ""
echo "==================== generate_report.R ===================="
echo "Start: $(date)"
Rscript generate_report.R
echo "Report regenerated: $(date)"

echo ""
echo "=== rerun4567 complete: $(date) ==="
