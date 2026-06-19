#!/usr/bin/env bash
#SBATCH --job-name=rew_0587
#SBATCH --output=cluster/logs/rewrite2_%j.out
#SBATCH --error=cluster/logs/rewrite2_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=defq

set -euo pipefail

# Batch 04 already done — run 05 → 08 → 07

if [ -f "${SLURM_SUBMIT_DIR}/methylation_pipeline/_config.R" ]; then
  PROJECT_DIR="$SLURM_SUBMIT_DIR"
elif [ -f "${SLURM_SUBMIT_DIR}/_config.R" ]; then
  PROJECT_DIR="$(dirname "$SLURM_SUBMIT_DIR")"
else
  echo "ERROR: Cannot find _config.R"
  exit 1
fi

cd "$PROJECT_DIR"
export PROJECT_DIR
mkdir -p cluster/logs

echo "========================================"
echo "Rewrite #2: 05 → 08 → 07 (batch 04 done)"
echo "Job: $SLURM_JOB_ID | Node: $(hostname) | $(date)"
echo "========================================"

for batch in 05 08 07; do
  script=$(ls methylation_pipeline/batch${batch}/code/*.R | head -1)
  echo ""
  echo "=== BATCH ${batch}: $(date) ==="
  Rscript "$script" || { echo "ERROR: Batch ${batch} failed at $(date)"; continue; }
  echo "=== Batch ${batch} done: $(date) ==="
done

echo ""
echo "=== All done: $(date) ==="
