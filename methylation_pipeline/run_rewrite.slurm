#!/usr/bin/env bash
#SBATCH --job-name=rewrite4
#SBATCH --output=cluster/logs/rewrite_%j.out
#SBATCH --error=cluster/logs/rewrite_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=defq

set -euo pipefail

# =============================================================================
# Run rewritten batches: 04 → 05 → 08 → 07
# Dependency order: 04 (BSseq) → 05 (uses BSseq) → 08 (WGCNA) → 07 (uses 08)
# =============================================================================

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
echo "Rewritten Batches: 04 → 05 → 08 → 07"
echo "Job ID:       $SLURM_JOB_ID"
echo "Node:         $(hostname)"
echo "RAM:          $(free -h | grep Mem | awk '{print $2}')"
echo "Project:      $PROJECT_DIR"
echo "Started:      $(date)"
echo "========================================"

for batch in 04 05 08 07; do
  script="methylation_pipeline/batch${batch}/code/$(ls methylation_pipeline/batch${batch}/code/*.R | head -1 | xargs basename)"
  echo ""
  echo "========================================"
  echo "=== BATCH ${batch}: ${script} ==="
  echo "=== Started: $(date) ==="
  echo "========================================"
  Rscript "$script" || {
    echo "ERROR: Batch ${batch} failed at $(date)"
    continue
  }
  echo "=== Batch ${batch} completed: $(date) ==="
done

echo ""
echo "========================================"
echo "All rewritten batches complete: $(date)"
echo "========================================"
