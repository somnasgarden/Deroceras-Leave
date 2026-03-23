#!/usr/bin/env bash
#SBATCH --job-name=meth_pipe
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G

set -euo pipefail

# =============================================================================
# Run the full methylation pipeline (batches 01-10)
# From repo root: sbatch methylation_pipeline/run_pipeline.slurm
# Or: cd methylation_pipeline && sbatch run_pipeline.slurm
# =============================================================================

# Find repo root
if [ -f "${SLURM_SUBMIT_DIR}/methylation_pipeline/_config.R" ]; then
  PROJECT_DIR="$SLURM_SUBMIT_DIR"
elif [ -f "${SLURM_SUBMIT_DIR}/_config.R" ]; then
  PROJECT_DIR="$(dirname "$SLURM_SUBMIT_DIR")"
else
  echo "ERROR: Cannot find _config.R. Submit from repo root or methylation_pipeline/"
  exit 1
fi

cd "$PROJECT_DIR"
export PROJECT_DIR

echo "========================================"
echo "Methylation Pipeline"
echo "Job ID:       $SLURM_JOB_ID"
echo "Node:         $(hostname)"
echo "RAM:          $(free -h | grep Mem | awk '{print $2}')"
echo "Project:      $PROJECT_DIR"
echo "Started:      $(date)"
echo "========================================"

# Run each batch in order
for batch in 01 02 03 04 05 06 07 08 09 10; do
  script=$(ls methylation_pipeline/batch${batch}/code/*.R 2>/dev/null | head -1)
  if [ -z "$script" ]; then
    echo "SKIP: No script found for batch${batch}"
    continue
  fi
  echo ""
  echo "========================================"
  echo "=== BATCH ${batch}: $(basename $script) ==="
  echo "=== Started: $(date) ==="
  echo "========================================"
  Rscript "$script" || {
    echo "ERROR: Batch ${batch} failed. Continuing to next batch."
    continue
  }
  echo "=== Batch ${batch} completed: $(date) ==="
done

echo ""
echo "========================================"
echo "Pipeline complete: $(date)"
echo "========================================"
