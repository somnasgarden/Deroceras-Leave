#!/usr/bin/env bash
#SBATCH --job-name=motif_annot
#SBATCH --output=cluster/logs/slurm_%j.out
#SBATCH --error=cluster/logs/slurm_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

# ============================================================
# Expanded Motif Annotation — JASPAR + HOMER
# Run from repo root: sbatch cluster/scripts/batch1.5_expanded_motif.slurm
# ============================================================

# Find repo root
if [ -f "${SLURM_SUBMIT_DIR}/methylation_pipeline/_config.R" ]; then
  PROJECT_DIR="$SLURM_SUBMIT_DIR"
elif [ -f "${SLURM_SUBMIT_DIR}/../methylation_pipeline/_config.R" ]; then
  PROJECT_DIR="$(cd "$SLURM_SUBMIT_DIR/.." && pwd)"
elif [ -f "${SLURM_SUBMIT_DIR}/../../methylation_pipeline/_config.R" ]; then
  PROJECT_DIR="$(cd "$SLURM_SUBMIT_DIR/../.." && pwd)"
else
  echo "ERROR: Cannot find _config.R. Submit from repo root."
  exit 1
fi
cd "$PROJECT_DIR"
export PROJECT_DIR
export XDG_CACHE_HOME="${PROJECT_DIR}/cluster/.cache"
mkdir -p "$XDG_CACHE_HOME"
mkdir -p cluster/logs cluster/results genome/cache

echo "========================================"
echo "Expanded Motif Annotation"
echo "Job ID:   $SLURM_JOB_ID"
echo "Node:     $(hostname)"
echo "RAM:      $(free -h | grep Mem | awk '{print $2}')"
echo "CPUs:     $SLURM_CPUS_PER_TASK"
echo "Project:  $PROJECT_DIR"
echo "Started:  $(date)"
echo "========================================"

# Step 1: Install HOMER (if not already done)
echo "=== Step 1: HOMER Install ==="
bash cluster/scripts/setup_homer.sh
export PATH="${PROJECT_DIR}/cluster/homer/bin:${PATH}"

# Step 2: Run expanded motif annotation (R)
echo "=== Step 2: Expanded Motif Annotation (R) ==="
Rscript cluster/scripts/batch1.5_expanded_motif.R

echo "========================================"
echo "=== Job Complete ==="
echo "Finished: $(date)"
echo "Output files:"
ls -lh cluster/results/motif_hits_extended.tsv.gz 2>/dev/null
ls -lh cluster/results/motif_summary_by_region.tsv 2>/dev/null
ls -lh cluster/results/homer_denovo/ 2>/dev/null
echo "========================================"
