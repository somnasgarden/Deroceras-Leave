#!/usr/bin/env bash
#SBATCH --job-name=dmltest
#SBATCH --output=cluster/logs/slurm_%j.out
#SBATCH --error=cluster/logs/slurm_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

set -euo pipefail

# Full-genome DMLtest with smoothing on 25M+ CpGs
# Run from repo root: sbatch cluster/scripts/dmltest_full.slurm

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

mkdir -p cluster/logs cluster/results genome/cache

echo "========================================"
echo "Full-Genome DMLtest"
echo "Job ID:   $SLURM_JOB_ID"
echo "Node:     $(hostname)"
echo "RAM:      $(free -h | grep Mem | awk '{print $2}')"
echo "Project:  $PROJECT_DIR"
echo "Started:  $(date)"
echo "========================================"

Rscript cluster/scripts/dmltest_full.R

echo "========================================"
echo "Finished: $(date)"
echo "Caches:"
ls -lh genome/cache/bsseq_tutorial.rds genome/cache/dmltest_full.rds 2>/dev/null
echo "Results:"
ls -lh cluster/results/dmps_strict.tsv cluster/results/dmrs.tsv 2>/dev/null
echo "========================================"
