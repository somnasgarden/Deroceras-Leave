#!/usr/bin/env bash
#SBATCH --job-name=dmltest
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G

set -euo pipefail

# Full-genome DMLtest with smoothing on 25M+ CpGs
# Run from: cd cluster && sbatch scripts/dmltest_full.slurm

if [ -f "${SLURM_SUBMIT_DIR}/scripts/dmltest_full.R" ]; then
  CLUSTER_ROOT="$SLURM_SUBMIT_DIR"
elif [ -f "${SLURM_SUBMIT_DIR}/dmltest_full.R" ]; then
  CLUSTER_ROOT="$(dirname "$SLURM_SUBMIT_DIR")"
else
  CLUSTER_ROOT="$(dirname "$SLURM_SUBMIT_DIR")"
fi
cd "$CLUSTER_ROOT"
export CLUSTER_ROOT

echo "========================================"
echo "Job ID:       $SLURM_JOB_ID"
echo "Node:         $(hostname)"
echo "RAM:          $(free -h | grep Mem | awk '{print $2}')"
echo "Cluster root: $CLUSTER_ROOT"
echo "Started:      $(date)"
echo "========================================"

mkdir -p results

Rscript scripts/dmltest_full.R

echo "========================================"
echo "Finished: $(date)"
echo "Results:"
ls -lh results/dmps_strict.tsv results/dmrs.tsv results/dmltest_summary.tsv 2>/dev/null
echo "========================================"
