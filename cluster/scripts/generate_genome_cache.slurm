#!/usr/bin/env bash
#SBATCH --job-name=gen_cache
#SBATCH --output=cluster/logs/slurm_%j.out
#SBATCH --error=cluster/logs/slurm_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

# Generate genome_chr1_31.rds + FASTA from BSgenome
# Run from repo root: sbatch cluster/scripts/generate_genome_cache.slurm

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

mkdir -p cluster/logs genome/cache

echo "=== Genome Cache Generation ==="
echo "Job ID: $SLURM_JOB_ID | Node: $(hostname) | Started: $(date)"
echo "Project: $PROJECT_DIR"

Rscript cluster/scripts/generate_genome_cache.R

echo "=== Done: $(date) ==="
