#!/usr/bin/env bash
#SBATCH --job-name=gen_cache
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

set -euo pipefail

# Generate genome_chr1_31.rds + FASTA from BSgenome
# Run FIRST before other cluster jobs
# cd cluster && sbatch scripts/generate_genome_cache.slurm

if [ -f "${SLURM_SUBMIT_DIR}/scripts/generate_genome_cache.R" ]; then
  CLUSTER_ROOT="$SLURM_SUBMIT_DIR"
else
  CLUSTER_ROOT="$(dirname "$SLURM_SUBMIT_DIR")"
fi
cd "$CLUSTER_ROOT"
export CLUSTER_ROOT

Rscript scripts/generate_genome_cache.R
