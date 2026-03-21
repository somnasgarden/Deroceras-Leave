#!/usr/bin/env bash
#SBATCH --job-name=motif_annot
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

# ============================================================
# Expanded Motif Annotation — JASPAR + HOMER
# Run from: STANDBY/cluster/
#   cd cluster && sbatch scripts/batch1.5_expanded_motif.slurm
# ============================================================

CLUSTER_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$CLUSTER_ROOT"
export CLUSTER_ROOT
export XDG_CACHE_HOME="${CLUSTER_ROOT}/.cache"
mkdir -p "$XDG_CACHE_HOME"
mkdir -p results

echo "========================================"
echo "Job ID:       $SLURM_JOB_ID"
echo "Node:         $(hostname)"
echo "RAM:          $(free -h | grep Mem | awk '{print $2}')"
echo "CPUs:         $SLURM_CPUS_PER_TASK"
echo "Cluster root: $CLUSTER_ROOT"
echo "Started:      $(date)"
echo "========================================"

# Step 1: Install HOMER (if not already done)
echo "=== Step 1: HOMER Install ==="
bash scripts/setup_homer.sh
export PATH="${CLUSTER_ROOT}/homer/bin:${PATH}"

# Step 2: Run expanded motif annotation (R)
echo "=== Step 2: Expanded Motif Annotation (R) ==="
Rscript scripts/batch1.5_expanded_motif.R

# Step 3: Re-run HOMER genome setup if FASTA was exported by R
if [ -f "${CLUSTER_ROOT}/genome/derLaeGenome_chr1_31.fasta" ]; then
  export GENOME_FASTA="${CLUSTER_ROOT}/genome/derLaeGenome_chr1_31.fasta"
  echo "=== Step 3: HOMER Genome Config (post-export) ==="
  bash scripts/setup_homer.sh
fi

echo "========================================"
echo "=== Job Complete ==="
echo "Finished: $(date)"
echo "Output files:"
ls -lh results/motif_hits_extended.tsv.gz 2>/dev/null
ls -lh results/motif_summary_by_region.tsv 2>/dev/null
ls -lh results/homer_denovo/ 2>/dev/null
echo "========================================"
