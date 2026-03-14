#!/bin/bash
#SBATCH --job-name=dlaeve_motif_pos
#SBATCH --output=results/motif_positions_%j.log
#SBATCH --error=results/motif_positions_%j.err
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --partition=general
# ============================================================================
# SLURM job: Batch 1.5 Section G — Motif Coordinate Export
# ============================================================================
# Runs matchMotifs(out="positions") on 1362 JASPAR motifs x 25K promoters.
# This OOM'd on local WSL (7.6 GB RAM). Requesting 32 GB here.
#
# The cluster/ folder is self-contained. Clone/copy it to the HPC.
# Required contents:
#   cluster/
#     genome/cache/genome_chr1_31.rds
#     genome/cache/gff_chr1_31.rds
#     genome/cache/promoters_2kb.rds
#     scripts/batch1.5_motif_positions.R
#     results/   (empty dir for output)
#
# BEFORE SUBMITTING:
#   1. cd into the cluster/ folder on the HPC
#   2. Adjust: --partition, module load line
#   3. Submit:  sbatch scripts/batch1.5_motif_positions.slurm
#
# AFTER COMPLETION:
#   Copy back to local STANDBY/:
#     cluster/results/motif_hits_annotated.tsv.gz   → results/batch1.5/
#     cluster/results/motif_annotation_summary.tsv   → results/batch1.5/
# ============================================================================

# --- CONFIGURE THESE ---
# MODULE_R="R/4.4"  # Uncomment and adjust for your cluster

# --- Load R ---
# Option A: module system (most HPC clusters)
# module load $MODULE_R

# Option B: conda/mamba (if R is installed via conda on the cluster)
# conda activate r_env

# --- Run ---
# Run from inside the cluster/ folder
CLUSTER_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$CLUSTER_ROOT"
export CLUSTER_ROOT
export XDG_CACHE_HOME="${CLUSTER_ROOT}/.cache"
mkdir -p "$XDG_CACHE_HOME"

echo "=== Job started: $(date) ==="
echo "Node: $(hostname)"
echo "RAM: $(free -h | grep Mem | awk '{print $2}')"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Cluster root: $CLUSTER_ROOT"

Rscript scripts/batch1.5_motif_positions.R

echo "=== Job finished: $(date) ==="
echo "Exit code: $?"
