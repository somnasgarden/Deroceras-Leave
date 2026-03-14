#!/bin/bash
#SBATCH --job-name=dlaeve_motif_pos
#SBATCH --output=cluster/results/motif_positions_%j.log
#SBATCH --error=cluster/results/motif_positions_%j.err
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
# BEFORE SUBMITTING:
#   1. Copy the STANDBY project to the cluster (or at minimum these files):
#        genome/cache/genome_chr1_31.rds
#        genome/cache/gff_chr1_31.rds
#        genome/cache/promoters_2kb.rds
#        cluster/scripts/batch1.5_motif_positions.R
#        cluster/results/   (empty dir for output)
#
#   2. Adjust these variables:
#        - SBATCH --partition  (your cluster's partition name)
#        - MODULE_R            (your cluster's R module name)
#        - PROJECT_ROOT        (path to STANDBY on the cluster)
#
#   3. Submit:  sbatch cluster/scripts/batch1.5_motif_positions.slurm
#
# AFTER COMPLETION:
#   Copy back to local:
#     cluster/results/motif_hits_annotated.tsv.gz   → results/batch1.5/
#     cluster/results/motif_annotation_summary.tsv   → results/batch1.5/
# ============================================================================

# --- CONFIGURE THESE ---
PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
# MODULE_R="R/4.4"  # Uncomment and adjust for your cluster

# --- Load R ---
# Option A: module system (most HPC clusters)
# module load $MODULE_R

# Option B: conda/mamba (if R is installed via conda on the cluster)
# conda activate r_env

# --- Run ---
cd "$PROJECT_ROOT"
export PROJECT_ROOT
export XDG_CACHE_HOME="${PROJECT_ROOT}/cluster/.cache"
mkdir -p "$XDG_CACHE_HOME"

echo "=== Job started: $(date) ==="
echo "Node: $(hostname)"
echo "RAM: $(free -h | grep Mem | awk '{print $2}')"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Project: $PROJECT_ROOT"

Rscript cluster/scripts/batch1.5_motif_positions.R

echo "=== Job finished: $(date) ==="
echo "Exit code: $?"
