#!/bin/bash
#SBATCH --job-name=dlaeve_motif_expanded
#SBATCH --output=results/motif_expanded_%j.log
#SBATCH --error=results/motif_expanded_%j.err
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --partition=general
# =============================================================================
# SLURM: Expanded Motif Annotation — JASPAR + HOMER
# =============================================================================
# Supersedes batch1.5_motif_positions.slurm.
#
# This job:
#   1. Installs HOMER (if not already done)
#   2. Runs JASPAR motif scanning on 10kb upstream + gene body + 10kb downstream
#   3. Runs HOMER de novo motif discovery on promoters
#   4. Annotates every motif hit with region type and distance to TSS
#
# BEFORE SUBMITTING:
#   1. cd into the cluster/ folder on the HPC
#   2. Ensure genome FASTA is accessible (set GENOME_FASTA below or place in
#      cluster/genome/derLaeGenome_chr1_31.fasta)
#   3. Ensure RDS caches are in cluster/genome/cache/:
#        - gff_chr1_31.rds (12 MB) — OR the GFF is at the cluster path below
#        - genome_chr1_31.rds (365 MB) — OR set GENOME_FASTA
#        - promoters_2kb.rds (optional, regions are recomputed)
#   4. Adjust --partition for your cluster
#   5. Submit: sbatch scripts/batch1.5_expanded_motif.slurm
#
# AFTER COMPLETION:
#   Copy cluster/results/* back to local results/batch1.5/
# =============================================================================

# ---- CONFIGURE THESE ----

# Genome: The R script loads BSgenome.Dlaeve.NCBI.dlgm (installed in cluster R env).
# Falls back to genome_chr1_31.rds cache or FASTA if BSgenome not available.
# You should NOT need to set this if BSgenome is installed.
export GENOME_FASTA=""

# R environment — uncomment ONE option:
# Option A: module system
# module load R/4.4

# Option B: conda/mamba
# source activate r_env
# OR: conda activate r_env

# REQUIRED PACKAGES (install once on your cluster R env if missing):
#   BiocManager::install(c("TFBSTools", "JASPAR2024", "motifmatchr", "universalmotif"))
# You already have: GenomicRanges, rtracklayer, Biostrings, data.table, dplyr,
#                   BSgenome.Dlaeve.NCBI.dlgm

# ---- SETUP ----
CLUSTER_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$CLUSTER_ROOT"
export CLUSTER_ROOT
export XDG_CACHE_HOME="${CLUSTER_ROOT}/.cache"
mkdir -p "$XDG_CACHE_HOME"
mkdir -p results

echo "========================================"
echo "=== Expanded Motif Annotation Job ==="
echo "========================================"
echo "Job ID:       $SLURM_JOB_ID"
echo "Node:         $(hostname)"
echo "RAM:          $(free -h | grep Mem | awk '{print $2}')"
echo "CPUs:         $SLURM_CPUS_PER_TASK"
echo "Cluster root: $CLUSTER_ROOT"
echo "Started:      $(date)"
echo "========================================"

# ---- STEP 1: Install HOMER base (may skip genome if FASTA not yet exported) ----
echo ""
echo "=== Step 1: HOMER Install ==="
bash scripts/setup_homer.sh
export PATH="${CLUSTER_ROOT}/homer/bin:${PATH}"

# ---- STEP 2: Run expanded motif annotation (R) ----
# This loads genome from BSgenome.Dlaeve.NCBI.dlgm, exports FASTA,
# runs JASPAR scanning, and calls HOMER de novo internally.
echo ""
echo "=== Step 2: Expanded Motif Annotation (R) ==="
Rscript scripts/batch1.5_expanded_motif.R

# ---- STEP 3: Re-run HOMER genome setup if FASTA was just exported ----
# The R script exports the genome FASTA from BSgenome. Now HOMER can configure it.
if [ -f "${CLUSTER_ROOT}/genome/derLaeGenome_chr1_31.fasta" ]; then
  export GENOME_FASTA="${CLUSTER_ROOT}/genome/derLaeGenome_chr1_31.fasta"
  echo ""
  echo "=== Step 3: HOMER Genome Config (post-export) ==="
  bash scripts/setup_homer.sh
fi

echo ""
echo "========================================"
echo "=== Job Complete ==="
echo "Finished: $(date)"
echo "Exit code: $?"
echo ""
echo "Output files:"
ls -lh results/motif_hits_extended.tsv.gz 2>/dev/null
ls -lh results/motif_summary_by_region.tsv 2>/dev/null
ls -lh results/motif_region_overview.tsv 2>/dev/null
ls -lh results/homer_denovo/ 2>/dev/null
echo "========================================"
