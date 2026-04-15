#!/bin/bash
#SBATCH --job-name=wg_motif
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=72:00:00
#SBATCH --output=/mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave/cluster/logs/wg_motif_%j.out
#SBATCH --error=/mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave/cluster/logs/wg_motif_%j.err

echo "=== Whole-Genome Motif Scan ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "RAM: $(free -h | grep Mem | awk '{print $2}')"
echo ""

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

export PROJECT_DIR=$(pwd)

Rscript cluster/scripts/whole_genome_motif_scan.R

echo ""
echo "=== Done: $(date) ==="
