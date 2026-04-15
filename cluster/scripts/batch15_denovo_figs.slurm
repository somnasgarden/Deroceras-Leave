#!/bin/bash
#SBATCH --job-name=denovo_figs
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=cluster/logs/denovo_figs_%j.out
#SBATCH --error=cluster/logs/denovo_figs_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

echo "=== Running batch 1.5 with de novo motif characterization ==="
echo "Start: $(date)"

# The script skips JASPAR scanning if motif_hits_extended.tsv.gz already exists
# and only runs the new de novo characterization section
Rscript methylation_pipeline/batch1.5/code/015_motif_annotation.R

echo ""
echo "=== Regenerating HTML report ==="
Rscript generate_report.R

echo "End: $(date)"
