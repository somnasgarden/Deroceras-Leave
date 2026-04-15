#!/bin/bash
#SBATCH --job-name=rerun_fixes
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=cluster/logs/rerun_fixes_%j.out
#SBATCH --error=cluster/logs/rerun_fixes_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

echo "=== Rerunning fixed batches ==="
echo "Start: $(date)"

# 1. batch03 (fix hardcoded subtitle) - fast
echo ""
echo "=== batch03: Promoter CpG Classification ==="
Rscript methylation_pipeline/batch03/code/03_promoter_cpg_classification.R

# 2. batch06 (fix duplicate annotation rows) - medium
echo ""
echo "=== batch06: Differential Methylation ==="
Rscript methylation_pipeline/batch06/code/06_differential_methylation.R

# 3. batch07 (depends on batch06 fix) - medium
echo ""
echo "=== batch07: Methylation-Expression Decoupling ==="
Rscript methylation_pipeline/batch07/code/07_methylation_expression_decoupling.R

# 4. batch09 (annot dedup fix) - medium
echo ""
echo "=== batch09: TF Motif Integration ==="
Rscript methylation_pipeline/batch09/code/09_tf_motif_integration.R

# 5. batch10 (new Fang framework analyses) - heavy
echo ""
echo "=== batch10: Entropy Analysis (with Fang framework) ==="
Rscript methylation_pipeline/batch10/code/10_methylation_entropy.R

# 6. Regenerate HTML report
echo ""
echo "=== Generating pipeline_report.html ==="
Rscript generate_report.R

echo ""
echo "=== All done ==="
echo "End: $(date)"
