#!/bin/bash
#SBATCH --job-name=batch10b
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output=cluster/logs/batch10b_%j.out
#SBATCH --error=cluster/logs/batch10b_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

echo "=== Batch 10b: Per-read NME (Fang/Xie method) ==="
echo "Start: $(date)"
echo ""

# --- Sort C1 BAM if not already sorted (the original is unsorted) ---
SORTED_DIR=/mnt/data/alfredvar/rlopezt/bams_sorted
SRC_C1=/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/08_deduplication/C1_paired_bismark_bt2_pe.deduplicated.bam
DST_C1=$SORTED_DIR/C1_paired_bismark_bt2_pe.deduplicated.sorted.bam

mkdir -p $SORTED_DIR
source /etc/profile.d/lmod.sh 2>/dev/null || true
module load samtools/1.22.1

if [ ! -f "$DST_C1" ] || [ ! -f "${DST_C1}.bai" ]; then
  echo "Sorting C1 BAM..."
  samtools sort -@ 8 -m 8G -T $SORTED_DIR/c1_tmp -o $DST_C1 $SRC_C1
  echo "Indexing C1 sorted BAM..."
  samtools index -@ 8 $DST_C1
  echo "C1 sort+index done: $(date)"
else
  echo "C1 sorted BAM already present, skipping sort."
fi

# Index other BAMs (these are already sorted, just need .bai)
for nm in C2 A1 A2; do
  SRC=/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/08_deduplication/${nm}_paired_bismark_bt2_pe.deduplicated.bam
  if [ ! -f "${SRC}.bai" ] && [ ! -f "${SRC%.bam}.bai" ]; then
    DST=$SORTED_DIR/${nm}.bai
    if [ ! -f "$DST" ]; then
      echo "Indexing $nm BAM (writing index alongside source if writable, else $SORTED_DIR)..."
      samtools index -@ 8 $SRC 2>/dev/null || {
        # Source dir not writable — copy index out via samtools quickcheck workaround
        echo "  source dir not writable; will let R/Rsamtools index into a writable location at runtime"
      }
    fi
  fi
done

echo ""
echo "Running batch10b R script..."
Rscript methylation_pipeline/batch10/code/10b_perread_entropy.R

echo ""
echo "=== Batch 10b complete: $(date) ==="
