#!/bin/bash
#SBATCH --job-name=sort_bams
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=cluster/logs/sort_bams_%j.out
#SBATCH --error=cluster/logs/sort_bams_%j.err

module load samtools 2>/dev/null || true

BAM_DIR="/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/08_deduplication"
OUT_DIR="/mnt/data/alfredvar/rlopezt/bams_sorted"
mkdir -p "$OUT_DIR"

for sample in A1 A2 C1 C2; do
  BAM="${BAM_DIR}/${sample}_paired_bismark_bt2_pe.deduplicated.bam"
  SORTED="${OUT_DIR}/${sample}_paired_bismark_bt2_pe.deduplicated.sorted.bam"

  if [ -f "${SORTED}.bai" ]; then
    echo "${sample}: already sorted and indexed"
    continue
  fi

  echo "Sorting ${sample}..."
  samtools sort -@ 8 -m 3G -o "$SORTED" "$BAM"
  echo "Indexing ${sample}..."
  samtools index -@ 8 "$SORTED"
  echo "${sample} done: $(ls -lh "$SORTED" | awk '{print $5}')"
done

echo "All BAMs sorted and indexed"
