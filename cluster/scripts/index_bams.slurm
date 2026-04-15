#!/bin/bash
#SBATCH --job-name=idx_bams
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=cluster/logs/index_bams_%j.out
#SBATCH --error=cluster/logs/index_bams_%j.err

module load samtools 2>/dev/null || true

BAM_DIR="/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/08_deduplication"

for sample in A1 A2 C2; do
  BAM="${BAM_DIR}/${sample}_paired_bismark_bt2_pe.deduplicated.bam"
  BAI="${BAM}.bai"
  if [ ! -f "$BAI" ]; then
    echo "Indexing ${sample}..."
    samtools index -@ 4 "$BAM"
    echo "${sample} done"
  else
    echo "${sample} already indexed"
  fi
done

echo "All indexing complete"
