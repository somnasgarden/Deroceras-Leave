#!/bin/bash
#SBATCH --job-name=fix_bams
#SBATCH --partition=defq
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=cluster/logs/fix_bams_%j.out
#SBATCH --error=cluster/logs/fix_bams_%j.err

cd /mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave

# Force temp files to shared storage (not home, which has 5GB quota)
export TMPDIR="/mnt/data/alfredvar/rlopezt/bams_sorted/tmp"
mkdir -p "$TMPDIR"

BAM_DIR="/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/08_deduplication"
OUT_DIR="/mnt/data/alfredvar/rlopezt/bams_sorted"

Rscript -e '
Sys.setenv(TMPDIR="/mnt/data/alfredvar/rlopezt/bams_sorted/tmp")
library(Rsamtools)

bam_dir <- Sys.getenv("BAM_DIR", "/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/08_deduplication")
out_dir <- "/mnt/data/alfredvar/rlopezt/bams_sorted"

# Sort C2
cat("Sorting C2 (TMPDIR =", Sys.getenv("TMPDIR"), ")...\n")
c2_in  <- file.path(bam_dir, "C2_paired_bismark_bt2_pe.deduplicated.bam")
c2_out <- file.path(out_dir, "C2_paired_bismark_bt2_pe.deduplicated.sorted")
sortBam(c2_in, c2_out, maxMemory = 4096)
cat("C2 sorted:", file.size(paste0(c2_out, ".bam")), "bytes\n")

cat("Indexing C2...\n")
indexBam(paste0(c2_out, ".bam"))

# Re-sort A2 from source (previous version was truncated)
cat("Re-sorting A2...\n")
a2_in  <- file.path(bam_dir, "A2_paired_bismark_bt2_pe.deduplicated.bam")
a2_out <- file.path(out_dir, "A2_paired_bismark_bt2_pe.deduplicated.sorted")
sortBam(a2_in, a2_out, maxMemory = 4096)
cat("A2 sorted:", file.size(paste0(a2_out, ".bam")), "bytes\n")

cat("Indexing A2...\n")
indexBam(paste0(a2_out, ".bam"))

# Verify
for (s in c("A1","A2","C1","C2")) {
  bam <- file.path(out_dir, paste0(s, "_paired_bismark_bt2_pe.deduplicated.sorted.bam"))
  bai <- paste0(bam, ".bai")
  if (file.exists(bam) && file.exists(bai)) {
    cat(sprintf("%s: OK (%.1f GB)\n", s, file.size(bam)/1e9))
  } else {
    cat(sprintf("%s: MISSING (bam=%s, bai=%s)\n", s, file.exists(bam), file.exists(bai)))
  }
}

# Cleanup tmp
unlink("/mnt/data/alfredvar/rlopezt/bams_sorted/tmp", recursive=TRUE)
cat("Done\n")
'
