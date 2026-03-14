# ==============================================================================
# Genome extraction pipeline — Deroceras laeve
# Extract chr1–chr31 from BSgenome and export FASTA to STANDBY/
# ==============================================================================

library(BSgenome.Dlaeve.NCBI.ASM5140357v1)
library(Biostrings)

# --- Load genome -------------------------------------------------------------
genome <- BSgenome.Dlaeve.NCBI.ASM5140357v1
genome

# --- Filter to chr1–chr31 (drop HiC scaffolds) -------------------------------
nm       <- seqnames(genome)
keep_chr <- intersect(paste0("chr", 1:31), nm)

genome_chr <- getSeq(genome, keep_chr)

cat("Chromosomes retained:", length(genome_chr), "\n")
cat("Names:", head(names(genome_chr)), "...\n")

# --- Export genome FASTA ------------------------------------------------------
out_fasta <- "/mnt/c/Users/rafae/Projects/STANDBY/derLaeGenome_chr1_31.fasta"
writeXStringSet(genome_chr, filepath = out_fasta)
cat("Genome FASTA written to:", out_fasta, "\n")
