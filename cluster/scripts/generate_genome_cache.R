#!/usr/bin/env Rscript
# =============================================================================
# Generate genome_chr1_31.rds cache on the cluster
# Run once: Rscript scripts/generate_genome_cache.R
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

CLUSTER_ROOT <- Sys.getenv("CLUSTER_ROOT", ".")
cache_dir <- file.path(CLUSTER_ROOT, "genome/cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(cache_dir, "genome_chr1_31.rds")

if (file.exists(out_rds)) {
  cat("genome_chr1_31.rds already exists at:", out_rds, "\n")
  cat("Size:", round(file.size(out_rds) / 1e6, 1), "MB\n")
  q("no")
}

cat("Loading BSgenome.Dlaeve.NCBI.ASM5140357v1...\n")
library(BSgenome.Dlaeve.NCBI.ASM5140357v1)
library(Biostrings)

genome <- BSgenome.Dlaeve.NCBI.ASM5140357v1
keep_chr <- paste0("chr", 1:31)
keep_chr <- intersect(keep_chr, seqnames(genome))

cat("Extracting", length(keep_chr), "chromosomes...\n")
seqs <- getSeq(genome, keep_chr)
names(seqs) <- keep_chr

cat("Saving to:", out_rds, "\n")
saveRDS(seqs, out_rds)
cat("Done! Size:", round(file.size(out_rds) / 1e6, 1), "MB\n")

# Also export FASTA for HOMER
fasta_out <- file.path(CLUSTER_ROOT, "genome/derLaeGenome_chr1_31.fasta")
if (!file.exists(fasta_out)) {
  cat("Exporting genome FASTA for HOMER...\n")
  writeXStringSet(seqs, filepath = fasta_out)
  cat("FASTA saved:", fasta_out, "\n")
}
