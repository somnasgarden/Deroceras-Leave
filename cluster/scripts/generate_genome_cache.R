#!/usr/bin/env Rscript
# =============================================================================
# Generate genome_chr1_31.rds cache on the cluster
# Sources _config.R for paths. Run from repo root.
# =============================================================================

source("methylation_pipeline/_config.R")

if (file.exists(CACHE$genome)) {
  cat("genome_chr1_31.rds already exists at:", CACHE$genome, "\n")
  cat("Size:", round(file.size(CACHE$genome) / 1e6, 1), "MB\n")
  q("no")
}

cat("Generating genome cache via load_genome()...\n")
genome <- load_genome()
cat("Done! Size:", round(file.size(CACHE$genome) / 1e6, 1), "MB\n")

# Also export FASTA for HOMER
fasta_out <- file.path(PROJECT_DIR, "genome/derLaeGenome_chr1_31.fasta")
if (!file.exists(fasta_out)) {
  cat("Exporting genome FASTA for HOMER...\n")
  library(Biostrings)
  writeXStringSet(genome, filepath = fasta_out)
  cat("FASTA saved:", fasta_out, "\n")
}
