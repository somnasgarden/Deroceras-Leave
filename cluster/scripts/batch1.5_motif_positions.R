#!/usr/bin/env Rscript
# =============================================================================
# Cluster job: Batch 1.5 Section G — Motif Coordinate Export
# =============================================================================
# This script runs the memory-heavy part of batch1.5 that OOM'd locally:
#   matchMotifs(..., out="positions") on 1362 JASPAR motifs x 25K promoters
#
# REQUIRES: ~16-32 GB RAM
#
# Input (must exist before running):
#   - genome/cache/genome_chr1_31.rds   (365 MB — cached genome)
#   - genome/cache/gff_chr1_31.rds      (12 MB  — cached GFF)
#   - genome/cache/promoters_2kb.rds    (from batch1.5 local run)
#
# Output (saved to cluster/results/):
#   - motif_hits_annotated.tsv.gz       (compressed coordinate table)
#   - motif_annotation_summary.tsv      (per-motif summary stats)
#
# After running, copy cluster/results/* back to results/batch1.5/
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

# Set cache dir BEFORE loading Bioconductor packages
cache_home <- Sys.getenv("XDG_CACHE_HOME", "~/.cache")
Sys.setenv(XDG_CACHE_HOME = cache_home)

library(data.table)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(GenomeInfoDb)
library(TFBSTools)
library(JASPAR2024)
library(motifmatchr)

# -- Paths --
# cluster/ is the self-contained root when running on the HPC.
# Set CLUSTER_ROOT to override (defaults to parent of this script's dir).
cluster_root <- Sys.getenv("CLUSTER_ROOT", getwd())
cache_dir    <- file.path(cluster_root, "genome/cache")
out_dir      <- file.path(cluster_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

keep_chr <- paste0("chr", 1:31)

cat("=== Motif Coordinate Export (Cluster) ===\n")
cat("Cluster root:", cluster_root, "\n")
cat("Cache dir:", cache_dir, "\n")
cat("Output dir:", out_dir, "\n")
cat("Available RAM:", system("free -h | grep Mem | awk '{print $2}'", intern = TRUE), "\n\n")


# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("--- Loading cached data ---\n")

genome <- readRDS(file.path(cache_dir, "genome_chr1_31.rds"))
cat("Genome:", length(genome), "chromosomes\n")

gff <- readRDS(file.path(cache_dir, "gff_chr1_31.rds"))
genes <- gff[gff$type == "gene"]
cat("Genes:", length(genes), "\n")

promoters <- readRDS(file.path(cache_dir, "promoters_2kb.rds"))
cat("Promoters:", length(promoters), "\n")


# =============================================================================
# 2. EXTRACT PROMOTER SEQUENCES
# =============================================================================
cat("\n--- Extracting promoter sequences ---\n")
prom_seqs <- getSeq(genome, promoters)
names(prom_seqs) <- promoters$gene_id
cat("Promoter sequences:", length(prom_seqs), "\n")

# Free genome (no longer needed)
rm(genome, gff, genes); gc()


# =============================================================================
# 3. LOAD JASPAR MOTIFS
# =============================================================================
cat("\n--- Loading JASPAR 2024 motifs ---\n")

# JASPAR2024 SQLite — let BiocFileCache handle download
jaspar_sqlite <- list.files(cache_home, pattern = "JASPAR.*sqlite",
                            recursive = TRUE, full.names = TRUE)
if (length(jaspar_sqlite) == 0) {
  cat("Downloading JASPAR2024 database...\n")
  jaspar_db <- JASPAR2024()
  jaspar_sqlite <- list.files(cache_home, pattern = "JASPAR.*sqlite",
                              recursive = TRUE, full.names = TRUE)
}
jaspar_sqlite <- jaspar_sqlite[1]
cat("JASPAR SQLite:", jaspar_sqlite, "\n")

# Get CORE metazoan motifs
all_pfms <- getMatrixSet(jaspar_sqlite, list(collection = "CORE", all_versions = FALSE))
metazoan_groups <- c("insects", "nematodes", "urochordates", "vertebrates")
is_metazoa <- sapply(all_pfms, function(x) {
  tg <- tags(x)$tax_group
  !is.null(tg) && tg %in% metazoan_groups
})
jaspar_pfms <- all_pfms[is_metazoa]
rm(all_pfms); gc()
cat("Metazoan motifs:", length(jaspar_pfms), "\n")

# Build motif metadata
motif_meta <- data.frame(
  motif_id   = sapply(jaspar_pfms, ID),
  motif_name = sapply(jaspar_pfms, name),
  tf_class   = vapply(jaspar_pfms, function(x) {
    tgs <- TFBSTools::tags(x)
    if ("family" %in% names(tgs) && length(tgs$family) > 0) tgs$family[1]
    else NA_character_
  }, character(1)),
  species = vapply(jaspar_pfms, function(x) {
    tgs <- TFBSTools::tags(x)
    if ("species" %in% names(tgs) && length(tgs$species) > 0)
      paste(tgs$species, collapse = ";")
    else NA_character_
  }, character(1)),
  stringsAsFactors = FALSE
)


# =============================================================================
# 4. MOTIF POSITION SCAN (THE MEMORY-HEAVY STEP)
# =============================================================================
cat("\n--- Scanning for motif positions ---\n")
cat("This is the memory-intensive step: matchMotifs(out='positions')\n")
cat(sprintf("  %d motifs x %d promoters\n", length(jaspar_pfms), length(prom_seqs)))

scan_start <- proc.time()
pos_result <- matchMotifs(
  jaspar_pfms,
  prom_seqs,
  out = "positions",
  p.cutoff = 5e-5
)
scan_time <- (proc.time() - scan_start)[3]
cat(sprintf("Position scan complete in %.1f seconds\n", scan_time))

# Extract positions — GRangesList of length = n_motifs
pos_list <- motifPositions(pos_result)
rm(pos_result); gc()


# =============================================================================
# 5. CONVERT TO GENOMIC COORDINATES
# =============================================================================
cat("\n--- Converting to genomic coordinates ---\n")

# Promoter lookup for coordinate conversion
prom_lookup <- data.frame(
  gene_id     = promoters$gene_id,
  prom_chr    = as.character(seqnames(promoters)),
  prom_start  = start(promoters),
  prom_end    = end(promoters),
  prom_strand = as.character(strand(promoters)),
  stringsAsFactors = FALSE
)

# Flatten GRangesList into annotated data.frame
cat("Flattening", length(pos_list), "motifs into coordinate table...\n")
all_hits <- vector("list", length(pos_list))

for (i in seq_along(pos_list)) {
  gr <- pos_list[[i]]
  if (length(gr) == 0) next

  hit_genes <- as.character(seqnames(gr))
  prom_idx  <- match(hit_genes, prom_lookup$gene_id)

  rel_start  <- start(gr)
  rel_end    <- end(gr)
  hit_strand <- as.character(strand(gr))

  ps      <- prom_lookup$prom_strand[prom_idx]
  p_start <- prom_lookup$prom_start[prom_idx]
  p_end   <- prom_lookup$prom_end[prom_idx]

  # Convert to genomic coordinates
  # + strand: genomic = prom_start + relative - 1
  # - strand: genomic_start = prom_end - rel_end + 1
  gen_start <- ifelse(ps == "+", p_start + rel_start - 1L, p_end - rel_end + 1L)
  gen_end   <- ifelse(ps == "+", p_start + rel_end - 1L,   p_end - rel_start + 1L)

  # Flip motif strand if promoter is on - strand
  gen_strand <- ifelse(hit_strand == "*", "*",
                       ifelse(ps == "+", hit_strand,
                              ifelse(hit_strand == "+", "-", "+")))

  score_vals <- if ("score" %in% colnames(mcols(gr))) {
    mcols(gr)$score
  } else {
    rep(NA_real_, length(gr))
  }

  all_hits[[i]] <- data.frame(
    chr        = prom_lookup$prom_chr[prom_idx],
    start      = gen_start,
    end        = gen_end,
    strand     = gen_strand,
    motif_id   = motif_meta$motif_id[i],
    motif_name = motif_meta$motif_name[i],
    tf_class   = motif_meta$tf_class[i],
    species    = motif_meta$species[i],
    gene_id    = hit_genes,
    score      = score_vals,
    stringsAsFactors = FALSE
  )

  if (i %% 200 == 0) cat(sprintf("  Processed %d / %d motifs\n", i, length(pos_list)))
}

hits_dt <- data.table::rbindlist(all_hits, fill = TRUE)
rm(all_hits, pos_list); gc()

cat(sprintf("Total motif hits: %s\n", format(nrow(hits_dt), big.mark = ",")))
cat(sprintf("Unique motifs with hits: %d / %d\n",
            length(unique(hits_dt$motif_id)), nrow(motif_meta)))
cat(sprintf("Unique genes with hits: %d / %d\n",
            length(unique(hits_dt$gene_id)), length(promoters)))


# =============================================================================
# 6. SAVE OUTPUT
# =============================================================================
cat("\n--- Saving output ---\n")

# Full coordinate table (compressed)
out_hits <- file.path(out_dir, "motif_hits_annotated.tsv.gz")
data.table::fwrite(hits_dt, out_hits, sep = "\t", compress = "gzip")
cat(sprintf("Coordinate table: %s (%.1f MB)\n", out_hits,
            file.size(out_hits) / 1e6))

# Per-motif summary
hit_summary <- hits_dt[, .(
  n_hits      = .N,
  n_promoters = uniqueN(gene_id)
), by = .(motif_id)]

motif_summary <- merge(
  as.data.table(motif_meta),
  hit_summary,
  by = "motif_id",
  all.x = TRUE
)
motif_summary[is.na(n_hits), c("n_hits", "n_promoters") := .(0L, 0L)]
motif_summary[, pct_promoters := round(100 * n_promoters / length(promoters), 2)]
motif_summary <- motif_summary[order(-n_hits)]

out_summary <- file.path(out_dir, "motif_annotation_summary.tsv")
data.table::fwrite(motif_summary, out_summary, sep = "\t")
cat(sprintf("Summary table: %d motifs → %s\n", nrow(motif_summary), out_summary))


# =============================================================================
# 7. VERIFICATION
# =============================================================================
cat("\n--- Verification ---\n")
cat(sprintf("  Rows in hits table: %s\n", format(nrow(hits_dt), big.mark = ",")))
cat(sprintf("  Columns: %s\n", paste(colnames(hits_dt), collapse = ", ")))

# Spot-check: sampled hits should overlap promoter regions
n_check <- min(1000, nrow(hits_dt))
check_idx <- sample(nrow(hits_dt), n_check)
check_gr <- GRanges(
  seqnames = hits_dt$chr[check_idx],
  ranges   = IRanges(hits_dt$start[check_idx], hits_dt$end[check_idx])
)
in_prom <- sum(overlapsAny(check_gr, promoters))
cat(sprintf("  Coordinate check: %d/%d sampled hits overlap promoters (%.1f%%)\n",
            in_prom, n_check, 100 * in_prom / n_check))

rm(hits_dt, motif_summary, check_gr); gc()

cat("\n=== Cluster job complete ===\n")
cat("Copy results back:\n")
cat("  cluster/results/motif_hits_annotated.tsv.gz  →  results/batch1.5/\n")
cat("  cluster/results/motif_annotation_summary.tsv  →  results/batch1.5/\n")
