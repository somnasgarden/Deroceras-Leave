#!/usr/bin/env Rscript
# =============================================================================
# Whole-Genome TFBS Motif Annotation — JASPAR 2024 CORE Metazoan
# =============================================================================
# Scans ALL of chr1-31 (not just gene-proximal regions).
# Per-chromosome × 50-motif chunks to stay within 64 GB RAM.
#
# Output: cluster/results/whole_genome_motif_hits.tsv.gz
#         + per-chromosome chunk files (merged at end)
#
# After completion, annotate with genomic context (genic, intergenic, TE,
# promoter, etc.) for enrichment testing.
#
# REQUIRES: ~64 GB RAM, ~24-48 hours
# =============================================================================

source("methylation_pipeline/_config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(Biostrings)
  library(TFBSTools)
  library(JASPAR2024)
  library(motifmatchr)
})

out_dir <- file.path(PROJECT_DIR, "cluster/results/whole_genome")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Whole-Genome Motif Scan ===\n")
cat("Project dir:", PROJECT_DIR, "\n")
cat("Output dir:", out_dir, "\n")
cat("Available RAM:", system("free -h | grep Mem | awk '{print $2}'", intern = TRUE), "\n")
cat("Start time:", format(Sys.time()), "\n\n")
t0_global <- proc.time()[3]

# =============================================================================
# 1. LOAD GENOME + JASPAR MOTIFS
# =============================================================================
cat("=== 1. Loading genome and motifs ===\n")

genome <- load_genome()
chr_lengths <- setNames(width(genome), names(genome))
cat(sprintf("Genome: %d chromosomes, %s bp total\n",
            length(genome), format(sum(chr_lengths), big.mark = ",")))

# JASPAR 2024 CORE metazoan
cache_home <- Sys.getenv("XDG_CACHE_HOME", file.path(PROJECT_DIR, "cluster/.cache"))
Sys.setenv(XDG_CACHE_HOME = cache_home)

jaspar_sqlite <- list.files(cache_home, pattern = "JASPAR.*sqlite",
                            recursive = TRUE, full.names = TRUE)
if (length(jaspar_sqlite) == 0) {
  cat("Downloading JASPAR2024 database...\n")
  jaspar_db <- JASPAR2024()
  jaspar_sqlite <- list.files(cache_home, pattern = "JASPAR.*sqlite",
                              recursive = TRUE, full.names = TRUE)
}

all_pfms <- getMatrixSet(jaspar_sqlite[1], list(collection = "CORE", all_versions = FALSE))
metazoan_groups <- c("insects", "nematodes", "urochordates", "vertebrates")
is_metazoa <- sapply(all_pfms, function(x) {
  tg <- tags(x)$tax_group
  !is.null(tg) && tg %in% metazoan_groups
})
jaspar_pfms <- all_pfms[is_metazoa]
rm(all_pfms); gc(verbose = FALSE)
cat(sprintf("JASPAR metazoan motifs: %d\n", length(jaspar_pfms)))

# Motif metadata
motif_meta <- data.frame(
  motif_id   = sapply(jaspar_pfms, ID),
  motif_name = sapply(jaspar_pfms, name),
  tf_class   = vapply(jaspar_pfms, function(x) {
    tgs <- TFBSTools::tags(x)
    if ("family" %in% names(tgs) && length(tgs$family) > 0) tgs$family[1]
    else NA_character_
  }, character(1)),
  stringsAsFactors = FALSE
)
fwrite(motif_meta, file.path(out_dir, "motif_metadata.tsv"), sep = "\t")

# =============================================================================
# 2. LOAD GFF FOR ANNOTATION
# =============================================================================
cat("\n=== 2. Loading GFF for post-scan annotation ===\n")

gff <- load_gff()
genes <- gff[gff$type == "gene"]
genes <- genes[as.character(seqnames(genes)) %in% keep_chr]
exons <- gff[gff$type == "exon"]
exons <- exons[as.character(seqnames(exons)) %in% keep_chr]

# Promoters (2kb upstream of TSS)
gene_strand <- as.character(strand(genes))
gene_tss <- ifelse(gene_strand == "+", start(genes), end(genes))
prom_start <- ifelse(gene_strand == "+", pmax(1L, gene_tss - 2000L), gene_tss + 1L)
prom_end   <- ifelse(gene_strand == "+", gene_tss - 1L,
                     pmin(as.integer(chr_lengths[as.character(seqnames(genes))]), gene_tss + 2000L))
promoters <- GRanges(seqnames = seqnames(genes),
                     ranges = IRanges(start = prom_start, end = prom_end),
                     gene_id = genes$ID)

# TSS lookup
tss_lookup <- data.table(
  gene_id = genes$ID,
  chr     = as.character(seqnames(genes)),
  tss     = gene_tss,
  strand  = gene_strand
)

# TE regions (for annotation)
te <- load_te()
te_gr <- te$te_gr
if (is.null(te_gr) || length(te_gr) == 0) {
  cat("  WARNING: TE data empty/NULL — TE annotation will be skipped\n")
  te_gr <- GRanges()  # empty placeholder so overlapsAny doesn't crash
}

cat(sprintf("  Genes: %d | Exons: %d | TEs: %d\n", length(genes), length(exons), length(te_gr)))

# =============================================================================
# 3. WHOLE-GENOME SCAN: per-chromosome × motif chunks
# =============================================================================
cat("\n=== 3. Scanning whole genome ===\n")

CHUNK_SIZE <- 50L  # motifs per chunk (same as batch 1.5)
n_motifs <- length(jaspar_pfms)
n_chunks <- ceiling(n_motifs / CHUNK_SIZE)
cat(sprintf("  %d motifs in %d chunks of %d\n", n_motifs, n_chunks, CHUNK_SIZE))

total_hits <- 0L
chunk_files <- c()

for (chr_name in keep_chr) {
  chr_t0 <- proc.time()[3]
  chr_seq <- genome[[chr_name]]
  chr_len <- length(chr_seq)
  cat(sprintf("\n  --- %s (%s bp) ---\n", chr_name, format(chr_len, big.mark = ",")))

  # Wrap in DNAStringSet with name (matchMotifs needs named set)
  chr_dna <- DNAStringSet(chr_seq)
  names(chr_dna) <- chr_name

  chr_hits_list <- vector("list", n_chunks)

  for (chunk_i in seq_len(n_chunks)) {
    chunk_start <- (chunk_i - 1L) * CHUNK_SIZE + 1L
    chunk_end   <- min(chunk_i * CHUNK_SIZE, n_motifs)
    pfm_chunk   <- jaspar_pfms[chunk_start:chunk_end]

    t1 <- proc.time()[3]
    pos_result <- tryCatch(
      matchMotifs(pfm_chunk, chr_dna, out = "positions", p.cutoff = 5e-5),
      error = function(e) {
        cat(sprintf("    ERROR chunk %d: %s\n", chunk_i, e$message))
        NULL
      }
    )
    scan_time <- proc.time()[3] - t1

    if (is.null(pos_result)) next

    # Extract hits from this chunk
    chunk_rows <- vector("list", length(pfm_chunk))
    for (mi in seq_along(pfm_chunk)) {
      x <- pos_result[[mi]]

      # Handle different return types
      if (is(x, "IRangesList")) {
        ir <- unlist(x)
        if (length(ir) == 0) next
        hit_start <- start(ir)
        hit_end   <- end(ir)
      } else if (is(x, "GRanges")) {
        if (length(x) == 0) next
        hit_start <- start(x)
        hit_end   <- end(x)
      } else if (is(x, "GRangesList")) {
        x <- unlist(x)
        if (length(x) == 0) next
        hit_start <- start(x)
        hit_end   <- end(x)
      } else next

      global_mi <- chunk_start + mi - 1L
      chunk_rows[[mi]] <- data.table(
        chr    = chr_name,
        start  = hit_start,
        end    = hit_end,
        motif_id   = motif_meta$motif_id[global_mi],
        motif_name = motif_meta$motif_name[global_mi],
        tf_class   = motif_meta$tf_class[global_mi]
      )
    }

    chunk_dt <- rbindlist(chunk_rows)
    chr_hits_list[[chunk_i]] <- chunk_dt
    n_chunk_hits <- nrow(chunk_dt)

    if (chunk_i %% 5 == 0 || chunk_i == n_chunks) {
      cat(sprintf("    Chunk %d/%d: motifs %d-%d, %s hits (%.0fs)\n",
                  chunk_i, n_chunks, chunk_start, chunk_end,
                  format(n_chunk_hits, big.mark = ","), scan_time))
      flush.console()
    }

    rm(pos_result, chunk_dt, chunk_rows); gc(verbose = FALSE)
  }

  # Combine all chunks for this chromosome
  chr_dt <- rbindlist(chr_hits_list)
  rm(chr_hits_list, chr_dna); gc(verbose = FALSE)

  if (nrow(chr_dt) > 0) {
    chr_file <- file.path(out_dir, sprintf("_chunk_%s.tsv", chr_name))
    fwrite(chr_dt, chr_file, sep = "\t")
    chunk_files <- c(chunk_files, chr_file)
    total_hits <- total_hits + nrow(chr_dt)
  }

  chr_time <- proc.time()[3] - chr_t0
  cat(sprintf("  %s DONE: %s hits (%.1f min) | Running total: %s\n",
              chr_name, format(nrow(chr_dt), big.mark = ","),
              chr_time / 60, format(total_hits, big.mark = ",")))
  flush.console()
  rm(chr_dt); gc(verbose = FALSE)
}

cat(sprintf("\nTotal scanning: %s hits in %.1f hours\n",
            format(total_hits, big.mark = ","),
            (proc.time()[3] - t0_global) / 3600))

# =============================================================================
# 4. COMBINE ALL CHROMOSOMES
# =============================================================================
cat("\n=== 4. Combining all chromosomes ===\n")

all_hits <- rbindlist(lapply(chunk_files, fread, sep = "\t"))
cat(sprintf("Combined: %s hits\n", format(nrow(all_hits), big.mark = ",")))

# =============================================================================
# 5. ANNOTATE WITH GENOMIC CONTEXT
# =============================================================================
cat("\n=== 5. Annotating with genomic context ===\n")

hit_gr <- GRanges(seqnames = all_hits$chr,
                  ranges = IRanges(start = all_hits$start, end = all_hits$end))

# Region annotation (priority: Promoter > Exon > Intron > Intergenic)
cat("  Overlapping with promoters...\n")
ov_prom <- overlapsAny(hit_gr, promoters)
cat("  Overlapping with exons...\n")
ov_exon <- overlapsAny(hit_gr, exons)
cat("  Overlapping with genes...\n")
ov_gene <- overlapsAny(hit_gr, genes)
cat("  Overlapping with TEs...\n")
ov_te   <- overlapsAny(hit_gr, te_gr)

all_hits[, region := fifelse(ov_prom, "Promoter",
                    fifelse(ov_exon, "Exon",
                    fifelse(ov_gene, "Intron", "Intergenic")))]
all_hits[, in_te := ov_te]

# Nearest gene + distance to TSS
cat("  Computing distance to nearest gene TSS...\n")
hit_mid <- (all_hits$start + all_hits$end) %/% 2L
hit_gr_mid <- GRanges(seqnames = all_hits$chr,
                      ranges = IRanges(start = hit_mid, width = 1))
gene_gr <- GRanges(seqnames = tss_lookup$chr,
                   ranges = IRanges(start = tss_lookup$tss, width = 1),
                   gene_id = tss_lookup$gene_id,
                   gene_strand = tss_lookup$strand)
nearest_idx <- nearest(hit_gr_mid, gene_gr)

all_hits[, nearest_gene := gene_gr$gene_id[nearest_idx]]
all_hits[, nearest_gene_strand := gene_gr$gene_strand[nearest_idx]]
all_hits[, nearest_tss := start(gene_gr)[nearest_idx]]
all_hits[, dist_to_nearest_tss := fifelse(nearest_gene_strand == "+",
                                          hit_mid - nearest_tss,
                                          nearest_tss - hit_mid)]

# Clean up helper columns
all_hits[, c("nearest_gene_strand", "nearest_tss") := NULL]

rm(hit_gr, hit_gr_mid, gene_gr, ov_prom, ov_exon, ov_gene, ov_te, nearest_idx, hit_mid)
gc(verbose = FALSE)

# =============================================================================
# 6. SAVE OUTPUT
# =============================================================================
cat("\n=== 6. Saving output ===\n")

# Full annotated hits
out_file <- file.path(out_dir, "whole_genome_motif_hits.tsv.gz")
fwrite(all_hits, out_file, sep = "\t", compress = "gzip")
cat(sprintf("Saved: %s (%s hits, %.1f MB)\n",
            out_file, format(nrow(all_hits), big.mark = ","),
            file.size(out_file) / 1e6))

# Region summary
region_summary <- all_hits[, .(
  total_hits    = .N,
  unique_motifs = uniqueN(motif_id),
  in_te_pct     = 100 * sum(in_te) / .N,
  mean_dist_tss = mean(abs(dist_to_nearest_tss), na.rm = TRUE)
), by = region]
cat("\nRegion summary:\n")
print(region_summary)
fwrite(region_summary, file.path(out_dir, "region_summary.tsv"), sep = "\t")

# Per-motif summary
motif_summary <- all_hits[, .(
  total_hits  = .N,
  n_promoter  = sum(region == "Promoter"),
  n_exon      = sum(region == "Exon"),
  n_intron    = sum(region == "Intron"),
  n_intergenic = sum(region == "Intergenic"),
  n_in_te     = sum(in_te),
  pct_intergenic = 100 * sum(region == "Intergenic") / .N
), by = .(motif_id, motif_name, tf_class)]
fwrite(motif_summary, file.path(out_dir, "motif_summary_whole_genome.tsv"), sep = "\t")

# TF class summary
class_summary <- all_hits[!is.na(tf_class), .(
  total_hits  = .N,
  unique_motifs = uniqueN(motif_id),
  pct_promoter = 100 * sum(region == "Promoter") / .N,
  pct_intergenic = 100 * sum(region == "Intergenic") / .N,
  pct_in_te = 100 * sum(in_te) / .N
), by = tf_class][order(-total_hits)]
fwrite(class_summary, file.path(out_dir, "tf_class_summary.tsv"), sep = "\t")

# Cleanup chunk files
cat("\nCleaning up chunk files...\n")
file.remove(chunk_files)

# =============================================================================
# 7. VERIFICATION
# =============================================================================
cat("\n=== 7. Verification ===\n")
cat(sprintf("  Total hits:         %s\n", format(nrow(all_hits), big.mark = ",")))
cat(sprintf("  Unique motifs:      %d / %d\n", uniqueN(all_hits$motif_id), nrow(motif_meta)))
cat(sprintf("  Regions:            %s\n", paste(names(table(all_hits$region)), collapse = ", ")))
cat(sprintf("  In TE:              %s (%.1f%%)\n",
            format(sum(all_hits$in_te), big.mark = ","),
            100 * sum(all_hits$in_te) / nrow(all_hits)))
cat(sprintf("  Intergenic:         %s (%.1f%%)\n",
            format(sum(all_hits$region == "Intergenic"), big.mark = ","),
            100 * sum(all_hits$region == "Intergenic") / nrow(all_hits)))

mem_info <- gc()
cat(sprintf("  Peak memory:        %.1f MB\n", sum(mem_info[, 6])))

elapsed_h <- (proc.time()[3] - t0_global) / 3600
cat(sprintf("\n=== Whole-genome scan complete. Total: %.1f hours ===\n", elapsed_h))
cat("End time:", format(Sys.time()), "\n")
