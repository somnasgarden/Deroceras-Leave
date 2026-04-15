#!/usr/bin/env Rscript
# =============================================================================
# Batch 10b per-chromosome worker: Per-read NME for one chromosome
# Called by batch10b_array.slurm with chromosome number as argument
# Output: cluster/results/nme/nme_chrNN.rds
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
chr_num <- as.integer(args[1])
chr_name <- paste0("chr", chr_num)

source("methylation_pipeline/_config.R")

library(data.table)
library(GenomicRanges)
library(Rsamtools)

out_dir <- file.path(PROJECT_DIR, "cluster/results/nme")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, sprintf("nme_%s.rds", chr_name))

cat(sprintf("=== NME worker: %s ===\n", chr_name))
t0 <- proc.time()[3]

# --- BAM files (use sorted versions from bams_sorted/) ---
sorted_dir <- "/mnt/data/alfredvar/rlopezt/bams_sorted"
bam_files <- c(
  C1 = file.path(sorted_dir, "C1_paired_bismark_bt2_pe.deduplicated.sorted.bam"),
  C2 = file.path(sorted_dir, "C2_paired_bismark_bt2_pe.deduplicated.sorted.bam"),
  A1 = file.path(sorted_dir, "A1_paired_bismark_bt2_pe.deduplicated.sorted.bam"),
  A2 = file.path(sorted_dir, "A2_paired_bismark_bt2_pe.deduplicated.sorted.bam")
)

for (nm in names(bam_files)) {
  stopifnot(file.exists(bam_files[nm]))
  # Ensure index exists
  bai <- paste0(bam_files[nm], ".bai")
  bai2 <- sub("\\.bam$", ".bam.bai", bam_files[nm])
  if (!file.exists(bai) && !file.exists(bai2)) {
    cat(sprintf("  Indexing %s...\n", nm))
    indexBam(bam_files[nm])
  }
}

# --- Load CpG positions for this chromosome ---
cat("Loading CpG positions...\n")
bs_obj <- readRDS(CACHE$bsseq)
cpg_gr <- granges(bs_obj)
cpg_dt <- data.table(chr = as.character(seqnames(cpg_gr)), pos = start(cpg_gr))
cpg_dt <- cpg_dt[chr == chr_name]
setkey(cpg_dt, pos)
cat(sprintf("  CpGs on %s: %s\n", chr_name, format(nrow(cpg_dt), big.mark = ",")))
rm(bs_obj, cpg_gr); gc(verbose = FALSE)

# --- Build 4-CpG sliding windows ---
N_CPG <- 4L
MAX_SPAN <- 500L

chr_pos <- cpg_dt$pos
if (length(chr_pos) < N_CPG) {
  cat("  Too few CpGs — skipping\n")
  saveRDS(data.table(), out_file)
  quit("no")
}

n_win <- length(chr_pos) - N_CPG + 1L
starts <- chr_pos[1:n_win]
ends   <- chr_pos[N_CPG:length(chr_pos)]
keep <- (ends - starts + 1L) <= MAX_SPAN

windows <- data.table(
  chr = chr_name,
  win_start = starts[keep],
  win_end   = ends[keep],
  cpg1 = chr_pos[which(keep)],
  cpg2 = chr_pos[which(keep) + 1L],
  cpg3 = chr_pos[which(keep) + 2L],
  cpg4 = chr_pos[which(keep) + 3L]
)
cat(sprintf("  Windows: %s\n", format(nrow(windows), big.mark = ",")))

# --- NME computation helper ---
compute_nme <- function(patterns, n_cpg = 4L) {
  if (length(patterns) < 4) return(NA_real_)
  freq <- table(patterns)
  p <- as.numeric(freq) / sum(freq)
  H <- -sum(p * log2(p))
  H / n_cpg
}

# --- Process each sample ---
for (sample_name in names(bam_files)) {
  cat(sprintf("\n  %s: scanning %d windows...\n", sample_name, nrow(windows)))
  bam <- bam_files[sample_name]
  t_s <- proc.time()[3]

  # Single scanBam call for whole chromosome
  chr_gr <- GRanges(seqnames = chr_name,
                     ranges = IRanges(start = windows$win_start, end = windows$win_end))
  param <- ScanBamParam(which = chr_gr, what = c("pos"), tag = "XM")
  all_reads <- scanBam(bam, param = param)

  nme_vec   <- rep(NA_real_, nrow(windows))
  nread_vec <- rep(0L, nrow(windows))

  for (k in seq_len(nrow(windows))) {
    rk <- all_reads[[k]]
    if (length(rk$pos) == 0) next

    xm_vec <- rk$tag$XM
    pos_vec <- rk$pos
    has_xm <- !is.na(xm_vec) & nchar(xm_vec) > 0
    if (!any(has_xm)) next

    xm_sub <- xm_vec[has_xm]; pos_sub <- pos_vec[has_xm]
    xm_len <- nchar(xm_sub)

    cpg_pos <- c(windows$cpg1[k], windows$cpg2[k], windows$cpg3[k], windows$cpg4[k])
    offset_mat <- outer(pos_sub, cpg_pos, function(p, c) c - p + 1L)
    valid_mask <- rowSums(offset_mat < 1L | offset_mat > xm_len) == 0L
    if (!any(valid_mask)) next

    xm_sub <- xm_sub[valid_mask]; offset_mat <- offset_mat[valid_mask, , drop = FALSE]
    result_mat <- matrix(NA_character_, nrow = nrow(offset_mat), ncol = N_CPG)
    for (ci in seq_len(N_CPG)) {
      result_mat[, ci] <- substr(xm_sub, offset_mat[, ci], offset_mat[, ci])
    }
    binary_mat <- ifelse(result_mat == "Z", "1",
                  ifelse(result_mat == "z", "0", NA_character_))
    complete <- rowSums(is.na(binary_mat)) == 0L
    if (!any(complete)) next

    patterns <- apply(binary_mat[complete, , drop = FALSE], 1, paste0, collapse = "")
    nme_vec[k] <- compute_nme(patterns)
    nread_vec[k] <- sum(complete)
  }

  rm(all_reads); gc(verbose = FALSE)

  windows[, paste0("nme_", sample_name) := nme_vec]
  windows[, paste0("reads_", sample_name) := nread_vec]

  elapsed <- proc.time()[3] - t_s
  n_valid <- sum(!is.na(nme_vec))
  cat(sprintf("  %s done: %d/%d valid (%.0f s)\n", sample_name, n_valid, nrow(windows), elapsed))
}

# --- Compute condition means ---
windows[, ctrl_nme := ifelse(!is.na(nme_C1) & !is.na(nme_C2), (nme_C1 + nme_C2) / 2, NA_real_)]
windows[, ampu_nme := ifelse(!is.na(nme_A1) & !is.na(nme_A2), (nme_A1 + nme_A2) / 2, NA_real_)]
windows[, delta_nme := ampu_nme - ctrl_nme]

n_valid <- sum(!is.na(windows$ctrl_nme) & !is.na(windows$ampu_nme))
cat(sprintf("\n  Valid windows: %d / %d\n", n_valid, nrow(windows)))

# Save
saveRDS(windows, out_file)
cat(sprintf("\nSaved: %s (%.1f MB)\n", out_file, file.size(out_file) / 1e6))
cat(sprintf("=== %s complete (%.1f min) ===\n", chr_name, (proc.time()[3] - t0) / 60))
