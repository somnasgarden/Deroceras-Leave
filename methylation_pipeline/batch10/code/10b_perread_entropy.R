#!/usr/bin/env Rscript
# =============================================================================
# Batch 10b: Per-read Methylation Entropy (Level 2 NME)
# Requires BAM files — run on cluster only (64+ GB RAM)
#
# Method: Xie et al. 2011 (NAR) 4-CpG sliding windows
#   For each window of 4 consecutive CpGs on the same chromosome:
#     1. Extract all reads spanning all 4 CpGs
#     2. Record per-read binary methylation pattern (e.g., 1010)
#     3. Count epiallele frequencies (up to 2^4 = 16 patterns)
#     4. Compute Shannon entropy: H = -sum(p_i * log2(p_i))
#     5. Normalize: NME = H / log2(2^4) = H / 4  (range [0, 1])
#
# Also computes 250bp window MML + NME (Fang et al. 2023 style)
#
# Input: Bismark deduplicated BAMs (XM tag for methylation calls)
# Output: Per-window NME, drift vs reprogramming test, figures
# =============================================================================

source("methylation_pipeline/_config.R")
t0 <- proc.time()

library(data.table)
library(GenomicRanges)
library(Rsamtools)
library(ggplot2)

BATCH_DIR <- file.path(PIPE_DIR, "batch10")

cat("=== Batch 10b: Per-read Methylation Entropy (Level 2 NME) ===\n\n")

# --- BAM file paths ---
# C1 is unsorted in source; the slurm wrapper sorts it into bams_sorted/
bam_dir    <- file.path(dirname(dirname(OG$cpg_C1)), "08_deduplication")
sorted_dir <- "/mnt/data/alfredvar/rlopezt/bams_sorted"
c1_sorted  <- file.path(sorted_dir, "C1_paired_bismark_bt2_pe.deduplicated.sorted.bam")
bam_files <- c(
  C1 = if (file.exists(c1_sorted)) c1_sorted else
       file.path(bam_dir, "C1_paired_bismark_bt2_pe.deduplicated.bam"),
  C2 = file.path(bam_dir, "C2_paired_bismark_bt2_pe.deduplicated.bam"),
  A1 = file.path(bam_dir, "A1_paired_bismark_bt2_pe.deduplicated.bam"),
  A2 = file.path(bam_dir, "A2_paired_bismark_bt2_pe.deduplicated.bam")
)

# Verify BAMs exist
for (nm in names(bam_files)) {
  if (!file.exists(bam_files[nm])) stop(sprintf("BAM not found: %s", bam_files[nm]))
  cat(sprintf("  %s: %s (%.1f GB)\n", nm, bam_files[nm],
              file.size(bam_files[nm]) / 1e9))
}

# --- Load CpG positions from BSseq ---
cat("\n[1/7] Loading CpG positions...\n")
bs_obj <- readRDS(CACHE$bsseq)
cpg_gr <- granges(bs_obj)
cpg_dt <- data.table(
  chr = as.character(seqnames(cpg_gr)),
  pos = start(cpg_gr)
)
cpg_dt <- cpg_dt[chr %in% keep_chr]
setkey(cpg_dt, chr, pos)
cat(sprintf("  CpG sites: %s\n", format(nrow(cpg_dt), big.mark = ",")))

# Also get Level 1 beta values for comparison
ctrl_beta <- rowMeans(bsseq::getMeth(bs_obj[, c("C1","C2")], type = "raw"), na.rm = TRUE)
ampu_beta <- rowMeans(bsseq::getMeth(bs_obj[, c("A1","A2")], type = "raw"), na.rm = TRUE)
rm(bs_obj); gc(verbose = FALSE)


# =============================================================================
# HELPER: Extract per-read methylation from Bismark BAM XM tag
# =============================================================================
# Bismark XM tag encodes methylation calls per base:
#   Z = methylated CpG, z = unmethylated CpG
#   H/h = CHH, X/x = CHG (we only care about Z/z)
#
# For each read overlapping a set of CpG positions, extract the binary
# methylation state at each CpG.

extract_read_patterns <- function(bam_file, window_gr, cpg_positions) {
  # window_gr: GRanges of one window
  # cpg_positions: integer vector of CpG positions within the window
  # Returns: character vector of binary patterns (e.g., "1010")

  param <- ScanBamParam(
    which = window_gr,
    what = c("pos", "seq", "qwidth"),
    tag = "XM"
  )

  reads <- scanBam(bam_file, param = param)[[1]]
  if (length(reads$pos) == 0) return(character(0))

  n_cpg <- length(cpg_positions)
  patterns <- character(length(reads$pos))
  valid <- logical(length(reads$pos))

  # Vectorized extraction across all reads
  xm_vec <- reads$tag$XM
  pos_vec <- reads$pos
  n_cpg <- length(cpg_positions)

  # Filter reads with valid XM tags
  has_xm <- !is.na(xm_vec) & nchar(xm_vec) > 0
  if (!any(has_xm)) return(character(0))

  xm_vec <- xm_vec[has_xm]
  pos_vec <- pos_vec[has_xm]
  xm_len <- nchar(xm_vec)

  # For each read, compute offsets for all CpG positions
  # cpg_offsets matrix: rows = reads, cols = CpG positions
  offset_mat <- outer(pos_vec, cpg_positions, function(p, c) c - p + 1L)

  # Valid reads: all offsets within XM string bounds
  valid_mask <- rowSums(offset_mat < 1L | offset_mat > xm_len) == 0L
  if (!any(valid_mask)) return(character(0))

  xm_vec <- xm_vec[valid_mask]
  offset_mat <- offset_mat[valid_mask, , drop = FALSE]

  # Extract characters at CpG positions using substr (vectorized per CpG)
  result_mat <- matrix(NA_character_, nrow = nrow(offset_mat), ncol = n_cpg)
  for (ci in seq_len(n_cpg)) {
    result_mat[, ci] <- substr(xm_vec, offset_mat[, ci], offset_mat[, ci])
  }

  # Convert to binary: Z=1, z=0, else NA
  binary_mat <- ifelse(result_mat == "Z", "1",
                ifelse(result_mat == "z", "0", NA_character_))

  # Keep only reads with no NAs
  complete <- rowSums(is.na(binary_mat)) == 0L
  if (!any(complete)) return(character(0))

  # Paste columns together for each row
  apply(binary_mat[complete, , drop = FALSE], 1, paste0, collapse = "")
}


# =============================================================================
# HELPER: Compute Shannon entropy from pattern vector
# =============================================================================
compute_nme <- function(patterns, n_cpg) {
  if (length(patterns) < 4) return(list(nme = NA, n_reads = length(patterns), n_patterns = NA))

  freq <- table(patterns)
  p <- as.numeric(freq) / sum(freq)
  H <- -sum(p * log2(p))
  # Normalize by max possible entropy: log2(2^n_cpg) = n_cpg
  nme <- H / n_cpg

  list(nme = nme, n_reads = length(patterns), n_patterns = length(freq))
}


# =============================================================================
# BUILD 4-CpG SLIDING WINDOWS (Xie et al. 2011)
# =============================================================================
cat("[2/7] Building 4-CpG sliding windows...\n")

N_CPG_WINDOW <- 4L
MAX_WINDOW_BP <- 500L  # max span of 4 CpGs (skip if too spread)
MIN_READS <- 4L        # minimum reads covering all 4 CpGs

# Build windows per chromosome
windows_list <- list()
for (chr_name in keep_chr) {
  chr_cpg <- cpg_dt[chr == chr_name]$pos
  if (length(chr_cpg) < N_CPG_WINDOW) next

  n_windows <- length(chr_cpg) - N_CPG_WINDOW + 1L
  starts <- chr_cpg[1:n_windows]
  ends   <- chr_cpg[N_CPG_WINDOW:length(chr_cpg)]
  spans  <- ends - starts + 1L

  # Filter: 4 CpGs must be within MAX_WINDOW_BP of each other
  keep <- spans <= MAX_WINDOW_BP
  if (sum(keep) == 0) next

  # Store window info
  win_dt <- data.table(
    chr = chr_name,
    win_start = starts[keep],
    win_end = ends[keep],
    cpg1 = chr_cpg[which(keep)],
    cpg2 = chr_cpg[which(keep) + 1L],
    cpg3 = chr_cpg[which(keep) + 2L],
    cpg4 = chr_cpg[which(keep) + 3L]
  )
  windows_list[[chr_name]] <- win_dt
}

windows <- rbindlist(windows_list)
cat(sprintf("  4-CpG windows: %s (max span %d bp)\n",
            format(nrow(windows), big.mark = ","), MAX_WINDOW_BP))


# =============================================================================
# COMPUTE PER-READ NME FOR EACH SAMPLE (per-chromosome for memory)
# =============================================================================
cat("[3/7] Computing per-read NME from BAMs...\n")
cat("  This is the slowest step — processing 4 BAMs x %s windows.\n",
    format(nrow(windows), big.mark = ","))

# Process in chunks per chromosome to manage memory
# Store results per sample
nme_results <- list()

for (sample_name in names(bam_files)) {
  cat(sprintf("\n  Processing %s...\n", sample_name))
  bam <- bam_files[sample_name]

  # Index BAM if needed
  bai <- paste0(bam, ".bai")
  if (!file.exists(bai)) {
    bai2 <- sub("\\.bam$", ".bam.bai", bam)
    if (!file.exists(bai2)) {
      cat("    Indexing BAM...\n")
      indexBam(bam)
    }
  }

  sample_nme <- rep(NA_real_, nrow(windows))
  sample_nreads <- rep(0L, nrow(windows))

  for (chr_name in keep_chr) {
    chr_idx <- which(windows$chr == chr_name)
    if (length(chr_idx) == 0) next

    cat(sprintf("    %s: %d windows...", chr_name, length(chr_idx)))
    t_chr <- proc.time()

    # Batch all windows on this chromosome into one scanBam call
    chr_gr <- GRanges(seqnames = windows$chr[chr_idx],
                       ranges = IRanges(start = windows$win_start[chr_idx],
                                        end = windows$win_end[chr_idx]))
    param <- ScanBamParam(which = chr_gr, what = c("pos", "seq", "qwidth"), tag = "XM")
    all_reads <- scanBam(bam, param = param)

    # Process each window's reads
    for (k in seq_along(chr_idx)) {
      j <- chr_idx[k]
      reads_k <- all_reads[[k]]
      if (length(reads_k$pos) == 0) next

      cpg_pos <- c(windows$cpg1[j], windows$cpg2[j], windows$cpg3[j], windows$cpg4[j])

      # Reuse vectorized extraction logic inline (avoid re-reading BAM)
      xm_vec <- reads_k$tag$XM
      pos_vec <- reads_k$pos
      has_xm <- !is.na(xm_vec) & nchar(xm_vec) > 0
      if (!any(has_xm)) next

      xm_sub <- xm_vec[has_xm]; pos_sub <- pos_vec[has_xm]
      xm_len <- nchar(xm_sub)
      offset_mat <- outer(pos_sub, cpg_pos, function(p, c) c - p + 1L)
      valid_mask <- rowSums(offset_mat < 1L | offset_mat > xm_len) == 0L
      if (!any(valid_mask)) next

      xm_sub <- xm_sub[valid_mask]; offset_mat <- offset_mat[valid_mask, , drop = FALSE]
      result_mat <- matrix(NA_character_, nrow = nrow(offset_mat), ncol = N_CPG_WINDOW)
      for (ci in seq_len(N_CPG_WINDOW)) {
        result_mat[, ci] <- substr(xm_sub, offset_mat[, ci], offset_mat[, ci])
      }
      binary_mat <- ifelse(result_mat == "Z", "1",
                    ifelse(result_mat == "z", "0", NA_character_))
      complete <- rowSums(is.na(binary_mat)) == 0L
      if (!any(complete)) next

      patterns <- apply(binary_mat[complete, , drop = FALSE], 1, paste0, collapse = "")
      result <- compute_nme(patterns, N_CPG_WINDOW)
      sample_nme[j] <- result$nme
      sample_nreads[j] <- result$n_reads
    }
    rm(all_reads); gc(verbose = FALSE)

    elapsed_chr <- (proc.time() - t_chr)[3]
    n_valid <- sum(!is.na(sample_nme[chr_idx]))
    cat(sprintf(" %.0fs (%d valid)\n", elapsed_chr, n_valid))
  }

  nme_results[[sample_name]] <- list(nme = sample_nme, n_reads = sample_nreads)
  cat(sprintf("  %s done: %d/%d windows with NME (%.1f%%)\n",
              sample_name, sum(!is.na(sample_nme)), nrow(windows),
              100 * mean(!is.na(sample_nme))))
}


# =============================================================================
# AGGREGATE AND COMPARE CONDITIONS
# =============================================================================
cat("\n[4/7] Aggregating NME across conditions...\n")

windows[, nme_C1 := nme_results$C1$nme]
windows[, nme_C2 := nme_results$C2$nme]
windows[, nme_A1 := nme_results$A1$nme]
windows[, nme_A2 := nme_results$A2$nme]
windows[, reads_C1 := nme_results$C1$n_reads]
windows[, reads_C2 := nme_results$C2$n_reads]
windows[, reads_A1 := nme_results$A1$n_reads]
windows[, reads_A2 := nme_results$A2$n_reads]

# Mean NME per condition (require both replicates)
windows[, ctrl_nme := ifelse(!is.na(nme_C1) & !is.na(nme_C2), (nme_C1 + nme_C2) / 2, NA_real_)]
windows[, ampu_nme := ifelse(!is.na(nme_A1) & !is.na(nme_A2), (nme_A1 + nme_A2) / 2, NA_real_)]
windows[, delta_nme := ampu_nme - ctrl_nme]

valid_windows <- windows[!is.na(ctrl_nme) & !is.na(ampu_nme)]
cat(sprintf("  Valid windows (both conditions): %s / %s\n",
            format(nrow(valid_windows), big.mark = ","),
            format(nrow(windows), big.mark = ",")))

cat(sprintf("  Mean NME — Control: %.4f, Amputated: %.4f\n",
            mean(valid_windows$ctrl_nme), mean(valid_windows$ampu_nme)))
cat(sprintf("  Mean delta NME: %.6f\n", mean(valid_windows$delta_nme)))

# Wilcoxon paired test
wt <- wilcox.test(valid_windows$ctrl_nme, valid_windows$ampu_nme, paired = TRUE)
cat(sprintf("  Wilcoxon paired: p = %s\n", format(wt$p.value, digits = 3)))

# Effect size (Cohen's d)
d <- mean(valid_windows$delta_nme) / sd(valid_windows$delta_nme)
cat(sprintf("  Cohen's d: %.4f\n", d))

save_data(valid_windows, BATCH_DIR, "perread_nme_windows")


# =============================================================================
# ANNOTATE WINDOWS BY REGION + DMP STATUS
# =============================================================================
cat("\n[5/7] Annotating windows...\n")

# Region annotation (midpoint of window)
gff <- load_gff(); genes <- gff[gff$type == "gene"]; exons <- gff[gff$type == "exon"]
promoters <- if (file.exists(CACHE$promoters)) readRDS(CACHE$promoters) else
  GenomicRanges::trim(GenomicRanges::promoters(genes, upstream = 2000, downstream = 0))

win_mid <- (valid_windows$win_start + valid_windows$win_end) %/% 2L
valid_windows[, region := annotate_regions(chr, win_mid, promoters, exons, genes)]

# DMP overlap: does this window contain a DMP?
dmp_file <- file.path(PIPE_DIR, "batch06/data/dmps_annotated.tsv")
if (file.exists(dmp_file)) {
  dmps <- fread(dmp_file)
  dmp_gr <- GRanges(seqnames = dmps$chr, ranges = IRanges(start = dmps$pos, width = 1))
  win_gr <- GRanges(seqnames = valid_windows$chr,
                     ranges = IRanges(start = valid_windows$win_start,
                                      end = valid_windows$win_end))
  valid_windows[, has_dmp := overlapsAny(win_gr, dmp_gr)]
  cat(sprintf("  Windows with DMP: %d (%.2f%%)\n",
              sum(valid_windows$has_dmp),
              100 * mean(valid_windows$has_dmp)))
} else {
  valid_windows[, has_dmp := FALSE]
  cat("  No DMP file found (run batch06 first).\n")
}

# CpG density (number of CpGs per kb in window)
valid_windows[, cpg_density := 4 / ((win_end - win_start + 1) / 1000)]

# NME by region
region_nme <- valid_windows[, .(
  mean_ctrl_nme = mean(ctrl_nme),
  mean_ampu_nme = mean(ampu_nme),
  mean_delta = mean(delta_nme),
  n_windows = .N
), by = region]
cat("\nPer-read NME by region:\n"); print(region_nme)
save_data(region_nme, BATCH_DIR, "perread_nme_by_region")


# =============================================================================
# DRIFT VS REPROGRAMMING TEST
# =============================================================================
cat("\n[6/7] Drift vs reprogramming test...\n")

if (any(valid_windows$has_dmp)) {
  dmp_nme <- valid_windows[has_dmp == TRUE, .(mean_ctrl = mean(ctrl_nme), n = .N)]
  non_dmp_nme <- valid_windows[has_dmp == FALSE, .(mean_ctrl = mean(ctrl_nme), n = .N)]

  cat(sprintf("  DMP windows — baseline NME: %.4f (n = %d)\n",
              dmp_nme$mean_ctrl, dmp_nme$n))
  cat(sprintf("  Non-DMP windows — baseline NME: %.4f (n = %d)\n",
              non_dmp_nme$mean_ctrl, non_dmp_nme$n))

  # Matched comparison: compare DMP windows to non-DMP windows with similar beta
  # Get mean beta for each window from Level 1 data
  cpg_beta_dt <- data.table(
    chr = cpg_dt$chr,
    pos = cpg_dt$pos,
    ctrl_beta = ctrl_beta[1:nrow(cpg_dt)],
    ampu_beta = ampu_beta[1:nrow(cpg_dt)]
  )

  # Mean beta per window (from the 4 CpGs) — vectorized via data.table join
  # Melt window CpG positions to long format for a single join
  win_cpgs <- rbindlist(lapply(seq_len(nrow(valid_windows)), function(i) {
    data.table(win_idx = i, chr = valid_windows$chr[i],
               pos = c(valid_windows$cpg1[i], valid_windows$cpg2[i],
                       valid_windows$cpg3[i], valid_windows$cpg4[i]))
  }))
  setkey(cpg_beta_dt, chr, pos)
  setkey(win_cpgs, chr, pos)
  win_cpgs <- cpg_beta_dt[win_cpgs, on = .(chr, pos), nomatch = NA]
  win_betas <- win_cpgs[, .(win_beta = mean(ctrl_beta, na.rm = TRUE)), by = win_idx]
  valid_windows[, win_beta := NA_real_]
  valid_windows[win_betas$win_idx, win_beta := win_betas$win_beta]

  # Beta-matched comparison: for each DMP window, find non-DMP windows within +/-0.05 beta
  cat("  Beta-matched comparison (within +/-0.05)...\n")
  dmp_win <- valid_windows[has_dmp == TRUE & !is.na(win_beta)]
  nondmp_win <- valid_windows[has_dmp == FALSE & !is.na(win_beta)]

  set.seed(42)
  matched_nondmp_nme <- sapply(seq_len(nrow(dmp_win)), function(i) {
    target_beta <- dmp_win$win_beta[i]
    candidates <- nondmp_win[abs(win_beta - target_beta) <= 0.05]
    if (nrow(candidates) < 5) return(NA)
    mean(sample(candidates$ctrl_nme, min(10, nrow(candidates))))
  })

  matched_valid <- !is.na(matched_nondmp_nme)
  if (sum(matched_valid) > 50) {
    wt_matched <- wilcox.test(dmp_win$ctrl_nme[matched_valid],
                               matched_nondmp_nme[matched_valid], paired = TRUE)
    cat(sprintf("  DMP baseline NME: %.4f\n", mean(dmp_win$ctrl_nme[matched_valid])))
    cat(sprintf("  Matched non-DMP NME: %.4f\n", mean(matched_nondmp_nme[matched_valid], na.rm = TRUE)))
    cat(sprintf("  Wilcoxon (matched): p = %s\n", format(wt_matched$p.value, digits = 3)))

    if (mean(dmp_win$ctrl_nme[matched_valid]) < mean(matched_nondmp_nme[matched_valid], na.rm = TRUE)) {
      cat("  --> DMPs at LOWER baseline entropy than matched non-DMPs\n")
      cat("  --> Suggests DIRECTED REPROGRAMMING (rejecting drift null)\n")
    } else {
      cat("  --> DMPs at HIGHER baseline entropy than matched non-DMPs\n")
      cat("  --> Consistent with STOCHASTIC DRIFT (failing to reject null)\n")
    }
  } else {
    cat("  Insufficient matched windows for comparison.\n")
  }

  save_data(data.frame(
    dmp_baseline_nme = mean(dmp_win$ctrl_nme, na.rm = TRUE),
    nondmp_baseline_nme = mean(nondmp_win$ctrl_nme, na.rm = TRUE),
    matched_nondmp_nme = mean(matched_nondmp_nme, na.rm = TRUE),
    matched_p = if (sum(matched_valid) > 50) wt_matched$p.value else NA,
    n_dmp_windows = nrow(dmp_win),
    n_matched = sum(matched_valid)
  ), BATCH_DIR, "drift_vs_reprogramming")
}


# =============================================================================
# FIGURES
# =============================================================================
cat("\n[7/7] Generating figures...\n")

# (A) NME distribution: control vs amputated
cat("  Fig 10e: per-read NME distribution...\n")
nme_long <- rbind(
  data.table(condition = "Control", nme = valid_windows$ctrl_nme),
  data.table(condition = "Amputated", nme = valid_windows$ampu_nme)
)
p10e <- ggplot(nme_long, aes(x = nme, fill = condition, color = condition)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values = COLORS$condition) +
  scale_color_manual(values = COLORS$condition) +
  labs(x = "Normalized methylation entropy (NME)", y = "Density",
       title = "Per-read methylation entropy: control vs amputated",
       subtitle = sprintf("4-CpG windows | Delta = %.5f | Wilcoxon p = %s | Cohen's d = %.4f",
                          mean(valid_windows$delta_nme), format(wt$p.value, digits = 3), d)) +
  theme_minimal(base_size = 12)
save_fig(p10e, BATCH_DIR, "fig10e_perread_nme_distribution", w = 9, h = 6)

# (B) NME by region
cat("  Fig 10f: per-read NME by region...\n")
p10f <- ggplot(region_nme, aes(x = region, y = mean_delta, fill = region)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = COLORS$region, guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "Mean delta NME (Amputated - Control)",
       title = "Per-read entropy change by genomic region",
       subtitle = "Positive = entropy increase, Negative = entropy decrease") +
  theme_minimal(base_size = 12)
save_fig(p10f, BATCH_DIR, "fig10f_perread_nme_by_region", w = 7, h = 6)

# (C) DMP vs non-DMP baseline NME (if available)
if (any(valid_windows$has_dmp)) {
  cat("  Fig 10g: DMP vs non-DMP baseline NME...\n")
  dmp_plot <- valid_windows[, .(chr, ctrl_nme, has_dmp, win_beta)]
  dmp_plot[, group := ifelse(has_dmp, "DMP window", "Non-DMP window")]

  p10g <- ggplot(dmp_plot[!is.na(ctrl_nme)], aes(x = group, y = ctrl_nme, fill = group)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.15, outlier.size = 0.3) +
    scale_fill_manual(values = c("DMP window" = "#C0392B", "Non-DMP window" = "#2471A3"), guide = "none") +
    labs(x = NULL, y = "Baseline NME (control)",
         title = "Baseline per-read entropy at DMP vs non-DMP windows",
         subtitle = "Drift test: are DMPs from high-entropy (disordered) sites?") +
    theme_minimal(base_size = 12)
  save_fig(p10g, BATCH_DIR, "fig10g_dmp_vs_nondmp_nme", w = 7, h = 6)

  # (D) NME vs CpG density (Fang et al. key finding: inverse relationship)
  cat("  Fig 10h: NME vs CpG density...\n")
  p10h <- ggplot(valid_windows[!is.na(cpg_density) & cpg_density < 100],
                  aes(x = cpg_density, y = ctrl_nme)) +
    geom_hex(bins = 60) +
    scale_fill_viridis_c(trans = "log10") +
    geom_smooth(method = "loess", color = "#C0392B", se = FALSE, linewidth = 1) +
    labs(x = "CpG density (CpGs per kb)", y = "Baseline NME (control)",
         title = "Per-read entropy vs CpG density",
         subtitle = "Fang et al. 2023: entropy inversely related to CpG density") +
    theme_minimal(base_size = 12)
  save_fig(p10h, BATCH_DIR, "fig10h_nme_vs_cpg_density", w = 8, h = 6)
}

# --- Summary ---
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\n=== Batch 10b complete (%.1f minutes) ===\n", elapsed / 60))
cat(sprintf("  Valid windows: %s\n", format(nrow(valid_windows), big.mark = ",")))
cat(sprintf("  Mean NME — Control: %.4f, Amputated: %.4f, Delta: %.6f\n",
            mean(valid_windows$ctrl_nme), mean(valid_windows$ampu_nme),
            mean(valid_windows$delta_nme)))
if (any(valid_windows$has_dmp)) {
  cat(sprintf("  DMP windows: %d\n", sum(valid_windows$has_dmp)))
}
