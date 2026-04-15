#!/usr/bin/env Rscript
# =============================================================================
# Batch 10b merge: Combine per-chromosome NME results + produce figures
# Run after all batch10b_array jobs complete
# =============================================================================

source("methylation_pipeline/_config.R")
t0 <- proc.time()[3]

library(data.table)
library(GenomicRanges)
library(ggplot2)

BATCH_DIR <- file.path(PIPE_DIR, "batch10")
nme_dir <- file.path(PROJECT_DIR, "cluster/results/nme")

cat("=== Batch 10b: Merging per-chromosome NME results ===\n\n")

# --- Load and combine all chromosome results ---
rds_files <- list.files(nme_dir, pattern = "^nme_chr.*\\.rds$", full.names = TRUE)
cat(sprintf("Found %d chromosome files\n", length(rds_files)))
stopifnot(length(rds_files) >= 20)  # sanity check

windows <- rbindlist(lapply(rds_files, readRDS))
cat(sprintf("Total windows: %s\n", format(nrow(windows), big.mark = ",")))

# Filter to valid windows (both conditions have NME)
valid <- windows[!is.na(ctrl_nme) & !is.na(ampu_nme)]
cat(sprintf("Valid windows: %s (%.1f%%)\n",
            format(nrow(valid), big.mark = ","),
            100 * nrow(valid) / nrow(windows)))

cat(sprintf("Mean NME — Control: %.4f, Amputated: %.4f, Delta: %.6f\n",
            mean(valid$ctrl_nme), mean(valid$ampu_nme), mean(valid$delta_nme)))

# Wilcoxon paired test
wt <- wilcox.test(valid$ctrl_nme, valid$ampu_nme, paired = TRUE)
d <- mean(valid$delta_nme) / sd(valid$delta_nme)
cat(sprintf("Wilcoxon p = %s, Cohen's d = %.4f\n", format(wt$p.value, digits = 3), d))

# --- Annotate ---
cat("\nAnnotating windows...\n")
gff <- load_gff(); genes <- gff[gff$type == "gene"]; exons <- gff[gff$type == "exon"]
promoters_gr <- if (file.exists(CACHE$promoters)) readRDS(CACHE$promoters) else
  trim(promoters(genes, upstream = 2000, downstream = 0))

win_mid <- (valid$win_start + valid$win_end) %/% 2L
valid[, region := annotate_regions(chr, win_mid, promoters_gr, exons, genes)]

# DMP overlap
dmp_file <- file.path(PIPE_DIR, "batch06/data/dmps_annotated.tsv")
if (file.exists(dmp_file)) {
  dmps <- unique(fread(dmp_file)[, .(chr, pos)])
  dmp_gr <- GRanges(seqnames = dmps$chr, ranges = IRanges(start = dmps$pos, width = 1))
  win_gr <- GRanges(seqnames = valid$chr,
                     ranges = IRanges(start = valid$win_start, end = valid$win_end))
  valid[, has_dmp := overlapsAny(win_gr, dmp_gr)]
  cat(sprintf("Windows with DMP: %d (%.2f%%)\n", sum(valid$has_dmp), 100 * mean(valid$has_dmp)))
}

valid[, cpg_density := 4 / ((win_end - win_start + 1) / 1000)]

# --- Save main output ---
save_data(valid, BATCH_DIR, "perread_nme_windows")

# Region NME
region_nme <- valid[, .(mean_ctrl = mean(ctrl_nme), mean_ampu = mean(ampu_nme),
                         mean_delta = mean(delta_nme), n = .N), by = region]
save_data(region_nme, BATCH_DIR, "perread_nme_by_region")

# --- Drift vs reprogramming ---
if ("has_dmp" %in% names(valid) && any(valid$has_dmp)) {
  cat("\nDrift vs reprogramming test...\n")
  dmp_base <- mean(valid[has_dmp == TRUE]$ctrl_nme)
  non_base <- mean(valid[has_dmp == FALSE]$ctrl_nme)
  wt_drift <- wilcox.test(valid[has_dmp == TRUE]$ctrl_nme,
                           valid[has_dmp == FALSE]$ctrl_nme)
  cat(sprintf("  DMP baseline NME: %.4f | Non-DMP: %.4f | p = %s\n",
              dmp_base, non_base, format(wt_drift$p.value, digits = 3)))
  if (dmp_base < non_base) {
    cat("  --> DIRECTED REPROGRAMMING (DMPs at lower entropy baseline)\n")
  } else {
    cat("  --> Consistent with STOCHASTIC DRIFT (DMPs at higher entropy)\n")
  }
  save_data(data.table(dmp_baseline_nme = dmp_base, nondmp_baseline_nme = non_base,
                        wilcox_p = wt_drift$p.value, n_dmp = sum(valid$has_dmp),
                        n_nondmp = sum(!valid$has_dmp)),
            BATCH_DIR, "drift_vs_reprogramming")
}

# --- Figures ---
cat("\nGenerating figures...\n")

# NME distribution
nme_long <- rbind(
  data.table(condition = "Control", nme = valid$ctrl_nme),
  data.table(condition = "Amputated", nme = valid$ampu_nme)
)
p1 <- ggplot(nme_long, aes(x = nme, fill = condition, color = condition)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values = COLORS$condition) +
  scale_color_manual(values = COLORS$condition) +
  labs(x = "Normalized methylation entropy (NME)", y = "Density",
       title = "Per-read methylation entropy: control vs amputated",
       subtitle = sprintf("4-CpG windows | Delta = %.5f | p = %s | d = %.4f",
                          mean(valid$delta_nme), format(wt$p.value, digits = 3), d)) +
  theme_minimal(base_size = 12)
save_fig(p1, BATCH_DIR, "fig10b_perread_nme_distribution", w = 9, h = 6)

# NME by region
p2 <- ggplot(region_nme, aes(x = region, y = mean_delta, fill = region)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = COLORS$region, guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "Mean delta NME",
       title = "Per-read entropy change by genomic region") +
  theme_minimal(base_size = 12)
save_fig(p2, BATCH_DIR, "fig10b_nme_by_region", w = 7, h = 6)

# DMP vs non-DMP
if ("has_dmp" %in% names(valid) && any(valid$has_dmp)) {
  valid[, group := ifelse(has_dmp, "DMP window", "Non-DMP window")]
  p3 <- ggplot(valid, aes(x = group, y = ctrl_nme, fill = group)) +
    geom_violin(alpha = 0.5) +
    geom_boxplot(width = 0.15, outlier.size = 0.3) +
    scale_fill_manual(values = c("DMP window" = "#C0392B", "Non-DMP window" = "#2471A3"), guide = "none") +
    labs(x = NULL, y = "Baseline NME (control)",
         title = "Baseline per-read entropy at DMP vs non-DMP windows") +
    theme_minimal(base_size = 12)
  save_fig(p3, BATCH_DIR, "fig10b_dmp_vs_nondmp_nme", w = 7, h = 6)
}

elapsed <- (proc.time()[3] - t0) / 60
cat(sprintf("\n=== Batch 10b merge complete (%.1f min) ===\n", elapsed))
