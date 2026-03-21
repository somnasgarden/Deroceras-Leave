#!/usr/bin/env Rscript
# =============================================================================
# Methylation Tutorial Part 3: From WGBS Methylation Calls to Biological Insights
# Adapted from NGS101 tutorial for D. laeve (non-model invertebrate)
#
# FULL RUN — all 25.4M CpGs, no subsampling, smoothing enabled
# Native Windows R — no WSL, no RAM limits
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

# --- 1. LIBRARIES ---
suppressPackageStartupMessages({
  library(bsseq)
  library(DSS)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(gridExtra)
})

cat("=== Part 3: WGBS Downstream Analysis for D. laeve ===\n")
cat("FULL RUN — all CpGs, smoothing enabled, native Windows R\n\n")

# --- 2. PATHS AND PARAMETERS (Windows native) ---
data_dir    <- "C:/Users/rafae/Projects/DATA"
project_dir <- "C:/Users/rafae/Projects/STANDBY"
cache_dir   <- file.path(project_dir, "genome/cache")
out_dir     <- file.path(project_dir, "methylation_tutorials/results")
fig_pdf_dir <- file.path(out_dir, "pdf")
fig_png_dir <- file.path(out_dir, "png")

dir.create(fig_pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_png_dir, recursive = TRUE, showWarnings = FALSE)

keep_chr <- paste0("chr", 1:31)

analysis_params <- list(
  min_coverage   = 5,
  min_diff       = 0.1,
  fdr_threshold  = 0.05,
  min_dmr_cpgs   = 3,
  min_dmr_length = 50
)

sample_info <- data.frame(
  sample_id = c("C1", "C2", "A1", "A2"),
  condition = c("Control", "Control", "Amputated", "Amputated"),
  file_path = file.path(data_dir, c(
    "C1.CpG_report.txt", "C2.CpG_report.txt",
    "A1.CpG_report.txt", "A2.CpG_report.txt"
  ))
)

cat("Samples:\n")
print(sample_info[, c("sample_id", "condition")])
cat("\n")

# Save helper (Windows native — no PDF lock workaround needed with cairo)
save_plot <- function(p, name, w = 8, h = 5) {
  png(file.path(fig_png_dir, paste0(name, ".png")),
      width = w, height = h, units = "in", res = 300)
  print(p); dev.off()
  tryCatch({
    cairo_pdf(file.path(fig_pdf_dir, paste0(name, ".pdf")),
              width = w, height = h)
    print(p); dev.off()
  }, error = function(e) {
    pdf(file.path(fig_pdf_dir, paste0(name, ".pdf")),
        width = w, height = h)
    print(p); dev.off()
  })
  cat("  Saved:", name, "\n")
}


# =============================================================================
# 3. LOADING METHYLATION DATA — FULL DATASET
# =============================================================================
cat("=== 3. Loading methylation data ===\n")

bsseq_cache <- file.path(cache_dir, "bsseq_tutorial.rds")

if (file.exists(bsseq_cache)) {
  cat("Loading FULL cached BSseq object...\n")
  bs_obj <- readRDS(bsseq_cache)
  pData(bs_obj)$condition <- sample_info$condition
  pData(bs_obj)$sample_id <- sample_info$sample_id
  cat(sprintf("BSseq loaded: %s sites x %d samples\n",
              format(nrow(bs_obj), big.mark = ","), ncol(bs_obj)))
} else {
  cat("No cache found. Building BSseq from CpG reports...\n\n")
  cat("Reading CpG reports with data.table (native Windows — no shell pipeline)...\n\n")

  # Windows-native reader: use data.table::fread directly on the 3 GB files
  # CpG_report.txt: chr, pos, strand, count_meth, count_unmeth, context, trinuc (no header)
  read_cpg_report <- function(file_path, sample_name, min_cov = 5) {
    cat(sprintf("  Reading %s (%s)...\n", sample_name, file_path))
    t0 <- proc.time()

    # Read only needed columns, filter as we go
    dt <- fread(file_path, header = FALSE, sep = "\t",
                select = c(1, 2, 3, 4, 5, 6),
                col.names = c("chr", "pos", "strand", "meth", "unmeth", "context"))

    cat(sprintf("    Raw rows: %s\n", format(nrow(dt), big.mark = ",")))

    # Filter: CG context only, chr1-31 only
    dt <- dt[context == "CG" & chr %in% keep_chr]
    cat(sprintf("    After CG + chr filter: %s\n", format(nrow(dt), big.mark = ",")))

    # Strand collapse: shift minus-strand positions back by 1
    dt[strand == "-", pos := pos - 1L]

    # Aggregate by chr+pos (collapse strands)
    dt <- dt[, .(meth = sum(meth), unmeth = sum(unmeth)), by = .(chr, pos)]

    # Coverage filter
    dt <- dt[meth + unmeth >= min_cov]
    cat(sprintf("    After strand collapse + cov>=%d: %s rows\n",
                min_cov, format(nrow(dt), big.mark = ",")))

    elapsed <- (proc.time() - t0)[3]
    cat(sprintf("    Done in %.0f seconds\n\n", elapsed))
    return(dt)
  }

  # PASS 1: Read each sample, save to disk, find common sites
  sample_names <- sample_info$sample_id
  common_dt <- NULL
  tmp_dir <- tempdir()

  for (i in seq_len(nrow(sample_info))) {
    dt <- read_cpg_report(sample_info$file_path[i],
                          sample_names[i],
                          analysis_params$min_coverage)
    tmp_rds <- file.path(tmp_dir, paste0("cpg_", sample_names[i], ".rds"))
    saveRDS(dt, tmp_rds)
    cat(sprintf("    Saved to disk: %s (%.0f MB)\n",
                tmp_rds, file.size(tmp_rds) / 1e6))

    site_dt <- dt[, .(chr, pos)]
    if (is.null(common_dt)) {
      common_dt <- site_dt
    } else {
      setkeyv(common_dt, c("chr", "pos"))
      setkeyv(site_dt, c("chr", "pos"))
      common_dt <- fintersect(common_dt, site_dt)
    }
    cat(sprintf("    Running common sites: %s\n\n",
                format(nrow(common_dt), big.mark = ",")))

    rm(dt, site_dt); gc(verbose = FALSE)
  }

  cat(sprintf("Common sites with coverage >= %d in all %d samples: %s\n\n",
              analysis_params$min_coverage, length(sample_names),
              format(nrow(common_dt), big.mark = ",")))

  # PASS 2: Build BSseq from common sites
  cat("Creating BSseq object from common sites...\n")
  setkeyv(common_dt, c("chr", "pos"))
  n_sites <- nrow(common_dt)
  M_matrix <- matrix(0L, nrow = n_sites, ncol = length(sample_names))
  Cov_matrix <- matrix(0L, nrow = n_sites, ncol = length(sample_names))
  colnames(M_matrix) <- colnames(Cov_matrix) <- sample_names
  gr <- NULL

  for (i in seq_along(sample_names)) {
    cat(sprintf("  Loading %s from disk...\n", sample_names[i]))
    tmp_rds <- file.path(tmp_dir, paste0("cpg_", sample_names[i], ".rds"))
    dt <- readRDS(tmp_rds)
    setkeyv(dt, c("chr", "pos"))
    dt <- dt[common_dt, nomatch = 0L]

    if (is.null(gr)) {
      gr <- GRanges(seqnames = dt$chr,
                    ranges = IRanges(start = dt$pos, width = 1))
    }

    M_matrix[, i] <- dt$meth
    Cov_matrix[, i] <- dt$meth + dt$unmeth

    rm(dt); gc(verbose = FALSE)
    file.remove(tmp_rds)
    cat("    Done, temp file removed.\n")
  }
  rm(common_dt); gc(verbose = FALSE)

  bs_obj <- BSseq(gr = gr, M = M_matrix, Cov = Cov_matrix)
  pData(bs_obj)$condition <- sample_info$condition
  pData(bs_obj)$sample_id <- sample_info$sample_id

  rm(gr, M_matrix, Cov_matrix); gc(verbose = FALSE)

  cat(sprintf("BSseq object created: %s CpG sites x %d samples\n",
              format(nrow(bs_obj), big.mark = ","), ncol(bs_obj)))

  cat("Caching BSseq object...\n")
  saveRDS(bs_obj, bsseq_cache)
  cat(sprintf("Cached to: %s (%.1f MB)\n\n",
              bsseq_cache, file.size(bsseq_cache) / 1e6))
}

cat(sprintf("\n*** FULL DATASET: %s CpG sites x %d samples — NO SUBSAMPLING ***\n\n",
            format(nrow(bs_obj), big.mark = ","), ncol(bs_obj)))


# =============================================================================
# 4. QUALITY CONTROL
# =============================================================================
cat("\n=== 4. Quality Control ===\n")

cat("Computing global methylation statistics...\n")
stats_list <- list()
for (i in 1:ncol(bs_obj)) {
  bv <- getMeth(bs_obj[, i], type = "raw")[, 1]
  bv <- bv[!is.nan(bv) & !is.na(bv)]
  stats_list[[i]] <- data.frame(
    sample    = colnames(bs_obj)[i],
    condition = pData(bs_obj)$condition[i],
    total_sites = length(bv),
    mean_methylation   = round(mean(bv), 4),
    median_methylation = round(median(bv), 4),
    sd_methylation     = round(sd(bv), 4),
    sites_high_meth   = sum(bv > 0.8),
    sites_low_meth    = sum(bv < 0.2),
    sites_intermediate = sum(bv >= 0.2 & bv <= 0.8)
  )
  rm(bv); gc(verbose = FALSE)
}
stats_df <- do.call(rbind, stats_list)
rm(stats_list)

cat("\nGlobal methylation statistics:\n")
print(stats_df)
write.table(stats_df, file.path(out_dir, "global_methylation_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Correlation — subsample 500K for speed
cat("\nComputing sample correlations (500K random sites)...\n")
set.seed(42)
n_cor <- min(500000, nrow(bs_obj))
cor_idx <- sort(sample(nrow(bs_obj), n_cor))
beta_sub <- getMeth(bs_obj[cor_idx, ], type = "raw")
complete_sites <- complete.cases(beta_sub)
beta_complete <- beta_sub[complete_sites, ]
rm(beta_sub, cor_idx); gc(verbose = FALSE)

cat(sprintf("Using %s complete sites for correlation\n",
            format(nrow(beta_complete), big.mark = ",")))
sample_cors <- cor(beta_complete)
cat("Sample correlations:\n")
print(round(sample_cors, 4))

# Correlation heatmap
annotation_df <- data.frame(
  Condition = pData(bs_obj)$condition,
  row.names = colnames(sample_cors)
)
ann_colors <- list(Condition = c("Control" = "#2471A3", "Amputated" = "#C0392B"))

p_cor <- pheatmap(sample_cors,
                  annotation_col = annotation_df,
                  annotation_row = annotation_df,
                  annotation_colors = ann_colors,
                  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                  main = "Sample-to-Sample Methylation Correlations",
                  fontsize = 12,
                  display_numbers = TRUE,
                  number_format = "%.4f",
                  silent = TRUE)

png(file.path(fig_png_dir, "sample_correlation_heatmap.png"),
    width = 8, height = 6, units = "in", res = 300)
grid::grid.draw(p_cor$gtable)
dev.off()
tryCatch({
  cairo_pdf(file.path(fig_pdf_dir, "sample_correlation_heatmap.pdf"), width = 8, height = 6)
  grid::grid.draw(p_cor$gtable); dev.off()
}, error = function(e) {
  pdf(file.path(fig_pdf_dir, "sample_correlation_heatmap.pdf"), width = 8, height = 6)
  grid::grid.draw(p_cor$gtable); dev.off()
})
cat("  Saved: sample_correlation_heatmap\n")

# PCA
cat("\nPerforming PCA...\n")
pca_result <- prcomp(t(beta_complete), center = TRUE, scale. = FALSE)
var_explained <- summary(pca_result)$importance[2, 1:min(4, ncol(pca_result$x))] * 100

pca_data <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Sample = rownames(pca_result$x),
  Condition = pData(bs_obj)$condition
)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 5) +
  geom_text(vjust = -1.2, hjust = 0.5, size = 4, show.legend = FALSE) +
  scale_color_manual(values = c("Control" = "#2471A3", "Amputated" = "#C0392B")) +
  labs(
    title = "PCA of DNA Methylation Profiles",
    subtitle = sprintf("D. laeve WGBS | %s CpG sites | PC1: %.1f%%, PC2: %.1f%%",
                       format(nrow(beta_complete), big.mark = ","),
                       var_explained[1], var_explained[2]),
    x = sprintf("PC1 (%.1f%%)", var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", var_explained[2])
  ) +
  theme_minimal(base_size = 13)

save_plot(p_pca, "pca_methylation", w = 8, h = 6)

rm(beta_complete, pca_result, sample_cors, p_cor, p_pca, pca_data, annotation_df, ann_colors)
gc(verbose = FALSE)
cat("Freed QC objects for DMLtest.\n")


# =============================================================================
# 5. DIFFERENTIAL METHYLATION — FULL GENOME, WITH SMOOTHING
# =============================================================================
cat("\n=== 5. Differential Methylation Analysis (DSS) ===\n")

group1 <- c("C1", "C2")
group2 <- c("A1", "A2")

cat(sprintf("Testing: %s vs %s\n", paste(group1, collapse = "+"), paste(group2, collapse = "+")))
cat(sprintf("BSseq object: %s sites — FULL GENOME, NO SUBSAMPLING\n",
            format(nrow(bs_obj), big.mark = ",")))
cat("Running DMLtest WITH smoothing (publication quality)...\n")
cat("Start time:", format(Sys.time()), "\n\n")

t0 <- proc.time()
dml_test <- DMLtest(bs_obj, group1 = group1, group2 = group2, smoothing = TRUE)
dml_time <- (proc.time() - t0)[3]
cat(sprintf("\nDMLtest completed in %.0f seconds (%.1f minutes)\n", dml_time, dml_time / 60))
cat(sprintf("Sites tested: %s\n", format(nrow(dml_test), big.mark = ",")))

# Cache DMLtest results
dml_cache <- file.path(cache_dir, "dmltest_full.rds")
saveRDS(dml_test, dml_cache)
cat(sprintf("DMLtest cached: %s (%.1f MB)\n", dml_cache, file.size(dml_cache) / 1e6))

# Free BSseq
rm(bs_obj); gc(verbose = FALSE)

# Call DMPs
cat("\nCalling DMPs (FDR < 0.05, |diff| > 10%)...\n")
dmp_results <- callDML(dml_test,
                       p.threshold = analysis_params$fdr_threshold,
                       delta = analysis_params$min_diff)

if (is.null(dmp_results) || nrow(dmp_results) == 0) {
  dmp_results <- data.frame(chr = character(), pos = integer(), mu1 = numeric(),
                            mu2 = numeric(), diff = numeric(), diff.se = numeric(),
                            stat = numeric(), phi1 = numeric(), phi2 = numeric(),
                            pval = numeric(), fdr = numeric())
  cat("No significant DMPs found.\n")
} else {
  cat(sprintf("Significant DMPs: %s\n", format(nrow(dmp_results), big.mark = ",")))
}

# Relaxed threshold for exploration
dmp_relaxed <- callDML(dml_test, p.threshold = 0.01)
cat(sprintf("DMPs at p < 0.01 (no delta filter): %s\n",
            format(nrow(dmp_relaxed), big.mark = ",")))
write.table(dmp_relaxed, file.path(out_dir, "dmps_relaxed_p01.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
rm(dmp_relaxed); gc(verbose = FALSE)

# Call DMRs
cat("\nCalling DMRs...\n")
dmr_results <- callDMR(dml_test,
                       p.threshold = analysis_params$fdr_threshold,
                       delta = analysis_params$min_diff,
                       minlen = analysis_params$min_dmr_length,
                       minCG = analysis_params$min_dmr_cpgs)

if (is.null(dmr_results)) {
  dmr_results <- data.frame(chr = character(), start = integer(), end = integer(),
                            length = integer(), nCG = integer(), areaStat = numeric())
  cat("No DMRs found.\n")
} else {
  cat(sprintf("Significant DMRs: %s\n", format(nrow(dmr_results), big.mark = ",")))
}

# Summary
n_hyper <- if (nrow(dmp_results) > 0) sum(dmp_results$diff > 0) else 0
n_hypo  <- if (nrow(dmp_results) > 0) sum(dmp_results$diff < 0) else 0
cat(sprintf("\nDMP breakdown:\n"))
cat(sprintf("  Hypermethylated (Amputated > Control): %s\n", format(n_hyper, big.mark = ",")))
cat(sprintf("  Hypomethylated  (Amputated < Control): %s\n", format(n_hypo, big.mark = ",")))

n_hyper_dmr <- 0; n_hypo_dmr <- 0
if (nrow(dmr_results) > 0) {
  n_hyper_dmr <- sum(dmr_results$areaStat > 0)
  n_hypo_dmr  <- sum(dmr_results$areaStat < 0)
  cat(sprintf("\nDMR breakdown:\n"))
  cat(sprintf("  Hypermethylated DMRs: %s\n", format(n_hyper_dmr, big.mark = ",")))
  cat(sprintf("  Hypomethylated DMRs:  %s\n", format(n_hypo_dmr, big.mark = ",")))
  cat(sprintf("  Mean DMR length: %.0f bp\n", mean(dmr_results$length)))
  cat(sprintf("  Mean CpGs per DMR: %.1f\n", mean(dmr_results$nCG)))
}

# Save results
write.table(dmp_results, file.path(out_dir, "differentially_methylated_positions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dmr_results, file.path(out_dir, "differentially_methylated_regions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_stats <- data.frame(
  Metric = c("Total CpG sites analyzed",
             "Significant DMPs (FDR < 0.05, |diff| > 10%)",
             "Significant DMRs",
             "Hypermethylated DMPs",
             "Hypomethylated DMPs",
             "Hypermethylated DMRs",
             "Hypomethylated DMRs",
             "Mean DMR length (bp)",
             "Mean CpGs per DMR",
             "DMLtest smoothing",
             "Analysis time (minutes)"),
  Value = c(
    format(nrow(dml_test), big.mark = ","),
    format(nrow(dmp_results), big.mark = ","),
    format(nrow(dmr_results), big.mark = ","),
    format(n_hyper, big.mark = ","),
    format(n_hypo, big.mark = ","),
    format(n_hyper_dmr, big.mark = ","),
    format(n_hypo_dmr, big.mark = ","),
    if (nrow(dmr_results) > 0) round(mean(dmr_results$length)) else "N/A",
    if (nrow(dmr_results) > 0) round(mean(dmr_results$nCG), 1) else "N/A",
    "TRUE",
    round(dml_time / 60, 1)
  )
)
write.table(summary_stats, file.path(out_dir, "analysis_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

rm(dml_test); gc(verbose = FALSE)

# Reload BSseq for visualization
cat("\nReloading BSseq for visualization...\n")
bs_obj <- readRDS(bsseq_cache)
pData(bs_obj)$condition <- sample_info$condition
pData(bs_obj)$sample_id <- sample_info$sample_id


# =============================================================================
# 6. VISUALIZATION
# =============================================================================
cat("\n=== 6. Visualization ===\n")

# 6.1 Volcano plot
if (nrow(dmp_results) > 0) {
  cat("Creating volcano plot...\n")
  volcano_data <- data.frame(
    diff   = dmp_results$diff,
    pvalue = dmp_results$pval,
    fdr    = dmp_results$fdr
  )
  volcano_data$category <- "Not Significant"
  volcano_data$category[volcano_data$fdr < analysis_params$fdr_threshold &
                         volcano_data$diff > analysis_params$min_diff] <- "Hypermethylated"
  volcano_data$category[volcano_data$fdr < analysis_params$fdr_threshold &
                         volcano_data$diff < -analysis_params$min_diff] <- "Hypomethylated"
  volcano_data$category <- factor(volcano_data$category,
                                  levels = c("Not Significant", "Hypermethylated", "Hypomethylated"))

  p_volcano <- ggplot(volcano_data, aes(x = diff, y = -log10(pvalue), color = category)) +
    geom_point(alpha = 0.5, size = 0.6) +
    scale_color_manual(
      values = c("Not Significant" = "gray70",
                 "Hypermethylated" = "#C0392B",
                 "Hypomethylated" = "#2471A3"),
      name = "Methylation Change"
    ) +
    geom_vline(xintercept = c(-analysis_params$min_diff, analysis_params$min_diff),
               linetype = "dashed", alpha = 0.7, color = "gray40") +
    geom_hline(yintercept = -log10(0.05),
               linetype = "dashed", alpha = 0.7, color = "gray40") +
    labs(
      title = "Volcano Plot -- Differentially Methylated Positions",
      subtitle = sprintf("D. laeve | Amputated vs Control | Hyper: %s | Hypo: %s",
                         format(n_hyper, big.mark = ","), format(n_hypo, big.mark = ",")),
      x = "Methylation Difference (Amputated - Control)",
      y = expression(-log[10](p-value))
    ) +
    theme_minimal(base_size = 12)

  save_plot(p_volcano, "volcano_plot_dmps", w = 10, h = 7)
  rm(volcano_data, p_volcano); gc(verbose = FALSE)
}

# 6.2 Genome-wide methylation profile
cat("Creating genome-wide methylation profile...\n")
set.seed(42)
plot_idx <- sort(sample(nrow(bs_obj), min(500000, nrow(bs_obj))))
window_size <- 1000000

genome_wide_plots <- list()
for (samp_i in 1:ncol(bs_obj)) {
  samp <- colnames(bs_obj)[samp_i]
  bv <- getMeth(bs_obj[plot_idx, samp_i], type = "raw")[, 1]

  meth_data <- data.table(
    chr = as.character(seqnames(bs_obj[plot_idx])),
    pos = start(bs_obj[plot_idx]),
    methylation = bv * 100
  )
  meth_data <- meth_data[!is.na(methylation) & chr %in% keep_chr]
  meth_data[, chr_num := as.integer(sub("chr", "", chr))]

  meth_binned <- meth_data[, .(
    mean_meth = mean(methylation, na.rm = TRUE),
    n_cpgs = .N
  ), by = .(chr, chr_num, bin = floor(pos / window_size) * window_size)]
  meth_binned <- meth_binned[n_cpgs >= 10]

  chr_lengths <- meth_binned[, .(max_pos = max(bin)), by = .(chr, chr_num)]
  chr_lengths <- chr_lengths[order(chr_num)]
  chr_lengths[, cumulative := cumsum(c(0, head(max_pos, -1)))]
  chr_lengths[, center := cumulative + max_pos / 2]

  meth_plot <- merge(meth_binned, chr_lengths[, .(chr, cumulative)], by = "chr")
  meth_plot[, genome_pos := cumulative + bin]

  genome_wide_plots[[samp]] <- ggplot(meth_plot,
    aes(x = genome_pos, y = mean_meth, color = factor(chr_num %% 2))) +
    geom_point(alpha = 0.5, size = 0.3) +
    scale_color_manual(values = c("0" = "#2471A3", "1" = "#C0392B"), guide = "none") +
    scale_x_continuous(breaks = chr_lengths$center,
                       labels = chr_lengths$chr_num, expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 100), expand = c(0.02, 0)) +
    labs(title = paste("Genome-wide Methylation:", samp),
         subtitle = "1 Mb windows | D. laeve chr1-31",
         x = "Chromosome", y = "Mean Methylation (%)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())

  rm(bv, meth_data, meth_binned, meth_plot, chr_lengths)
}
gc(verbose = FALSE)

p_genome <- gridExtra::arrangeGrob(grobs = genome_wide_plots, ncol = 1)
png(file.path(fig_png_dir, "genome_wide_methylation.png"),
    width = 16, height = 12, units = "in", res = 300)
grid::grid.draw(p_genome)
dev.off()
tryCatch({
  cairo_pdf(file.path(fig_pdf_dir, "genome_wide_methylation.pdf"), width = 16, height = 12)
  grid::grid.draw(p_genome); dev.off()
}, error = function(e) {
  pdf(file.path(fig_pdf_dir, "genome_wide_methylation.pdf"), width = 16, height = 12)
  grid::grid.draw(p_genome); dev.off()
})
cat("  Saved: genome_wide_methylation\n")
rm(genome_wide_plots, p_genome); gc(verbose = FALSE)

# 6.3 Manhattan plot
if (nrow(dmp_results) > 0) {
  cat("Creating Manhattan plot...\n")
  manhattan_data <- data.table(
    chr = dmp_results$chr,
    pos = dmp_results$pos,
    diff = dmp_results$diff,
    fdr = dmp_results$fdr
  )
  manhattan_data[, chr_num := as.integer(sub("chr", "", chr))]
  manhattan_data <- manhattan_data[!is.na(chr_num)]
  manhattan_data <- manhattan_data[order(chr_num, pos)]

  chr_info <- manhattan_data[, .(max_pos = max(pos)), by = .(chr, chr_num)]
  chr_info <- chr_info[order(chr_num)]
  chr_info[, cumulative := cumsum(c(0, head(max_pos, -1)))]
  chr_info[, center := cumulative + max_pos / 2]

  manhattan_data <- merge(manhattan_data, chr_info[, .(chr, cumulative)], by = "chr")
  manhattan_data[, genome_pos := cumulative + pos]

  p_manhattan <- ggplot(manhattan_data,
    aes(x = genome_pos, y = diff, color = factor(chr_num %% 2))) +
    geom_point(alpha = 0.4, size = 0.4) +
    scale_color_manual(values = c("0" = "#2471A3", "1" = "#C0392B"), guide = "none") +
    scale_x_continuous(breaks = chr_info$center, labels = chr_info$chr_num,
                       expand = c(0.01, 0)) +
    geom_hline(yintercept = c(-analysis_params$min_diff, analysis_params$min_diff),
               linetype = "dashed", alpha = 0.5, color = "gray40") +
    labs(
      title = "Manhattan Plot -- Methylation Differences Across Genome",
      subtitle = sprintf("D. laeve | %s DMPs | Amputated vs Control",
                         format(nrow(dmp_results), big.mark = ",")),
      x = "Chromosome", y = "Methylation Difference"
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank())

  save_plot(p_manhattan, "manhattan_plot_dmps", w = 14, h = 6)
  rm(manhattan_data, chr_info, p_manhattan); gc(verbose = FALSE)
}

# Beta distribution
cat("Creating beta distribution plot...\n")
set.seed(42)
beta_sample_idx <- sort(sample(nrow(bs_obj), min(500000, nrow(bs_obj))))
beta_long <- data.frame()
for (i in 1:ncol(bs_obj)) {
  bv <- getMeth(bs_obj[beta_sample_idx, i], type = "raw")[, 1]
  bv <- bv[!is.na(bv) & !is.nan(bv)]
  beta_long <- rbind(beta_long, data.frame(
    sample = colnames(bs_obj)[i],
    condition = pData(bs_obj)$condition[i],
    beta = bv
  ))
  rm(bv)
}

p_beta <- ggplot(beta_long, aes(x = beta, fill = condition, color = condition)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values = c("Control" = "#2471A3", "Amputated" = "#C0392B")) +
  scale_color_manual(values = c("Control" = "#2471A3", "Amputated" = "#C0392B")) +
  facet_wrap(~sample, ncol = 2) +
  labs(
    title = "CpG Methylation Level Distribution",
    subtitle = sprintf("D. laeve WGBS | %s CpG sites sampled",
                       format(length(beta_sample_idx), big.mark = ",")),
    x = "Beta value (methylation level)", y = "Density"
  ) +
  theme_minimal(base_size = 12)

save_plot(p_beta, "beta_distribution", w = 10, h = 7)
rm(beta_long, p_beta, beta_sample_idx); gc(verbose = FALSE)


# =============================================================================
# 7. GENOMIC ANNOTATION
# =============================================================================
cat("\n=== 7. Genomic Annotation ===\n")

gff <- readRDS(file.path(cache_dir, "gff_chr1_31.rds"))
genes <- gff[gff$type == "gene"]
exons <- gff[gff$type == "exon"]

te_rds <- file.path(cache_dir, "te_chr1_31.rds")
has_te <- file.exists(te_rds)
if (has_te) {
  te_data <- readRDS(te_rds)
  if (is(te_data, "list") && "te_gr" %in% names(te_data)) {
    te_gr <- te_data$te_gr
  } else if (is(te_data, "GRanges")) {
    te_gr <- te_data
  } else {
    te_gr <- te_data[[1]]
  }
  cat("TE annotations loaded:", length(te_gr), "elements\n")
  rm(te_data); gc(verbose = FALSE)
}

promoters_rds <- file.path(cache_dir, "promoters_2kb.rds")
if (file.exists(promoters_rds)) {
  promoters <- readRDS(promoters_rds)
} else {
  chr_lens <- seqlengths(genes)
  gene_strand <- as.character(strand(genes))
  gene_chr <- as.character(seqnames(genes))
  prom_start <- ifelse(gene_strand == "+",
                       pmax(1L, start(genes) - 2000L), end(genes) + 1L)
  prom_end <- ifelse(gene_strand == "+",
                     start(genes) - 1L,
                     pmin(chr_lens[gene_chr], end(genes) + 2000L))
  valid <- prom_start < prom_end
  promoters <- GRanges(seqnames = gene_chr[valid],
                       ranges = IRanges(start = prom_start[valid], end = prom_end[valid]),
                       gene_id = genes$ID[valid])
}
cat("Promoters:", length(promoters), "\n")

annotate_positions <- function(chr, pos) {
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1))
  in_promoter <- overlapsAny(gr, promoters)
  in_exon     <- overlapsAny(gr, exons)
  in_gene     <- overlapsAny(gr, genes)
  in_te       <- if (has_te) overlapsAny(gr, te_gr) else rep(FALSE, length(gr))
  annotation <- rep("Intergenic", length(gr))
  annotation[in_te] <- "TE"
  annotation[in_gene & !in_exon] <- "Intron"
  annotation[in_exon] <- "Exon"
  annotation[in_promoter] <- "Promoter"
  return(annotation)
}

dmp_anno_summary <- NULL
if (nrow(dmp_results) > 0) {
  cat("Annotating DMPs...\n")
  dmp_results$annotation <- annotate_positions(dmp_results$chr, dmp_results$pos)

  dmp_anno_summary <- as.data.frame(table(dmp_results$annotation))
  colnames(dmp_anno_summary) <- c("Region", "Count")
  dmp_anno_summary$Percent <- round(100 * dmp_anno_summary$Count / sum(dmp_anno_summary$Count), 1)
  dmp_anno_summary <- dmp_anno_summary[order(-dmp_anno_summary$Count), ]
  cat("DMP annotation summary:\n")
  print(dmp_anno_summary)

  region_colors <- c("Promoter" = "#8E44AD", "Exon" = "#2471A3",
                     "Intron" = "#1ABC9C", "TE" = "#F39C12", "Intergenic" = "#C0392B")

  p_anno <- ggplot(dmp_anno_summary, aes(x = "", y = Count, fill = Region)) +
    geom_col(width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = region_colors) +
    labs(title = "DMP Genomic Distribution",
         subtitle = sprintf("D. laeve | %s significant DMPs",
                            format(nrow(dmp_results), big.mark = ","))) +
    theme_void(base_size = 12) +
    theme(legend.position = "right")

  save_plot(p_anno, "dmp_annotation_pie", w = 8, h = 6)

  dmp_results$direction <- ifelse(dmp_results$diff > 0, "Hyper", "Hypo")
  dir_anno <- as.data.frame(table(dmp_results$annotation, dmp_results$direction))
  colnames(dir_anno) <- c("Region", "Direction", "Count")

  p_dir <- ggplot(dir_anno, aes(x = Region, y = Count, fill = Direction)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Hyper" = "#C0392B", "Hypo" = "#2471A3")) +
    labs(title = "DMP Direction by Genomic Region",
         subtitle = "D. laeve | Amputated vs Control",
         x = NULL, y = "Number of DMPs") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  save_plot(p_dir, "dmp_direction_by_region", w = 8, h = 6)

  cat("Finding nearest genes for DMPs...\n")
  dmp_gr <- GRanges(seqnames = dmp_results$chr,
                    ranges = IRanges(start = dmp_results$pos, width = 1))
  nearest_idx <- nearest(dmp_gr, genes)
  dmp_results$nearest_gene <- ifelse(!is.na(nearest_idx), genes$ID[nearest_idx], NA)
  dmp_results$dist_to_gene <- ifelse(!is.na(nearest_idx),
                                     mcols(distanceToNearest(dmp_gr, genes))$distance, NA)

  write.table(dmp_results, file.path(out_dir, "dmps_annotated.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  rm(dmp_gr, p_anno, p_dir, dir_anno); gc(verbose = FALSE)
}

dmr_anno_summary <- NULL
if (nrow(dmr_results) > 0) {
  cat("\nAnnotating DMRs...\n")
  dmr_mid <- (dmr_results$start + dmr_results$end) %/% 2
  dmr_results$annotation <- annotate_positions(dmr_results$chr, dmr_mid)

  dmr_anno_summary <- as.data.frame(table(dmr_results$annotation))
  colnames(dmr_anno_summary) <- c("Region", "Count")
  dmr_anno_summary$Percent <- round(100 * dmr_anno_summary$Count / sum(dmr_anno_summary$Count), 1)
  cat("DMR annotation summary:\n")
  print(dmr_anno_summary)

  dmr_gr <- GRanges(seqnames = dmr_results$chr,
                    ranges = IRanges(start = dmr_results$start, end = dmr_results$end))
  nearest_idx <- nearest(dmr_gr, genes)
  dmr_results$nearest_gene <- ifelse(!is.na(nearest_idx), genes$ID[nearest_idx], NA)

  write.table(dmr_results, file.path(out_dir, "dmrs_annotated.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

rm(gff, genes, exons); gc(verbose = FALSE)


# =============================================================================
# 8. METHYLATION HEATMAP
# =============================================================================
cat("\n=== 8. Methylation Heatmaps ===\n")

ann_colors <- list(Condition = c("Control" = "#2471A3", "Amputated" = "#C0392B"))

if (nrow(dmp_results) > 0) {
  cat("Creating DMP heatmap...\n")
  n_top <- min(100, nrow(dmp_results))
  top_dmps <- dmp_results[order(abs(dmp_results$diff), decreasing = TRUE), ][1:n_top, ]

  dmp_ids <- paste(top_dmps$chr, top_dmps$pos, sep = ":")
  bs_ids  <- paste(as.character(seqnames(bs_obj)), start(bs_obj), sep = ":")
  match_idx <- which(bs_ids %in% dmp_ids)

  if (length(match_idx) >= 2) {
    hm_beta <- getMeth(bs_obj[match_idx, ], type = "raw")
    rownames(hm_beta) <- paste0(as.character(seqnames(bs_obj[match_idx])), ":",
                                start(bs_obj[match_idx]))
    colnames(hm_beta) <- colnames(bs_obj)

    valid <- complete.cases(hm_beta) & apply(hm_beta, 1, var, na.rm = TRUE) > 0
    hm_beta <- hm_beta[valid, , drop = FALSE]

    if (nrow(hm_beta) >= 2) {
      col_anno <- data.frame(Condition = pData(bs_obj)$condition,
                             row.names = colnames(hm_beta))

      p_hm <- pheatmap(hm_beta,
                       annotation_col = col_anno,
                       annotation_colors = ann_colors,
                       color = colorRampPalette(c("#2471A3", "white", "#C0392B"))(100),
                       main = sprintf("Top %d DMPs -- Methylation Levels", nrow(hm_beta)),
                       fontsize = 10,
                       show_rownames = (nrow(hm_beta) <= 50),
                       clustering_method = "ward.D2",
                       silent = TRUE)

      png(file.path(fig_png_dir, "dmp_heatmap.png"),
          width = 8, height = 10, units = "in", res = 300)
      grid::grid.draw(p_hm$gtable)
      dev.off()
      tryCatch({
        cairo_pdf(file.path(fig_pdf_dir, "dmp_heatmap.pdf"), width = 8, height = 10)
        grid::grid.draw(p_hm$gtable); dev.off()
      }, error = function(e) {
        pdf(file.path(fig_pdf_dir, "dmp_heatmap.pdf"), width = 8, height = 10)
        grid::grid.draw(p_hm$gtable); dev.off()
      })
      cat("  Saved: dmp_heatmap\n")
    }
  }
}

rm(bs_obj); gc(verbose = FALSE)


# =============================================================================
# 9. HTML REPORT
# =============================================================================
cat("\n=== 9. Generating HTML Report ===\n")

report_path <- file.path(out_dir, "part3_report.html")

sink(report_path)
cat("<!DOCTYPE html>\n<html><head><meta charset='utf-8'>\n")
cat("<title>Methylation Tutorial Part 3: D. laeve Full-Genome Analysis</title>\n")
cat("<style>\n")
cat("body { font-family: 'Segoe UI', Arial, sans-serif; max-width: 1200px; margin: 40px auto; padding: 0 20px; color: #333; }\n")
cat("h1 { color: #2c3e50; border-bottom: 2px solid #2471A3; padding-bottom: 10px; }\n")
cat("h2 { color: #2471A3; margin-top: 30px; }\n")
cat("table { border-collapse: collapse; width: auto; margin: 15px 0; }\n")
cat("th { background: #2471A3; color: white; padding: 8px 12px; text-align: left; }\n")
cat("td { padding: 6px 12px; border-bottom: 1px solid #ddd; }\n")
cat("tr:nth-child(even) { background: #f8f9fa; }\n")
cat("img { max-width: 100%; border: 1px solid #ddd; margin: 10px 0; }\n")
cat(".key-finding { background: #eaf2f8; padding: 12px; border-left: 4px solid #2471A3; margin: 15px 0; }\n")
cat(".metric { font-size: 1.2em; font-weight: bold; color: #2471A3; }\n")
cat("</style>\n</head><body>\n")

cat("<h1>Methylation Tutorial Part 3: Full-Genome WGBS Analysis</h1>\n")
cat("<p>Adapted from <a href='https://ngs101.com'>NGS101</a> for <em>D. laeve</em> tail regeneration | Generated:", format(Sys.Date()), "</p>\n")
cat("<p><strong>Design:</strong> 2 Control (C1, C2) vs 2 Amputated (A1, A2) tail samples</p>\n")
cat("<p><strong>FULL GENOME RUN</strong> -- all CpG sites, DSS smoothing enabled, native Windows R</p>\n\n")

cat("<h2>1. Analysis Summary</h2>\n")
cat("<table>\n")
for (i in seq_len(nrow(summary_stats))) {
  cat(sprintf("<tr><td>%s</td><td class='metric'>%s</td></tr>\n",
              summary_stats$Metric[i], summary_stats$Value[i]))
}
cat("</table>\n\n")

cat("<h2>2. Global Methylation Statistics</h2>\n")
cat("<table>\n")
cat("<tr><th>Sample</th><th>Condition</th><th>Sites</th><th>Mean Beta</th><th>Median Beta</th><th>High (&gt;0.8)</th><th>Low (&lt;0.2)</th><th>Intermediate</th></tr>\n")
for (i in seq_len(nrow(stats_df))) {
  r <- stats_df[i, ]
  cat(sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%.4f</td><td>%.4f</td><td>%s</td><td>%s</td><td>%s</td></tr>\n",
              r$sample, r$condition, format(r$total_sites, big.mark = ","),
              r$mean_methylation, r$median_methylation,
              format(r$sites_high_meth, big.mark = ","),
              format(r$sites_low_meth, big.mark = ","),
              format(r$sites_intermediate, big.mark = ",")))
}
cat("</table>\n\n")

cat("<h2>3. Beta Value Distribution</h2>\n")
cat("<img src='png/beta_distribution.png'>\n\n")

cat("<h2>4. Sample Correlation</h2>\n")
cat("<img src='png/sample_correlation_heatmap.png'>\n\n")

cat("<h2>5. PCA of Methylation Profiles</h2>\n")
cat("<img src='png/pca_methylation.png'>\n\n")

if (nrow(dmp_results) > 0) {
  cat("<h2>6. Volcano Plot -- DMPs</h2>\n")
  cat("<img src='png/volcano_plot_dmps.png'>\n\n")

  cat("<h2>7. Manhattan Plot -- DMPs</h2>\n")
  cat("<img src='png/manhattan_plot_dmps.png'>\n\n")
}

cat("<h2>8. Genome-Wide Methylation Profiles</h2>\n")
cat("<img src='png/genome_wide_methylation.png'>\n\n")

if (nrow(dmp_results) > 0) {
  cat("<h2>9. DMP Genomic Distribution</h2>\n")
  cat("<img src='png/dmp_annotation_pie.png'>\n\n")

  cat("<h2>10. DMP Direction by Region</h2>\n")
  cat("<img src='png/dmp_direction_by_region.png'>\n\n")

  if (!is.null(dmp_anno_summary)) {
    cat("<h3>DMP Annotation Summary</h3>\n")
    cat("<table><tr><th>Region</th><th>Count</th><th>%</th></tr>\n")
    for (i in seq_len(nrow(dmp_anno_summary))) {
      cat(sprintf("<tr><td>%s</td><td>%s</td><td>%.1f</td></tr>\n",
                  dmp_anno_summary$Region[i],
                  format(dmp_anno_summary$Count[i], big.mark = ","),
                  dmp_anno_summary$Percent[i]))
    }
    cat("</table>\n\n")
  }

  cat("<h2>11. DMP Heatmap (Top DMPs)</h2>\n")
  cat("<img src='png/dmp_heatmap.png'>\n\n")
}

if (nrow(dmr_results) > 0 && !is.null(dmr_anno_summary)) {
  cat("<h2>12. DMR Summary</h2>\n")
  cat("<table><tr><th>Region</th><th>Count</th><th>%</th></tr>\n")
  for (i in seq_len(nrow(dmr_anno_summary))) {
    cat(sprintf("<tr><td>%s</td><td>%s</td><td>%.1f</td></tr>\n",
                dmr_anno_summary$Region[i],
                format(dmr_anno_summary$Count[i], big.mark = ","),
                dmr_anno_summary$Percent[i]))
  }
  cat("</table>\n\n")
}

cat("<h2>Methods</h2>\n")
cat("<div class='key-finding'>\n")
cat("<p><strong>Full-genome analysis (native Windows R):</strong></p>\n")
cat("<ul>\n")
cat("<li><em>D. laeve</em> tail regeneration, custom genome (1.78 Gb, 31 chromosomes)</li>\n")
cat("<li>Input: Bismark CpG_report.txt files -- strand collapsed, coverage >= 5</li>\n")
cat("<li>Differential methylation: DSS with smoothing enabled (full genome)</li>\n")
cat(sprintf("<li>Significance: FDR < %.2f, |methylation difference| > %.0f%%</li>\n",
            analysis_params$fdr_threshold, analysis_params$min_diff * 100))
cat("<li>Annotation: custom GFF (EviAnn) -- promoter/exon/intron/TE/intergenic</li>\n")
cat("</ul>\n</div>\n")

cat("\n</body></html>\n")
sink()
cat("Report saved:", report_path, "\n")

cat("\n=== Full-Genome Analysis Complete ===\n")
cat("Output directory:", out_dir, "\n")
cat("Total time:", format(Sys.time()), "\n")
