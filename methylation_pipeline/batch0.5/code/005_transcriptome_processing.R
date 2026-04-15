#!/usr/bin/env Rscript
# =============================================================================
# Batch 0.5: Transcriptome Processing
# Question: What is the clean, filtered expression dataset for all downstream analysis?
# Output: data/ (sample metadata, filtered counts, DE results) + figures/ (QC plots)
#
# This batch processes ALL transcriptome data once. Downstream batches (02, 04,
# 07, 08) read from here instead of running DESeq2 independently.
#
# Samples: 44 total across 6 tissues × multiple conditions
#   Tail control:        C1S1, C2S2, C3S3, C4S4
#   Tail amputated:      T1S5, T2S6, T3S7, T4S8
#   Eye control:         C1, C4, C5, C8
#   Eye amputated:       R1, R3, R4, R6
#   Bodywall control:    dcrep1-4, dcrep6, fungicide_l0a-d
#   Bodywall irradiated: irrep1, irrep4-7
#   Bodywall fungicide:  fungicide_l30a-e
#   Head control:        fpDLHead1-3
#   Juvenile control:    fpDLJuv1-3
#   Ovotestis control:   fpDLOvo1-3
# NOTE: filtering (>=10 counts in >=min_group_size) is within-group (per tissue
# DESeq2 model), so groups don't need to be balanced across tissues.
#
# Outlier detection: Cook's distance (DESeq2 built-in) + PCA visual inspection
# Gene filtering: rowSums(counts >= 10) >= min_group_size (Bioconductor standard)
# LFC shrinkage: apeglm (recommended over normal)
# VST blind=TRUE for QC, subset per tissue for testing (GTEx/Broad standard)
# =============================================================================

source("methylation_pipeline/_config.R")
t0 <- proc.time()

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
})

BATCH_DIR <- file.path(PIPE_DIR, "batch0.5")

# Clean old output
unlink(list.files(file.path(BATCH_DIR, "figures"), full.names = TRUE))
unlink(list.files(file.path(BATCH_DIR, "data"), full.names = TRUE))
dir.create(file.path(BATCH_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BATCH_DIR, "data"), showWarnings = FALSE, recursive = TRUE)

cat("=== Batch 0.5: Transcriptome Processing ===\n\n")

# =============================================================================
# 1. SAMPLE METADATA — assign tissue + condition from filenames
# =============================================================================
cat("[1/6] Building sample metadata...\n")
all_files <- list.files(OG$counts_dir, pattern = "\\.txt$")

meta <- data.frame(
  sampleName = all_files,
  fileName   = all_files,
  tissue     = NA_character_,
  condition  = NA_character_,
  stringsAsFactors = FALSE
)

# Tail control: C1S1, C2S2, C3S3, C4S4
idx <- grepl("^C\\d+S\\d+_", all_files)
meta$tissue[idx] <- "tail"; meta$condition[idx] <- "control"

# Tail amputated: T1S5, T2S6, T3S7, T4S8
idx <- grepl("^T\\d+S\\d+_", all_files)
meta$tissue[idx] <- "tail"; meta$condition[idx] <- "amputated"

# Eye control: C1, C4, C5, C8 (no S suffix)
idx <- grepl("^C\\d+_htseq", all_files)
meta$tissue[idx] <- "eye"; meta$condition[idx] <- "control"

# Bodywall control: dcrep1-4, dcrep6
idx <- grepl("^dcrep", all_files)
meta$tissue[idx] <- "bodywall"; meta$condition[idx] <- "control"

# Eye amputated: R1, R3, R4, R6
idx <- grepl("^R\\d+_", all_files)
meta$tissue[idx] <- "eye"; meta$condition[idx] <- "amputated"

# Bodywall irradiated: irrep1, irrep4-7
idx <- grepl("^irrep", all_files)
meta$tissue[idx] <- "bodywall"; meta$condition[idx] <- "irradiated"

# Bodywall control (solvent): fungicide_l0*
idx <- grepl("^fungicide_l0", all_files)
meta$tissue[idx] <- "bodywall"; meta$condition[idx] <- "control"

# Bodywall fungicide: fungicide_l30*
idx <- grepl("^fungicide_l30", all_files)
meta$tissue[idx] <- "bodywall"; meta$condition[idx] <- "fungicide"

# Head control
idx <- grepl("fpDLHead", all_files)
meta$tissue[idx] <- "head"; meta$condition[idx] <- "control"

# Juvenile control
idx <- grepl("fpDLJuv", all_files)
meta$tissue[idx] <- "juvenile"; meta$condition[idx] <- "control"

# Ovotestis control
idx <- grepl("fpDLOvo", all_files)
meta$tissue[idx] <- "ovotestis"; meta$condition[idx] <- "control"

# Drop unassigned
meta <- meta[!is.na(meta$tissue), ]
meta$group <- paste(meta$tissue, meta$condition, sep = "_")

# Readable labels
tissue_counts <- list()
meta$label <- vapply(seq_len(nrow(meta)), function(i) {
  key <- meta$group[i]
  tissue_counts[[key]] <<- (tissue_counts[[key]] %||% 0L) + 1L
  sprintf("%s_%s_%d", tools::toTitleCase(meta$tissue[i]),
          tools::toTitleCase(meta$condition[i]), tissue_counts[[key]])
}, character(1))

cat(sprintf("  Total samples: %d\n", nrow(meta)))
cat("  Breakdown:\n")
for (g in sort(unique(meta$group))) {
  cat(sprintf("    %-25s %d samples\n", g, sum(meta$group == g)))
}

save_data(meta, BATCH_DIR, "sample_metadata")

# =============================================================================
# 2. BUILD DESeq2 OBJECT (all samples)
# =============================================================================
cat("\n[2/6] Building DESeq2 object (all samples)...\n")

dds_all <- DESeqDataSetFromHTSeqCount(
  sampleTable = meta,
  directory = OG$counts_dir,
  design = ~ 1  # no design — just for normalization + QC
)
# Gene filter: require >= 10 counts in at least min_group_size samples
# (Bioconductor standard, stricter than rowSums > 10)
min_group <- min(table(meta$group))
keep_genes <- rowSums(counts(dds_all) >= 10) >= min_group
dds_all <- dds_all[keep_genes, ]
cat(sprintf("  Genes after filter (>= 10 counts in >= %d samples): %d\n",
            min_group, nrow(dds_all)))

# VST for QC (blind=TRUE for unbiased QC, as per DESeq2 vignette)
vsd_all <- varianceStabilizingTransformation(dds_all, blind = TRUE)
vst_mat <- assay(vsd_all)
colnames(vst_mat) <- meta$label[match(colnames(vst_mat), meta$sampleName)]

cat(sprintf("  VST matrix: %d genes x %d samples\n", nrow(vst_mat), ncol(vst_mat)))

# =============================================================================
# 3. WITHIN-GROUP OUTLIER DETECTION (data-driven)
# =============================================================================
cat("\n[3/6] Within-group outlier detection...\n")

# Two-step outlier detection:
# Step A: Within-group Euclidean distance to centroid (VST space)
#   Flag samples with distance > 2x median within-group distance
# Step B: Visual PCA inspection (saved as figures for manual review)
# Cook's distance handles per-gene outliers automatically in DESeq2
outliers <- c()

for (g in unique(meta$group)) {
  g_idx <- which(meta$group == g)
  if (length(g_idx) < 3) next

  g_labels <- meta$label[g_idx]
  g_mat <- t(vst_mat[, g_labels, drop = FALSE])
  centroid <- colMeans(g_mat)
  dists <- apply(g_mat, 1, function(x) sqrt(sum((x - centroid)^2)))
  med_dist <- median(dists)

  bad <- which(dists > 2 * med_dist)
  if (length(bad) > 0) {
    bad_labels <- rownames(g_mat)[bad]
    cat(sprintf("  %s: outlier(s) = %s (dist %.0f vs median %.0f)\n",
                g, paste(bad_labels, collapse = ", "),
                max(dists[bad]), med_dist))
    outliers <- c(outliers, bad_labels)
  }
}
cat("  NOTE: Cook's distance (per-gene) handled automatically by DESeq2.\n")
cat("  NOTE: Check PCA figures for visual confirmation of outliers.\n")

if (length(outliers) > 0) {
  cat(sprintf("\n  Removing %d outlier(s): %s\n", length(outliers), paste(outliers, collapse = ", ")))
  meta$outlier <- meta$label %in% outliers
} else {
  cat("  No within-group outliers detected.\n")
  meta$outlier <- FALSE
}

meta_clean <- meta[!meta$outlier, ]
cat(sprintf("  Clean samples: %d / %d\n", nrow(meta_clean), nrow(meta)))
save_data(meta, BATCH_DIR, "sample_metadata_with_outliers")

# =============================================================================
# 4. DESeq2: TAIL CONTROL vs AMPUTATED (clean samples)
# =============================================================================
cat("\n[4/6] DESeq2: tail control vs amputated...\n")

tail_meta <- meta_clean[meta_clean$tissue == "tail" &
                          meta_clean$condition %in% c("control", "amputated"), ]
cat(sprintf("  Tail samples: %d control + %d amputated\n",
            sum(tail_meta$condition == "control"), sum(tail_meta$condition == "amputated")))

# Subset to tail only — per-tissue testing is the standard approach
# (GTEx/Broad: run DESeq2 separately per tissue, not one big model)
dds_tail <- DESeqDataSetFromHTSeqCount(
  sampleTable = tail_meta[, c("sampleName", "fileName", "condition")],
  directory = OG$counts_dir,
  design = ~ condition
)
dds_tail$condition <- relevel(dds_tail$condition, ref = "control")

# Proper gene filter: >= 10 counts in >= smallest group size
min_tail_group <- min(table(tail_meta$condition))
keep_tail <- rowSums(counts(dds_tail) >= 10) >= min_tail_group
dds_tail <- dds_tail[keep_tail, ]
cat(sprintf("  Genes after filter: %d\n", sum(keep_tail)))

dds_tail <- DESeq(dds_tail)

# LFC shrinkage with apeglm (recommended, replaces type="normal")
if (requireNamespace("apeglm", quietly = TRUE)) {
  res_tail <- as.data.frame(lfcShrink(dds_tail, coef = "condition_amputated_vs_control",
                                        type = "apeglm"))
  cat("  LFC shrinkage: apeglm\n")
} else {
  res_tail <- as.data.frame(results(dds_tail, alpha = 0.05))
  cat("  LFC shrinkage: none (apeglm not installed)\n")
}
res_tail$gene_id <- rownames(res_tail)
res_tail <- res_tail[!is.na(res_tail$padj), ]

n_de <- sum(res_tail$padj < 0.05)
cat(sprintf("  DE genes (FDR < 0.05): %d\n", n_de))
cat(sprintf("  Up: %d | Down: %d\n",
            sum(res_tail$padj < 0.05 & res_tail$log2FoldChange > 0),
            sum(res_tail$padj < 0.05 & res_tail$log2FoldChange < 0)))

save_data(res_tail, BATCH_DIR, "deseq2_tail_results")

# =============================================================================
# 5. MULTI-TISSUE VST (clean samples, for batch02 + batch08)
# =============================================================================
cat("\n[5/6] Multi-tissue VST (clean samples)...\n")

# Multi-tissue VST: design = ~1 for blind transformation (QC + WGCNA input)
# No DE testing here — per-tissue testing done above
dds_clean <- DESeqDataSetFromHTSeqCount(
  sampleTable = meta_clean[, c("sampleName", "fileName", "condition")],
  directory = OG$counts_dir,
  design = ~ 1
)
min_group_clean <- min(table(meta_clean$group))
keep_clean <- rowSums(counts(dds_clean) >= 10) >= min_group_clean
dds_clean <- dds_clean[keep_clean, ]
vsd_clean <- varianceStabilizingTransformation(dds_clean, blind = TRUE)
vst_clean <- assay(vsd_clean)

# Rename columns to readable labels
col_match <- match(colnames(vst_clean), meta_clean$sampleName)
colnames(vst_clean) <- meta_clean$label[col_match]

cat(sprintf("  Clean VST: %d genes x %d samples\n", nrow(vst_clean), ncol(vst_clean)))

save_data(as.data.frame(vst_clean, check.names = FALSE), BATCH_DIR, "vst_all_tissues")
save_data(meta_clean, BATCH_DIR, "sample_metadata_clean")

# Cache for downstream batches
cache_file <- CACHE$transcriptome
saveRDS(list(
  meta_all     = meta,
  meta_clean   = meta_clean,
  outliers     = outliers,
  res_tail     = res_tail,
  vst_clean    = vst_clean
), cache_file)
cat(sprintf("  Cached: %s\n", cache_file))

# =============================================================================
# 6. QC FIGURES
# =============================================================================
cat("\n[6/6] QC figures...\n")

# Fig 0.5A: Sample dendrogram (all samples, outliers marked)
sample_tree <- hclust(dist(t(vst_mat)), method = "average")
outlier_colors <- ifelse(colnames(vst_mat) %in% outliers, "red", "black")

fig_dir <- file.path(BATCH_DIR, "figures")
p05a_fn <- function() {
  par(mar = c(12, 4, 2, 1))
  dend <- as.dendrogram(sample_tree)
  plot(dend, main = "All samples — hierarchical clustering",
       xlab = "", sub = "", cex = 0.5, hang = -1)
  # Mark outliers in title
  if (length(outliers) > 0) {
    mtext(paste("Outliers (red):", paste(outliers, collapse = ", ")),
          side = 1, line = 10, cex = 0.7, col = "red")
  }
}
png(file.path(fig_dir, "fig05a_sample_dendrogram.png"),
    width = 16, height = 8, units = "in", res = 300)
p05a_fn(); dev.off()
tryCatch({
  cairo_pdf(file.path(fig_dir, "fig05a_sample_dendrogram.pdf"), width = 16, height = 8)
  p05a_fn(); dev.off()
}, error = function(e) {
  pdf(file.path(fig_dir, "fig05a_sample_dendrogram.pdf"), width = 16, height = 8)
  p05a_fn(); dev.off()
})
cat("  Saved: fig05a_sample_dendrogram\n")

# Fig 0.5B: PCA (all samples, colored by tissue)
pca <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
ve <- summary(pca)$importance[2, 1:2] * 100
pca_df <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  label = colnames(vst_mat),
  tissue = meta$tissue[match(colnames(vst_mat), meta$label)],
  condition = meta$condition[match(colnames(vst_mat), meta$label)],
  outlier = colnames(vst_mat) %in% outliers
)

tissue_colors <- c(tail = "#C0392B", eye = "#2471A3", bodywall = "#27AE60",
                   head = "#8E44AD", juvenile = "#F39C12", ovotestis = "#1ABC9C")

p05b <- ggplot(pca_df, aes(x = PC1, y = PC2, color = tissue, shape = outlier)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = tissue_colors, name = "Tissue") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), name = "Outlier") +
  labs(title = "PCA — all samples",
       subtitle = sprintf("PC1 (%.1f%%) vs PC2 (%.1f%%)", ve[1], ve[2]),
       x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
  theme_minimal(base_size = 12)
save_fig(p05b, BATCH_DIR, "fig05b_pca_all_samples", w = 10, h = 8)

# Fig 0.5C: PCA tail only (control vs amputated)
tail_labels <- meta_clean$label[meta_clean$tissue == "tail" &
                                  meta_clean$condition %in% c("control", "amputated")]
tail_vst <- vst_clean[, colnames(vst_clean) %in% tail_labels, drop = FALSE]
if (ncol(tail_vst) >= 4) {
  pca_tail <- prcomp(t(tail_vst), center = TRUE, scale. = FALSE)
  ve_tail <- summary(pca_tail)$importance[2, 1:2] * 100
  pca_tail_df <- data.frame(
    PC1 = pca_tail$x[, 1], PC2 = pca_tail$x[, 2],
    label = colnames(tail_vst),
    condition = meta_clean$condition[match(colnames(tail_vst), meta_clean$label)]
  )
  p05c <- ggplot(pca_tail_df, aes(x = PC1, y = PC2, color = condition, label = label)) +
    geom_point(size = 5) +
    geom_text(vjust = -1.2, size = 3.5, show.legend = FALSE) +
    scale_color_manual(values = COLORS$condition) +
    labs(title = "PCA — tail control vs amputated (clean)",
         x = sprintf("PC1 (%.1f%%)", ve_tail[1]),
         y = sprintf("PC2 (%.1f%%)", ve_tail[2])) +
    theme_minimal(base_size = 12)
  save_fig(p05c, BATCH_DIR, "fig05c_pca_tail", w = 8, h = 6)
}

# Fig 0.5D: Sample-to-sample correlation heatmap (all clean)
cat("  Sample correlation heatmap...\n")
sample_cors <- cor(vst_clean, method = "spearman")

ann_df <- data.frame(
  Tissue = meta_clean$tissue[match(colnames(vst_clean), meta_clean$label)],
  Condition = meta_clean$condition[match(colnames(vst_clean), meta_clean$label)],
  row.names = colnames(vst_clean)
)
ann_colors <- list(
  Tissue = tissue_colors,
  Condition = c(control = "#2471A3", amputated = "#C0392B",
                irradiated = "#F39C12", fungicide = "#27AE60")
)

png(file.path(fig_dir, "fig05d_correlation_heatmap.png"),
    width = 14, height = 12, units = "in", res = 300)
pheatmap(sample_cors, annotation_col = ann_df, annotation_row = ann_df,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main = "Sample-to-sample correlation (Spearman, all clean samples)",
         fontsize = 7, show_rownames = TRUE, show_colnames = FALSE)
dev.off()
tryCatch({
  cairo_pdf(file.path(fig_dir, "fig05d_correlation_heatmap.pdf"), width = 14, height = 12)
  pheatmap(sample_cors, annotation_col = ann_df, annotation_row = ann_df,
           annotation_colors = ann_colors,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = "Sample-to-sample correlation (Spearman, all clean samples)",
           fontsize = 7, show_rownames = TRUE, show_colnames = FALSE)
  dev.off()
}, error = function(e) NULL)
cat("  Saved: fig05d_correlation_heatmap\n")

# =============================================================================
# SUMMARY
# =============================================================================
elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("\n=== Batch 0.5 complete (%.1f min) ===\n", elapsed / 60))
cat(sprintf("Samples: %d total, %d clean (%d outliers)\n",
            nrow(meta), nrow(meta_clean), sum(meta$outlier)))
cat(sprintf("Tail DE genes: %d (FDR < 0.05)\n", n_de))
cat(sprintf("Cached to: %s\n", cache_file))
n_figs <- length(list.files(fig_dir, pattern = "\\.png$"))
n_data <- length(list.files(file.path(BATCH_DIR, "data"), pattern = "\\.tsv$"))
cat(sprintf("Figures: %d | Data files: %d\n", n_figs, n_data))
