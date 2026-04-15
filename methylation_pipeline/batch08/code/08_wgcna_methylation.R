#!/usr/bin/env Rscript
# =============================================================================
# Batch 08: WGCNA Modules + Methylation Enrichment
# Question: Are methylation changes concentrated in specific expression modules?
# Output: data/ (module assignments, enrichment, hub genes) + figures/ (10 plots)
# Requires: Batch 06 DMPs, DESeq2 multi-tissue counts, STRING enrichment terms
# =============================================================================

source("methylation_pipeline/_config.R")
t0 <- proc.time()

suppressPackageStartupMessages({
  library(WGCNA)
  library(DESeq2)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(clusterProfiler)
  library(patchwork)
})

BATCH_DIR <- file.path(PIPE_DIR, "batch08")

# Clean old output
unlink(list.files(file.path(BATCH_DIR, "figures"), full.names = TRUE))
unlink(list.files(file.path(BATCH_DIR, "data"), full.names = TRUE))
dir.create(file.path(BATCH_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BATCH_DIR, "data"), showWarnings = FALSE, recursive = TRUE)

cat("=== Batch 08: WGCNA + Methylation ===\n\n")
allowWGCNAThreads(4)

# =============================================================================
# 1. BUILD MULTI-TISSUE EXPRESSION MATRIX
# =============================================================================
cat("[1/9] Building multi-tissue expression matrix...\n")

if (file.exists(CACHE$transcriptome)) {
  cat("  Loading from CACHE$transcriptome...\n")
  tc <- readRDS(CACHE$transcriptome)
  vst_mat_raw <- tc$vst_clean
  sample_table <- tc$meta_clean
  cat(sprintf("  VST matrix: %d genes x %d samples (outliers from batch0.5 already removed)\n",
              nrow(vst_mat_raw), ncol(vst_mat_raw)))
  expr_mat <- t(vst_mat_raw)  # samples x genes
  rm(tc, vst_mat_raw); gc(verbose = FALSE)
} else {
  cat("  Cache not found — running DESeq2...\n")
  all_files <- list.files(OG$counts_dir, pattern = "\\.txt$")

  sample_table <- data.frame(sampleName = all_files, fileName = all_files, stringsAsFactors = FALSE)
  sample_table$condition <- NA_character_
  sample_table$tissue    <- NA_character_

  idx <- grepl("^C\\d+S\\d+_", all_files)
  sample_table$condition[idx] <- "control"; sample_table$tissue[idx] <- "tail"
  idx <- grepl("^T\\d+S\\d+_", all_files)
  sample_table$condition[idx] <- "amputated"; sample_table$tissue[idx] <- "tail"
  idx <- grepl("^C\\d+_htseq", all_files)
  sample_table$condition[idx] <- "control"; sample_table$tissue[idx] <- "eye"
  idx <- grepl("^dcrep", all_files)
  sample_table$condition[idx] <- "control"; sample_table$tissue[idx] <- "bodywall"
  idx <- grepl("^R\\d+_", all_files)
  sample_table$condition[idx] <- "amputated"; sample_table$tissue[idx] <- "eye"
  idx <- grepl("^irrep", all_files)
  sample_table$condition[idx] <- "irradiated"; sample_table$tissue[idx] <- "bodywall"
  idx <- grepl("^fungicide_l0", all_files)
  sample_table$condition[idx] <- "control"; sample_table$tissue[idx] <- "bodywall"
  idx <- grepl("^fungicide_l30", all_files)
  sample_table$condition[idx] <- "fungicide"; sample_table$tissue[idx] <- "bodywall"
  idx <- grepl("fpDLHead", all_files)
  sample_table$condition[idx] <- "control"; sample_table$tissue[idx] <- "head"
  idx <- grepl("fpDLJuv", all_files)
  sample_table$condition[idx] <- "control"; sample_table$tissue[idx] <- "juvenile"
  idx <- grepl("fpDLOvo", all_files)
  sample_table$condition[idx] <- "control"; sample_table$tissue[idx] <- "ovotestis"

  sample_table <- sample_table[!is.na(sample_table$condition), ]

  dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table,
                                     directory = OG$counts_dir, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 10, ]
  vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
  expr_mat <- t(assay(vsd))  # samples x genes
  rm(dds, vsd); gc(verbose = FALSE)
}
cat(sprintf("  Samples: %d total\n", nrow(expr_mat)))
cat(sprintf("  Expression matrix: %d samples x %d genes\n", nrow(expr_mat), ncol(expr_mat)))

# =============================================================================
# 2. VARIANCE FILTER + OUTLIER DETECTION
# =============================================================================
cat("[2/9] Filtering low-variance genes + outlier detection...\n")
gene_var <- apply(expr_mat, 2, var)
keep_genes <- gene_var > quantile(gene_var, 0.25)
expr_mat <- expr_mat[, keep_genes]
cat(sprintf("  After variance filter (top 75%%): %d genes\n", ncol(expr_mat)))

# Rename samples: file names → Tissue_N (unique, readable)
# Match by sampleName (DESeq2 path) or label (cache path)
match_col <- {
  if ("sampleName" %in% names(sample_table) &&
      any(rownames(expr_mat) %in% sample_table$sampleName)) "sampleName"
  else if ("label" %in% names(sample_table)) "label"
  else "sampleName"
}
matched_pre <- sample_table[match(rownames(expr_mat), sample_table[[match_col]]), ]

# If already labeled from cache, skip relabeling
if (match_col == "label") {
  sample_labels <- rownames(expr_mat)
  if (!"label" %in% names(sample_table)) sample_table$label <- sample_table[[match_col]]
} else {
  tissue_counts <- list()
  sample_labels <- character(nrow(matched_pre))
  for (i in seq_len(nrow(matched_pre))) {
    tis <- matched_pre$tissue[i]
    cond <- matched_pre$condition[i]
    key <- paste0(tis, "_", cond)
    tissue_counts[[key]] <- (tissue_counts[[key]] %||% 0L) + 1L
    sample_labels[i] <- sprintf("%s_%s_%d",
      tools::toTitleCase(tis), tools::toTitleCase(cond), tissue_counts[[key]])
  }
  rownames(expr_mat) <- sample_labels
  sample_table$label <- sample_labels
}

cat("  Tissue composition:\n")
for (tis in sort(unique(matched_pre$tissue))) {
  n <- sum(matched_pre$tissue == tis)
  cat(sprintf("    %-12s %d samples\n", tis, n))
}

# Sample dendrogram (with readable labels)
sample_tree <- hclust(dist(expr_mat), method = "average")
fig_dir <- file.path(BATCH_DIR, "figures")

# Fig 8A: Sample dendrogram
p8a_fn <- function() {
  par(mar = c(10, 4, 2, 1))
  plot(sample_tree, main = "Sample clustering dendrogram",
       xlab = "", sub = "", cex = 0.6, hang = -1)
}
png(file.path(fig_dir, "fig8a_sample_dendrogram.png"), width = 14, height = 7, units = "in", res = 300)
p8a_fn(); dev.off()
tryCatch({
  cairo_pdf(file.path(fig_dir, "fig8a_sample_dendrogram.pdf"), width = 14, height = 7)
  p8a_fn(); dev.off()
}, error = function(e) {
  pdf(file.path(fig_dir, "fig8a_sample_dendrogram.pdf"), width = 14, height = 7)
  p8a_fn(); dev.off()
})
cat("  Saved: fig8a_sample_dendrogram\n")

# Outlier detection: within each tissue group, flag samples that are far from
# their tissue peers. A sample is an outlier if its distance to its tissue
# centroid is > 2x the median within-tissue distance.
cat("  Within-tissue outlier detection...\n")
outliers <- c()
for (tis in unique(matched_pre$tissue)) {
  tis_idx <- which(matched_pre$tissue == tis)
  if (length(tis_idx) < 3) next  # can't detect outliers with < 3 samples
  tis_mat <- expr_mat[tis_idx, , drop = FALSE]
  centroid <- colMeans(tis_mat)
  dists <- apply(tis_mat, 1, function(x) sqrt(sum((x - centroid)^2)))
  med_dist <- median(dists)
  bad <- which(dists > 2 * med_dist)
  if (length(bad) > 0) {
    cat(sprintf("    %s: outlier(s) = %s (dist %.0f vs median %.0f)\n",
                tis, paste(rownames(tis_mat)[bad], collapse = ", "),
                max(dists[bad]), med_dist))
    outliers <- c(outliers, rownames(tis_mat)[bad])
  }
}

if (length(outliers) > 0) {
  cat(sprintf("  Removing %d within-tissue outlier(s): %s\n",
              length(outliers), paste(outliers, collapse = ", ")))
  keep <- !rownames(expr_mat) %in% outliers
  expr_mat <- expr_mat[keep, ]
  sample_table <- sample_table[sample_table$label %in% rownames(expr_mat), ]
} else {
  cat("  No within-tissue outliers detected.\n")
}
cat(sprintf("  Final matrix: %d samples x %d genes\n", nrow(expr_mat), ncol(expr_mat)))

# =============================================================================
# 3. WGCNA MODULE DETECTION
# =============================================================================
power <- 12

# Fig 8B: Scale-free topology fit
cat(sprintf("[3/9] Scale-free topology + module detection (power = %d)...\n", power))
sft <- pickSoftThreshold(expr_mat, powerVector = 1:20, verbose = 0)
sft_df <- data.frame(Power = sft$fitIndices$Power,
                     SFT = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
                     MeanK = sft$fitIndices$mean.k.)

p8b_sft <- ggplot(sft_df, aes(x = Power, y = SFT)) +
  geom_text(aes(label = Power), size = 3) +
  geom_hline(yintercept = 0.85, color = "red", linetype = "dashed") +
  geom_vline(xintercept = power, color = "blue", linetype = "dashed") +
  labs(x = "Soft threshold power", y = "Scale-free topology fit (R²)",
       title = "Scale-free topology fit") +
  theme_minimal(base_size = 12)

p8b_conn <- ggplot(sft_df, aes(x = Power, y = MeanK)) +
  geom_text(aes(label = Power), size = 3) +
  geom_vline(xintercept = power, color = "blue", linetype = "dashed") +
  labs(x = "Soft threshold power", y = "Mean connectivity",
       title = "Mean connectivity") +
  theme_minimal(base_size = 12)

p8b <- p8b_sft + p8b_conn + patchwork::plot_layout(ncol = 2)
save_fig(p8b, BATCH_DIR, "fig8b_soft_threshold", w = 12, h = 5)

# Build network
cor <- WGCNA::cor
net <- blockwiseModules(expr_mat, power = power, TOMType = "signed",
                        minModuleSize = 50, mergeCutHeight = 0.25,
                        numericLabels = FALSE, verbose = 0)
cor <- stats::cor

module_colors <- net$colors
n_modules <- length(unique(module_colors)) - 1
cat(sprintf("  Modules: %d (+ %d grey genes)\n", n_modules, sum(module_colors == "grey")))

module_df <- data.frame(gene_id = colnames(expr_mat), module = module_colors, stringsAsFactors = FALSE)
module_sizes <- module_df %>% filter(module != "grey") %>% count(module, name = "size") %>% arrange(desc(size))
cat("  Module sizes:\n"); print(as.data.frame(module_sizes))
save_data(module_df, BATCH_DIR, "wgcna_module_assignments")

# Fig 8C: Module dendrogram
p8c_fn <- function() {
  plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
                      "Module", dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Gene dendrogram and module assignments")
}
png(file.path(fig_dir, "fig8c_module_dendrogram.png"), width = 12, height = 6, units = "in", res = 300)
p8c_fn(); dev.off()
tryCatch({
  cairo_pdf(file.path(fig_dir, "fig8c_module_dendrogram.pdf"), width = 12, height = 6)
  p8c_fn(); dev.off()
}, error = function(e) {
  pdf(file.path(fig_dir, "fig8c_module_dendrogram.pdf"), width = 12, height = 6)
  p8c_fn(); dev.off()
})
cat("  Saved: fig8c_module_dendrogram\n")

# =============================================================================
# 4. MODULE EIGENGENES + MODULE-TRAIT CORRELATION
# =============================================================================
cat("[4/9] Module eigengenes + trait correlations...\n")

MEs <- moduleEigengenes(expr_mat, module_colors)$eigengenes
MEs <- orderMEs(MEs)

# Build binary trait matrix: tissue × condition interaction terms
# Each column = 1 specific tissue+condition combination (no mixing!)
trait_mat <- data.frame(row.names = rownames(expr_mat))
matched_st <- sample_table[match(rownames(expr_mat), sample_table$label), ]

tc_combos <- unique(paste0(matched_st$tissue, "_", matched_st$condition))
tc_combos <- sort(tc_combos)
for (tc in tc_combos) {
  parts <- strsplit(tc, "_")[[1]]
  tis <- parts[1]; cond <- paste(parts[-1], collapse = "_")
  trait_mat[[tc]] <- as.integer(matched_st$tissue == tis & matched_st$condition == cond)
}
trait_mat <- as.data.frame(trait_mat)

# Correlate eigengenes with traits
n_samples <- nrow(expr_mat)
module_trait_cor <- WGCNA::cor(MEs, trait_mat, use = "p")
module_trait_p   <- corPvalueStudent(module_trait_cor, n_samples)

save_data(as.data.frame(module_trait_cor), BATCH_DIR, "module_trait_correlations")
save_data(as.data.frame(module_trait_p), BATCH_DIR, "module_trait_pvalues")

# Fig 8D: Module-trait heatmap
text_matrix <- paste0(signif(module_trait_cor, 2), "\n(",
                       signif(module_trait_p, 1), ")")
dim(text_matrix) <- dim(module_trait_cor)

p8d_fn <- function() {
  par(mar = c(8, 10, 3, 3))
  labeledHeatmap(Matrix = module_trait_cor,
                 xLabels = colnames(trait_mat),
                 yLabels = gsub("^ME", "", colnames(MEs)),
                 ySymbols = gsub("^ME", "", colnames(MEs)),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = text_matrix,
                 setStdMargins = FALSE,
                 cex.text = 0.5,
                 zlim = c(-1, 1),
                 main = "Module-trait relationships")
}
png(file.path(fig_dir, "fig8d_module_trait_heatmap.png"), width = 12, height = 10, units = "in", res = 300)
p8d_fn(); dev.off()
tryCatch({
  cairo_pdf(file.path(fig_dir, "fig8d_module_trait_heatmap.pdf"), width = 12, height = 10)
  p8d_fn(); dev.off()
}, error = function(e) {
  pdf(file.path(fig_dir, "fig8d_module_trait_heatmap.pdf"), width = 12, height = 10)
  p8d_fn(); dev.off()
})
cat("  Saved: fig8d_module_trait_heatmap\n")

# --- Fig 8d2: Module eigengene adjacency/correlation heatmap ---
cat("  Module eigengene correlation heatmap...\n")
me_cor <- WGCNA::cor(MEs, use = "p")
me_dist <- as.dist(1 - me_cor)
me_tree <- hclust(me_dist, method = "average")

p8d2_fn <- function() {
  par(mar = c(1, 6, 3, 1))
  plotEigengeneNetworks(MEs, "Module eigengene relationships",
                        marDendro = c(1, 6, 3, 1), marHeatmap = c(6, 6, 1, 2),
                        plotDendrograms = TRUE, xLabelsAngle = 90)
}
png(file.path(fig_dir, "fig8d2_eigengene_network.png"), width = 10, height = 10, units = "in", res = 300)
p8d2_fn(); dev.off()
tryCatch({
  cairo_pdf(file.path(fig_dir, "fig8d2_eigengene_network.pdf"), width = 10, height = 10)
  p8d2_fn(); dev.off()
}, error = function(e) {
  pdf(file.path(fig_dir, "fig8d2_eigengene_network.pdf"), width = 10, height = 10)
  p8d2_fn(); dev.off()
})
cat("  Saved: fig8d2_eigengene_network\n")

# --- Fig 8d3: Module membership vs gene significance (tail regeneration) ---
# Gene significance = correlation of gene expression with tail_amputated trait
cat("  Module membership vs gene significance...\n")
tail_ampu_col <- grep("tail_amputated", colnames(trait_mat), value = TRUE)
if (length(tail_ampu_col) > 0) {
  gene_sig <- WGCNA::cor(expr_mat, trait_mat[, tail_ampu_col[1], drop = FALSE], use = "p")
  gene_sig_p <- corPvalueStudent(gene_sig, n_samples)
  colnames(gene_sig) <- "GS_tail_amputated"

  # Module membership = correlation with module eigengene
  mm_mat <- WGCNA::cor(expr_mat, MEs, use = "p")
  colnames(mm_mat) <- gsub("^ME", "", colnames(mm_mat))

  # Plot for each non-grey module (top 6 by size)
  top_mods <- head(module_sizes$module, 6)
  mm_gs_plots <- list()
  for (mod in top_mods) {
    mod_genes_idx <- which(module_colors == mod)
    if (length(mod_genes_idx) < 10) next
    df_mmgs <- data.frame(
      MM = mm_mat[mod_genes_idx, mod],
      GS = abs(gene_sig[mod_genes_idx, 1]))
    ct <- cor.test(df_mmgs$MM, df_mmgs$GS, method = "pearson")
    mm_gs_plots[[mod]] <- ggplot(df_mmgs, aes(x = MM, y = GS)) +
      geom_point(alpha = 0.3, size = 1, color = mod) +
      geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
      labs(x = paste0("Module membership (", mod, ")"),
           y = "|GS: tail amputated|",
           title = mod,
           subtitle = sprintf("r = %.2f, p = %s", ct$estimate,
                              format(ct$p.value, digits = 2))) +
      theme_minimal(base_size = 10)
  }
  if (length(mm_gs_plots) > 0) {
    p8d3 <- patchwork::wrap_plots(mm_gs_plots, ncol = 3)
    save_fig(p8d3, BATCH_DIR, "fig8d3_mm_vs_gs", w = 14, h = 10)
    cat("  Saved: fig8d3_mm_vs_gs\n")
  }
} else {
  cat("  WARNING: tail_amputated trait not found, skipping MM vs GS plot.\n")
}

# =============================================================================
# 5. WILCOXON: MODULE EIGENGENES CONTROL vs AMPUTATED (TAIL ONLY)
# =============================================================================
cat("[5/9] Wilcoxon test: eigengenes control vs amputated (tail)...\n")

tail_idx <- which(matched_st$tissue == "tail" &
                  matched_st$condition %in% c("control", "amputated"))
tail_cond <- matched_st$condition[tail_idx]

wilcox_results <- data.frame(module = gsub("^ME", "", colnames(MEs)),
                             stringsAsFactors = FALSE)
wilcox_results$wilcox_p <- NA_real_
wilcox_results$cohens_d <- NA_real_
wilcox_results$direction <- NA_character_

for (i in seq_len(ncol(MEs))) {
  me_vals <- MEs[tail_idx, i]
  ctrl_vals <- me_vals[tail_cond == "control"]
  ampu_vals <- me_vals[tail_cond == "amputated"]
  if (length(ctrl_vals) >= 2 && length(ampu_vals) >= 2) {
    wt <- wilcox.test(ctrl_vals, ampu_vals)
    wilcox_results$wilcox_p[i] <- wt$p.value
    pooled_sd <- sqrt((var(ctrl_vals) + var(ampu_vals)) / 2)
    if (pooled_sd > 0) {
      wilcox_results$cohens_d[i] <- (mean(ampu_vals) - mean(ctrl_vals)) / pooled_sd
    }
    wilcox_results$direction[i] <- ifelse(mean(ampu_vals) > mean(ctrl_vals), "Up", "Down")
  }
}
wilcox_results$wilcox_padj <- p.adjust(wilcox_results$wilcox_p, method = "BH")
wilcox_results <- wilcox_results[order(wilcox_results$wilcox_p), ]

cat("  Wilcoxon results (top 5):\n")
print(head(wilcox_results[, c("module", "wilcox_p", "wilcox_padj", "cohens_d", "direction")], 5))
save_data(wilcox_results, BATCH_DIR, "module_eigengene_wilcoxon")

# Fig 8E: Module eigengene bar (tail control vs amputated)
me_tail <- data.frame(MEs[tail_idx, ])
me_tail$condition <- tail_cond
me_long <- me_tail %>%
  pivot_longer(-condition, names_to = "module", values_to = "eigengene") %>%
  mutate(module = gsub("^ME", "", module))

me_summary <- me_long %>%
  group_by(module, condition) %>%
  summarize(mean_me = mean(eigengene), se = sd(eigengene) / sqrt(n()), .groups = "drop")

sig_mods <- wilcox_results$module[!is.na(wilcox_results$wilcox_padj) &
                                   wilcox_results$wilcox_padj < 0.2]
if (length(sig_mods) < 3) sig_mods <- head(wilcox_results$module, 6)

p8e <- ggplot(me_summary %>% filter(module %in% sig_mods),
              aes(x = module, y = mean_me, fill = condition)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = mean_me - se, ymax = mean_me + se),
                position = position_dodge(width = 0.7), width = 0.2) +
  scale_fill_manual(values = c(control = "#2471A3", amputated = "#C0392B")) +
  labs(x = "Module", y = "Mean eigengene", fill = "Condition",
       title = "Module eigengenes: tail control vs amputated",
       subtitle = "Modules with Wilcoxon padj < 0.2 or top 6") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p8e, BATCH_DIR, "fig8e_eigengene_control_vs_amputated", w = 10, h = 6)

# =============================================================================
# 6. DMP ENRICHMENT PER MODULE (FISHER'S EXACT)
# =============================================================================
cat("[6/9] Fisher's exact test: DMP enrichment per module...\n")

dmp_file <- file.path(PIPE_DIR, "batch06/data/dmps_annotated.tsv")
if (!file.exists(dmp_file)) stop("Run batch06 first: dmps_annotated.tsv not found.")
dmp <- fread(dmp_file)

gene_dmp <- dmp[, .(n_dmp = .N, mean_diff = mean(diff),
                     n_hyper = sum(diff > 0), n_hypo = sum(diff < 0)),
                 by = .(gene_id = nearest_gene)]
gene_dmp <- gene_dmp[!is.na(gene_id) & gene_id != ""]

module_meth <- merge(module_df, gene_dmp, by = "gene_id", all.x = TRUE)
module_meth$has_dmp <- !is.na(module_meth$n_dmp) & module_meth$n_dmp > 0
module_meth$n_dmp[is.na(module_meth$n_dmp)] <- 0

real_modules <- setdiff(unique(module_df$module), "grey")
total_in_network <- sum(module_df$module != "grey")
total_dmp_genes  <- sum(module_meth$has_dmp & module_meth$module != "grey")

fisher_list <- list()
for (mod in real_modules) {
  mod_genes <- module_meth[module_meth$module == mod, ]
  mod_size  <- nrow(mod_genes)
  mod_dmp   <- sum(mod_genes$has_dmp)

  a <- mod_dmp
  b <- mod_size - mod_dmp
  c_val <- total_dmp_genes - mod_dmp
  d_val <- total_in_network - mod_size - c_val

  ft <- tryCatch(
    fisher.test(matrix(c(a, c_val, b, d_val), nrow = 2), alternative = "greater"),
    error = function(e) list(p.value = NA, estimate = NA))

  fisher_list[[mod]] <- data.frame(
    module = mod, size = mod_size, genes_with_dmp = mod_dmp,
    pct_with_dmp = round(100 * mod_dmp / mod_size, 1),
    total_dmps = sum(mod_genes$n_dmp),
    n_hyper = sum(mod_genes$n_hyper, na.rm = TRUE),
    n_hypo = sum(mod_genes$n_hypo, na.rm = TRUE),
    fisher_p = ft$p.value, odds_ratio = as.numeric(ft$estimate),
    stringsAsFactors = FALSE)
}
fisher_df <- do.call(rbind, fisher_list)
fisher_df$fisher_padj <- p.adjust(fisher_df$fisher_p, method = "BH")
fisher_df <- fisher_df[order(fisher_df$fisher_p), ]

fisher_df <- merge(fisher_df,
                   wilcox_results[, c("module", "wilcox_p", "wilcox_padj", "cohens_d", "direction")],
                   by = "module", all.x = TRUE)

cat("  Fisher enrichment (all modules):\n")
print(fisher_df[, c("module", "size", "genes_with_dmp", "pct_with_dmp",
                     "fisher_padj", "odds_ratio", "wilcox_padj")])
save_data(fisher_df, BATCH_DIR, "module_fisher_dmp_enrichment")

# --- DMR enrichment per module (Fisher's exact) ---
cat("  Fisher's exact: DMR enrichment per module...\n")
dmr_file <- file.path(PIPE_DIR, "batch06/data/dmrs_annotated.tsv")
if (file.exists(dmr_file)) {
  dmr <- fread(dmr_file)
  gene_dmr <- dmr[, .(n_dmr = .N, dmr_mean_diff = mean(diff.Methy, na.rm = TRUE),
                       dmr_hyper = sum(diff.Methy > 0, na.rm = TRUE),
                       dmr_hypo = sum(diff.Methy < 0, na.rm = TRUE)),
                   by = .(gene_id = nearest_gene)]
  gene_dmr <- gene_dmr[!is.na(gene_id) & gene_id != ""]

  module_dmr <- merge(module_df, gene_dmr, by = "gene_id", all.x = TRUE)
  module_dmr$has_dmr <- !is.na(module_dmr$n_dmr) & module_dmr$n_dmr > 0
  module_dmr$n_dmr[is.na(module_dmr$n_dmr)] <- 0

  total_dmr_genes <- sum(module_dmr$has_dmr & module_dmr$module != "grey")

  fisher_dmr_list <- list()
  for (mod in real_modules) {
    mod_genes <- module_dmr[module_dmr$module == mod, ]
    mod_size  <- nrow(mod_genes)
    mod_dmr   <- sum(mod_genes$has_dmr)
    a <- mod_dmr; b <- mod_size - mod_dmr
    c_val <- total_dmr_genes - mod_dmr
    d_val <- total_in_network - mod_size - c_val
    ft <- tryCatch(
      fisher.test(matrix(c(a, c_val, b, d_val), nrow = 2), alternative = "greater"),
      error = function(e) list(p.value = NA, estimate = NA))
    fisher_dmr_list[[mod]] <- data.frame(
      module = mod, size = mod_size, genes_with_dmr = mod_dmr,
      pct_with_dmr = round(100 * mod_dmr / mod_size, 1),
      total_dmrs = sum(mod_genes$n_dmr),
      dmr_hyper = sum(mod_genes$dmr_hyper, na.rm = TRUE),
      dmr_hypo = sum(mod_genes$dmr_hypo, na.rm = TRUE),
      fisher_p = ft$p.value, odds_ratio = as.numeric(ft$estimate),
      stringsAsFactors = FALSE)
  }
  fisher_dmr_df <- do.call(rbind, fisher_dmr_list)
  fisher_dmr_df$fisher_padj <- p.adjust(fisher_dmr_df$fisher_p, method = "BH")
  fisher_dmr_df <- fisher_dmr_df[order(fisher_dmr_df$fisher_p), ]
  cat("  DMR Fisher enrichment (all modules):\n")
  print(fisher_dmr_df[, c("module", "size", "genes_with_dmr", "pct_with_dmr",
                           "fisher_padj", "odds_ratio")])
  save_data(fisher_dmr_df, BATCH_DIR, "module_fisher_dmr_enrichment")

  # Merge DMP + DMR into combined enrichment for plotting
  fisher_df <- merge(fisher_df,
    fisher_dmr_df[, c("module", "genes_with_dmr", "pct_with_dmr", "total_dmrs",
                       "dmr_hyper", "dmr_hypo")],
    by = "module", all.x = TRUE, suffixes = c("", "_dmr"))
} else {
  cat("  WARNING: dmrs_annotated.tsv not found, skipping DMR enrichment.\n")
}

# Fig 8F: DMP burden bar
p8f <- ggplot(fisher_df, aes(x = reorder(module, -pct_with_dmp),
                              y = pct_with_dmp, fill = module)) +
  geom_col(width = 0.7) + scale_fill_identity() +
  geom_text(aes(label = sprintf("%.1f%%", pct_with_dmp)), vjust = -0.3, size = 3) +
  labs(x = "Module", y = "% genes with DMPs",
       title = "DMP burden by WGCNA module",
       subtitle = sprintf("Network: %s genes | DMPs: %s sites",
                          format(total_in_network, big.mark = ","),
                          format(nrow(dmp), big.mark = ","))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p8f, BATCH_DIR, "fig8f_module_dmp_burden", w = 10, h = 6)

# Fig 8G: Fisher enrichment dot plot
p8g <- ggplot(fisher_df, aes(x = reorder(module, -log10(fisher_padj)),
                              y = -log10(fisher_padj),
                              size = genes_with_dmp, color = module)) +
  geom_point(alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_identity() +
  scale_size_continuous(name = "Genes with DMP", range = c(2, 10)) +
  coord_flip() +
  labs(x = "Module", y = "-log10(adjusted p-value)",
       title = "DMP enrichment per WGCNA module (Fisher's exact)",
       subtitle = "Red dashed: padj = 0.05") +
  theme_minimal(base_size = 12)
save_fig(p8g, BATCH_DIR, "fig8g_fisher_dmp_enrichment", w = 10, h = 7)

# Fig 8H: Bubble — Wilcoxon x DMP burden
p8h <- ggplot(fisher_df, aes(x = -log10(wilcox_p), y = pct_with_dmp,
                              size = total_dmps, color = module)) +
  geom_point(alpha = 0.8) +
  scale_color_identity() +
  scale_size_continuous(name = "Total DMPs", range = c(2, 12)) +
  ggrepel::geom_text_repel(aes(label = module), size = 3.5, max.overlaps = 15) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
  labs(x = "-log10(Wilcoxon p) expression change in tail regeneration",
       y = "% module genes with DMPs",
       title = "Module expression change vs DMP burden") +
  theme_minimal(base_size = 12)
save_fig(p8h, BATCH_DIR, "fig8h_wilcoxon_vs_dmp_burden", w = 11, h = 8)

# =============================================================================
# 7. HUB GENE IDENTIFICATION
# =============================================================================
cat("[7/9] Identifying hub genes...\n")

kME <- WGCNA::cor(expr_mat, MEs, use = "p")
colnames(kME) <- gsub("^ME", "", colnames(kME))

adj <- adjacency(expr_mat, power = power, type = "signed")
connectivity <- intramodularConnectivity(adj, module_colors)
rm(adj); gc(verbose = FALSE)

hub_df <- data.frame(
  gene_id = colnames(expr_mat),
  module = module_colors,
  kTotal = connectivity$kTotal,
  kWithin = connectivity$kWithin,
  stringsAsFactors = FALSE
)

hub_df$kME <- sapply(seq_len(nrow(hub_df)), function(i) {
  mod <- hub_df$module[i]
  if (mod == "grey" || !(mod %in% colnames(kME))) return(NA_real_)
  kME[i, mod]
})

hub_df <- hub_df %>%
  group_by(module) %>%
  mutate(kWithin_rank = percent_rank(kWithin),
         is_hub = kWithin_rank >= 0.9 & abs(kME) > 0.7 & module != "grey") %>%
  ungroup()

hub_df <- merge(hub_df, gene_dmp[, c("gene_id", "n_dmp", "mean_diff")],
                by = "gene_id", all.x = TRUE)
hub_df$has_dmp <- !is.na(hub_df$n_dmp)

# Add gene names if annotation available
if (file.exists(OG$annot)) {
  annot <- fread(OG$annot, header = TRUE)
  name_col <- intersect(c("Name", "name", "gene_name", "display_name"), colnames(annot))
  id_col   <- intersect(c("ID", "gene_id", "#gene_id"), colnames(annot))
  if (length(name_col) > 0 && length(id_col) > 0) {
    hub_df <- merge(hub_df, annot[, c(id_col[1], name_col[1]), with = FALSE],
                    by.x = "gene_id", by.y = id_col[1], all.x = TRUE)
  }
}

hubs_only <- hub_df[hub_df$is_hub == TRUE, ]
cat(sprintf("  Hub genes: %d / %d (%.1f%%)\n",
            nrow(hubs_only), sum(module_colors != "grey"),
            100 * nrow(hubs_only) / sum(module_colors != "grey")))
cat(sprintf("  Hub genes with DMPs: %d (%.1f%%)\n",
            sum(hubs_only$has_dmp), 100 * mean(hubs_only$has_dmp)))

save_data(hub_df, BATCH_DIR, "hub_genes_all")
save_data(hubs_only, BATCH_DIR, "hub_genes_identified")

hub_summary <- hubs_only %>%
  group_by(module) %>%
  summarize(n_hubs = n(), hubs_with_dmp = sum(has_dmp),
            pct_hub_dmp = round(100 * mean(has_dmp), 1),
            mean_kME = round(mean(kME, na.rm = TRUE), 3), .groups = "drop") %>%
  arrange(desc(hubs_with_dmp))
save_data(hub_summary, BATCH_DIR, "hub_gene_summary_per_module")

# =============================================================================
# 8. MODULE GO ENRICHMENT (STRING-BASED)
# =============================================================================
cat("[8/9] Module GO enrichment (STRING)...\n")

if (file.exists(OG$string_enrich)) {
  string_raw <- fread(OG$string_enrich, sep = "\t", header = TRUE, quote = "")
  string_raw[, gene_id := sub("^.*\\.(LOC_\\d+)$", "\\1", `#string_protein_id`)]
  string_raw[!grepl("^LOC_", gene_id), gene_id := sub("^.*\\.(XLOC_\\S+)$", "\\1", `#string_protein_id`)]
  string_raw <- string_raw[grepl("^LOC_|^XLOC_", gene_id)]

  string_df <- as.data.frame(string_raw)
  go_bp <- string_df[grepl("Biological Process", string_df$category), ]
  T2G_bp <- unique(go_bp[, c("term", "gene_id")])
  T2N_bp <- unique(go_bp[, c("term", "description")])
  rm(string_raw, string_df, go_bp); gc(verbose = FALSE)

  cat(sprintf("  GO BP terms: %d\n", nrow(T2N_bp)))

  all_genes_universe <- unique(module_df$gene_id)
  go_results_list <- list()

  for (mod in real_modules) {
    mod_genes <- module_df$gene_id[module_df$module == mod]
    if (length(mod_genes) < 20) next

    res <- tryCatch(
      enricher(gene = mod_genes, TERM2GENE = T2G_bp, TERM2NAME = T2N_bp,
               pvalueCutoff = 0.05, pAdjustMethod = "BH",
               minGSSize = 5, maxGSSize = 500, universe = all_genes_universe),
      error = function(e) NULL)

    if (!is.null(res) && nrow(res@result[res@result$p.adjust < 0.05, ]) > 0) {
      sig <- res@result[res@result$p.adjust < 0.05, ]
      sig$module <- mod
      go_results_list[[mod]] <- sig
      cat(sprintf("    %s: %d significant GO BP terms\n", mod, nrow(sig)))
    }
  }

  if (length(go_results_list) > 0) {
    go_all <- do.call(rbind, go_results_list)
    save_data(go_all, BATCH_DIR, "module_go_bp_enrichment")

    # Fig 8I: Top GO terms per DMP-enriched modules
    enriched_mods <- fisher_df$module[fisher_df$fisher_padj < 0.1]
    if (length(enriched_mods) == 0) enriched_mods <- head(fisher_df$module[order(fisher_df$fisher_p)], 3)

    go_top <- go_all %>%
      filter(module %in% enriched_mods) %>%
      group_by(module) %>%
      arrange(p.adjust) %>%
      slice_head(n = 5) %>%
      ungroup()

    if (nrow(go_top) > 0) {
      p8i <- ggplot(go_top, aes(x = Count, y = reorder(Description, Count),
                                 fill = module)) +
        geom_col() + scale_fill_identity() +
        facet_wrap(~ module, scales = "free_y", ncol = 1) +
        labs(x = "Gene count", y = NULL,
             title = "GO Biological Process: DMP-enriched modules",
             subtitle = "Top 5 terms per module") +
        theme_minimal(base_size = 11) +
        theme(axis.text.y = element_text(size = 8), strip.text = element_text(face = "bold"))
      save_fig(p8i, BATCH_DIR, "fig8i_module_go_enrichment",
               w = 12, h = max(5, nrow(go_top) * 0.5))
    }
  }
} else {
  cat("  STRING enrichment file not found, skipping GO.\n")
}

# =============================================================================
# 9. IDENTIFY TISSUE-SPECIFIC MODULES + SAVE
# =============================================================================
cat("[9/9] Identifying tissue-specific modules...\n")

ovo_col <- grep("ovotestis", colnames(trait_mat), value = TRUE)
if (length(ovo_col) > 0) {
  ovo_cors <- module_trait_cor[, ovo_col[1]]
  names(ovo_cors) <- gsub("^ME", "", rownames(module_trait_cor))
  best_ovo <- names(which.max(abs(ovo_cors)))
  cat(sprintf("  Module most correlated with %s: %s (r = %.3f, p = %s)\n",
              ovo_col[1], best_ovo, ovo_cors[best_ovo],
              signif(module_trait_p[paste0("ME", best_ovo), ovo_col[1]], 3)))

  ovo_row <- fisher_df[fisher_df$module == best_ovo, ]
  if (nrow(ovo_row) > 0) {
    cat(sprintf("  %s module: %d genes, %d with DMP (%.1f%%), Fisher padj = %s\n",
                best_ovo, ovo_row$size, ovo_row$genes_with_dmp,
                ovo_row$pct_with_dmp, signif(ovo_row$fisher_padj, 3)))
  }
}

# Save module eigengenes
me_df <- as.data.frame(MEs)
me_df$sample <- rownames(MEs)
me_df$tissue <- matched_st$tissue[match(me_df$sample, matched_st$label)]
me_df$condition <- matched_st$condition[match(me_df$sample, matched_st$label)]
save_data(me_df, BATCH_DIR, "module_eigengenes")

# Fig 8J: Module size bar
p8j <- ggplot(module_sizes, aes(x = reorder(module, -size), y = size, fill = module)) +
  geom_col() + scale_fill_identity() +
  geom_text(aes(label = size), vjust = -0.3, size = 3) +
  labs(x = "Module", y = "Number of genes",
       title = sprintf("WGCNA module sizes (%d modules, power = %d)", n_modules, power)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p8j, BATCH_DIR, "fig8j_module_sizes", w = 10, h = 6)

# Save key objects for batch07
saveRDS(list(module_df = module_df, fisher_df = fisher_df,
             wilcox_results = wilcox_results, hub_df = hub_df,
             MEs = MEs, module_trait_cor = module_trait_cor,
             module_meth = module_meth),
        file.path(BATCH_DIR, "data", "wgcna_objects.rds"))

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("\n=== Batch 08 complete (%.1f min) ===\n", elapsed / 60))
n_figs <- length(list.files(file.path(BATCH_DIR, "figures"), pattern = "\\.png$"))
n_data <- length(list.files(file.path(BATCH_DIR, "data"), pattern = "\\.tsv$"))
cat(sprintf("Figures: %d | Data files: %d\n", n_figs, n_data))
