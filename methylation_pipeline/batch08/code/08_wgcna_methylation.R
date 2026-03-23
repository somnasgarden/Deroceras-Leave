#!/usr/bin/env Rscript
# =============================================================================
# Batch 08: WGCNA Modules + Methylation Enrichment
# Question: Are methylation changes concentrated in specific expression modules?
# Output: data/ (module assignments, enrichment) + figures/ (3 plots)
# Requires: Batch 06 DMPs, DESeq2 multi-tissue counts
# =============================================================================

source("methylation_pipeline/_config.R")

library(WGCNA)
library(DESeq2)
library(data.table)
library(ggplot2)
library(dplyr)

BATCH_DIR <- file.path(PIPE_DIR, "batch08")
cat("=== Batch 08: WGCNA + Methylation ===\n\n")
allowWGCNAThreads(4)

# --- Build multi-tissue expression matrix ---
cat("Building multi-tissue expression matrix...\n")
all_files <- list.files(OG$counts_dir, pattern = "\\.txt$")

sample_table <- data.frame(sampleName = all_files, fileName = all_files, stringsAsFactors = FALSE)
sample_table$condition <- "unknown"
sample_table$condition[grepl("control", sample_table$sampleName, ignore.case = TRUE)] <- "control"
sample_table$condition[grepl("amput", sample_table$sampleName, ignore.case = TRUE)] <- "amputated"
sample_table$condition[grepl("irrad", sample_table$sampleName, ignore.case = TRUE)] <- "irradiated"
sample_table$condition[grepl("fungi", sample_table$sampleName, ignore.case = TRUE)] <- "fungicide"
sample_table <- sample_table[sample_table$condition != "unknown", ]

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table, directory = OG$counts_dir, design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 10, ]
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
expr_mat <- t(assay(vsd))  # samples x genes

cat(sprintf("Expression matrix: %d samples x %d genes\n", nrow(expr_mat), ncol(expr_mat)))

# --- WGCNA ---
cat("Running WGCNA...\n")
powers <- c(1:20)
sft <- pickSoftThreshold(expr_mat, powerVector = powers, verbose = 0)
power <- sft$powerEstimate
if (is.na(power)) power <- 6
cat(sprintf("Soft threshold power: %d\n", power))

net <- blockwiseModules(expr_mat, power = power, TOMType = "signed",
                        minModuleSize = 50, mergeCutHeight = 0.25,
                        numericLabels = FALSE, verbose = 0)

module_colors <- net$colors
n_modules <- length(unique(module_colors)) - 1  # exclude grey
cat(sprintf("Modules detected: %d (+ grey)\n", n_modules))

module_df <- data.frame(gene_id = colnames(expr_mat), module = module_colors, stringsAsFactors = FALSE)
module_sizes <- module_df %>% filter(module != "grey") %>% count(module, name = "size") %>% arrange(desc(size))
cat("Module sizes:\n"); print(as.data.frame(module_sizes))

save_data(module_df, BATCH_DIR, "wgcna_module_assignments")

# --- Merge with DMP data ---
dmp_file <- file.path(PIPE_DIR, "batch06/data/dmps_annotated.tsv")
if (file.exists(dmp_file)) {
  dmp <- fread(dmp_file)
  gene_dmp <- dmp %>%
    group_by(nearest_gene) %>%
    summarize(n_dmp = n(), mean_diff = mean(diff), .groups = "drop") %>%
    rename(gene_id = nearest_gene)

  module_meth <- merge(module_df, gene_dmp, by = "gene_id", all.x = TRUE)
  module_meth$has_dmp <- !is.na(module_meth$n_dmp)
  module_meth$n_dmp[is.na(module_meth$n_dmp)] <- 0

  # Per-module DMP burden
  module_burden <- module_meth %>%
    filter(module != "grey") %>%
    group_by(module) %>%
    summarize(size = n(), genes_with_dmp = sum(has_dmp),
              pct_with_dmp = round(100 * mean(has_dmp), 1),
              total_dmps = sum(n_dmp), .groups = "drop") %>%
    arrange(desc(pct_with_dmp))

  cat("\nModule DMP burden:\n"); print(as.data.frame(module_burden))
  save_data(module_burden, BATCH_DIR, "module_dmp_burden")

  # Fig 8A: Module DMP burden
  p8a <- ggplot(module_burden, aes(x = reorder(module, -pct_with_dmp), y = pct_with_dmp, fill = module)) +
    geom_col(width = 0.7) + scale_fill_identity() +
    labs(x = "Module", y = "% genes with DMPs", title = "DMP burden by WGCNA module") +
    theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p8a, BATCH_DIR, "fig8a_module_dmp_burden", w = 9, h = 6)
}

cat("\n=== Batch 08 complete ===\n")
