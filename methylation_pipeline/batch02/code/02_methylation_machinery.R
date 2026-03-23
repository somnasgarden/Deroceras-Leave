#!/usr/bin/env Rscript
# =============================================================================
# Batch 02: Methylation Machinery
# Question: Does D. laeve have a functional methylation toolkit? Where expressed?
# Output: data/ (gene table) + figures/ (3 plots)
# =============================================================================

source("methylation_pipeline/_config.R")

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(rtracklayer)
library(RColorBrewer)
library(scales)
library(grid)
library(pheatmap)

BATCH_DIR <- file.path(PIPE_DIR, "batch02")

cat("=== Batch 02: Methylation Machinery ===\n\n")

# --- Load data ---
gff <- load_gff()
annot <- data.table::fread(OG$annot)

# --- Curated gene list (18 methylation machinery genes) ---
meth_genes <- data.frame(
  Gene = c("DNMT1", "TET2", "TET3", "CFP1", "MBD2", "MBD4",
           "UHRF1_1", "UHRF1_2", "APEX1", "TDG_1", "TDG_2",
           "GADD45B", "GADD45G_1", "GADD45G_2", "GADD45GIP1",
           "HELLS_1", "HELLS_2", "HDAC2"),
  Gene_ID = c("LOC_00009210", "LOC_00008803", "LOC_00017436",
              "LOC_00024130", "LOC_00005931", "LOC_00023929",
              "LOC_00011698", "LOC_00017302", "LOC_00000327",
              "LOC_00018573", "LOC_00024530", "LOC_00012796",
              "LOC_00002628", "LOC_00012923", "LOC_00002413",
              "LOC_00010497", "LOC_00013941", "LOC_00003192"),
  Category = c("Writer", "Eraser", "Eraser", "Reader", "Reader", "Reader",
               "Recruiter", "Recruiter", "BER", "BER", "BER",
               "Demethylation", "Demethylation", "Demethylation", "Demethylation",
               "Remodeler", "Remodeler", "Deacetylase"),
  stringsAsFactors = FALSE
)

# Verify in GFF
genes_gff <- gff[gff$type == "gene"]
found <- meth_genes$Gene_ID %in% genes_gff$ID
cat(sprintf("Genes found in GFF: %d / %d\n\n", sum(found), nrow(meth_genes)))

# --- DESeq2 for expression ---
sample_table <- data.frame(
  sampleName = list.files(OG$counts_dir, pattern = "\\.txt$"),
  stringsAsFactors = FALSE
)
sample_table$condition <- ifelse(grepl("control|ctrl|C\\d+S", sample_table$sampleName, ignore.case = TRUE),
                                  "control", NA)
sample_table$condition <- ifelse(grepl("irradiat", sample_table$sampleName, ignore.case = TRUE),
                                  "irradiated", sample_table$condition)
sample_table$condition <- ifelse(grepl("amput|T\\d+S", sample_table$sampleName, ignore.case = TRUE),
                                  "amputated", sample_table$condition)
sample_table$condition <- ifelse(grepl("fungi", sample_table$sampleName, ignore.case = TRUE),
                                  "fungicide", sample_table$condition)
sample_table$tissue <- sub("_.*", "", sub("^(tail|eye|bodywall|head|juvenile|ovotestis).*", "\\1",
                           tolower(sample_table$sampleName)))
sample_table$fileName <- sample_table$sampleName
sample_table <- sample_table[!is.na(sample_table$condition), ]

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table,
                                   directory = OG$counts_dir,
                                   design = ~ condition)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vst_mat <- assay(vsd)

cat(sprintf("DESeq2 done. Genes retained: %d\n\n", nrow(vst_mat)))

# --- Extract expression for methylation genes ---
found_ids <- meth_genes$Gene_ID[meth_genes$Gene_ID %in% rownames(vst_mat)]
expr_mat <- vst_mat[found_ids, ]

# Build tissue means (controls only)
ctrl_samples <- sample_table$sampleName[sample_table$condition == "control"]
tissue_means <- do.call(cbind, lapply(unique(sample_table$tissue[sample_table$condition == "control"]), function(t) {
  samps <- sample_table$sampleName[sample_table$tissue == t & sample_table$condition == "control"]
  if (length(samps) > 1) rowMeans(expr_mat[, samps, drop = FALSE])
  else expr_mat[, samps, drop = TRUE]
}))
colnames(tissue_means) <- tools::toTitleCase(unique(sample_table$tissue[sample_table$condition == "control"]))

expr_df <- as.data.frame(tissue_means[found_ids, , drop = FALSE])
expr_df$Gene_ID <- rownames(expr_df)
expr_df <- merge(meth_genes, expr_df, by = "Gene_ID")

save_data(expr_df, BATCH_DIR, "methylation_machinery_expression")

# --- Fig 2A: Presence/absence ---
presence_df <- meth_genes
presence_df$Present <- meth_genes$Gene_ID %in% genes_gff$ID
presence_df$Gene <- factor(presence_df$Gene, levels = rev(meth_genes$Gene))

p2a <- ggplot(presence_df, aes(x = 1, y = Gene, fill = Present)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_manual(values = c("TRUE" = "#27AE60", "FALSE" = "#2C3E50"), guide = "none") +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Methylation toolkit in D. laeve genome", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), strip.text.y = element_text(angle = 0))
save_fig(p2a, BATCH_DIR, "fig2a_methylation_toolkit_presence", w = 6, h = 8)

# --- Fig 2B: Expression boxplot ---
expr_long <- expr_df %>%
  pivot_longer(cols = -c(Gene_ID, Gene, Category), names_to = "Tissue", values_to = "VST")

p2b <- ggplot(expr_long, aes(x = Gene, y = VST, fill = Category)) +
  geom_boxplot(outlier.size = 0.5) +
  coord_flip() +
  labs(title = "Methylation machinery expression by tissue",
       subtitle = "Controls only, VST-normalized", y = "VST expression", x = NULL) +
  theme_minimal() + theme(legend.position = "bottom")
save_fig(p2b, BATCH_DIR, "fig2b_methylation_machinery_expression", w = 10, h = 7)

# --- Fig 2C: Heatmap ---
hm_mat <- as.matrix(expr_df[, colnames(tissue_means)])
rownames(hm_mat) <- expr_df$Gene
hm_z <- t(scale(t(hm_mat)))

ann_row <- data.frame(Category = expr_df$Category, row.names = expr_df$Gene)
p2c <- pheatmap(hm_z, annotation_row = ann_row,
                color = colorRampPalette(c("#2471A3", "white", "#C0392B"))(100),
                main = "Methylation toolkit expression (z-scored)",
                clustering_method = "ward.D2", fontsize = 10, silent = TRUE)

png(file.path(BATCH_DIR, "figures", "fig2c_methylation_toolkit_heatmap.png"),
    width = 8, height = 10, units = "in", res = 300)
grid::grid.draw(p2c$gtable); dev.off()
tryCatch({
  cairo_pdf(file.path(BATCH_DIR, "figures", "fig2c_methylation_toolkit_heatmap.pdf"), width = 8, height = 10)
  grid::grid.draw(p2c$gtable); dev.off()
}, error = function(e) {
  pdf(file.path(BATCH_DIR, "figures", "fig2c_methylation_toolkit_heatmap.pdf"), width = 8, height = 10)
  grid::grid.draw(p2c$gtable); dev.off()
})
cat("  Saved: fig2c_methylation_toolkit_heatmap\n")

cat("\n=== Batch 02 complete ===\n")
