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

# --- Curated gene list: present + absent methylation machinery ---
# Present genes (18 found in D. laeve genome)
present_genes <- data.frame(
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
# Absent genes — key absences that define the methylation system
absent_genes <- data.frame(
  Gene = c("DNMT3A", "DNMT3B", "DNMT3L", "DNMT2", "TET1", "MBD1", "MBD3", "MECP2"),
  Gene_ID = NA_character_,
  Category = c("Writer (de novo)", "Writer (de novo)", "Writer (de novo)", "Writer (tRNA)",
               "Eraser", "Reader", "Reader", "Reader"),
  stringsAsFactors = FALSE
)
meth_genes <- rbind(present_genes, absent_genes)

# Verify in GFF
genes_gff <- gff[gff$type == "gene"]
found <- !is.na(meth_genes$Gene_ID) & meth_genes$Gene_ID %in% genes_gff$ID
meth_genes$Present <- found
cat(sprintf("Genes present: %d / %d (absent: %d)\n\n",
            sum(found), nrow(meth_genes), sum(!found)))

# --- Expression data (prefer CACHE$transcriptome from batch0.5) ---
if (file.exists(CACHE$transcriptome)) {
  cat("Loading expression from CACHE$transcriptome...\n")
  tc <- readRDS(CACHE$transcriptome)
  vst_mat <- tc$vst_clean
  sample_table <- tc$meta_clean
  cat(sprintf("  VST matrix: %d genes x %d samples (outliers removed)\n",
              nrow(vst_mat), ncol(vst_mat)))
  rm(tc); gc(verbose = FALSE)
} else {
  cat("Cache not found — running DESeq2...\n")
  all_files <- list.files(OG$counts_dir, pattern = "\\.txt$")

  cond <- rep(NA_character_, length(all_files))
  cond[grepl("^C\\d+S\\d+_", all_files)]    <- "control"
  cond[grepl("^T\\d+S\\d+_", all_files)]    <- "amputated"
  cond[grepl("^C\\d+_htseq", all_files)]    <- "control"
  cond[grepl("^dcrep", all_files)]           <- "control"
  cond[grepl("^R\\d+_", all_files)]          <- "amputated"
  cond[grepl("^irrep", all_files)]           <- "irradiated"
  cond[grepl("^fungicide_l0", all_files)]    <- "control"
  cond[grepl("^fungicide_l30", all_files)]   <- "fungicide"
  cond[grepl("fpDLHead", all_files)]         <- "control"
  cond[grepl("fpDLJuv", all_files)]          <- "control"
  cond[grepl("fpDLOvo", all_files)]          <- "control"

  tissue <- rep(NA_character_, length(all_files))
  tissue[grepl("^C\\d+S\\d+_|^T\\d+S\\d+_", all_files)] <- "tail"
  tissue[grepl("^C\\d+_htseq", all_files)]   <- "eye"
  tissue[grepl("^R\\d+_", all_files)]        <- "eye"
  tissue[grepl("^dcrep", all_files)]          <- "bodywall"
  tissue[grepl("^irrep", all_files)]         <- "bodywall"
  tissue[grepl("^fungicide", all_files)]      <- "bodywall"
  tissue[grepl("fpDLHead", all_files)]        <- "head"
  tissue[grepl("fpDLJuv", all_files)]         <- "juvenile"
  tissue[grepl("fpDLOvo", all_files)]         <- "ovotestis"

  sample_table <- data.frame(
    sampleName = all_files, fileName = all_files,
    condition = cond, tissue = tissue, stringsAsFactors = FALSE)
  sample_table <- sample_table[!is.na(sample_table$condition), ]

  dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table,
                                     directory = OG$counts_dir, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  vst_mat <- assay(vsd)
  rm(dds, vsd); gc(verbose = FALSE)
}
cat(sprintf("Genes with expression: %d\n\n", nrow(vst_mat)))

# --- Extract expression for present genes only ---
found_ids <- meth_genes$Gene_ID[meth_genes$Present]
found_ids <- found_ids[found_ids %in% rownames(vst_mat)]
expr_mat <- vst_mat[found_ids, ]

# Build tissue × condition means (show both control and amputated separately)
sample_id_col <- if ("label" %in% names(sample_table)) "label" else "fileName"
tissue_cond <- unique(sample_table[, c("tissue", "condition")])
tissue_cond <- tissue_cond[order(tissue_cond$tissue, tissue_cond$condition), ]
tissue_means_list <- lapply(seq_len(nrow(tissue_cond)), function(i) {
  t <- tissue_cond$tissue[i]; cond <- tissue_cond$condition[i]
  samps <- sample_table[[sample_id_col]][sample_table$tissue == t & sample_table$condition == cond]
  samps <- intersect(samps, colnames(expr_mat))
  if (length(samps) == 0) return(NULL)
  if (length(samps) > 1) rowMeans(expr_mat[, samps, drop = FALSE])
  else expr_mat[, samps]
})
keep <- !sapply(tissue_means_list, is.null)
tissue_means <- do.call(cbind, tissue_means_list[keep])
colnames(tissue_means) <- paste0(tools::toTitleCase(tissue_cond$tissue[keep]),
                                  " (", tissue_cond$condition[keep], ")")

expr_df <- as.data.frame(tissue_means[found_ids, , drop = FALSE])
expr_df$Gene_ID <- rownames(expr_df)
expr_df <- merge(meth_genes, expr_df, by = "Gene_ID")

save_data(expr_df, BATCH_DIR, "methylation_machinery_expression")

# --- Fig 2A: Presence/absence (including absent genes) ---
presence_df <- meth_genes
presence_df$Gene <- factor(presence_df$Gene, levels = rev(meth_genes$Gene))
presence_df$Status <- ifelse(presence_df$Present, "Present", "Absent")

p2a <- ggplot(presence_df, aes(x = 1, y = Gene, fill = Status)) +
  geom_tile(color = "white", linewidth = 1) +
  scale_fill_manual(values = c(Present = "#27AE60", Absent = "#C0392B"),
                    name = NULL) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Methylation toolkit in D. laeve genome",
       subtitle = "Notable: DNMT3 (de novo) absent, DNMT1 (maintenance) present",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), strip.text.y = element_text(angle = 0),
        legend.position = "bottom")
save_fig(p2a, BATCH_DIR, "fig2a_methylation_toolkit_presence", w = 6, h = 10)

# --- Fig 2C: Heatmap (tissue × condition) ---
hm_mat <- as.matrix(expr_df[, colnames(tissue_means)])
rownames(hm_mat) <- expr_df$Gene
hm_z <- t(scale(t(hm_mat)))

ann_row <- data.frame(Category = expr_df$Category, row.names = expr_df$Gene)
# Column annotation: tissue and condition parsed from "Tissue (condition)" names
col_tissue <- sub(" \\(.*", "", colnames(tissue_means))
col_condition <- sub(".*\\((.*)\\)", "\\1", colnames(tissue_means))
ann_col <- data.frame(Tissue = col_tissue, Condition = col_condition,
                       row.names = colnames(tissue_means))
ann_colors <- list(
  Condition = c(control = "#2471A3", amputated = "#C0392B",
                irradiated = "#F39C12", fungicide = "#27AE60"),
  Category = setNames(RColorBrewer::brewer.pal(
    length(unique(expr_df$Category)), "Set3"),
    unique(expr_df$Category))
)
p2c <- pheatmap(hm_z, annotation_row = ann_row, annotation_col = ann_col,
                annotation_colors = ann_colors,
                color = colorRampPalette(c("#2471A3", "white", "#C0392B"))(100),
                main = "Methylation toolkit expression (z-scored, all conditions)",
                clustering_method = "ward.D2", fontsize = 10, silent = TRUE)

png(file.path(BATCH_DIR, "figures", "fig2c_methylation_toolkit_heatmap.png"),
    width = 12, height = 10, units = "in", res = 300)
grid::grid.draw(p2c$gtable); dev.off()
tryCatch({
  cairo_pdf(file.path(BATCH_DIR, "figures", "fig2c_methylation_toolkit_heatmap.pdf"), width = 12, height = 10)
  grid::grid.draw(p2c$gtable); dev.off()
}, error = function(e) {
  pdf(file.path(BATCH_DIR, "figures", "fig2c_methylation_toolkit_heatmap.pdf"), width = 12, height = 10)
  grid::grid.draw(p2c$gtable); dev.off()
})
cat("  Saved: fig2c_methylation_toolkit_heatmap\n")

cat("\n=== Batch 02 complete ===\n")
