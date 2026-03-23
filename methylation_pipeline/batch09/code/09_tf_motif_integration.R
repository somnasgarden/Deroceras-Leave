#!/usr/bin/env Rscript
# =============================================================================
# Batch 09: TF + GENIE3 + Motif Integration
# Question: Do TF binding sites near DMPs explain regulatory connections?
# Output: data/ (TF-DMP overlap, GENIE3 network) + figures/ (4 plots)
# Requires: Batch 06 DMPs, Batch 1.5 motif hits, GENIE3 data
# =============================================================================

source("methylation_pipeline/_config.R")

library(data.table)
library(ggplot2)
library(GenomicRanges)
library(dplyr)

BATCH_DIR <- file.path(PIPE_DIR, "batch09")
cat("=== Batch 09: TF + Motif + GENIE3 Integration ===\n\n")

# --- Load DMPs ---
dmp <- fread(file.path(PIPE_DIR, "batch06/data/dmps_annotated.tsv"))
cat(sprintf("DMPs: %s\n", format(nrow(dmp), big.mark = ",")))

# --- Load TF predictions ---
tf_pred <- fread(OG$deeptf)
# prediction column is logical TRUE/FALSE
tf_genes <- unique(sub("\\..*", "", tf_pred$sequence_ID[tf_pred$prediction == TRUE]))
cat(sprintf("Predicted TFs: %d\n", length(tf_genes)))

# --- TF enrichment for DMPs ---
gff <- load_gff(); genes <- gff[gff$type == "gene"]
all_genes <- genes$ID

tf_with_dmp <- sum(tf_genes %in% unique(dmp$nearest_gene))
tf_without_dmp <- length(tf_genes) - tf_with_dmp
nontf_with_dmp <- length(setdiff(unique(dmp$nearest_gene), tf_genes))
nontf_without_dmp <- length(all_genes) - length(tf_genes) - nontf_with_dmp

ft <- fisher.test(matrix(c(tf_with_dmp, tf_without_dmp, nontf_with_dmp, nontf_without_dmp), nrow = 2))
cat(sprintf("\nTF DMP enrichment: OR = %.2f, p = %s\n", ft$estimate, format(ft$p.value, digits = 3)))
cat(sprintf("  TFs with DMP: %d/%d (%.1f%%)\n", tf_with_dmp, length(tf_genes),
            100 * tf_with_dmp / length(tf_genes)))

save_data(data.frame(OR = ft$estimate, pvalue = ft$p.value, tf_with_dmp = tf_with_dmp,
                     tf_total = length(tf_genes)), BATCH_DIR, "tf_dmp_enrichment")

# --- Load motif hits (if available from batch 1.5) ---
motif_file <- file.path(PIPE_DIR, "../results/batch1.5/motif_hits_annotated.tsv.gz")
if (!file.exists(motif_file)) {
  motif_file <- file.path(PIPE_DIR, "../local/scripts/batch1.5/motif_hits_annotated.tsv.gz")
}

if (file.exists(motif_file)) {
  cat("Loading motif hits...\n")
  motif_hits <- fread(motif_file)
  cat(sprintf("Motif hits: %s\n", format(nrow(motif_hits), big.mark = ",")))

  # DMP-motif overlap: which DMPs fall within a TFBS?
  dmp_gr <- GRanges(seqnames = dmp$chr, ranges = IRanges(start = dmp$pos, width = 1))
  motif_gr <- GRanges(seqnames = motif_hits$chr,
                      ranges = IRanges(start = motif_hits$start, end = motif_hits$end))
  dmp$in_tfbs <- overlapsAny(dmp_gr, motif_gr)
  cat(sprintf("DMPs in TFBS: %d / %d (%.1f%%)\n",
              sum(dmp$in_tfbs), nrow(dmp), 100 * mean(dmp$in_tfbs)))
} else {
  cat("No motif hits file found. Run batch 1.5 section G or cluster job.\n")
  dmp$in_tfbs <- NA
}

# --- GENIE3 integration ---
if (file.exists(OG$genie3)) {
  cat("\nLoading GENIE3 network...\n")
  genie3 <- fread(OG$genie3, nrows = 500000)
  colnames(genie3) <- c("regulator", "target", "weight")
  cat(sprintf("GENIE3 edges: %s\n", format(nrow(genie3), big.mark = ",")))

  # TFs with DMPs that are also GENIE3 regulators
  tf_dmp_genes <- intersect(tf_genes, unique(dmp$nearest_gene))
  tf_in_network <- intersect(tf_dmp_genes, genie3$regulator)
  cat(sprintf("TFs with DMPs in GENIE3 network: %d\n", length(tf_in_network)))

  save_data(data.frame(tf_gene = tf_in_network), BATCH_DIR, "tf_dmp_in_genie3")
}

# --- Fig 9A: TF family census (from batch1.5 data) ---
# This reuses batch1.5 plots. Copy or reference them.
cat("\nFigures for batch09 depend on batch1.5 motif data.\n")
cat("Plots A-C come from batch1.5. Plot D (Sox19a) requires motif positions.\n")

save_data(dmp[, c("chr", "pos", "diff", "direction", "annotation", "nearest_gene",
                   "dist_to_gene", "in_te", "in_tfbs")],
          BATCH_DIR, "dmps_with_tfbs_annotation")

cat("\n=== Batch 09 complete ===\n")
