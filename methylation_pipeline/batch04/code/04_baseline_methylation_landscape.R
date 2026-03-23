#!/usr/bin/env Rscript
# =============================================================================
# Batch 04: Baseline Methylation Landscape + Gene Body Meth vs Expression
# Question: Where is methylation? Does gene body meth correlate with expression?
# Output: data/ + figures/ (4 plots)
# Requires: BSseq cache, DESeq2 counts
# =============================================================================

source("methylation_pipeline/_config.R")

library(bsseq)
library(GenomicRanges)
library(rtracklayer)
library(DESeq2)
library(data.table)
library(ggplot2)
library(dplyr)
library(scales)

BATCH_DIR <- file.path(PIPE_DIR, "batch04")

cat("=== Batch 04: Baseline Methylation Landscape ===\n\n")

# --- Load BSseq + GFF ---
bs_obj <- readRDS(CACHE$bsseq)
cat(sprintf("BSseq: %s sites x %d samples\n", format(nrow(bs_obj), big.mark = ","), ncol(bs_obj)))

gff <- load_gff()
genes <- gff[gff$type == "gene"]
exons <- gff[gff$type == "exon"]
te_data <- load_te()

# Promoters (2kb upstream)
promoters <- if (file.exists(CACHE$promoters)) readRDS(CACHE$promoters) else {
  trim(promoters(genes, upstream = 2000, downstream = 0))
}

# --- Per-sample global stats ---
cat("\nPer-sample methylation:\n")
sample_stats <- list()
for (i in 1:ncol(bs_obj)) {
  bv <- getMeth(bs_obj[, i], type = "raw")[, 1]
  bv <- bv[!is.nan(bv) & !is.na(bv)]
  sample_stats[[i]] <- data.frame(
    sample = colnames(bs_obj)[i], total_sites = length(bv),
    mean_beta = round(mean(bv), 4), median_beta = round(median(bv), 4),
    high = sum(bv > 0.8), low = sum(bv < 0.2), intermediate = sum(bv >= 0.2 & bv <= 0.8))
  rm(bv); gc(verbose = FALSE)
}
sample_stats <- do.call(rbind, sample_stats)
print(sample_stats)
save_data(sample_stats, BATCH_DIR, "per_sample_methylation_stats")

# --- Region-wise methylation (control samples, mean of C1+C2) ---
cat("\nComputing region methylation (controls)...\n")
ctrl_idx <- which(colnames(bs_obj) %in% c("C1", "C2"))
ctrl_beta <- rowMeans(getMeth(bs_obj[, ctrl_idx], type = "raw"), na.rm = TRUE)

bs_chr <- as.character(seqnames(bs_obj))
bs_pos <- start(bs_obj)

region_meth <- function(region_gr, label) {
  gr <- GRanges(seqnames = bs_chr, ranges = IRanges(start = bs_pos, width = 1))
  in_reg <- overlapsAny(gr, region_gr)
  bv <- ctrl_beta[in_reg]
  bv <- bv[!is.na(bv) & !is.nan(bv)]
  data.frame(region = label, n_cpgs = length(bv), mean_beta = round(mean(bv), 4),
             median_beta = round(median(bv), 4))
}

reg_stats <- rbind(
  region_meth(promoters, "Promoter"),
  region_meth(exons, "Exon"),
  region_meth(GenomicRanges::setdiff(reduce(genes), reduce(exons)), "Intron"),
  region_meth(reduce(genes), "Gene body")
)
# Intergenic = everything not in genes or promoters
genic_prom <- reduce(c(reduce(genes), promoters))
strand(genic_prom) <- "*"; genic_prom <- reduce(genic_prom)
cat("Region methylation:\n"); print(reg_stats)
save_data(reg_stats, BATCH_DIR, "region_methylation_stats")

# --- Fig 4D: Global methylation density per sample ---
set.seed(42)
sub_idx <- sort(sample(nrow(bs_obj), min(500000, nrow(bs_obj))))
beta_long <- data.frame()
for (i in 1:ncol(bs_obj)) {
  bv <- getMeth(bs_obj[sub_idx, i], type = "raw")[, 1]
  bv <- bv[!is.na(bv) & !is.nan(bv)]
  beta_long <- rbind(beta_long, data.frame(sample = colnames(bs_obj)[i], beta = bv))
}

p4d <- ggplot(beta_long, aes(x = beta, color = sample)) +
  geom_density(linewidth = 0.8) +
  labs(x = "Beta value", y = "Density", title = "CpG methylation distribution per sample",
       subtitle = sprintf("%s CpGs sampled", format(length(sub_idx), big.mark = ","))) +
  theme_minimal(base_size = 12)
save_fig(p4d, BATCH_DIR, "fig4d_methylation_density_per_sample", w = 9, h = 6)
rm(beta_long); gc(verbose = FALSE)

# --- Gene body methylation vs expression ---
cat("\nLoading expression data (DESeq2)...\n")

# Build DESeq2 from counts (tail controls only for baseline)
sample_table <- data.frame(
  sampleName = list.files(OG$counts_dir, pattern = "tail.*control|control.*tail", ignore.case = TRUE),
  stringsAsFactors = FALSE
)
# If no tail-control specific files, use all tail files
if (nrow(sample_table) == 0) {
  sample_table <- data.frame(
    sampleName = list.files(OG$counts_dir, pattern = "tail", ignore.case = TRUE),
    stringsAsFactors = FALSE
  )
}
sample_table$condition <- "control"
sample_table$fileName <- sample_table$sampleName

if (nrow(sample_table) > 0) {
  dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table,
                                     directory = OG$counts_dir,
                                     design = ~ 1)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  expr_vec <- rowMeans(assay(vsd))

  # Gene body methylation (mean beta in gene body, controls)
  gene_meth <- data.frame(gene_id = genes$ID, stringsAsFactors = FALSE)
  gr_all <- GRanges(seqnames = bs_chr, ranges = IRanges(start = bs_pos, width = 1))

  gene_meth$mean_beta <- sapply(seq_along(genes), function(i) {
    ov <- overlapsAny(gr_all, genes[i])
    bv <- ctrl_beta[ov]
    bv <- bv[!is.na(bv) & !is.nan(bv)]
    if (length(bv) == 0) return(NA)
    mean(bv)
  })

  gene_meth$expression <- expr_vec[match(gene_meth$gene_id, names(expr_vec))]
  gene_meth <- gene_meth[!is.na(gene_meth$mean_beta) & !is.na(gene_meth$expression), ]

  cor_test <- cor.test(gene_meth$mean_beta, gene_meth$expression, method = "spearman")
  cat(sprintf("\nGene body meth vs expression: rho = %.3f, p = %s, n = %d\n",
              cor_test$estimate, format(cor_test$p.value, digits = 3), nrow(gene_meth)))

  save_data(gene_meth, BATCH_DIR, "gene_body_meth_vs_expression")

  # Fig 4B: Scatter
  p4b <- ggplot(gene_meth, aes(x = mean_beta, y = expression)) +
    geom_point(alpha = 0.1, size = 0.5, color = "#2471A3") +
    geom_smooth(method = "loess", color = "#C0392B", se = FALSE) +
    labs(x = "Mean gene body methylation (beta)", y = "Expression (VST)",
         title = "Gene body methylation vs expression",
         subtitle = sprintf("Spearman rho = %.3f, p %s, n = %s genes",
                            cor_test$estimate,
                            ifelse(cor_test$p.value < 2.2e-16, "< 2.2e-16",
                                   format(cor_test$p.value, digits = 3)),
                            format(nrow(gene_meth), big.mark = ","))) +
    theme_minimal(base_size = 12)
  save_fig(p4b, BATCH_DIR, "fig4b_genebody_meth_vs_expression", w = 8, h = 7)

  # Fig 4C: Expression by methylation quartile
  gene_meth$meth_quartile <- cut(gene_meth$mean_beta,
                                  breaks = quantile(gene_meth$mean_beta, probs = c(0, 0.25, 0.5, 0.75, 1)),
                                  labels = c("Q1 (low)", "Q2", "Q3", "Q4 (high)"),
                                  include.lowest = TRUE)

  p4c <- ggplot(gene_meth, aes(x = meth_quartile, y = expression, fill = meth_quartile)) +
    geom_boxplot(outlier.size = 0.3) +
    scale_fill_manual(values = c("#2471A3", "#1ABC9C", "#F39C12", "#C0392B"), guide = "none") +
    labs(x = "Methylation quartile", y = "Expression (VST)",
         title = "Expression by gene body methylation level") +
    theme_minimal(base_size = 12)
  save_fig(p4c, BATCH_DIR, "fig4c_expression_by_meth_quartile", w = 7, h = 6)
} else {
  cat("WARNING: No tail control samples found for expression analysis.\n")
}

# --- Fig 4A: Metagene methylation profile ---
# TODO: This requires computing methylation in windows around TSS
# (upstream -> gene body -> downstream) which is memory-intensive.
# Placeholder for now — will be implemented with per-chromosome approach.
cat("\nFig 4A (metagene profile): requires per-gene window computation. Skipped for now.\n")

cat("\n=== Batch 04 complete ===\n")
