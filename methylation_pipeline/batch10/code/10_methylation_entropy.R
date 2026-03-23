#!/usr/bin/env Rscript
# =============================================================================
# Batch 10: Methylation Entropy Analysis
# Question: Does regeneration increase or decrease methylation entropy?
# Output: data/ (entropy stats) + figures/ (3 plots)
# Requires: BSseq cache
# =============================================================================

source("methylation_pipeline/_config.R")

library(bsseq)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(dplyr)

BATCH_DIR <- file.path(PIPE_DIR, "batch10")
cat("=== Batch 10: Methylation Entropy ===\n\n")

# --- Load BSseq ---
bs_obj <- readRDS(CACHE$bsseq)
cat(sprintf("BSseq: %s sites\n", format(nrow(bs_obj), big.mark = ",")))

# --- Per-CpG Shannon entropy ---
# For each CpG, compute entropy across samples
# H = -sum(p * log2(p)) where p = proportion of reads methylated
# Binary entropy: H(p) = -p*log2(p) - (1-p)*log2(1-p)
cat("Computing per-CpG entropy...\n")

binary_entropy <- function(p) {
  p[p == 0] <- 1e-10; p[p == 1] <- 1 - 1e-10
  -p * log2(p) - (1 - p) * log2(1 - p)
}

ctrl_beta <- rowMeans(getMeth(bs_obj[, c("C1", "C2")], type = "raw"), na.rm = TRUE)
ampu_beta <- rowMeans(getMeth(bs_obj[, c("A1", "A2")], type = "raw"), na.rm = TRUE)

ctrl_entropy <- binary_entropy(ctrl_beta)
ampu_entropy <- binary_entropy(ampu_beta)

valid <- !is.na(ctrl_entropy) & !is.na(ampu_entropy) & !is.nan(ctrl_entropy) & !is.nan(ampu_entropy)

entropy_df <- data.frame(
  chr = as.character(seqnames(bs_obj))[valid],
  pos = start(bs_obj)[valid],
  ctrl_beta = ctrl_beta[valid],
  ampu_beta = ampu_beta[valid],
  ctrl_entropy = ctrl_entropy[valid],
  ampu_entropy = ampu_entropy[valid],
  delta_entropy = ampu_entropy[valid] - ctrl_entropy[valid]
)

cat(sprintf("Valid CpGs for entropy: %s\n", format(nrow(entropy_df), big.mark = ",")))
cat(sprintf("Mean entropy - Control: %.4f, Amputated: %.4f\n",
            mean(entropy_df$ctrl_entropy), mean(entropy_df$ampu_entropy)))
cat(sprintf("Mean delta entropy: %.6f\n", mean(entropy_df$delta_entropy)))

# Wilcoxon test (subsample for speed)
set.seed(42)
sub_n <- min(100000, nrow(entropy_df))
sub_idx <- sample(nrow(entropy_df), sub_n)
wt <- wilcox.test(entropy_df$ctrl_entropy[sub_idx], entropy_df$ampu_entropy[sub_idx], paired = TRUE)
cat(sprintf("Wilcoxon paired test (n=%s): p = %s\n", format(sub_n, big.mark = ","),
            format(wt$p.value, digits = 3)))

save_data(data.frame(
  mean_ctrl_entropy = mean(entropy_df$ctrl_entropy),
  mean_ampu_entropy = mean(entropy_df$ampu_entropy),
  mean_delta = mean(entropy_df$delta_entropy),
  wilcox_p = wt$p.value,
  n_cpgs = nrow(entropy_df)
), BATCH_DIR, "entropy_summary")

# --- Annotate by region ---
gff <- load_gff(); genes <- gff[gff$type == "gene"]; exons <- gff[gff$type == "exon"]
promoters <- if(file.exists(CACHE$promoters)) readRDS(CACHE$promoters) else trim(promoters(genes, upstream=2000, downstream=0))
entropy_df$region <- annotate_regions(entropy_df$chr, entropy_df$pos, promoters, exons, genes)

region_entropy <- entropy_df %>%
  group_by(region) %>%
  summarize(mean_ctrl = mean(ctrl_entropy), mean_ampu = mean(ampu_entropy),
            mean_delta = mean(delta_entropy), n = n(), .groups = "drop")
cat("\nEntropy by region:\n"); print(as.data.frame(region_entropy))
save_data(region_entropy, BATCH_DIR, "entropy_by_region")

rm(bs_obj); gc(verbose = FALSE)

# --- Figures ---
# (A) Entropy distribution: control vs amputated
set.seed(42)
plot_idx <- sample(nrow(entropy_df), min(500000, nrow(entropy_df)))
ent_long <- rbind(
  data.frame(condition = "Control", entropy = entropy_df$ctrl_entropy[plot_idx]),
  data.frame(condition = "Amputated", entropy = entropy_df$ampu_entropy[plot_idx])
)

p10a <- ggplot(ent_long, aes(x = entropy, fill = condition, color = condition)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values = COLORS$condition) +
  scale_color_manual(values = COLORS$condition) +
  labs(x = "Binary Shannon entropy", y = "Density",
       title = "CpG methylation entropy: control vs amputated",
       subtitle = sprintf("Delta = %.5f | Wilcoxon p = %s",
                          mean(entropy_df$delta_entropy), format(wt$p.value, digits = 3))) +
  theme_minimal(base_size = 12)
save_fig(p10a, BATCH_DIR, "fig10a_entropy_distribution", w = 9, h = 6)

# (B) Per-region entropy change
p10b <- ggplot(region_entropy, aes(x = region, y = mean_delta, fill = region)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = COLORS$region, guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "Mean delta entropy (Amputated - Control)",
       title = "Entropy change by genomic region") +
  theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p10b, BATCH_DIR, "fig10b_entropy_by_region", w = 7, h = 6)

# (C) Entropy vs expression (if batch07 data exists)
cat("\nFig 10C (entropy vs expression) depends on batch07 data.\n")

cat("\n=== Batch 10 complete ===\n")
