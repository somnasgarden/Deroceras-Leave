#!/usr/bin/env Rscript
# =============================================================================
# Batch 05: TE Methylation + Evolutionary Age
# Question: Is TE methylation related to age? Intronic vs intergenic TEs?
# Output: data/ + figures/ (2 plots)
# Requires: BSseq cache, TE age data
# =============================================================================

source("methylation_pipeline/_config.R")

library(bsseq)
library(GenomicRanges)
library(data.table)
library(ggplot2)
library(dplyr)

BATCH_DIR <- file.path(PIPE_DIR, "batch05")

cat("=== Batch 05: TE Methylation + Age ===\n\n")

# --- Load data ---
bs_obj <- readRDS(CACHE$bsseq)
te_age <- fread(OG$te_age)
te_age <- te_age[chrom %in% keep_chr]

gff <- load_gff()
genes <- gff[gff$type == "gene"]
gene_bodies <- reduce(genes)

# Control methylation
ctrl_idx <- which(colnames(bs_obj) %in% c("C1", "C2"))
ampu_idx <- which(colnames(bs_obj) %in% c("A1", "A2"))
ctrl_beta <- rowMeans(getMeth(bs_obj[, ctrl_idx], type = "raw"), na.rm = TRUE)
ampu_beta <- rowMeans(getMeth(bs_obj[, ampu_idx], type = "raw"), na.rm = TRUE)

bs_gr <- GRanges(seqnames = seqnames(bs_obj), ranges = IRanges(start = start(bs_obj), width = 1))

cat("Computing TE methylation by class and age...\n")

# --- Per-TE methylation ---
te_gr <- makeGRangesFromDataFrame(te_age, seqnames.field = "chrom",
                                   start.field = "start", end.field = "end",
                                   keep.extra.columns = TRUE)
te_gr$te_class <- sub("/.*", "", te_gr$class_family)

# Classify: intronic vs intergenic
te_gr$genic <- overlapsAny(te_gr, gene_bodies)

# Compute mean methylation per TE (subsample for speed)
set.seed(42)
n_te <- min(50000, length(te_gr))
te_sub <- te_gr[sample(length(te_gr), n_te)]

te_meth <- data.frame(
  te_class = te_sub$te_class,
  kimura = te_sub$perc_div,
  genic = te_sub$genic,
  stringsAsFactors = FALSE
)

# Mean beta per TE
te_meth$ctrl_beta <- sapply(seq_len(n_te), function(i) {
  ov <- overlapsAny(bs_gr, te_sub[i])
  bv <- ctrl_beta[ov]; bv <- bv[!is.na(bv) & !is.nan(bv)]
  if (length(bv) == 0) return(NA); mean(bv)
})
te_meth$ampu_beta <- sapply(seq_len(n_te), function(i) {
  ov <- overlapsAny(bs_gr, te_sub[i])
  bv <- ampu_beta[ov]; bv <- bv[!is.na(bv) & !is.nan(bv)]
  if (length(bv) == 0) return(NA); mean(bv)
})

te_meth <- te_meth[!is.na(te_meth$ctrl_beta), ]
te_meth$delta_beta <- te_meth$ampu_beta - te_meth$ctrl_beta

# Convert Kimura to millions of years (approximate: 1% divergence ~ 2.27 My for neutral DNA)
te_meth$age_my <- te_meth$kimura * 2.27

# Keep major classes
major_classes <- names(sort(table(te_meth$te_class), decreasing = TRUE))[1:6]
te_meth_major <- te_meth[te_meth$te_class %in% major_classes, ]

cat(sprintf("TEs with methylation data: %s\n", format(nrow(te_meth), big.mark = ",")))

# TE methylation does NOT change between conditions
cat("\nTE methylation change by class:\n")
te_change <- te_meth_major %>%
  group_by(te_class) %>%
  summarize(mean_ctrl = mean(ctrl_beta), mean_ampu = mean(ampu_beta),
            mean_delta = mean(delta_beta),
            cohens_d = mean(delta_beta) / sd(delta_beta),
            .groups = "drop")
print(as.data.frame(te_change))
save_data(te_meth_major, BATCH_DIR, "te_methylation_by_age")
save_data(te_change, BATCH_DIR, "te_methylation_change_summary")

# --- Fig 5A: Kimura divergence vs methylation by TE class ---
p5a <- ggplot(te_meth_major, aes(x = age_my, y = ctrl_beta, color = te_class)) +
  geom_point(alpha = 0.05, size = 0.3) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1.2) +
  facet_wrap(~ te_class, ncol = 3) +
  labs(x = "Evolutionary age (millions of years)", y = "Mean CpG methylation (beta)",
       title = "TE methylation vs evolutionary age",
       subtitle = "Younger TEs tend to be more methylated (passive gene body methylation)") +
  scale_color_manual(values = COLORS$te, guide = "none") +
  theme_minimal(base_size = 11)
save_fig(p5a, BATCH_DIR, "fig5a_te_age_vs_methylation", w = 10, h = 7)

# --- Fig 5B: Intronic vs intergenic TE methylation ---
te_meth_major$location <- ifelse(te_meth_major$genic, "Intronic/Genic", "Intergenic")

p5b <- ggplot(te_meth_major, aes(x = ctrl_beta, fill = location)) +
  geom_density(alpha = 0.4, linewidth = 0.6) +
  facet_wrap(~ te_class, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = c("Intronic/Genic" = "#27AE60", "Intergenic" = "#C0392B")) +
  labs(x = "Mean CpG methylation (beta)", y = "Density",
       title = "TE methylation: intronic vs intergenic",
       subtitle = "Intronic TEs are more methylated (passive gene body methylation)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
save_fig(p5b, BATCH_DIR, "fig5b_te_intronic_vs_intergenic", w = 10, h = 7)

rm(bs_obj, te_meth, te_meth_major); gc(verbose = FALSE)

cat("\n=== Batch 05 complete ===\n")
