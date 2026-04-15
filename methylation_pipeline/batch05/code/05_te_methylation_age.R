#!/usr/bin/env Rscript
# =============================================================================
# Batch 05: TE Methylation + Evolutionary Age
# Question: Is TE methylation related to age? Intronic vs intergenic? DMP/DMR overlap?
# Output: data/ (6 tables) + figures/ (10+ plots)
# Requires: BSseq cache, TE age data, Batch 06 DMPs/DMRs
# =============================================================================

source("methylation_pipeline/_config.R")
t0 <- proc.time()

suppressPackageStartupMessages({
  library(bsseq)
  library(GenomicRanges)
  library(data.table)
  library(ggplot2)
  library(ggridges)
  library(dplyr)
  library(tidyr)
  library(scales)
})

BATCH_DIR <- file.path(PIPE_DIR, "batch05")

# Clean old output
unlink(list.files(file.path(BATCH_DIR, "figures"), full.names = TRUE))
unlink(list.files(file.path(BATCH_DIR, "data"), full.names = TRUE))
dir.create(file.path(BATCH_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BATCH_DIR, "data"), showWarnings = FALSE, recursive = TRUE)

cat("=== Batch 05: TE Methylation + Age ===\n\n")

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("[1/7] Loading data...\n")
bs_obj <- readRDS(CACHE$bsseq)
cat(sprintf("  BSseq: %s sites x %d samples\n", format(nrow(bs_obj), big.mark = ","), ncol(bs_obj)))

te_age <- fread(OG$te_age)
te_age <- te_age[chrom %in% keep_chr]
cat(sprintf("  TEs: %s\n", format(nrow(te_age), big.mark = ",")))

gff <- load_gff()
genes <- gff[gff$type == "gene"]
gene_bodies <- GenomicRanges::reduce(genes)

# Beta values
ctrl_idx <- which(colnames(bs_obj) %in% c("C1", "C2"))
ampu_idx <- which(colnames(bs_obj) %in% c("A1", "A2"))
ctrl_beta <- rowMeans(getMeth(bs_obj[, ctrl_idx], type = "raw"), na.rm = TRUE)
ampu_beta <- rowMeans(getMeth(bs_obj[, ampu_idx], type = "raw"), na.rm = TRUE)

bs_gr <- GRanges(seqnames = seqnames(bs_obj), ranges = IRanges(start = start(bs_obj), width = 1))

# =============================================================================
# 2. PER-TE METHYLATION
# =============================================================================
cat("[2/7] Computing per-TE methylation...\n")

te_gr <- makeGRangesFromDataFrame(te_age, seqnames.field = "chrom",
                                   start.field = "start", end.field = "end",
                                   keep.extra.columns = TRUE)
te_gr$te_class <- sub("/.*", "", te_gr$class_family)

# Genic vs intergenic
te_gr$genic <- overlapsAny(te_gr, gene_bodies)

# Vectorized methylation per TE
hits <- findOverlaps(bs_gr, te_gr)
hit_dt <- data.table(te_idx = subjectHits(hits),
                     ctrl_beta = ctrl_beta[queryHits(hits)],
                     ampu_beta = ampu_beta[queryHits(hits)])
hit_dt <- hit_dt[!is.na(ctrl_beta) & !is.nan(ctrl_beta)]

te_means <- hit_dt[, .(ctrl_beta = mean(ctrl_beta), ampu_beta = mean(ampu_beta),
                        n_cpgs = .N), by = te_idx]

te_meth <- data.frame(
  te_idx = seq_along(te_gr),
  te_class = te_gr$te_class,
  kimura = te_gr$kimura_div,
  genic = te_gr$genic,
  stringsAsFactors = FALSE
)
te_meth$ctrl_beta <- NA_real_; te_meth$ampu_beta <- NA_real_; te_meth$n_cpgs <- 0L
te_meth$ctrl_beta[te_means$te_idx] <- te_means$ctrl_beta
te_meth$ampu_beta[te_means$te_idx] <- te_means$ampu_beta
te_meth$n_cpgs[te_means$te_idx] <- te_means$n_cpgs
te_meth <- te_meth[!is.na(te_meth$ctrl_beta), ]
te_meth$delta_beta <- te_meth$ampu_beta - te_meth$ctrl_beta

# Convert Kimura to millions of years (1% divergence ~ 2.27 My for neutral DNA)
te_meth$age_my <- te_meth$kimura * 2.27

# Keep major classes
class_counts <- sort(table(te_meth$te_class), decreasing = TRUE)
major_classes <- names(class_counts)[1:min(6, length(class_counts))]
te_major <- te_meth[te_meth$te_class %in% major_classes, ]
te_major$location <- ifelse(te_major$genic, "Overlapping genes", "Intergenic")

cat(sprintf("  TEs with methylation: %s (%d major classes)\n",
            format(nrow(te_meth), big.mark = ","), length(major_classes)))

# =============================================================================
# 3. TE METHYLATION STATISTICS
# =============================================================================
cat("[3/7] Computing TE methylation statistics...\n")

# By class
te_class_stats <- te_major %>%
  group_by(te_class) %>%
  summarize(n_tes = n(), total_cpgs = sum(n_cpgs),
            mean_ctrl = mean(ctrl_beta), mean_ampu = mean(ampu_beta),
            mean_delta = mean(delta_beta),
            cohens_d = mean(delta_beta) / sd(delta_beta),
            sd_ctrl = sd(ctrl_beta),
            .groups = "drop")
cat("\n  TE methylation by class:\n"); print(as.data.frame(te_class_stats))
save_data(te_class_stats, BATCH_DIR, "te_methylation_by_class")
save_data(te_major, BATCH_DIR, "te_methylation_by_age")

# =============================================================================
# 4. DMP/DMR OVERLAP WITH TEs
# =============================================================================
cat("[4/7] DMP/DMR overlap with TEs...\n")

dmp_file <- file.path(PIPE_DIR, "batch06/data/dmps_annotated.tsv")
dmr_file <- file.path(PIPE_DIR, "batch06/data/dmrs_annotated.tsv")

has_dmp <- file.exists(dmp_file)
has_dmr <- file.exists(dmr_file)

if (has_dmp) {
  dmp <- fread(dmp_file)
  dmp_gr <- GRanges(seqnames = dmp$chr, ranges = IRanges(start = dmp$pos, width = 1))
  dmp_gr$diff <- dmp$diff
  dmp_gr$direction <- ifelse(dmp$diff > 0, "Hyper", "Hypo")

  dmp_te_hits <- findOverlaps(dmp_gr, te_gr)
  n_dmp_in_te <- length(unique(queryHits(dmp_te_hits)))
  cat(sprintf("  DMPs in TEs: %d / %d (%.1f%%)\n",
              n_dmp_in_te, length(dmp_gr), 100 * n_dmp_in_te / length(dmp_gr)))

  # DMP-TE mapping
  dmp_te_df <- data.frame(
    dmp_idx = queryHits(dmp_te_hits),
    te_idx = subjectHits(dmp_te_hits),
    diff = dmp_gr$diff[queryHits(dmp_te_hits)],
    direction = dmp_gr$direction[queryHits(dmp_te_hits)],
    te_class = te_gr$te_class[subjectHits(dmp_te_hits)],
    stringsAsFactors = FALSE
  )

  # DMP counts by TE class
  dmp_by_class <- dmp_te_df %>%
    group_by(te_class, direction) %>%
    summarize(n_dmps = n(), .groups = "drop")
  save_data(dmp_by_class, BATCH_DIR, "dmp_count_by_te_class")

  # Fisher's exact test: are DMPs over/under-represented per TE class
  # relative to background CpG count in each class?
  cat("  Fisher's exact test — DMP enrichment per TE class:\n")
  bg_cpgs_in_te <- countOverlaps(te_gr, bs_gr)
  bg_by_class <- data.frame(te_class = te_gr$te_class, n_cpg = bg_cpgs_in_te) %>%
    group_by(te_class) %>% summarize(bg_cpgs = sum(n_cpg), .groups = "drop")
  dmp_total_by_class <- dmp_te_df %>%
    group_by(te_class) %>% summarize(dmp_cpgs = n_distinct(dmp_idx), .groups = "drop")
  fisher_df <- merge(dmp_total_by_class, bg_by_class, by = "te_class")
  total_dmp <- length(dmp_gr); total_bg <- length(bs_gr)
  fisher_df$fisher_p <- sapply(seq_len(nrow(fisher_df)), function(i) {
    fisher.test(matrix(c(fisher_df$dmp_cpgs[i], total_dmp - fisher_df$dmp_cpgs[i],
                          fisher_df$bg_cpgs[i], total_bg - fisher_df$bg_cpgs[i]), nrow = 2))$p.value
  })
  fisher_df$fold <- (fisher_df$dmp_cpgs / total_dmp) / (fisher_df$bg_cpgs / total_bg)
  fisher_df$sig <- ifelse(fisher_df$fisher_p < 0.001, "***",
                    ifelse(fisher_df$fisher_p < 0.01, "**",
                    ifelse(fisher_df$fisher_p < 0.05, "*", "ns")))
  cat("  TE class DMP enrichment:\n")
  print(as.data.frame(fisher_df[order(fisher_df$fisher_p), ]))
  save_data(fisher_df, BATCH_DIR, "dmp_te_class_fisher_enrichment")
}

if (has_dmr) {
  dmr <- fread(dmr_file)
  dmr_gr <- GRanges(seqnames = dmr$chr,
                     ranges = IRanges(start = dmr$start, end = dmr$end))
  dmr_gr$diff_methy <- dmr$diff.Methy
  dmr_gr$direction <- ifelse(dmr$diff.Methy > 0, "Hyper", "Hypo")

  dmr_te_hits <- findOverlaps(dmr_gr, te_gr)
  n_dmr_in_te <- length(unique(queryHits(dmr_te_hits)))
  cat(sprintf("  DMRs overlapping TEs: %d / %d (%.1f%%)\n",
              n_dmr_in_te, length(dmr_gr), 100 * n_dmr_in_te / length(dmr_gr)))

  dmr_te_df <- data.frame(
    dmr_idx = queryHits(dmr_te_hits),
    te_idx = subjectHits(dmr_te_hits),
    diff_methy = dmr_gr$diff_methy[queryHits(dmr_te_hits)],
    direction = dmr_gr$direction[queryHits(dmr_te_hits)],
    te_class = te_gr$te_class[subjectHits(dmr_te_hits)],
    stringsAsFactors = FALSE
  )

  dmr_by_class <- dmr_te_df %>%
    group_by(te_class, direction) %>%
    summarize(n_dmrs = n(), .groups = "drop")
  save_data(dmr_by_class, BATCH_DIR, "dmr_count_by_te_class")
}

# =============================================================================
# 5. FIGURES — RIDGE PLOTS
# =============================================================================
cat("[5/7] Generating ridge plots...\n")

# Fig 5A: TE methylation ridge by class — control and amputated SEPARATE panels
te_ridge <- te_major %>%
  select(te_class, ctrl_beta, ampu_beta) %>%
  pivot_longer(c(ctrl_beta, ampu_beta), names_to = "condition", values_to = "meth") %>%
  mutate(condition = ifelse(condition == "ctrl_beta", "Control", "Amputated"))

p5a <- ggplot(te_ridge, aes(x = meth * 100, y = te_class, fill = te_class)) +
  geom_density_ridges(alpha = 0.6, scale = 0.9) +
  scale_fill_manual(values = COLORS$te, guide = "none") +
  facet_wrap(~ condition) +
  labs(title = "TE methylation distribution by class",
       subtitle = "Bimodal distribution expected — panels separated by condition",
       x = "Mean CpG methylation (%)", y = "TE class") +
  theme_minimal(base_size = 12)
save_fig(p5a, BATCH_DIR, "fig5a_te_ridge_by_class", w = 14, h = 7)

# Fig 5B: Intronic vs intergenic ridge
p5b <- ggplot(te_major, aes(x = ctrl_beta * 100, y = te_class, fill = location)) +
  geom_density_ridges(alpha = 0.5, scale = 0.9) +
  scale_fill_manual(values = c("Overlapping genes" = "#27AE60", "Intergenic" = "#C0392B")) +
  labs(title = "TE methylation: overlapping genes vs intergenic",
       subtitle = "Gene-overlapping TEs more methylated (passive gene body methylation)",
       x = "Mean CpG methylation (%)", y = "TE class", fill = "Location") +
  theme_minimal(base_size = 12)
save_fig(p5b, BATCH_DIR, "fig5b_te_intronic_vs_intergenic_ridge", w = 10, h = 7)

# fig5c (dmp_diff_ridge_by_te_class) and fig5d (violin_dmp_diff_in_tes) removed
# at user request — uninformative when meth diffs cluster around the 10% DMP threshold.

# =============================================================================
# 6. FIGURES — VIOLIN, BAR, SCATTER
# =============================================================================
cat("[6/7] Generating bar and scatter plots...\n")

# Filter for age analysis
te_plot <- te_major[!is.na(te_major$kimura) & te_major$kimura < 60, ]
cat(sprintf("  TEs for age plot: %s\n", format(nrow(te_plot), big.mark = ",")))

# Fig 5G: Binned bar plot of methylation by TE age, split genic vs intergenic
# All TEs with reads (n_cpgs >= 1) — no CpG count filter to avoid length-age confound
cat("  Fig 5G: TE age bins (bar plot, genic vs intergenic)...\n")
te_plot$location <- ifelse(te_plot$genic, "Overlapping genes", "Intergenic")
age_breaks <- c(-Inf, 5, 15, 30, 45, Inf)
age_labels <- c("Very Young\n(0-5%)", "Young\n(5-15%)", "Intermediate\n(15-30%)",
                "Old\n(30-45%)", "Very Old\n(>45%)")
te_plot$age_bin <- cut(te_plot$kimura, breaks = age_breaks, labels = age_labels)
te_age_bar <- as.data.frame(
  data.table::as.data.table(te_plot)[, .(
    mean_meth = mean(ctrl_beta) * 100,
    se = sd(ctrl_beta) * 100 / sqrt(.N),
    n = .N
  ), by = .(te_class, age_bin, location)]
)
te_age_bar <- te_age_bar[!is.na(te_age_bar$age_bin), ]
save_data(te_age_bar, BATCH_DIR, "te_methylation_by_age_bin")
p5g2 <- ggplot(te_age_bar, aes(x = age_bin, y = mean_meth, fill = location)) +
  geom_col(position = "dodge", width = 0.75) +
  geom_errorbar(aes(ymin = mean_meth - se, ymax = mean_meth + se),
                position = position_dodge(0.75), width = 0.25, linewidth = 0.4) +
  facet_wrap(~ te_class, ncol = 3) +
  scale_fill_manual(values = c("Overlapping genes" = "#2471A3", "Intergenic" = "#E67E22"), name = "Location") +
  labs(title = "Mean CpG Methylation by TE Age, Class, and Genomic Context",
       subtitle = "All TEs with reads | Error bars = \u00b1 1 SE | Kimura divergence bins",
       x = "Kimura Divergence Bin", y = "Mean Methylation (%)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 8))
save_fig(p5g2, BATCH_DIR, "fig5g_te_age_bar", w = 12, h = 7)

# fig5h (te_age_vs_meth_diff) and fig5i (te_intronic_vs_intergenic_density)
# removed at user request — redundant with fig5g (age vs methylation) and
# fig5b (intronic vs intergenic ridge), respectively.

# =============================================================================
# 7. SUMMARY
# =============================================================================
cat("\n[7/7] Summary:\n")
cat(sprintf("  TEs with methylation data: %s\n", format(nrow(te_meth), big.mark = ",")))
cat(sprintf("  Mean TE methylation (ctrl): %.1f%%\n", mean(te_meth$ctrl_beta) * 100))
cat(sprintf("  Mean TE methylation (ampu): %.1f%%\n", mean(te_meth$ampu_beta) * 100))
if (has_dmp) cat(sprintf("  DMPs in TEs: %d (%.1f%%)\n",
                          n_dmp_in_te, 100 * n_dmp_in_te / length(dmp_gr)))
if (has_dmr) cat(sprintf("  DMRs in TEs: %d (%.1f%%)\n",
                          n_dmr_in_te, 100 * n_dmr_in_te / length(dmr_gr)))

rm(bs_obj); gc(verbose = FALSE)

elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("\n=== Batch 05 complete (%.1f min) ===\n", elapsed / 60))
cat(sprintf("Figures: 10 | Data files: 5\n"))
