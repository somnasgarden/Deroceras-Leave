#!/usr/bin/env Rscript
# =============================================================================
# Batch 03: Promoter CpG Classification (Weber criteria)
# Question: What does the D. laeve promoter CpG landscape look like?
#           How does it compare to mammals?
# Output: data/ (promoter classification) + figures/ (4 plots)
# Requires: Batch 01 (genome cache)
# =============================================================================

source("methylation_pipeline/_config.R")

library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(data.table)

BATCH_DIR <- file.path(PIPE_DIR, "batch03")

cat("=== Batch 03: Promoter CpG Classification ===\n\n")

# --- Load genome + GFF ---
genome <- load_genome()
gff <- load_gff()
genes <- gff[gff$type == "gene"]
cat(sprintf("Genes: %d\n", length(genes)))

# --- Define promoter windows: TSS +/- 1 kb (Weber et al. 2007) ---
# +strand: TSS = start(gene), window = [TSS-1000, TSS+1000]
# -strand: TSS = end(gene),   window = [TSS-1000, TSS+1000]
gene_strand <- as.character(strand(genes))
gene_chr <- as.character(seqnames(genes))

tss <- ifelse(gene_strand == "+", start(genes), end(genes))
prom_start <- pmax(1L, tss - 1000L)
prom_end <- tss + 1000L

# Get chromosome lengths to clip
chr_lengths <- width(genome)
names(chr_lengths) <- names(genome)
prom_end <- pmin(prom_end, chr_lengths[gene_chr])

valid <- prom_start < prom_end & gene_chr %in% keep_chr
prom_gr <- GRanges(seqnames = gene_chr[valid],
                   ranges = IRanges(start = prom_start[valid], end = prom_end[valid]),
                   gene_id = genes$ID[valid],
                   gene_strand = gene_strand[valid])
cat(sprintf("Promoter windows (TSS +/- 1kb): %d\n\n", length(prom_gr)))

# --- Compute CpG O/E and GC% for each promoter ---
cat("Computing per-promoter CpG stats...\n")

prom_stats <- do.call(rbind, lapply(unique(as.character(seqnames(prom_gr))), function(chr_name) {
  chr_prom <- prom_gr[seqnames(prom_gr) == chr_name]
  if (length(chr_prom) == 0) return(NULL)

  v <- Views(genome[[chr_name]], start = start(chr_prom), end = end(chr_prom))
  seqs <- as(v, "DNAStringSet")

  # Per-promoter stats
  n_seq <- length(seqs)
  result <- data.frame(
    gene_id = mcols(chr_prom)$gene_id,
    chr = chr_name,
    tss = ifelse(mcols(chr_prom)$gene_strand == "+", start(chr_prom) + 1000L, end(chr_prom) - 1000L),
    width = width(seqs),
    stringsAsFactors = FALSE
  )

  bf <- letterFrequency(seqs, letters = c("A", "C", "G", "T"))
  dn <- dinucleotideFrequency(seqs)

  result$n_C <- bf[, "C"]
  result$n_G <- bf[, "G"]
  result$n_CpG <- dn[, "CG"]
  result$gc_pct <- 100 * (bf[, "C"] + bf[, "G"]) / rowSums(bf)
  result$cpg_oe <- ifelse(result$n_C > 0 & result$n_G > 0,
                          (result$n_CpG * result$width) / (result$n_C * result$n_G), NA)

  rm(seqs, v); gc(verbose = FALSE)
  result
}))

cat(sprintf("Promoters computed: %d\n", nrow(prom_stats)))

# --- Weber classification ---
# HCP: CpG O/E >= 0.75 AND GC% >= 55%
# LCP: CpG O/E < 0.48
# ICP: everything else
prom_stats$weber_class <- "ICP"
prom_stats$weber_class[prom_stats$cpg_oe >= 0.75 & prom_stats$gc_pct >= 55] <- "HCP"
prom_stats$weber_class[prom_stats$cpg_oe < 0.48] <- "LCP"
prom_stats$weber_class <- factor(prom_stats$weber_class, levels = c("HCP", "ICP", "LCP"))

class_counts <- table(prom_stats$weber_class)
class_pct <- round(100 * class_counts / sum(class_counts), 1)

cat("\nWeber promoter classification:\n")
for (cl in names(class_counts)) {
  cat(sprintf("  %s: %s (%.1f%%)\n", cl, format(class_counts[cl], big.mark = ","), class_pct[cl]))
}

save_data(prom_stats, BATCH_DIR, "promoter_cpg_classification")
save_data(as.data.frame(class_counts), BATCH_DIR, "promoter_class_summary")

# --- Plots ---
weber_colors <- c(HCP = "#C0392B", ICP = "#F39C12", LCP = "#2471A3")

# Fig 3A: CpG O/E histogram of all promoters
p3a <- ggplot(prom_stats, aes(x = cpg_oe)) +
  geom_histogram(bins = 80, fill = "#2471A3", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "#C0392B", linewidth = 0.8) +
  geom_vline(xintercept = 0.48, linetype = "dashed", color = "#2471A3", linewidth = 0.8) +
  annotate("text", x = 0.78, y = Inf, label = "HCP threshold (0.75)", color = "#C0392B",
           vjust = 2, hjust = 0, size = 3) +
  annotate("text", x = 0.35, y = Inf, label = "LCP threshold (0.48)", color = "#2471A3",
           vjust = 2, hjust = 0, size = 3) +
  labs(x = "CpG O/E ratio", y = "Number of promoters",
       title = "Promoter CpG O/E distribution",
       subtitle = sprintf("D. laeve | %s promoters (TSS +/- 1kb) | Weber thresholds marked",
                          format(nrow(prom_stats), big.mark = ","))) +
  theme_minimal(base_size = 12)
save_fig(p3a, BATCH_DIR, "fig3a_promoter_cpg_oe_histogram", w = 9, h = 6)

# Fig 3B: CpG O/E vs GC% scatter (Weber classification)
p3b <- ggplot(prom_stats, aes(x = gc_pct, y = cpg_oe, color = weber_class)) +
  geom_point(alpha = 0.15, size = 0.5) +
  scale_color_manual(values = weber_colors, name = "Class") +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = 0.48, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 55, linetype = "dashed", color = "gray40") +
  labs(x = "GC content (%)", y = "CpG O/E ratio",
       title = "Promoter CpG landscape (Weber classification)",
       subtitle = sprintf("HCP: %s (%.1f%%) | ICP: %s (%.1f%%) | LCP: %s (%.1f%%)",
                          format(class_counts["HCP"], big.mark = ","), class_pct["HCP"],
                          format(class_counts["ICP"], big.mark = ","), class_pct["ICP"],
                          format(class_counts["LCP"], big.mark = ","), class_pct["LCP"])) +
  theme_minimal(base_size = 12)
save_fig(p3b, BATCH_DIR, "fig3b_promoter_oe_vs_gc_scatter", w = 9, h = 7)

# Fig 3C: Class proportions bar
class_df <- data.frame(Class = names(class_counts), Count = as.integer(class_counts),
                       Pct = as.numeric(class_pct))
p3c <- ggplot(class_df, aes(x = Class, y = Pct, fill = Class)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(Pct, "%\n(", format(Count, big.mark = ","), ")")),
            vjust = -0.3, size = 4) +
  scale_fill_manual(values = weber_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "% of promoters",
       title = "Promoter class proportions",
       subtitle = "D. laeve vs mammals: 91.6% ICP (mammals: ~16%)") +
  theme_minimal(base_size = 12)
save_fig(p3c, BATCH_DIR, "fig3c_promoter_class_proportions", w = 7, h = 6)

# Fig 3D: Faceted histograms by class
p3d <- ggplot(prom_stats, aes(x = cpg_oe, fill = weber_class)) +
  geom_histogram(bins = 50, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = weber_colors, guide = "none") +
  facet_wrap(~ weber_class, scales = "free_y", ncol = 1) +
  labs(x = "CpG O/E ratio", y = "Count",
       title = "CpG O/E distribution by promoter class") +
  theme_minimal(base_size = 12)
save_fig(p3d, BATCH_DIR, "fig3d_promoter_oe_by_class", w = 8, h = 9)

cat("\n=== Batch 03 complete ===\n")
cat("Figures:", length(list.files(file.path(BATCH_DIR, "figures"), pattern = "\\.png$")), "\n")
cat("Data files:", length(list.files(file.path(BATCH_DIR, "data"), pattern = "\\.tsv$")), "\n")
