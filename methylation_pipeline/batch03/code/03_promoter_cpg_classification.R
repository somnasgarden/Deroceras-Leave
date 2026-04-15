#!/usr/bin/env Rscript
# =============================================================================
# Batch 03: Promoter CpG Classification (Weber criteria)
# Question: What does the D. laeve promoter CpG landscape look like?
#           How does it compare to mammals?
# Output: data/ (promoter classification) + figures/ (5 plots)
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
weber_colors <- c(HCP = "#E31A1C", ICP = "#B8B8B8", LCP = "#1F78B4")

# Fig 3A: CpG O/E histogram of all promoters
p3a <- ggplot(prom_stats, aes(x = cpg_oe)) +
  geom_histogram(bins = 80, fill = "#2471A3", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "#E31A1C", linewidth = 0.8) +
  geom_vline(xintercept = 0.48, linetype = "dashed", color = "#1F78B4", linewidth = 0.8) +
  annotate("text", x = 0.78, y = Inf, label = "HCP threshold (0.75)", color = "#E31A1C",
           vjust = 2, hjust = 0, size = 3) +
  annotate("text", x = 0.35, y = Inf, label = "LCP threshold (0.48)", color = "#1F78B4",
           vjust = 2, hjust = 0, size = 3) +
  labs(x = "CpG O/E ratio", y = "Number of promoters",
       title = "Promoter CpG O/E distribution",
       subtitle = sprintf("D. laeve | %s promoters (TSS +/- 1kb) | Weber thresholds marked",
                          format(nrow(prom_stats), big.mark = ","))) +
  theme_minimal(base_size = 12)
save_fig(p3a, BATCH_DIR, "fig3a_promoter_cpg_oe_histogram", w = 9, h = 6)

# Fig 3B: CpG O/E vs GC% scatter (Weber classification)
# Plot ICP/LCP first (background), HCP last on top with larger points
prom_stats$plot_order <- match(prom_stats$weber_class, c("ICP", "LCP", "HCP"))
prom_sorted <- prom_stats[order(prom_stats$plot_order), ]
p3b <- ggplot(prom_sorted, aes(x = gc_pct, y = cpg_oe, color = weber_class, size = weber_class)) +
  geom_point(alpha = 0.2, shape = 16) +
  scale_color_manual(values = weber_colors, name = "Class") +
  scale_size_manual(values = c(HCP = 3, ICP = 0.3, LCP = 0.3), guide = "none") +
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
       subtitle = paste0("D. laeve vs mammals: ", class_pct["ICP"], "% ICP (mammals: ~16%)")) +
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

# --- HCP gene list (only ~15-20 exist, biologically significant) ---
hcp_genes <- prom_stats[prom_stats$weber_class == "HCP", ]
cat(sprintf("\nHCP genes: %d\n", nrow(hcp_genes)))
if (file.exists(OG$annot)) {
  annot_dt <- data.table::fread(OG$annot, header = FALSE,
                                 col.names = c("gene_id", "gene_name", "description"))
  hcp_genes <- merge(hcp_genes, annot_dt, by = "gene_id", all.x = TRUE)
  cat("HCP genes with names:\n")
  for (i in seq_len(nrow(hcp_genes))) {
    cat(sprintf("  %s  %-20s  O/E=%.3f  GC=%.1f%%  %s\n",
                hcp_genes$gene_id[i],
                ifelse(!is.na(hcp_genes$gene_name[i]) & hcp_genes$gene_name[i] != "",
                       hcp_genes$gene_name[i], "(unnamed)"),
                hcp_genes$cpg_oe[i], hcp_genes$gc_pct[i],
                ifelse(!is.na(hcp_genes$description[i]), substr(hcp_genes$description[i], 1, 60), "")))
  }
}
save_data(hcp_genes, BATCH_DIR, "hcp_gene_list")

# Fig 3E: HCP gene summary (O/E, methylation, DMP/DMR status)
hcp_plot <- hcp_genes[!duplicated(hcp_genes$gene_id), ]
hcp_plot$label <- ifelse(!is.na(hcp_plot$gene_name) & hcp_plot$gene_name != "",
                         hcp_plot$gene_name, hcp_plot$gene_id)
hcp_plot$desc_short <- ifelse(!is.na(hcp_plot$description) & hcp_plot$description != "",
                              substr(hcp_plot$description, 1, 50), "Unknown function")

# Add methylation + DMP/DMR info if available
hcp_plot$ctrl_meth <- NA_real_
hcp_plot$has_dmp <- FALSE
hcp_plot$has_dmr <- FALSE

if (file.exists(CACHE$bsseq)) {
  cat("  Adding methylation levels to HCP plot...\n")
  suppressPackageStartupMessages(library(bsseq))
  bs <- readRDS(CACHE$bsseq)
  bs_gr <- GRanges(seqnames = seqnames(bs), ranges = IRanges(start = start(bs), width = 1))
  ctrl_beta <- rowMeans(getMeth(bs[, c("C1","C2")], type = "raw"), na.rm = TRUE)

  for (i in seq_len(nrow(hcp_plot))) {
    pg <- GRanges(seqnames = hcp_plot$chr[i],
                  ranges = IRanges(start = hcp_plot$tss[i] - 1000, end = hcp_plot$tss[i] + 1000))
    idx <- queryHits(findOverlaps(bs_gr, pg))
    if (length(idx) > 0) hcp_plot$ctrl_meth[i] <- mean(ctrl_beta[idx], na.rm = TRUE) * 100
  }
  rm(bs, bs_gr, ctrl_beta); gc(verbose = FALSE)
}

# Check for DMPs/DMRs in HCP promoter regions
dmp_file <- file.path(PIPE_DIR, "batch06", "data", "dmps_annotated.tsv")
dmr_file <- file.path(PIPE_DIR, "batch06", "data", "dmrs_annotated.tsv")
if (file.exists(dmp_file)) {
  dmps <- data.table::fread(dmp_file)
  dmp_gr <- GRanges(seqnames = dmps$chr, ranges = IRanges(start = dmps$pos, width = 1))
  for (i in seq_len(nrow(hcp_plot))) {
    pg <- GRanges(seqnames = hcp_plot$chr[i],
                  ranges = IRanges(start = hcp_plot$tss[i] - 1000, end = hcp_plot$tss[i] + 1000))
    hcp_plot$has_dmp[i] <- any(overlapsAny(dmp_gr, pg))
  }
}
if (file.exists(dmr_file)) {
  dmrs <- data.table::fread(dmr_file)
  dmr_gr <- GRanges(seqnames = dmrs$chr, ranges = IRanges(start = dmrs$start, end = dmrs$end))
  for (i in seq_len(nrow(hcp_plot))) {
    pg <- GRanges(seqnames = hcp_plot$chr[i],
                  ranges = IRanges(start = hcp_plot$tss[i] - 1000, end = hcp_plot$tss[i] + 1000))
    hcp_plot$has_dmr[i] <- any(overlapsAny(dmr_gr, pg))
  }
}

# Build display label with methylation and DM status
hcp_plot$meth_label <- ifelse(!is.na(hcp_plot$ctrl_meth),
                              sprintf("%.1f%%", hcp_plot$ctrl_meth), "no cov")
hcp_plot$dm_label <- ifelse(hcp_plot$has_dmp | hcp_plot$has_dmr,
                            paste0(ifelse(hcp_plot$has_dmp, "DMP", ""),
                                   ifelse(hcp_plot$has_dmp & hcp_plot$has_dmr, "+", ""),
                                   ifelse(hcp_plot$has_dmr, "DMR", "")), "")
hcp_plot$display <- paste0(hcp_plot$label, " — ", hcp_plot$desc_short)
hcp_plot <- hcp_plot[order(hcp_plot$cpg_oe), ]
hcp_plot$display <- factor(hcp_plot$display, levels = hcp_plot$display)

# Determine if all HCP are unmethylated (CpG island-like)
hcp_mean_meth <- mean(hcp_plot$ctrl_meth, na.rm = TRUE)
cat(sprintf("  HCP mean promoter methylation: %.2f%% (%s)\n",
            hcp_mean_meth, ifelse(hcp_mean_meth < 5, "unmethylated — CpG island-like", "methylated")))
cat(sprintf("  HCP with DMPs: %d/%d | HCP with DMRs: %d/%d\n",
            sum(hcp_plot$has_dmp), nrow(hcp_plot), sum(hcp_plot$has_dmr), nrow(hcp_plot)))

p3e <- ggplot(hcp_plot, aes(x = cpg_oe, y = display)) +
  geom_segment(aes(xend = 0.75, yend = display), color = "gray60", linewidth = 0.5) +
  geom_point(aes(fill = ctrl_meth), size = 4, shape = 21, color = "black", stroke = 0.5) +
  geom_text(aes(label = meth_label), hjust = -0.3, size = 3, color = "gray30") +
  geom_point(data = hcp_plot[hcp_plot$has_dmp | hcp_plot$has_dmr, ],
             aes(x = cpg_oe), shape = 8, size = 2, color = "#E31A1C") +
  geom_vline(xintercept = 0.75, linetype = "dashed", color = "gray40") +
  scale_fill_gradient(low = "white", high = "#E74C3C", na.value = "gray80",
                      name = "Meth (%)", limits = c(0, 10)) +
  labs(x = "CpG O/E ratio", y = NULL,
       title = sprintf("HCP promoters in D. laeve (n = %d) — CpG island candidates",
                       nrow(hcp_plot)),
       subtitle = sprintf("Fill: control methylation. Star: has DMP/DMR. Mean meth = %.1f%%",
                          hcp_mean_meth)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 9))
save_fig(p3e, BATCH_DIR, "fig3e_hcp_gene_list", w = 11, h = 6)
save_data(hcp_plot[, c("gene_id", "label", "desc_short", "cpg_oe", "gc_pct",
                        "ctrl_meth", "has_dmp", "has_dmr")],
          BATCH_DIR, "hcp_gene_summary")

# =============================================================================
# NEW ANALYSES: Promoter class × methylation × gene structure × DM
# =============================================================================

# --- Fig 3F: Weber-style histogram colored by class (single panel, overlaid) ---
cat("\n[NEW] Weber-style CpG O/E histogram colored by class...\n")
p3f <- ggplot(prom_stats, aes(x = cpg_oe, fill = weber_class)) +
  geom_histogram(bins = 80, color = "white", linewidth = 0.15, position = "stack") +
  scale_fill_manual(values = weber_colors, name = "Promoter class") +
  geom_vline(xintercept = c(0.48, 0.75), linetype = "dashed", color = "gray30", linewidth = 0.6) +
  labs(x = "CpG O/E ratio", y = "Number of promoters",
       title = "Promoter CpG O/E distribution by class (Weber et al. 2007)",
       subtitle = sprintf("D. laeve: %s promoters | HCP %.1f%% | ICP %.1f%% | LCP %.1f%%",
                          format(nrow(prom_stats), big.mark = ","),
                          class_pct["HCP"], class_pct["ICP"], class_pct["LCP"])) +
  theme_minimal(base_size = 12)
save_fig(p3f, BATCH_DIR, "fig3f_weber_histogram_by_class", w = 10, h = 6)

# --- Fig 3G: Mean methylation level by promoter class (boxplot) ---
cat("[NEW] Methylation by promoter class...\n")
if (file.exists(CACHE$bsseq)) {
  suppressPackageStartupMessages(library(bsseq))
  bs <- readRDS(CACHE$bsseq)
  bs_gr <- GRanges(seqnames = seqnames(bs), ranges = IRanges(start = start(bs), width = 1))
  ctrl_beta <- rowMeans(getMeth(bs[, c("C1","C2")], type = "raw"), na.rm = TRUE)
  ampu_beta <- rowMeans(getMeth(bs[, c("A1","A2")], type = "raw"), na.rm = TRUE)

  # Build promoter GRanges for all genes (TSS +/- 1kb, same as classification)
  prom_class_gr <- GRanges(seqnames = prom_stats$chr,
                            ranges = IRanges(start = pmax(1L, prom_stats$tss - 1000L),
                                             end = prom_stats$tss + 1000L))
  prom_class_gr$gene_id <- prom_stats$gene_id
  prom_class_gr$weber_class <- prom_stats$weber_class

  hits_pm <- findOverlaps(bs_gr, prom_class_gr)
  prom_meth_dt <- data.table(
    gene_id = prom_class_gr$gene_id[subjectHits(hits_pm)],
    weber_class = as.character(prom_class_gr$weber_class[subjectHits(hits_pm)]),
    ctrl_beta = ctrl_beta[queryHits(hits_pm)],
    ampu_beta = ampu_beta[queryHits(hits_pm)]
  )
  prom_meth_dt <- prom_meth_dt[!is.na(ctrl_beta)]

  prom_gene_meth <- prom_meth_dt[, .(ctrl_meth = mean(ctrl_beta) * 100,
                                      ampu_meth = mean(ampu_beta) * 100,
                                      n_cpgs = .N), by = .(gene_id, weber_class)]

  cat("  Promoter methylation by class:\n")
  class_meth_summary <- prom_gene_meth[, .(n = .N,
    mean_ctrl = round(mean(ctrl_meth), 2),
    median_ctrl = round(median(ctrl_meth), 2),
    mean_ampu = round(mean(ampu_meth), 2)), by = weber_class]
  print(class_meth_summary)
  save_data(prom_gene_meth, BATCH_DIR, "promoter_methylation_by_class")

  prom_gene_meth[, weber_class := factor(weber_class, levels = c("HCP", "ICP", "LCP"))]
  p3g <- ggplot(prom_gene_meth, aes(x = weber_class, y = ctrl_meth, fill = weber_class)) +
    geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
    scale_fill_manual(values = weber_colors, guide = "none") +
    labs(x = "Promoter class", y = "Mean promoter methylation (%)",
         title = "Promoter methylation by Weber class (control)",
         subtitle = sprintf("HCP: %.1f%% | ICP: %.1f%% | LCP: %.1f%% (median)",
                            class_meth_summary[weber_class == "HCP", median_ctrl],
                            class_meth_summary[weber_class == "ICP", median_ctrl],
                            class_meth_summary[weber_class == "LCP", median_ctrl])) +
    theme_minimal(base_size = 12)
  save_fig(p3g, BATCH_DIR, "fig3g_methylation_by_promoter_class", w = 8, h = 6)

  # --- Fig 3H: DMP and DMR proportion by promoter class (100% stacked) ---
  cat("[NEW] DMP/DMR proportion by promoter class...\n")
  dmp_file <- file.path(PIPE_DIR, "batch06", "data", "dmps_annotated.tsv")
  dmr_file <- file.path(PIPE_DIR, "batch06", "data", "dmrs_annotated.tsv")

  if (file.exists(dmp_file) && file.exists(dmr_file)) {
    dmps <- fread(dmp_file)
    dmrs <- fread(dmr_file)
    dmp_gr_all <- GRanges(seqnames = dmps$chr, ranges = IRanges(start = dmps$pos, width = 1))
    dmr_gr_all <- GRanges(seqnames = dmrs$chr, ranges = IRanges(start = dmrs$start, end = dmrs$end))

    # Assign DMPs to promoter class
    dmp_prom_hits <- findOverlaps(dmp_gr_all, prom_class_gr)
    dmp_class <- data.table(
      dmp_idx = queryHits(dmp_prom_hits),
      weber_class = as.character(prom_class_gr$weber_class[subjectHits(dmp_prom_hits)])
    )
    dmp_class <- dmp_class[!duplicated(dmp_idx)]  # each DMP counted once
    dmp_class_counts <- dmp_class[, .N, by = weber_class]
    # DMPs NOT in any promoter
    dmp_not_prom <- nrow(dmps) - nrow(dmp_class)
    dmp_class_counts <- rbind(dmp_class_counts,
                               data.table(weber_class = "Non-promoter", N = dmp_not_prom))

    # Assign DMRs to promoter class
    dmr_prom_hits <- findOverlaps(dmr_gr_all, prom_class_gr)
    dmr_class <- data.table(
      dmr_idx = queryHits(dmr_prom_hits),
      weber_class = as.character(prom_class_gr$weber_class[subjectHits(dmr_prom_hits)])
    )
    dmr_class <- dmr_class[!duplicated(dmr_idx)]
    dmr_class_counts <- dmr_class[, .N, by = weber_class]
    dmr_not_prom <- nrow(dmrs) - nrow(dmr_class)
    dmr_class_counts <- rbind(dmr_class_counts,
                               data.table(weber_class = "Non-promoter", N = dmr_not_prom))

    # Combine for 100% stacked bar
    dmp_class_counts[, type := "DMPs"]
    dmr_class_counts[, type := "DMRs"]
    dm_bar <- rbind(dmp_class_counts, dmr_class_counts)
    dm_bar[, pct := 100 * N / sum(N), by = type]
    dm_bar[, weber_class := factor(weber_class, levels = c("HCP", "ICP", "LCP", "Non-promoter"))]

    bar_colors <- c(weber_colors, "Non-promoter" = "#95A5A6")
    p3h <- ggplot(dm_bar, aes(x = type, y = pct, fill = weber_class)) +
      geom_col(width = 0.6, color = "white", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.1f%%", pct)),
                position = position_stack(vjust = 0.5), size = 3.5) +
      scale_fill_manual(values = bar_colors, name = "Promoter class") +
      labs(x = NULL, y = "Proportion (%)",
           title = "DMP and DMR distribution by promoter class",
           subtitle = "100% stacked — most DM events are outside promoters") +
      theme_minimal(base_size = 12)
    save_fig(p3h, BATCH_DIR, "fig3h_dmp_dmr_by_promoter_class", w = 8, h = 6)
    save_data(dm_bar, BATCH_DIR, "dmp_dmr_by_promoter_class")
    cat("  DMP/DMR by class:\n"); print(dm_bar)
  } else {
    cat("  Skipping DMP/DMR plot — batch06 output not found.\n")
  }

  # --- Fig 3I: Metagene methylation by promoter class ---
  cat("[NEW] Metagene methylation by promoter class...\n")
  exons <- gff[gff$type == "exon"]
  N_BINS <- 20; FLANK <- 5000

  gene_strand_v <- as.character(strand(genes))
  gene_starts <- start(genes); gene_ends <- end(genes)
  gene_tss <- ifelse(gene_strand_v == "+", gene_starts, gene_ends)
  gene_tts <- ifelse(gene_strand_v == "+", gene_ends, gene_starts)
  gene_chr_v <- as.character(seqnames(genes))

  # Map promoter class to gene index
  prom_class_map <- prom_stats[, c("gene_id", "weber_class")]
  gene_class <- prom_class_map$weber_class[match(genes$ID, prom_class_map$gene_id)]

  # Build segments: distal upstream, promoter, gene body, downstream
  seg_order <- c("Distal upstream\n(5-2 kb)", "Promoter\n(2 kb)", "Gene body", "TTS downstream\n(5 kb)")

  seg_list <- vector("list", length(genes) * 4L)
  k <- 0L
  valid_idx <- which(!is.na(gene_class) & genes$ID %in% keep_chr)
  # Use all genes with a class assignment
  valid_idx <- which(!is.na(gene_class))

  for (gi in valid_idx) {
    g_chr <- gene_chr_v[gi]; g_strand <- gene_strand_v[gi]
    g_tss <- gene_tss[gi]; g_tts <- gene_tts[gi]
    g_start <- gene_starts[gi]; g_end <- gene_ends[gi]
    gc <- as.character(gene_class[gi])

    segs <- list(
      `Distal upstream\n(5-2 kb)` = c(max(1L, g_tss - FLANK), g_tss - 2001L),
      `Promoter\n(2 kb)` = c(max(1L, g_tss - 2000L), g_tss - 1L),
      `Gene body` = c(g_start, g_end),
      `TTS downstream\n(5 kb)` = c(g_tts + 1L, g_tts + FLANK))

    for (sn in names(segs)) {
      s <- segs[[sn]]
      if (is.na(s[1]) || is.na(s[2]) || s[1] >= s[2] || s[1] < 1L) next
      k <- k + 1L
      seg_list[[k]] <- data.table(
        chr = g_chr, start = s[1], end = s[2], strand = g_strand,
        gene_id = genes$ID[gi], segment = sn, weber_class = gc)
    }
  }
  seg_dt <- rbindlist(seg_list[1:k])
  cat(sprintf("  Segments built: %s for %d genes\n", format(nrow(seg_dt), big.mark = ","),
              length(valid_idx)))

  seg_gr <- GRanges(seqnames = seg_dt$chr,
                     ranges = IRanges(start = seg_dt$start, end = seg_dt$end))
  bs_chr <- as.character(seqnames(bs_gr))
  bs_pos <- start(bs_gr)
  ov <- findOverlaps(bs_gr, seg_gr)

  ov_dt <- data.table(
    cpg_idx = queryHits(ov), seg_idx = subjectHits(ov),
    beta = ctrl_beta[queryHits(ov)], pos = bs_pos[queryHits(ov)])
  ov_dt <- ov_dt[!is.na(beta)]
  ov_dt[, `:=`(seg_start = seg_dt$start[seg_idx], seg_end = seg_dt$end[seg_idx],
               seg_strand = seg_dt$strand[seg_idx], segment = seg_dt$segment[seg_idx],
               weber_class = seg_dt$weber_class[seg_idx],
               gene_id = seg_dt$gene_id[seg_idx])]

  ov_dt[, rel_pos := (pos - seg_start) / (seg_end - seg_start)]
  ov_dt[seg_strand == "-" & segment != "Gene body", rel_pos := 1 - rel_pos]
  ov_dt[, bin := pmin(N_BINS, pmax(1L, ceiling(rel_pos * N_BINS)))]

  gene_bins <- ov_dt[, .(mean_beta = mean(beta) * 100),
                     by = .(gene_id, segment, bin, weber_class)]
  mg_agg <- gene_bins[, .(mean_meth = mean(mean_beta), se = sd(mean_beta) / sqrt(.N)),
                       by = .(segment, bin, weber_class)]
  mg_agg[, segment := factor(segment, levels = seg_order)]
  mg_agg[, x_pos := (as.integer(segment) - 1) * N_BINS + bin]
  mg_agg[, weber_class := factor(weber_class, levels = c("HCP", "ICP", "LCP"))]

  n_segs <- length(seg_order)
  seg_breaks <- (1:(n_segs - 1)) * N_BINS + 0.5
  seg_labels_pos <- (0:(n_segs - 1)) * N_BINS + N_BINS / 2

  p3i <- ggplot(mg_agg, aes(x = x_pos, y = mean_meth, color = weber_class)) +
    geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 1.2) +
    geom_vline(xintercept = seg_breaks, linetype = "dashed", color = "gray70", linewidth = 0.3) +
    scale_color_manual(values = weber_colors, name = "Promoter class") +
    scale_x_continuous(breaks = seg_labels_pos, labels = seg_order) +
    labs(title = "Metagene methylation profile by promoter class",
         subtitle = sprintf("HCP n=%d | ICP n=%d | LCP n=%d",
                            sum(gene_class == "HCP", na.rm = TRUE),
                            sum(gene_class == "ICP", na.rm = TRUE),
                            sum(gene_class == "LCP", na.rm = TRUE)),
         x = NULL, y = "Methylation level (%)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(size = 9), panel.grid.minor = element_blank())
  save_fig(p3i, BATCH_DIR, "fig3i_metagene_by_promoter_class", w = 11, h = 6)

  rm(ov_dt, gene_bins, mg_agg, seg_dt, seg_gr, ov); gc(verbose = FALSE)
  rm(bs, bs_gr, ctrl_beta, ampu_beta); gc(verbose = FALSE)
} else {
  cat("  Skipping methylation-based figures — BSseq cache not found.\n")
}

# --- Fig 3J: Gene structure by promoter class (intronless, 2-exon, 3+ exon) ---
cat("[NEW] Gene structure × promoter class...\n")
exons_all <- gff[gff$type == "exon"]
exon_hits <- findOverlaps(exons_all, genes)
exon_per_gene <- data.table(gene_idx = subjectHits(exon_hits),
                             exon_idx = queryHits(exon_hits))
exon_counts <- exon_per_gene[, .(n_exons = uniqueN(exon_idx)), by = gene_idx]
exon_counts[, gene_id := genes$ID[gene_idx]]
exon_counts[, structure := fifelse(n_exons == 1, "Intronless (1 exon)",
                             fifelse(n_exons == 2, "2-exon", "3+ exons"))]

# Merge with promoter class
struct_class <- merge(
  data.table(gene_id = prom_stats$gene_id, weber_class = prom_stats$weber_class),
  exon_counts[, .(gene_id, structure, n_exons)],
  by = "gene_id", all.x = TRUE)
struct_class[is.na(structure), structure := "Intronless (1 exon)"]
struct_class[, weber_class := factor(weber_class, levels = c("HCP", "ICP", "LCP"))]
struct_class[, structure := factor(structure, levels = c("Intronless (1 exon)", "2-exon", "3+ exons"))]

# Cross-tab
cross_tab <- struct_class[, .N, by = .(weber_class, structure)]
cross_tab[, pct := round(100 * N / sum(N), 1), by = weber_class]
cat("\n  Gene structure × promoter class:\n")
print(dcast(cross_tab, weber_class ~ structure, value.var = "N", fill = 0))
cat("\n  Percentages within each class:\n")
print(dcast(cross_tab, weber_class ~ structure, value.var = "pct", fill = 0))
save_data(cross_tab, BATCH_DIR, "gene_structure_by_promoter_class")

# Plot: grouped bar (percentage within each class)
struct_colors <- c("Intronless (1 exon)" = "#E74C3C", "2-exon" = "#F39C12", "3+ exons" = "#2471A3")
p3j <- ggplot(cross_tab, aes(x = weber_class, y = pct, fill = structure)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", N, pct)),
            position = position_dodge(width = 0.7), vjust = -0.2, size = 3) +
  scale_fill_manual(values = struct_colors, name = "Gene structure") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Promoter class", y = "% of genes in class",
       title = "Gene structure composition by promoter class",
       subtitle = "Do intronless genes cluster in a specific promoter class?") +
  theme_minimal(base_size = 12)
save_fig(p3j, BATCH_DIR, "fig3j_gene_structure_by_class", w = 10, h = 7)

# --- Reverse view: what promoter class are intronless vs multi-exon genes? ---
cross_tab_rev <- struct_class[, .N, by = .(structure, weber_class)]
cross_tab_rev[, pct := round(100 * N / sum(N), 1), by = structure]

p3j2 <- ggplot(cross_tab_rev, aes(x = structure, y = pct, fill = weber_class)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", N, pct)),
            position = position_dodge(width = 0.7), vjust = -0.2, size = 3) +
  scale_fill_manual(values = weber_colors, name = "Promoter class") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Gene structure", y = "% of genes in structure group",
       title = "Promoter class composition by gene structure",
       subtitle = "Are intronless gene promoters enriched for LCP?") +
  theme_minimal(base_size = 12)
save_fig(p3j2, BATCH_DIR, "fig3j2_promoter_class_by_structure", w = 10, h = 7)
save_data(cross_tab_rev, BATCH_DIR, "promoter_class_by_gene_structure")

# --- Print interpretation ---
cat("\n=== PATTERN ANALYSIS ===\n")
n_hcp <- sum(prom_stats$weber_class == "HCP")
n_icp <- sum(prom_stats$weber_class == "ICP")
n_lcp <- sum(prom_stats$weber_class == "LCP")
cat(sprintf("Promoter classes: HCP=%d (%.1f%%) | ICP=%d (%.1f%%) | LCP=%d (%.1f%%)\n",
            n_hcp, 100*n_hcp/nrow(prom_stats), n_icp, 100*n_icp/nrow(prom_stats),
            n_lcp, 100*n_lcp/nrow(prom_stats)))

intronless_lcp <- cross_tab[structure == "Intronless (1 exon)" & weber_class == "LCP", N]
intronless_icp <- cross_tab[structure == "Intronless (1 exon)" & weber_class == "ICP", N]
multi_icp <- cross_tab[structure == "3+ exons" & weber_class == "ICP", N]
cat(sprintf("Intronless genes: LCP=%d, ICP=%d\n", intronless_lcp, intronless_icp))
cat(sprintf("3+ exon genes in ICP: %d\n", multi_icp))

cat("\nSpeculations to investigate:\n")
cat("  1. If LCP is enriched in intronless genes: these may be evolutionary\n")
cat("     newcomers (retrogenes, TF-derived) that lack CpG enrichment at promoters.\n")
cat("     Without introns, they also lack gene body methylation signal.\n")
cat("  2. If ICP 3+ exon genes have higher gene body meth: the classic invertebrate\n")
cat("     pattern. Gene body meth scales with gene length and intron content.\n")
cat("     Weber's ICP = the 'default' promoter in a DNMT3-absent genome.\n")
cat("  3. HCPs (if any) should be unmethylated and housekeeping — like mammalian\n")
cat("     CpG islands. In Weber, HCPs are protected from de novo methylation.\n")
cat("     Without DNMT3, the protection mechanism is irrelevant, which may explain\n")
cat("     why D. laeve has almost none: no selective pressure to maintain CpG islands\n")
cat("     when there is no de novo methylation threat.\n")
cat("  4. DMP/DMR enrichment: if DMPs cluster in ICP promoters, regeneration\n")
cat("     may remodel the 'active' promoter class. If they avoid promoters entirely\n")
cat("     (non-promoter dominant), the gene body methylation story holds:\n")
cat("     regeneration acts on gene bodies, not promoters.\n")
cat("  5. Metagene by class: if ICP genes show the classic bell-shaped gene body\n")
cat("     methylation curve while LCP genes are flat, this confirms that promoter\n")
cat("     CpG architecture predicts the whole gene's methylation pattern.\n")

cat("\n=== Batch 03 complete ===\n")
cat("Figures:", length(list.files(file.path(BATCH_DIR, "figures"), pattern = "\\.png$")), "\n")
cat("Data files:", length(list.files(file.path(BATCH_DIR, "data"), pattern = "\\.tsv$")), "\n")
