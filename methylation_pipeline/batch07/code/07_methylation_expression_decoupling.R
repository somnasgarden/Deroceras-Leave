#!/usr/bin/env Rscript
# =============================================================================
# Batch 07: Differential Methylation vs Expression
# Question: Does differential methylation predict expression changes?
# Output: data/ + figures/
# Requires: Batch 06 DMPs/DMRs, CACHE$transcriptome or HTSeq counts
# NOTE: DMPs and DMRs are treated SEPARATELY (DMPs can be inside DMRs)
#
# Gene assignment: DMP/DMR assigned to ALL genes within reach:
#   promoter (2kb upstream) + distal upstream (5kb) + gene body + downstream (5kb)
#   A DMP between two genes can be assigned to BOTH.
# =============================================================================

source("methylation_pipeline/_config.R")
t0 <- proc.time()

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(ggplot2)
  library(ggrepel)
})

BATCH_DIR <- file.path(PIPE_DIR, "batch07")

# Clean old output
unlink(file.path(BATCH_DIR, "figures"), recursive = TRUE)
unlink(file.path(BATCH_DIR, "data"), recursive = TRUE)
dir.create(file.path(BATCH_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BATCH_DIR, "data"), showWarnings = FALSE, recursive = TRUE)

cat("=== Batch 07: Differential Methylation vs Expression ===\n\n")

# =============================================================================
# 1. LOAD DMPs + DMRs FROM BATCH 06
# =============================================================================
cat("[1/7] Loading DMPs and DMRs from batch06...\n")
dmp_file <- file.path(PIPE_DIR, "batch06/data/dmps_annotated.tsv")
dmr_file <- file.path(PIPE_DIR, "batch06/data/dmrs_annotated.tsv")
if (!file.exists(dmp_file)) stop("Run batch06 first.")

dmp <- fread(dmp_file)
# Deduplicate (both strands produce duplicate rows)
dmp <- unique(dmp, by = c("chr", "pos"))
cat(sprintf("  DMPs: %s (unique positions)\n", format(nrow(dmp), big.mark = ",")))

has_dmr <- file.exists(dmr_file)
if (has_dmr) {
  dmr <- fread(dmr_file)
  dmr <- unique(dmr, by = c("chr", "start", "end"))
  has_dmr <- nrow(dmr) > 0
  if (has_dmr) cat(sprintf("  DMRs: %s (unique regions)\n", format(nrow(dmr), big.mark = ",")))
}

# =============================================================================
# 2. EXPRESSION DATA
# =============================================================================
cat("[2/7] Loading expression data...\n")
if (file.exists(CACHE$transcriptome)) {
  cat("  Loading from CACHE$transcriptome...\n")
  tc <- readRDS(CACHE$transcriptome)
  de_results <- as.data.table(as.data.frame(tc$res_tail), keep.rownames = "gene_id")
  de_results <- de_results[!is.na(padj)]
  cat(sprintf("  Outliers removed: %s\n",
              if (length(tc$outliers) > 0) paste(tc$outliers, collapse = ", ") else "none"))
  rm(tc); gc(verbose = FALSE)
} else {
  cat("  Cache not found — running DESeq2...\n")
  suppressPackageStartupMessages(library(DESeq2))
  all_files <- list.files(OG$counts_dir, pattern = "\\.txt$")
  tail_ctrl <- all_files[grepl("^C\\d+S\\d+_", all_files)]
  tail_ampu <- all_files[grepl("^T\\d+S\\d+_", all_files)]
  tail_ampu <- tail_ampu[!grepl("^T1S5_", tail_ampu)]
  tail_files <- c(tail_ctrl, tail_ampu)
  sample_table <- data.frame(sampleName = tail_files, fileName = tail_files,
                             condition = ifelse(grepl("^C\\d+S", tail_files), "control", "amputated"),
                             stringsAsFactors = FALSE)
  dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table,
                                     directory = OG$counts_dir, design = ~ condition)
  dds$condition <- relevel(dds$condition, ref = "control")
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  de_results <- as.data.table(as.data.frame(results(dds, alpha = 0.05)), keep.rownames = "gene_id")
  de_results <- de_results[!is.na(padj)]
  rm(dds); gc(verbose = FALSE)
}
n_de <- sum(de_results$padj < 0.05)
cat(sprintf("  DE genes (FDR < 0.05): %d / %d\n", n_de, nrow(de_results)))

# Gene names
gene_annot <- if (file.exists(OG$annot)) {
  fread(OG$annot, header = FALSE, col.names = c("gene_id", "gene_name", "description"))
} else { data.table(gene_id = character(), gene_name = character(), description = character()) }
# Deduplicate: keep first name per gene
gene_annot <- gene_annot[!duplicated(gene_id)]

# =============================================================================
# 3. BUILD GENE WINDOWS: promoter 2kb + distal 5kb + gene body + downstream 5kb
# =============================================================================
cat("[3/7] Building gene windows (promoter 2kb + distal 5kb + body + downstream 5kb)...\n")
gff <- load_gff()
genes <- gff[gff$type == "gene"]

# Get chromosome lengths for trimming
genome <- load_genome()
chr_lengths <- setNames(as.integer(width(genome)), names(genome))
rm(genome); gc(verbose = FALSE)

# Build per-gene windows: extend 5kb upstream of promoter + 5kb downstream
gene_strand <- as.character(strand(genes))
gene_chr <- as.character(seqnames(genes))
gene_start <- start(genes)
gene_end <- end(genes)

# Upstream: 5kb distal + 2kb promoter = 7kb upstream of TSS
# Downstream: 5kb downstream of gene end
tss <- ifelse(gene_strand == "+", gene_start, gene_end)
tes <- ifelse(gene_strand == "+", gene_end, gene_start)  # transcription end site

win_start <- pmax(1L, ifelse(gene_strand == "+", tss - 7000L, tes - 5000L))
win_end   <- ifelse(gene_strand == "+", tes + 5000L, tss + 7000L)
win_end   <- pmin(win_end, chr_lengths[gene_chr])

gene_windows <- GRanges(seqnames = gene_chr,
                        ranges = IRanges(start = win_start, end = win_end),
                        gene_id = genes$ID)

cat(sprintf("  Gene windows built: %d genes (body + 7kb up + 5kb down)\n", length(gene_windows)))

# Also build expanded annotation for per-region analysis
expanded <- build_expanded_regions(gff)

# =============================================================================
# 4. ASSIGN DMPs/DMRs TO GENES VIA WINDOW OVERLAP
# =============================================================================
cat("[4/7] Assigning DMPs/DMRs to genes via window overlap...\n")

# DMP → gene assignment (one DMP can map to multiple genes)
dmp_gr <- GRanges(seqnames = dmp$chr, ranges = IRanges(start = dmp$pos, width = 1))
dmp_hits <- findOverlaps(dmp_gr, gene_windows)
dmp_gene_map <- data.table(
  dmp_idx = queryHits(dmp_hits),
  gene_id = mcols(gene_windows)$gene_id[subjectHits(dmp_hits)]
)
# Add DMP info
dmp_gene_map[, `:=`(chr = dmp$chr[dmp_idx], pos = dmp$pos[dmp_idx],
                     diff = dmp$diff[dmp_idx], direction = dmp$direction[dmp_idx],
                     annotation = dmp$annotation[dmp_idx])]

n_dmp_assigned <- length(unique(dmp_gene_map$dmp_idx))
n_dmp_genes <- length(unique(dmp_gene_map$gene_id))
cat(sprintf("  DMPs assigned to genes: %d / %d (%.1f%%)\n",
            n_dmp_assigned, nrow(dmp), 100 * n_dmp_assigned / nrow(dmp)))
cat(sprintf("  Genes with DMPs: %d\n", n_dmp_genes))

# Expanded annotation for each DMP
dmp[, annotation_expanded := annotate_regions_expanded(chr, pos, expanded)]

# DMR → gene assignment
if (has_dmr) {
  dmr_gr <- GRanges(seqnames = dmr$chr, ranges = IRanges(start = dmr$start, end = dmr$end))
  dmr_hits <- findOverlaps(dmr_gr, gene_windows)
  dmr_gene_map <- data.table(
    dmr_idx = queryHits(dmr_hits),
    gene_id = mcols(gene_windows)$gene_id[subjectHits(dmr_hits)]
  )
  dmr_gene_map[, `:=`(chr = dmr$chr[dmr_idx], start = dmr$start[dmr_idx],
                        end = dmr$end[dmr_idx], diff.Methy = dmr$diff.Methy[dmr_idx],
                        nCG = dmr$nCG[dmr_idx], direction = dmr$direction[dmr_idx],
                        annotation = dmr$annotation[dmr_idx])]
  n_dmr_assigned <- length(unique(dmr_gene_map$dmr_idx))
  n_dmr_genes <- length(unique(dmr_gene_map$gene_id))
  cat(sprintf("  DMRs assigned to genes: %d / %d (%.1f%%)\n",
              n_dmr_assigned, nrow(dmr), 100 * n_dmr_assigned / nrow(dmr)))
  cat(sprintf("  Genes with DMRs: %d\n", n_dmr_genes))
}

# =============================================================================
# 5. DMP ANALYSIS
# =============================================================================
cat("\n[5/7] === DMP ANALYSIS ===\n")

# Gene-level DMP table (using window assignment)
gene_dmp <- dmp_gene_map[, .(n_dmp = .N, mean_diff = mean(diff),
                              n_hyper = sum(diff > 0), n_hypo = sum(diff < 0)),
                          by = gene_id]
gene_dmp <- merge(gene_dmp, de_results[, .(gene_id, log2FoldChange, padj, baseMean)],
                  by = "gene_id")
gene_dmp[, de := padj < 0.05]
gene_dmp[, meth_direction := ifelse(mean_diff > 0, "Hyper", "Hypo")]
gene_dmp <- merge(gene_dmp, gene_annot[, .(gene_id, gene_name)],
                  by = "gene_id", all.x = TRUE)
gene_dmp[, dmp_bin := cut(n_dmp, breaks = c(0, 1, 3, 5, 10, Inf),
                           labels = c("1", "2-3", "4-5", "6-10", ">10"))]

cat(sprintf("  Genes with DMP + expression: %d\n", nrow(gene_dmp)))
cat(sprintf("  Genes with DMP + DE: %d (%.1f%%)\n",
            sum(gene_dmp$de), 100 * mean(gene_dmp$de)))

# DMP correlation by basic region
gene_dmp_region <- dmp_gene_map[, .(n_dmp = .N, mean_diff = mean(diff)),
                                 by = .(gene_id, region = annotation)]
dmp_region_cors <- merge(gene_dmp_region, de_results[, .(gene_id, log2FoldChange, padj)],
                         by = "gene_id")[, {
  if (.N >= 20) {
    ct <- tryCatch(cor.test(mean_diff, log2FoldChange, method = "spearman"),
                   error = function(e) list(estimate = NA, p.value = NA))
    .(n_genes = .N, spearman_rho = ct$estimate, spearman_p = ct$p.value,
      pct_de = round(100 * mean(padj < 0.05), 1))
  }
}, by = region][order(spearman_p)]
cat("\n  DMP correlation by region:\n"); print(as.data.frame(dmp_region_cors))

# DMP overall correlation
cor_dmp <- cor.test(gene_dmp$mean_diff, gene_dmp$log2FoldChange, method = "spearman")
cat(sprintf("  DMP overall: rho = %.4f, p = %.2e\n", cor_dmp$estimate, cor_dmp$p.value))

# DMP DE rate by region
dmp_de_rate <- merge(gene_dmp_region, de_results[, .(gene_id, padj)],
                     by = "gene_id")[, .(n_genes = uniqueN(gene_id),
                                         n_de = uniqueN(gene_id[padj < 0.05])),
                                     by = region]
dmp_de_rate[, de_rate := round(100 * n_de / n_genes, 1)]

# DMP dose-response
dmp_dose <- gene_dmp[, .(genes = .N, de_genes = sum(de), de_rate = round(100 * mean(de), 1)),
                      by = dmp_bin][order(dmp_bin)]
cat("\n  DMP dose-response:\n"); print(as.data.frame(dmp_dose))

# DMP expanded annotation per-region correlation
dmp_expanded_map <- merge(
  dmp[, .(chr, pos, diff, annotation_expanded)],
  dmp_gene_map[, .(chr, pos, gene_id)],
  by = c("chr", "pos"), allow.cartesian = TRUE
)
gene_dmp_expanded <- dmp_expanded_map[, .(n_dmp = .N, mean_diff = mean(diff)),
                                       by = .(gene_id, region = annotation_expanded)]
dmp_expanded_cors <- merge(gene_dmp_expanded, de_results[, .(gene_id, log2FoldChange, padj)],
                           by = "gene_id")[, {
  if (.N >= 20) {
    ct <- tryCatch(cor.test(mean_diff, log2FoldChange, method = "spearman"),
                   error = function(e) list(estimate = NA, p.value = NA))
    .(n_genes = .N, spearman_rho = ct$estimate, spearman_p = ct$p.value,
      pct_de = round(100 * mean(padj < 0.05), 1))
  }
}, by = region][order(spearman_p)]

# DMP location distributions
dmp_loc <- dmp[, .N, by = annotation]
dmp_loc[, pct := round(100 * N / sum(N), 1)]

save_data(gene_dmp, BATCH_DIR, "dmp_gene_expression")
save_data(dmp_region_cors, BATCH_DIR, "dmp_correlation_by_region")
save_data(dmp_expanded_cors, BATCH_DIR, "dmp_correlation_by_expanded_region")
save_data(dmp_de_rate, BATCH_DIR, "dmp_de_rate_by_region")
save_data(dmp_dose, BATCH_DIR, "dmp_dose_response")

# =============================================================================
# 6. DMR ANALYSIS
# =============================================================================
if (has_dmr) {
  cat("\n[6/7] === DMR ANALYSIS ===\n")

  gene_dmr <- dmr_gene_map[, .(n_dmr = .N, mean_dmr_diff = mean(diff.Methy, na.rm = TRUE),
                                 total_cpgs = sum(nCG, na.rm = TRUE)),
                             by = gene_id]
  gene_dmr <- merge(gene_dmr, de_results[, .(gene_id, log2FoldChange, padj, baseMean)],
                    by = "gene_id")
  gene_dmr[, de := padj < 0.05]
  gene_dmr[, meth_direction := ifelse(mean_dmr_diff > 0, "Hyper", "Hypo")]
  gene_dmr <- merge(gene_dmr, gene_annot[, .(gene_id, gene_name)],
                    by = "gene_id", all.x = TRUE)

  cat(sprintf("  Genes with DMR + expression: %d\n", nrow(gene_dmr)))
  cat(sprintf("  Genes with DMR + DE: %d (%.1f%%)\n",
              sum(gene_dmr$de), 100 * mean(gene_dmr$de)))

  # DMR correlation by region
  gene_dmr_region <- dmr_gene_map[, .(n_dmr = .N, mean_dmr_diff = mean(diff.Methy, na.rm = TRUE)),
                                   by = .(gene_id, region = annotation)]
  dmr_region_cors <- merge(gene_dmr_region, de_results[, .(gene_id, log2FoldChange, padj)],
                           by = "gene_id")[, {
    if (.N >= 10) {
      ct <- tryCatch(cor.test(mean_dmr_diff, log2FoldChange, method = "spearman"),
                     error = function(e) list(estimate = NA, p.value = NA))
      .(n_genes = .N, spearman_rho = ct$estimate, spearman_p = ct$p.value,
        pct_de = round(100 * mean(padj < 0.05), 1))
    }
  }, by = region][order(spearman_p)]
  cat("\n  DMR correlation by region:\n"); print(as.data.frame(dmr_region_cors))

  cor_dmr <- cor.test(gene_dmr$mean_dmr_diff, gene_dmr$log2FoldChange, method = "spearman")
  cat(sprintf("  DMR overall: rho = %.4f, p = %.2e\n", cor_dmr$estimate, cor_dmr$p.value))

  # DMR DE rate
  dmr_de_rate <- merge(gene_dmr_region, de_results[, .(gene_id, padj)],
                       by = "gene_id")[, .(n_genes = uniqueN(gene_id),
                                           n_de = uniqueN(gene_id[padj < 0.05])),
                                       by = region]
  dmr_de_rate[, de_rate := round(100 * n_de / n_genes, 1)]

  # DMR expanded
  dmr_mid_dt <- dmr[, .(chr, mid = (start + end) %/% 2L)]
  dmr[, annotation_expanded := annotate_regions_expanded(chr, (start + end) %/% 2L, expanded)]
  dmr_expanded_map <- merge(
    dmr[, .(chr, start, end, diff.Methy, annotation_expanded)],
    dmr_gene_map[, .(chr, start, end, gene_id)],
    by = c("chr", "start", "end"), allow.cartesian = TRUE
  )
  gene_dmr_expanded <- dmr_expanded_map[, .(n_dmr = .N, mean_dmr_diff = mean(diff.Methy, na.rm = TRUE)),
                                         by = .(gene_id, region = annotation_expanded)]
  dmr_expanded_cors <- merge(gene_dmr_expanded, de_results[, .(gene_id, log2FoldChange, padj)],
                             by = "gene_id")[, {
    if (.N >= 10) {
      ct <- tryCatch(cor.test(mean_dmr_diff, log2FoldChange, method = "spearman"),
                     error = function(e) list(estimate = NA, p.value = NA))
      .(n_genes = .N, spearman_rho = ct$estimate, spearman_p = ct$p.value,
        pct_de = round(100 * mean(padj < 0.05), 1))
    }
  }, by = region][order(spearman_p)]

  dmr_loc <- dmr[, .N, by = annotation]
  dmr_loc[, pct := round(100 * N / sum(N), 1)]

  save_data(gene_dmr, BATCH_DIR, "dmr_gene_expression")
  save_data(dmr_region_cors, BATCH_DIR, "dmr_correlation_by_region")
  save_data(dmr_expanded_cors, BATCH_DIR, "dmr_correlation_by_expanded_region")
  save_data(dmr_de_rate, BATCH_DIR, "dmr_de_rate_by_region")
} else {
  cat("\n[6/7] No DMRs — skipping.\n")
}

# =============================================================================
# 7. CONCORDANT GENE LISTS + FIGURES
# =============================================================================
cat("\n[7/7] Concordant gene lists + figures...\n")

# --- Concordant DMP genes: >=3 DMPs AND DE ---
concordant_dmp <- gene_dmp[n_dmp >= 3 & de == TRUE]
concordant_dmp[, display_name := fifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id)]
concordant_dmp[, concordant := (meth_direction == "Hyper" & log2FoldChange > 0) |
                               (meth_direction == "Hypo" & log2FoldChange < 0)]
concordant_dmp[, quadrant := paste0(meth_direction, "_",
                                     ifelse(log2FoldChange > 0, "Up", "Down"))]
cat(sprintf("  Concordant DMP genes (>=3 DMPs + DE): %d\n", nrow(concordant_dmp)))
cat(sprintf("    Hyper+Up: %d | Hyper+Down: %d | Hypo+Up: %d | Hypo+Down: %d\n",
            sum(concordant_dmp$quadrant == "Hyper_Up"),
            sum(concordant_dmp$quadrant == "Hyper_Down"),
            sum(concordant_dmp$quadrant == "Hypo_Up"),
            sum(concordant_dmp$quadrant == "Hypo_Down")))
save_data(concordant_dmp[order(-abs(log2FoldChange))], BATCH_DIR, "concordant_dmp_genes")

# --- Concordant DMR genes: DMR + DE ---
if (has_dmr) {
  concordant_dmr <- gene_dmr[de == TRUE]
  concordant_dmr[, display_name := fifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id)]
  concordant_dmr[, concordant := (meth_direction == "Hyper" & log2FoldChange > 0) |
                                 (meth_direction == "Hypo" & log2FoldChange < 0)]
  concordant_dmr[, quadrant := paste0(meth_direction, "_",
                                       ifelse(log2FoldChange > 0, "Up", "Down"))]
  cat(sprintf("  Concordant DMR genes (DMR + DE): %d\n", nrow(concordant_dmr)))
  cat(sprintf("    Hyper+Up: %d | Hyper+Down: %d | Hypo+Up: %d | Hypo+Down: %d\n",
              sum(concordant_dmr$quadrant == "Hyper_Up"),
              sum(concordant_dmr$quadrant == "Hyper_Down"),
              sum(concordant_dmr$quadrant == "Hypo_Up"),
              sum(concordant_dmr$quadrant == "Hypo_Down")))
  save_data(concordant_dmr[order(-abs(log2FoldChange))], BATCH_DIR, "concordant_dmr_genes")
}

# --- Summary statistics ---
decoupling_stats <- data.frame(
  metric = c("DMP_overall_rho", "DMP_overall_p", "DMP_n_genes",
             "DMP_genes_with_DE", "DMP_pct_DE",
             "DMP_concordant_3plus", "DMP_concordant_3plus_DE"),
  value = c(cor_dmp$estimate, cor_dmp$p.value, nrow(gene_dmp),
            sum(gene_dmp$de), round(100 * mean(gene_dmp$de), 2),
            sum(gene_dmp$n_dmp >= 3), nrow(concordant_dmp))
)
if (has_dmr) {
  decoupling_stats <- rbind(decoupling_stats, data.frame(
    metric = c("DMR_overall_rho", "DMR_overall_p", "DMR_n_genes",
               "DMR_genes_with_DE", "DMR_pct_DE",
               "DMR_concordant", "DMR_concordant_DE"),
    value = c(cor_dmr$estimate, cor_dmr$p.value, nrow(gene_dmr),
              sum(gene_dmr$de), round(100 * mean(gene_dmr$de), 2),
              nrow(gene_dmr), nrow(concordant_dmr))
  ))
}
save_data(decoupling_stats, BATCH_DIR, "decoupling_statistics")

# =============================================================================
# FIGURES
# =============================================================================
cat("  Generating figures...\n")

# Helper: rho-bar plot with labels placed INSIDE bars (large rhos) or OUTSIDE
# (small rhos) so nothing gets clipped at the plot edge regardless of rho size.
plot_rho_bars <- function(dt, title, subtitle, ylab, fname) {
  dt <- copy(dt)
  dt[, lab := sprintf("rho=%.3f\np=%.1e\nn=%d", spearman_rho, spearman_p, n_genes)]
  dt[, lab_hj  := ifelse(spearman_rho >= 0, 1.1, -0.1)]
  dt[, lab_col := ifelse(abs(spearman_rho) >= 0.10, "white", "black")]
  rho_max <- max(abs(dt$spearman_rho), na.rm = TRUE)
  p <- ggplot(dt, aes(x = reorder(region, spearman_rho),
                      y = spearman_rho, fill = spearman_rho)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = lab, hjust = lab_hj, color = lab_col),
              size = 3.2, lineheight = 0.85) +
    scale_color_identity() +
    scale_fill_gradient2(low = "#2471A3", mid = "white", high = "#C0392B",
                         midpoint = 0, guide = "none") +
    scale_y_continuous(limits = c(-rho_max * 1.45, rho_max * 1.45),
                       expand = expansion(mult = 0.02)) +
    coord_flip() +
    labs(title = title, subtitle = subtitle, x = NULL, y = ylab) +
    theme_minimal(base_size = 12)
  save_fig(p, BATCH_DIR, fname, w = 11, h = 6)
}

# --- Fig 7A: DMP correlation by region ---
if (nrow(dmp_region_cors) > 0) {
  plot_rho_bars(
    dmp_region_cors,
    "DMP: methylation-expression correlation by region",
    sprintf("Overall rho = %.4f (p = %.2e) | Gene windows: body + 7kb up + 5kb down",
            cor_dmp$estimate, cor_dmp$p.value),
    "Spearman rho (mean DMP diff vs log2FC)",
    "fig7a_dmp_correlation_by_region")
}

# --- Fig 7B: DMR correlation by region ---
if (has_dmr && exists("dmr_region_cors") && nrow(dmr_region_cors) > 0) {
  plot_rho_bars(
    dmr_region_cors,
    "DMR: methylation-expression correlation by region",
    sprintf("Overall rho = %.4f (p = %.2e)", cor_dmr$estimate, cor_dmr$p.value),
    "Spearman rho (mean DMR diff vs log2FC)",
    "fig7b_dmr_correlation_by_region")
}

# --- Fig 7C/D: DE rate by region (DMP + DMR) ---
p7c <- ggplot(dmp_de_rate, aes(x = reorder(region, de_rate), y = de_rate, fill = region)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", de_rate, n_de, n_genes)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = COLORS$region, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "DMP: DE rate by region", x = NULL, y = "% differentially expressed") +
  theme_minimal(base_size = 12)
save_fig(p7c, BATCH_DIR, "fig7c_dmp_de_rate_by_region", w = 8, h = 6)

if (has_dmr && exists("dmr_de_rate")) {
  p7d <- ggplot(dmr_de_rate, aes(x = reorder(region, de_rate), y = de_rate, fill = region)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", de_rate, n_de, n_genes)),
              vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = COLORS$region, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(title = "DMR: DE rate by region", x = NULL, y = "% differentially expressed") +
    theme_minimal(base_size = 12)
  save_fig(p7d, BATCH_DIR, "fig7d_dmr_de_rate_by_region", w = 8, h = 6)
}

# --- Fig 7E: DMP dose-response ---
p7e <- ggplot(dmp_dose, aes(x = dmp_bin, y = de_rate)) +
  geom_col(fill = "#2471A3", width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", de_rate, de_genes, genes)),
            vjust = -0.3, size = 3) +
  labs(x = "DMPs per gene", y = "% differentially expressed",
       title = "DMP dose-response: more DMPs -> more likely DE?") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_minimal(base_size = 12)
save_fig(p7e, BATCH_DIR, "fig7e_dmp_dose_response", w = 8, h = 6)

# --- Fig 7F/G: Location distribution bars ---
p7f <- ggplot(dmp_loc, aes(x = reorder(annotation, -N), y = N, fill = annotation)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%s\n(%.1f%%)", format(N, big.mark = ","), pct)),
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = COLORS$region, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "DMP genomic distribution", x = NULL, y = "Count") +
  theme_minimal(base_size = 12)
save_fig(p7f, BATCH_DIR, "fig7f_dmp_location", w = 8, h = 6)

if (has_dmr) {
  p7g <- ggplot(dmr_loc, aes(x = reorder(annotation, -N), y = N, fill = annotation)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%s\n(%.1f%%)", format(N, big.mark = ","), pct)),
              vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = COLORS$region, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "DMR genomic distribution", x = NULL, y = "Count") +
    theme_minimal(base_size = 12)
  save_fig(p7g, BATCH_DIR, "fig7g_dmr_location", w = 8, h = 6)
}

# --- Fig 7H/I: Expanded correlation by region ---
if (nrow(dmp_expanded_cors) > 0) {
  p7h <- ggplot(dmp_expanded_cors, aes(x = reorder(region, spearman_rho),
                                        y = spearman_rho, fill = spearman_rho)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("rho=%.3f\nn=%d", spearman_rho, n_genes)),
              hjust = ifelse(dmp_expanded_cors$spearman_rho > 0, -0.1, 1.1), size = 2.8) +
    scale_fill_gradient2(low = "#2471A3", mid = "white", high = "#C0392B",
                         midpoint = 0, guide = "none") +
    coord_flip() +
    labs(title = "DMP: meth-expression correlation by expanded region",
         x = NULL, y = "Spearman rho") +
    theme_minimal(base_size = 12)
  save_fig(p7h, BATCH_DIR, "fig7h_dmp_expanded_correlation", w = 10, h = 8)
}

if (has_dmr && exists("dmr_expanded_cors") && nrow(dmr_expanded_cors) > 0) {
  p7i <- ggplot(dmr_expanded_cors, aes(x = reorder(region, spearman_rho),
                                        y = spearman_rho, fill = spearman_rho)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("rho=%.3f\nn=%d", spearman_rho, n_genes)),
              hjust = ifelse(dmr_expanded_cors$spearman_rho > 0, -0.1, 1.1), size = 2.8) +
    scale_fill_gradient2(low = "#2471A3", mid = "white", high = "#C0392B",
                         midpoint = 0, guide = "none") +
    coord_flip() +
    labs(title = "DMR: meth-expression correlation by expanded region",
         x = NULL, y = "Spearman rho") +
    theme_minimal(base_size = 12)
  save_fig(p7i, BATCH_DIR, "fig7i_dmr_expanded_correlation", w = 10, h = 8)
}

# --- Fig 7J: DMP methylation vs expression scatter ---
cat("  Fig 7J: DMP meth vs expression scatter...\n")
gene_dmp[, label_col := fifelse(de, "DE (padj < 0.05)", "Not DE")]
gene_dmp[, display_name := fifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id)]
top_dmp_label <- gene_dmp[de == TRUE][order(-abs(mean_diff))][1:min(20, sum(gene_dmp$de))]

p7j <- ggplot(gene_dmp, aes(x = mean_diff, y = log2FoldChange)) +
  geom_point(aes(color = label_col), alpha = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_text_repel(data = top_dmp_label, aes(label = display_name),
                  size = 2.8, max.overlaps = 30, segment.alpha = 0.5,
                  min.segment.length = 0) +
  scale_color_manual(values = c("DE (padj < 0.05)" = "#E74C3C", "Not DE" = "#AEB6BF"),
                     name = NULL) +
  labs(title = "Differential methylation vs expression change (DMPs)",
       subtitle = sprintf("Spearman rho = %.4f (p = %.2e) | %d genes (window: body + 7kb up + 5kb down)",
                          cor_dmp$estimate, cor_dmp$p.value, nrow(gene_dmp)),
       x = "Mean methylation difference (DMP)", y = "log2 Fold Change (amputated / control)") +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
save_fig(p7j, BATCH_DIR, "fig7j_dmp_meth_vs_expression_scatter", w = 10, h = 8)

# --- Fig 7K: DMR methylation vs expression scatter ---
if (has_dmr && exists("gene_dmr") && nrow(gene_dmr) > 0) {
  cat("  Fig 7K: DMR meth vs expression scatter...\n")
  gene_dmr[, label_col := fifelse(de, "DE (padj < 0.05)", "Not DE")]
  gene_dmr[, display_name := fifelse(!is.na(gene_name) & gene_name != "", gene_name, gene_id)]
  top_dmr_label <- gene_dmr[de == TRUE][order(-abs(mean_dmr_diff))][1:min(20, sum(gene_dmr$de))]

  p7k <- ggplot(gene_dmr, aes(x = mean_dmr_diff, y = log2FoldChange)) +
    geom_point(aes(color = label_col), alpha = 0.4, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_text_repel(data = top_dmr_label, aes(label = display_name),
                    size = 2.8, max.overlaps = 30, segment.alpha = 0.5,
                    min.segment.length = 0) +
    scale_color_manual(values = c("DE (padj < 0.05)" = "#E74C3C", "Not DE" = "#AEB6BF"),
                       name = NULL) +
    labs(title = "Differential methylation vs expression change (DMRs)",
         subtitle = sprintf("Spearman rho = %.4f (p = %.2e) | %d genes",
                            cor_dmr$estimate, cor_dmr$p.value, nrow(gene_dmr)),
         x = "Mean methylation difference (DMR)", y = "log2 Fold Change (amputated / control)") +
    theme_minimal(base_size = 12) + theme(legend.position = "top")
  save_fig(p7k, BATCH_DIR, "fig7k_dmr_meth_vs_expression_scatter", w = 10, h = 8)
}

# =============================================================================
# DONE
# =============================================================================
elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("\n=== Batch 07 complete (%.1f min) ===\n", elapsed / 60))
n_figs <- length(list.files(file.path(BATCH_DIR, "figures"), pattern = "\\.png$"))
n_data <- length(list.files(file.path(BATCH_DIR, "data"), pattern = "\\.tsv$"))
cat(sprintf("Figures: %d | Data files: %d\n", n_figs, n_data))
