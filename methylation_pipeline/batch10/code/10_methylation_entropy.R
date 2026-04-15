#!/usr/bin/env Rscript
# =============================================================================
# Batch 10: Methylation Entropy Analysis
# Question: Does regeneration increase or decrease methylation entropy?
# Output: data/ (entropy stats) + figures/ (8 plots: A-H)
# Requires: BSseq cache, batch06 DMPs, batch07 DMP-expression
# =============================================================================

source("methylation_pipeline/_config.R")
t0 <- proc.time()

library(bsseq)
library(GenomicRanges)
library(data.table)
library(ggplot2)

BATCH_DIR <- file.path(PIPE_DIR, "batch10")
cat("=== Batch 10: Methylation Entropy ===\n\n")

# Clean old output
unlink(file.path(BATCH_DIR, "figures"), recursive = TRUE)
unlink(file.path(BATCH_DIR, "data"), recursive = TRUE)

# =============================================================================
# 1. Load BSseq
# =============================================================================
cat("Loading BSseq object...\n")
bs_obj <- readRDS(CACHE$bsseq)
cat(sprintf("BSseq: %s sites\n", format(nrow(bs_obj), big.mark = ",")))

# =============================================================================
# 2. Compute per-CpG binary Shannon entropy
# =============================================================================
cat("Computing per-CpG entropy...\n")

binary_entropy <- function(p) {
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  -p * log2(p) - (1 - p) * log2(1 - p)
}

ctrl_beta <- rowMeans(getMeth(bs_obj[, c("C1", "C2")], type = "raw"), na.rm = TRUE)
ampu_beta <- rowMeans(getMeth(bs_obj[, c("A1", "A2")], type = "raw"), na.rm = TRUE)

ctrl_entropy <- binary_entropy(ctrl_beta)
ampu_entropy <- binary_entropy(ampu_beta)

valid <- !is.na(ctrl_entropy) & !is.na(ampu_entropy) &
         !is.nan(ctrl_entropy) & !is.nan(ampu_entropy)

entropy_dt <- data.table(
  chr           = as.character(seqnames(bs_obj))[valid],
  pos           = start(bs_obj)[valid],
  ctrl_beta     = ctrl_beta[valid],
  ampu_beta     = ampu_beta[valid],
  ctrl_entropy  = ctrl_entropy[valid],
  ampu_entropy  = ampu_entropy[valid],
  delta_entropy = ampu_entropy[valid] - ctrl_entropy[valid]
)

cat(sprintf("Valid CpGs for entropy: %s\n", format(nrow(entropy_dt), big.mark = ",")))
cat(sprintf("Mean entropy - Control: %.4f, Amputated: %.4f\n",
            mean(entropy_dt$ctrl_entropy), mean(entropy_dt$ampu_entropy)))
cat(sprintf("Mean delta entropy: %.6f\n", mean(entropy_dt$delta_entropy)))

# Wilcoxon test — all data, no subsampling
wt <- wilcox.test(entropy_dt$ctrl_entropy, entropy_dt$ampu_entropy, paired = TRUE)
cat(sprintf("Wilcoxon paired test (n=%s): p = %s\n",
            format(nrow(entropy_dt), big.mark = ","), format(wt$p.value, digits = 3)))

save_data(data.table(
  mean_ctrl_entropy = mean(entropy_dt$ctrl_entropy),
  mean_ampu_entropy = mean(entropy_dt$ampu_entropy),
  mean_delta        = mean(entropy_dt$delta_entropy),
  wilcox_p          = wt$p.value,
  n_cpgs            = nrow(entropy_dt)
), BATCH_DIR, "entropy_summary")

# Free BSseq
rm(bs_obj, ctrl_beta, ampu_beta, ctrl_entropy, ampu_entropy, valid)
gc(verbose = FALSE)

# =============================================================================
# 3. Annotate by genomic region
# =============================================================================
cat("\nAnnotating CpGs by genomic region...\n")
gff <- load_gff()
genes <- gff[gff$type == "gene"]
exons <- gff[gff$type == "exon"]
promoters_gr <- if (file.exists(CACHE$promoters)) {
  readRDS(CACHE$promoters)
} else {
  trim(promoters(genes, upstream = 2000, downstream = 0))
}

entropy_dt[, region := annotate_regions(chr, pos, promoters_gr, exons, genes)]

region_entropy <- entropy_dt[, .(
  mean_ctrl  = mean(ctrl_entropy),
  mean_ampu  = mean(ampu_entropy),
  mean_delta = mean(delta_entropy),
  n          = .N
), by = region]
cat("\nEntropy by region:\n")
print(region_entropy)
save_data(region_entropy, BATCH_DIR, "entropy_by_region")

# =============================================================================
# 4. Load DMPs from batch06
# =============================================================================
dmp_path <- file.path(PIPE_DIR, "batch06/data/dmps_annotated.tsv")
has_dmps <- file.exists(dmp_path)
if (has_dmps) {
  cat("\nLoading DMPs from batch06...\n")
  dmps <- fread(dmp_path)
  cat(sprintf("DMPs loaded: %s\n", format(nrow(dmps), big.mark = ",")))

  # Merge DMPs with entropy data
  setkey(entropy_dt, chr, pos)
  setkey(dmps, chr, pos)
  dmp_entropy <- entropy_dt[dmps, nomatch = 0L]
  cat(sprintf("DMPs matched to entropy: %s\n", format(nrow(dmp_entropy), big.mark = ",")))

  # Mark DMP status in entropy_dt
  entropy_dt[, is_dmp := FALSE]
  entropy_dt[dmp_entropy[, .(chr, pos)], is_dmp := TRUE, on = .(chr, pos)]
} else {
  cat("\nWARNING: batch06 DMPs not found at", dmp_path, "\n")
  cat("Figures C-F will be skipped.\n")
}

# =============================================================================
# 5. Beta-matched non-DMP sampling (for figures C, D, E)
# =============================================================================
if (has_dmps && nrow(dmp_entropy) > 0) {
  cat("\nMatching non-DMPs by control beta...\n")
  set.seed(42)

  # Bin all CpGs by control beta (width 0.05)
  entropy_dt[, beta_bin := floor(ctrl_beta / 0.05)]
  dmp_entropy[, beta_bin := floor(ctrl_beta / 0.05)]

  # Non-DMP pool
  non_dmp_pool <- entropy_dt[is_dmp == FALSE]

  # For each DMP, sample one non-DMP from same beta bin
  dmp_bins <- dmp_entropy[, .(n_needed = .N), by = beta_bin]
  matched_list <- vector("list", nrow(dmp_bins))

  for (i in seq_len(nrow(dmp_bins))) {
    bb <- dmp_bins$beta_bin[i]
    nn <- dmp_bins$n_needed[i]
    pool_i <- non_dmp_pool[beta_bin == bb]
    if (nrow(pool_i) == 0) next
    samp_idx <- sample(nrow(pool_i), min(nn, nrow(pool_i)), replace = nrow(pool_i) < nn)
    matched_list[[i]] <- pool_i[samp_idx]
  }
  matched_nondmp <- rbindlist(matched_list[!sapply(matched_list, is.null)])
  cat(sprintf("Matched non-DMPs: %s\n", format(nrow(matched_nondmp), big.mark = ",")))

  # Clean up beta_bin columns
  entropy_dt[, beta_bin := NULL]
  dmp_entropy[, beta_bin := NULL]
}

# =============================================================================
# FIGURES
# =============================================================================
cat("\nGenerating figures...\n")

# --- (A) Entropy distribution: control vs amputated ---
set.seed(42)
ent_long <- rbindlist(list(
  data.table(condition = "Control",   entropy = entropy_dt$ctrl_entropy),
  data.table(condition = "Amputated", entropy = entropy_dt$ampu_entropy)
))

p10a <- ggplot(ent_long, aes(x = entropy, fill = condition, color = condition)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  scale_fill_manual(values = COLORS$condition) +
  scale_color_manual(values = COLORS$condition) +
  labs(x = "Binary Shannon entropy", y = "Density",
       title = "CpG methylation entropy: control vs amputated",
       subtitle = sprintf("Delta = %.5f | Wilcoxon p = %s",
                          mean(entropy_dt$delta_entropy),
                          format(wt$p.value, digits = 3))) +
  theme_minimal(base_size = 12)
save_fig(p10a, BATCH_DIR, "fig10a_entropy_distribution", w = 9, h = 6)
rm(ent_long); gc(verbose = FALSE)

# --- (B) Per-region entropy change ---
p10b <- ggplot(region_entropy, aes(x = region, y = mean_delta, fill = region)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = COLORS$region, guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "Mean delta entropy (Amputated - Control)",
       title = "Entropy change by genomic region") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p10b, BATCH_DIR, "fig10b_entropy_by_region", w = 7, h = 6)

# --- (C) DMP baseline entropy vs beta-matched non-DMPs ---
if (has_dmps && nrow(dmp_entropy) > 0 && exists("matched_nondmp")) {
  cat("\nFig 10C: DMP baseline entropy vs beta-matched non-DMPs\n")

  plot_c_dt <- rbindlist(list(
    dmp_entropy[direction == "Hyper", .(group = "DMP (Hyper)", baseline_entropy = ctrl_entropy)],
    dmp_entropy[direction == "Hypo",  .(group = "DMP (Hypo)",  baseline_entropy = ctrl_entropy)],
    matched_nondmp[, .(group = "Matched non-DMP", baseline_entropy = ctrl_entropy)]
  ))
  plot_c_dt[, group := factor(group, levels = c("DMP (Hyper)", "DMP (Hypo)", "Matched non-DMP"))]

  # Wilcoxon: all DMPs vs matched
  all_dmp_entropy <- dmp_entropy$ctrl_entropy
  wt_c <- wilcox.test(all_dmp_entropy, matched_nondmp$ctrl_entropy)
  cohens_d_c <- (mean(all_dmp_entropy) - mean(matched_nondmp$ctrl_entropy)) /
    sqrt((var(all_dmp_entropy) + var(matched_nondmp$ctrl_entropy)) / 2)

  colors_c <- c("DMP (Hyper)" = unname(COLORS$direction["Hyper"]),
                "DMP (Hypo)"  = unname(COLORS$direction["Hypo"]),
                "Matched non-DMP" = "gray60")

  p10c <- ggplot(plot_c_dt, aes(x = group, y = baseline_entropy, fill = group)) +
    geom_violin(alpha = 0.6, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
    scale_fill_manual(values = colors_c, guide = "none") +
    labs(x = NULL, y = "Baseline (control) binary entropy",
         title = "Baseline entropy: DMPs vs beta-matched non-DMPs",
         subtitle = sprintf("Wilcoxon p = %s | Cohen's d = %.3f",
                            format(wt_c$p.value, digits = 3), cohens_d_c)) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  save_fig(p10c, BATCH_DIR, "fig10c_dmp_vs_matched_entropy", w = 8, h = 6)

  save_data(data.table(
    test = "DMP_vs_matched_baseline_entropy",
    wilcox_p = wt_c$p.value,
    cohens_d = cohens_d_c,
    mean_dmp = mean(all_dmp_entropy),
    mean_matched = mean(matched_nondmp$ctrl_entropy),
    n_dmp = length(all_dmp_entropy),
    n_matched = nrow(matched_nondmp)
  ), BATCH_DIR, "dmp_vs_matched_entropy_test")

  rm(plot_c_dt); gc(verbose = FALSE)
} else {
  cat("Skipping fig10c (no DMPs)\n")
}

# --- (D) Entropy change (delta) at DMPs vs non-DMPs ---
if (has_dmps && nrow(dmp_entropy) > 0 && exists("matched_nondmp")) {
  cat("\nFig 10D: Delta entropy at DMPs vs matched non-DMPs\n")

  plot_d_dt <- rbindlist(list(
    dmp_entropy[, .(group = "DMP sites", delta = delta_entropy)],
    matched_nondmp[, .(group = "Matched non-DMP", delta = delta_entropy)]
  ))
  plot_d_dt[, group := factor(group, levels = c("DMP sites", "Matched non-DMP"))]

  wt_d <- wilcox.test(dmp_entropy$delta_entropy, matched_nondmp$delta_entropy)

  p10d <- ggplot(plot_d_dt, aes(x = group, y = delta, fill = group)) +
    geom_violin(alpha = 0.6, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = c("DMP sites" = "#C0392B", "Matched non-DMP" = "gray60"),
                      guide = "none") +
    labs(x = NULL, y = "Delta entropy (Amputated - Control)",
         title = "Entropy change at DMP sites vs matched non-DMPs",
         subtitle = sprintf("Wilcoxon p = %s | DMP mean delta = %.5f | Matched mean delta = %.5f",
                            format(wt_d$p.value, digits = 3),
                            mean(dmp_entropy$delta_entropy),
                            mean(matched_nondmp$delta_entropy))) +
    theme_minimal(base_size = 12)
  save_fig(p10d, BATCH_DIR, "fig10d_delta_entropy_dmp_vs_nondmp", w = 7, h = 6)

  rm(plot_d_dt); gc(verbose = FALSE)
} else {
  cat("Skipping fig10d (no DMPs)\n")
}

# --- (E) Entropy by DMP direction ---
if (has_dmps && nrow(dmp_entropy) > 0 && exists("matched_nondmp")) {
  cat("\nFig 10E: Baseline entropy by DMP direction\n")

  plot_e_dt <- rbindlist(list(
    dmp_entropy[direction == "Hyper", .(group = "Hyper DMPs", baseline_entropy = ctrl_entropy)],
    dmp_entropy[direction == "Hypo",  .(group = "Hypo DMPs",  baseline_entropy = ctrl_entropy)],
    matched_nondmp[, .(group = "Matched non-DMP", baseline_entropy = ctrl_entropy)]
  ))
  plot_e_dt[, group := factor(group, levels = c("Hyper DMPs", "Hypo DMPs", "Matched non-DMP"))]

  # Pairwise Wilcoxon tests
  wt_hyper_hypo <- wilcox.test(
    dmp_entropy[direction == "Hyper"]$ctrl_entropy,
    dmp_entropy[direction == "Hypo"]$ctrl_entropy
  )

  colors_e <- c("Hyper DMPs" = unname(COLORS$direction["Hyper"]),
                "Hypo DMPs"  = unname(COLORS$direction["Hypo"]),
                "Matched non-DMP" = "gray60")

  p10e <- ggplot(plot_e_dt, aes(x = group, y = baseline_entropy, fill = group)) +
    geom_violin(alpha = 0.6, scale = "width") +
    geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.7) +
    scale_fill_manual(values = colors_e, guide = "none") +
    labs(x = NULL, y = "Baseline (control) binary entropy",
         title = "Baseline entropy by DMP direction",
         subtitle = sprintf("Hyper vs Hypo: Wilcoxon p = %s | Hyper median = %.4f | Hypo median = %.4f",
                            format(wt_hyper_hypo$p.value, digits = 3),
                            median(dmp_entropy[direction == "Hyper"]$ctrl_entropy),
                            median(dmp_entropy[direction == "Hypo"]$ctrl_entropy))) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  save_fig(p10e, BATCH_DIR, "fig10e_entropy_by_dmp_direction", w = 8, h = 6)

  rm(plot_e_dt); gc(verbose = FALSE)
} else {
  cat("Skipping fig10e (no DMPs)\n")
}

# --- (F) Regional entropy at DMP sites vs genome-wide ---
if (has_dmps && nrow(dmp_entropy) > 0) {
  cat("\nFig 10F: Regional delta entropy at DMP sites vs genome-wide\n")

  # DMP delta by region (use annotation from dmp_entropy which inherited from entropy_dt)
  dmp_entropy[, region := annotate_regions(chr, pos, promoters_gr, exons, genes)]

  dmp_region <- dmp_entropy[, .(mean_delta = mean(delta_entropy), source = "DMP sites"), by = region]
  bg_region  <- region_entropy[, .(region, mean_delta, source = "Genome-wide")]

  plot_f_dt <- rbindlist(list(dmp_region, bg_region))
  plot_f_dt[, region := factor(region, levels = c("Promoter", "Exon", "Intron", "Intergenic"))]
  plot_f_dt[, source := factor(source, levels = c("DMP sites", "Genome-wide"))]

  p10f <- ggplot(plot_f_dt, aes(x = region, y = mean_delta, fill = source)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = c("DMP sites" = "#C0392B", "Genome-wide" = "gray60")) +
    labs(x = NULL, y = "Mean delta entropy (Amputated - Control)",
         title = "Regional entropy change: DMP sites vs genome-wide",
         fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top")
  save_fig(p10f, BATCH_DIR, "fig10f_regional_entropy_dmp_vs_genomewide", w = 8, h = 6)

  save_data(plot_f_dt, BATCH_DIR, "regional_entropy_dmp_vs_genomewide")
} else {
  cat("Skipping fig10f (no DMPs)\n")
}

# --- (G) Mutual Information: methylation change vs expression change ---
mi_path <- file.path(PIPE_DIR, "batch07/data/dmp_gene_expression.tsv")
if (file.exists(mi_path)) {
  cat("\nFig 10G: Mutual Information (methylation vs expression)\n")

  mi_dt <- fread(mi_path)
  mi_dt <- mi_dt[!is.na(mean_diff) & !is.na(log2FoldChange)]
  cat(sprintf("Genes with both meth diff and LFC: %d\n", nrow(mi_dt)))

  if (nrow(mi_dt) >= 50) {
    # Discretize into 10 equal-frequency bins
    mi_dt[, meth_bin := as.integer(cut(mean_diff, breaks = quantile(mean_diff, probs = seq(0, 1, 0.1)),
                                       include.lowest = TRUE, labels = FALSE))]
    mi_dt[, expr_bin := as.integer(cut(log2FoldChange, breaks = quantile(log2FoldChange, probs = seq(0, 1, 0.1)),
                                       include.lowest = TRUE, labels = FALSE))]

    # Compute MI
    compute_mi <- function(x_bin, y_bin) {
      joint <- table(x_bin, y_bin)
      p_xy <- joint / sum(joint)
      p_x <- rowSums(p_xy)
      p_y <- colSums(p_xy)
      mi <- 0
      for (i in seq_along(p_x)) {
        for (j in seq_along(p_y)) {
          if (p_xy[i, j] > 0 && p_x[i] > 0 && p_y[j] > 0) {
            mi <- mi + p_xy[i, j] * log2(p_xy[i, j] / (p_x[i] * p_y[j]))
          }
        }
      }
      mi
    }

    observed_mi <- compute_mi(mi_dt$meth_bin, mi_dt$expr_bin)
    cat(sprintf("Observed MI: %.6f bits\n", observed_mi))

    # Permutation test (1000 shuffles)
    set.seed(42)
    n_perm <- 1000L
    perm_mi <- vapply(seq_len(n_perm), function(i) {
      compute_mi(mi_dt$meth_bin, sample(mi_dt$expr_bin))
    }, numeric(1))

    perm_p <- mean(perm_mi >= observed_mi)
    cat(sprintf("Permutation p-value: %.4f\n", perm_p))

    perm_dt <- data.table(perm_mi = perm_mi)

    p10g <- ggplot(perm_dt, aes(x = perm_mi)) +
      geom_histogram(bins = 40, fill = "gray70", color = "gray40", alpha = 0.8) +
      geom_vline(xintercept = observed_mi, color = "red", linewidth = 1.2, linetype = "solid") +
      labs(x = "Mutual Information (bits)", y = "Count (permutations)",
           title = "Mutual Information: methylation change vs expression change",
           subtitle = sprintf("MI = %.6f bits, permutation p = %.3f (n = %d genes)",
                              observed_mi, perm_p, nrow(mi_dt))) +
      theme_minimal(base_size = 12)
    save_fig(p10g, BATCH_DIR, "fig10g_mutual_information", w = 8, h = 6)

    save_data(data.table(
      observed_mi = observed_mi,
      perm_p      = perm_p,
      perm_mean   = mean(perm_mi),
      perm_sd     = sd(perm_mi),
      n_genes     = nrow(mi_dt)
    ), BATCH_DIR, "mutual_information_test")
  } else {
    cat("Too few genes for MI analysis (need >= 50).\n")
  }
} else {
  cat("\nSkipping fig10g: batch07 data not found at", mi_path, "\n")
}

# --- (H) Entropy vs expression variability (Fang framework) ---
cat("\nFig 10H: Entropy vs expression variability\n")
tryCatch({
  # Try loading transcriptome cache
  has_transcriptome <- FALSE
  norm_counts <- NULL

  if (file.exists(CACHE$transcriptome)) {
    cat("Loading CACHE$transcriptome...\n")
    txome <- readRDS(CACHE$transcriptome)

    # The transcriptome cache is a list; extract normalized counts for tail controls
    # Try common structures
    if (is.list(txome) && "counts_norm" %in% names(txome)) {
      norm_counts <- txome$counts_norm
      has_transcriptome <- TRUE
    } else if (is.list(txome) && "dds" %in% names(txome)) {
      # DESeq2 object inside
      norm_counts <- DESeq2::counts(txome$dds, normalized = TRUE)
      has_transcriptome <- TRUE
    } else if (is(txome, "DESeqDataSet")) {
      norm_counts <- DESeq2::counts(txome, normalized = TRUE)
      has_transcriptome <- TRUE
    } else if (is.list(txome) && "vsd" %in% names(txome)) {
      norm_counts <- SummarizedExperiment::assay(txome$vsd)
      has_transcriptome <- TRUE
    } else if (is.list(txome)) {
      # Try to find anything with counts in the name
      count_names <- grep("count|norm|expr", names(txome), value = TRUE, ignore.case = TRUE)
      if (length(count_names) > 0) {
        candidate <- txome[[count_names[1]]]
        if (is.matrix(candidate) || is.data.frame(candidate)) {
          norm_counts <- as.matrix(candidate)
          has_transcriptome <- TRUE
        }
      }
    }
    rm(txome); gc(verbose = FALSE)
  }

  if (!has_transcriptome) {
    # Fallback: try loading raw counts from counts_dir
    cat("Transcriptome cache not usable. Trying raw counts...\n")
    count_files <- list.files(OG$counts_dir, pattern = "tail.*control|T.*C.*\\.txt",
                              full.names = TRUE, ignore.case = TRUE)
    if (length(count_files) == 0) {
      count_files <- list.files(OG$counts_dir, pattern = "\\.txt$", full.names = TRUE)
    }
    if (length(count_files) >= 2) {
      count_list <- lapply(count_files, function(f) {
        dt <- fread(f, header = FALSE, col.names = c("gene_id", "count"))
        dt <- dt[!grepl("^__", gene_id)]
        dt
      })
      # Merge all by gene_id
      merged <- count_list[[1]]
      setnames(merged, "count", basename(count_files[1]))
      for (i in 2:length(count_list)) {
        setnames(count_list[[i]], "count", basename(count_files[i]))
        merged <- merge(merged, count_list[[i]], by = "gene_id", all = FALSE)
      }
      gene_ids <- merged$gene_id
      count_mat <- as.matrix(merged[, -1, with = FALSE])
      rownames(count_mat) <- gene_ids
      # Simple size-factor normalization
      lib_sizes <- colSums(count_mat)
      norm_counts <- sweep(count_mat, 2, lib_sizes / median(lib_sizes), "/")
      has_transcriptome <- TRUE
      cat(sprintf("Loaded %d count files, %d genes\n", length(count_files), nrow(norm_counts)))
    }
  }

  if (has_transcriptome && !is.null(norm_counts) && ncol(norm_counts) >= 2) {
    # Identify tail control columns (heuristic: columns with "C" or "control" or "tail")
    ctrl_cols <- grep("tail.*control|T.*C[0-9]|^C[0-9]", colnames(norm_counts),
                      value = TRUE, ignore.case = TRUE)
    if (length(ctrl_cols) < 2) {
      # Use all columns as fallback
      ctrl_cols <- colnames(norm_counts)
    }
    cat(sprintf("Using %d columns for CV: %s\n", length(ctrl_cols),
                paste(head(ctrl_cols, 5), collapse = ", ")))

    if (length(ctrl_cols) >= 2) {
      ctrl_mat <- norm_counts[, ctrl_cols, drop = FALSE]
      gene_mean <- rowMeans(ctrl_mat, na.rm = TRUE)
      gene_sd   <- apply(ctrl_mat, 1, sd, na.rm = TRUE)
      gene_cv   <- gene_sd / gene_mean
      gene_cv[!is.finite(gene_cv)] <- NA

      expr_cv_dt <- data.table(
        gene_id = rownames(ctrl_mat),
        mean_expr = gene_mean,
        cv = gene_cv
      )
      expr_cv_dt <- expr_cv_dt[!is.na(cv) & is.finite(cv) & mean_expr > 1]

      # Compute mean gene-body entropy per gene
      # Annotate CpGs to genes using findOverlaps
      cpg_gr <- GRanges(seqnames = entropy_dt$chr,
                         ranges = IRanges(start = entropy_dt$pos, width = 1))
      hits <- findOverlaps(cpg_gr, genes)
      gene_body_dt <- data.table(
        cpg_idx = queryHits(hits),
        gene_id = genes$ID[subjectHits(hits)]
      )
      gene_body_dt[, ctrl_entropy := entropy_dt$ctrl_entropy[cpg_idx]]
      gene_entropy <- gene_body_dt[, .(mean_gb_entropy = mean(ctrl_entropy), n_cpgs = .N),
                                   by = gene_id]

      # Merge
      merged_h <- merge(gene_entropy, expr_cv_dt, by = "gene_id")
      merged_h <- merged_h[n_cpgs >= 5]  # require at least 5 CpGs
      cat(sprintf("Genes for entropy-CV plot: %d\n", nrow(merged_h)))

      if (nrow(merged_h) >= 30) {
        rho_h <- cor.test(merged_h$mean_gb_entropy, merged_h$cv, method = "spearman")

        # NEVER subsample — use all genes. Use geom_hex for visual density instead.
        plot_h <- merged_h

        p10h <- ggplot(plot_h, aes(x = mean_gb_entropy, y = cv)) +
          geom_hex(bins = 80) +
          scale_fill_viridis_c(trans = "log10", name = "Genes (log10)") +
          geom_smooth(method = "loess", se = TRUE, color = "#C0392B", linewidth = 1) +
          labs(x = "Mean gene body entropy (control)",
               y = "Expression CV (control samples)",
               title = "Entropy vs expression variability (Fang framework)",
               subtitle = sprintf("Spearman rho = %.3f, p = %s (n = %d genes)",
                                  rho_h$estimate, format(rho_h$p.value, digits = 3),
                                  nrow(merged_h))) +
          theme_minimal(base_size = 12)
        save_fig(p10h, BATCH_DIR, "fig10h_entropy_vs_expression_cv", w = 8, h = 7)

        save_data(data.table(
          spearman_rho = as.numeric(rho_h$estimate),
          spearman_p   = rho_h$p.value,
          n_genes      = nrow(merged_h)
        ), BATCH_DIR, "entropy_vs_expression_cv")
      } else {
        cat("Too few genes for entropy-CV plot.\n")
      }

      rm(cpg_gr, hits, gene_body_dt, gene_entropy); gc(verbose = FALSE)
    }
  } else {
    cat("Could not load expression data for fig10h.\n")
  }
}, error = function(e) {
  cat("Error in fig10h (skipping gracefully):", conditionMessage(e), "\n")
})

# =============================================================================
# 6. FANG FRAMEWORK: Sliding-window entropy-gaining regions (1kb windows)
# =============================================================================
cat("\n[Fang] Sliding-window entropy analysis (1kb)...\n")

# Tile genome in 1kb windows, assign CpGs
chr_lens <- entropy_dt[, .(max_pos = max(pos)), by = chr]
window_results <- lapply(seq_len(nrow(chr_lens)), function(i) {
  this_chr <- chr_lens$chr[i]
  chr_cpgs <- entropy_dt[chr == this_chr]
  chr_cpgs[, win_id := paste0(this_chr, ":", floor(pos / 1000L) * 1000L)]
  # Per-window stats: require >= 5 CpGs
  chr_cpgs[, .(
    chr        = this_chr,
    win_start  = floor(pos[1] / 1000L) * 1000L,
    n_cpgs     = .N,
    mean_delta = mean(delta_entropy),
    sd_delta   = sd(delta_entropy),
    mean_ctrl_entropy = mean(ctrl_entropy),
    mean_ctrl_beta    = mean(ctrl_beta)
  ), by = win_id][n_cpgs >= 5]
})
win_dt <- rbindlist(window_results)
cat(sprintf("  Windows with >=5 CpGs: %s\n", format(nrow(win_dt), big.mark = ",")))

# One-sample Wilcoxon per window is too slow for millions of windows.
# Use z-test: z = mean_delta / (sd_delta / sqrt(n))
win_dt[, z_score := mean_delta / (sd_delta / sqrt(n_cpgs))]
win_dt[, p_value := 2 * pnorm(-abs(z_score))]
win_dt[, fdr := p.adjust(p_value, method = "BH")]

# Entropy-gaining windows (significant + positive delta)
gain_windows <- win_dt[fdr < 0.05 & mean_delta > 0]
loss_windows <- win_dt[fdr < 0.05 & mean_delta < 0]
cat(sprintf("  Significant entropy-GAINING windows: %s\n", format(nrow(gain_windows), big.mark = ",")))
cat(sprintf("  Significant entropy-LOSING windows:  %s\n", format(nrow(loss_windows), big.mark = ",")))

# Annotate gaining AND losing windows by region
if (nrow(gain_windows) > 0) {
  gain_gr <- GRanges(seqnames = gain_windows$chr,
                     ranges = IRanges(start = gain_windows$win_start, width = 1000L))
  gain_windows[, region := annotate_regions(chr, win_start + 500L, promoters_gr, exons, genes)]

  # Nearest gene for gaining windows
  gain_nearest <- nearest(gain_gr, genes)
  gain_windows[, nearest_gene := ifelse(!is.na(gain_nearest), genes$ID[gain_nearest], NA_character_)]
}
if (nrow(loss_windows) > 0) {
  loss_windows[, region := annotate_regions(chr, win_start + 500L, promoters_gr, exons, genes)]
}

save_data(win_dt[fdr < 0.05], BATCH_DIR, "significant_entropy_windows")
save_data(gain_windows, BATCH_DIR, "entropy_gaining_windows")

# Fig 10I: Distribution of window-level entropy change
p10i <- ggplot(win_dt, aes(x = mean_delta)) +
  geom_histogram(bins = 200, fill = "gray60", color = "gray40", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  annotate("text", x = max(win_dt$mean_delta) * 0.7, y = Inf, vjust = 2,
           label = sprintf("Gaining: %s windows\nLosing: %s windows",
                           format(nrow(gain_windows), big.mark = ","),
                           format(nrow(loss_windows), big.mark = ",")),
           hjust = 1, size = 4) +
  labs(x = "Mean delta entropy per 1kb window",
       y = "Number of windows",
       title = "Window-level entropy change distribution (1kb, >=5 CpGs)",
       subtitle = sprintf("FDR < 0.05: %s gaining, %s losing (of %s tested)",
                          format(nrow(gain_windows), big.mark = ","),
                          format(nrow(loss_windows), big.mark = ","),
                          format(nrow(win_dt), big.mark = ","))) +
  theme_minimal(base_size = 12)
save_fig(p10i, BATCH_DIR, "fig10i_window_entropy_distribution", w = 9, h = 6)

# Fig 10J: Regional breakdown of entropy-gaining windows
if (nrow(gain_windows) > 0) {
  gain_region_counts <- gain_windows[, .N, by = region]
  loss_region_counts <- loss_windows[, .N, by = region]
  gain_region_counts[, direction := "Gaining"]
  loss_region_counts[, direction := "Losing"]
  region_change <- rbindlist(list(gain_region_counts, loss_region_counts))

  p10j <- ggplot(region_change, aes(x = region, y = N, fill = direction)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    scale_fill_manual(values = c("Gaining" = "#E74C3C", "Losing" = "#3498DB")) +
    labs(x = NULL, y = "Number of significant 1kb windows",
         title = "Entropy-gaining vs -losing windows by genomic region",
         fill = "Direction") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p10j, BATCH_DIR, "fig10j_gaining_losing_by_region", w = 8, h = 6)
}

# =============================================================================
# 7. FANG FRAMEWORK: "Entropy-only" hidden layer
# =============================================================================
cat("\n[Fang] Identifying entropy-only hidden layer...\n")

# Sites where mean methylation barely changes but entropy does
entropy_dt[, abs_diff := abs(ampu_beta - ctrl_beta)]
entropy_dt[, abs_delta_ent := abs(delta_entropy)]

# Threshold: |delta_beta| < 0.05 (not a DMP candidate) + top 5% delta_entropy
beta_stable <- entropy_dt[abs_diff < 0.05]
ent_threshold <- quantile(beta_stable$abs_delta_ent, 0.95)
hidden_layer <- beta_stable[abs_delta_ent >= ent_threshold]
cat(sprintf("  Beta-stable CpGs (|diff| < 0.05): %s\n",
            format(nrow(beta_stable), big.mark = ",")))
cat(sprintf("  Entropy-only hidden layer (top 5%%): %s CpGs\n",
            format(nrow(hidden_layer), big.mark = ",")))
cat(sprintf("  Entropy threshold: %.4f\n", ent_threshold))

# Annotate hidden layer
hidden_layer[, region := annotate_regions(chr, pos, promoters_gr, exons, genes)]
hidden_region <- hidden_layer[, .N, by = region]
hidden_region[, pct := round(100 * N / sum(N), 1)]
cat("  Hidden layer by region:\n")
print(hidden_region)

# Direction of entropy change in hidden layer
hidden_layer[, ent_direction := ifelse(delta_entropy > 0, "Gaining", "Losing")]
hidden_dir <- hidden_layer[, .N, by = ent_direction]
cat("  Hidden layer direction:\n")
print(hidden_dir)

save_data(hidden_layer[, .(chr, pos, ctrl_beta, ampu_beta, ctrl_entropy, ampu_entropy,
                           delta_entropy, region, ent_direction)],
          BATCH_DIR, "entropy_only_hidden_layer")

# Fig 10K: Hidden layer — entropy change vs beta change scatter
set.seed(42)
# Use all beta_stable (no subsample), hex binning for density
p10k <- ggplot(beta_stable, aes(x = abs_diff, y = delta_entropy)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(trans = "log10", name = "CpGs (log10)") +
  geom_hline(yintercept = c(-ent_threshold, ent_threshold),
             linetype = "dashed", color = "red", linewidth = 0.7) +
  labs(x = "|Delta beta| (mean methylation change)",
       y = "Delta entropy (Amputated - Control)",
       title = "Entropy-only hidden layer: sites invisible to DMP analysis",
       subtitle = sprintf("Beta-stable (|diff|<0.05): %s CpGs | Hidden layer (top 5%% entropy): %s CpGs",
                          format(nrow(beta_stable), big.mark = ","),
                          format(nrow(hidden_layer), big.mark = ","))) +
  theme_minimal(base_size = 12)
save_fig(p10k, BATCH_DIR, "fig10k_entropy_only_hidden_layer", w = 9, h = 7)

# =============================================================================
# 8. FANG FRAMEWORK: Per-gene delta entropy + pathway analysis
# =============================================================================
cat("\n[Fang] Computing per-gene delta entropy...\n")

cpg_gr2 <- GRanges(seqnames = entropy_dt$chr,
                    ranges = IRanges(start = entropy_dt$pos, width = 1))
hits2 <- findOverlaps(cpg_gr2, genes)
gene_ent_dt <- data.table(
  cpg_idx = queryHits(hits2),
  gene_id = genes$ID[subjectHits(hits2)]
)
gene_ent_dt[, `:=`(
  ctrl_entropy  = entropy_dt$ctrl_entropy[cpg_idx],
  ampu_entropy  = entropy_dt$ampu_entropy[cpg_idx],
  delta_entropy = entropy_dt$delta_entropy[cpg_idx],
  ctrl_beta     = entropy_dt$ctrl_beta[cpg_idx],
  ampu_beta     = entropy_dt$ampu_beta[cpg_idx]
)]
per_gene_entropy <- gene_ent_dt[, .(
  n_cpgs            = .N,
  mean_ctrl_entropy = mean(ctrl_entropy),
  mean_ampu_entropy = mean(ampu_entropy),
  mean_delta        = mean(delta_entropy),
  sd_delta          = sd(delta_entropy),
  mean_ctrl_beta    = mean(ctrl_beta)
), by = gene_id][n_cpgs >= 5]

# Z-test per gene
per_gene_entropy[, z_score := mean_delta / (sd_delta / sqrt(n_cpgs))]
per_gene_entropy[, p_value := 2 * pnorm(-abs(z_score))]
per_gene_entropy[, fdr := p.adjust(p_value, method = "BH")]
per_gene_entropy[, direction := ifelse(mean_delta > 0, "Gaining", "Losing")]

sig_genes <- per_gene_entropy[fdr < 0.05]
cat(sprintf("  Genes with >=5 CpGs: %s\n", format(nrow(per_gene_entropy), big.mark = ",")))
cat(sprintf("  Genes with significant entropy change (FDR<0.05): %s\n",
            format(nrow(sig_genes), big.mark = ",")))
cat(sprintf("    Gaining: %d | Losing: %d\n",
            sum(sig_genes$direction == "Gaining"),
            sum(sig_genes$direction == "Losing")))

save_data(per_gene_entropy, BATCH_DIR, "per_gene_entropy")

# Fig 10L: Per-gene delta entropy volcano
p10l <- ggplot(per_gene_entropy, aes(x = mean_delta, y = -log10(fdr),
                                      color = ifelse(fdr < 0.05 & mean_delta > 0, "Gaining",
                                              ifelse(fdr < 0.05 & mean_delta < 0, "Losing", "NS")))) +
  geom_point(alpha = 0.3, size = 0.8) +
  scale_color_manual(values = c("Gaining" = "#E74C3C", "Losing" = "#3498DB", "NS" = "gray70"),
                     name = "Direction") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  labs(x = "Mean delta entropy per gene",
       y = "-log10(FDR)",
       title = "Per-gene entropy change volcano",
       subtitle = sprintf("Gaining: %d | Losing: %d (FDR < 0.05, >=5 CpGs)",
                          sum(sig_genes$direction == "Gaining"),
                          sum(sig_genes$direction == "Losing"))) +
  theme_minimal(base_size = 12)
save_fig(p10l, BATCH_DIR, "fig10l_per_gene_entropy_volcano", w = 9, h = 7)

rm(cpg_gr2, hits2, gene_ent_dt); gc(verbose = FALSE)

# =============================================================================
# 9. FANG FRAMEWORK: CpG density stratification
# =============================================================================
cat("\n[Fang] CpG density stratification...\n")

# Compute local CpG density: count CpGs within 500bp neighborhood
setkey(entropy_dt, chr, pos)
entropy_dt[, cpg_density := {
  counts <- integer(.N)
  for (ci in unique(chr)) {
    idx <- which(chr == ci)
    p <- pos[idx]
    # Vectorized: use findInterval for efficiency
    left  <- findInterval(p - 500L, p, left.open = TRUE) + 1L
    right <- findInterval(p + 500L, p)
    counts[idx] <- right - left + 1L
  }
  counts
}, by = chr]

# Bin into density quartiles
entropy_dt[, density_q := cut(cpg_density, breaks = quantile(cpg_density, probs = c(0, 0.25, 0.5, 0.75, 1)),
                               include.lowest = TRUE, labels = c("Q1 (sparse)", "Q2", "Q3", "Q4 (dense)"))]

density_strat <- entropy_dt[, .(
  mean_delta    = mean(delta_entropy),
  mean_ctrl_ent = mean(ctrl_entropy),
  n_cpgs        = .N
), by = density_q]
cat("  Entropy change by CpG density quartile:\n")
print(density_strat)

save_data(density_strat, BATCH_DIR, "entropy_by_cpg_density")

# Fig 10M: Delta entropy by CpG density quartile
p10m <- ggplot(density_strat, aes(x = density_q, y = mean_delta, fill = density_q)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_brewer(palette = "YlOrRd", guide = "none") +
  labs(x = "Local CpG density quartile (500bp window)",
       y = "Mean delta entropy (Amputated - Control)",
       title = "Entropy change stratified by local CpG density",
       subtitle = "Fang 2023: entropy inversely related to CpG density") +
  theme_minimal(base_size = 12)
save_fig(p10m, BATCH_DIR, "fig10m_entropy_by_cpg_density", w = 8, h = 6)

# Also stratify DMPs by density if available
if (has_dmps && nrow(dmp_entropy) > 0) {
  dmp_with_density <- merge(dmp_entropy, entropy_dt[, .(chr, pos, density_q)],
                            by = c("chr", "pos"), all.x = TRUE)
  dmp_density <- dmp_with_density[!is.na(density_q), .(
    mean_delta_ent = mean(delta_entropy),
    n_dmps = .N
  ), by = density_q]
  cat("  DMP entropy change by density:\n")
  print(dmp_density)
}

# =============================================================================
# 10. Direction-specific entropy (Hyper vs Hypo DMPs)
# =============================================================================
if (has_dmps && nrow(dmp_entropy) > 0) {
  cat("\n[Fang] Direction-specific entropy analysis...\n")

  # For each DMP direction: what happens to entropy?
  dir_entropy <- dmp_entropy[, .(
    mean_ctrl_entropy  = mean(ctrl_entropy),
    mean_ampu_entropy  = mean(ampu_entropy),
    mean_delta_entropy = mean(delta_entropy),
    mean_ctrl_beta     = mean(ctrl_beta),
    mean_ampu_beta     = mean(ampu_beta),
    pct_entropy_gain   = round(100 * mean(delta_entropy > 0), 1),
    n = .N
  ), by = direction]
  cat("  Direction-specific entropy:\n")
  print(dir_entropy)

  save_data(dir_entropy, BATCH_DIR, "direction_specific_entropy")

  # Fig 10N: Delta entropy density by DMP direction
  p10n <- ggplot(dmp_entropy, aes(x = delta_entropy, fill = direction, color = direction)) +
    geom_density(alpha = 0.3, linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = COLORS$direction) +
    scale_color_manual(values = COLORS$direction) +
    labs(x = "Delta entropy (Amputated - Control)",
         y = "Density",
         title = "Entropy change distribution by DMP direction",
         subtitle = sprintf("Hyper: %.1f%% gain entropy | Hypo: %.1f%% gain entropy",
                            dir_entropy[direction == "Hyper"]$pct_entropy_gain,
                            dir_entropy[direction == "Hypo"]$pct_entropy_gain)) +
    theme_minimal(base_size = 12)
  save_fig(p10n, BATCH_DIR, "fig10n_delta_entropy_by_direction", w = 9, h = 6)

  # Drift vs reprogramming interpretation
  # Low baseline entropy + significant change = directed reprogramming
  # High baseline entropy + change = could be drift
  cat("\n  Drift vs reprogramming summary:\n")
  dmp_entropy[, baseline_cat := ifelse(ctrl_entropy < 0.3, "Low (<0.3)",
                                ifelse(ctrl_entropy < 0.7, "Mid (0.3-0.7)", "High (>0.7)"))]
  drift_table <- dmp_entropy[, .(
    n = .N,
    mean_delta_ent = round(mean(delta_entropy), 5),
    pct_of_dmps = round(100 * .N / nrow(dmp_entropy), 1)
  ), by = .(direction, baseline_cat)]
  setorder(drift_table, direction, baseline_cat)
  cat("  DMPs by direction x baseline entropy:\n")
  print(drift_table)
  save_data(drift_table, BATCH_DIR, "drift_vs_reprogramming")
}

# Clean up temp columns
entropy_dt[, c("abs_diff", "abs_delta_ent", "cpg_density", "density_q") :=
             list(NULL, NULL, NULL, NULL)]

# =============================================================================
# Done
# =============================================================================
elapsed <- (proc.time() - t0)["elapsed"]
cat(sprintf("\n=== Batch 10 complete (%.1f s) ===\n", elapsed))
