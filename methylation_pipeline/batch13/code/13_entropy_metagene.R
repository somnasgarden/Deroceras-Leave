#!/usr/bin/env Rscript
# =============================================================================
# Batch 13: Entropy Metagene Profiles + Differential Entropy Gene Lists
# Question: Does normalized entropy vary across the gene structure by expression?
#           Which genes gain/lose entropy during regeneration? Do motifs explain it?
# Requires: CACHE$bsseq, CACHE$transcriptome, batch06 DMPs, batch1.5 motif hits
# =============================================================================

source("methylation_pipeline/_config.R")
t0 <- proc.time()

BATCH_DIR <- file.path(PIPE_DIR, "batch13")
dir.create(file.path(BATCH_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BATCH_DIR, "data"),    showWarnings = FALSE, recursive = TRUE)

# Clean old output
old_figs <- list.files(file.path(BATCH_DIR, "figures"), full.names = TRUE)
old_data <- list.files(file.path(BATCH_DIR, "data"), full.names = TRUE)
if (length(old_figs) > 0) { unlink(old_figs); cat(sprintf("  Cleaned %d old figures\n", length(old_figs))) }
if (length(old_data) > 0) { unlink(old_data); cat(sprintf("  Cleaned %d old data files\n", length(old_data))) }

cat("\n=== Batch 13: Entropy Metagene Profiles ===\n\n")

library(data.table)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(IRanges)
library(bsseq)

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("[1/6] Loading data...\n")

# BSseq
bs_obj <- readRDS(CACHE$bsseq)
bs_chr <- as.character(seqnames(bs_obj))
bs_pos <- start(bs_obj)

# Beta values (mean across biological replicates)
ctrl_beta <- rowMeans(getMeth(bs_obj, type = "raw")[, c("C1", "C2")], na.rm = TRUE)
ampu_beta <- rowMeans(getMeth(bs_obj, type = "raw")[, c("A1", "A2")], na.rm = TRUE)

# Per-site binary Shannon entropy
binary_entropy <- function(p) {
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  -p * log2(p) - (1 - p) * log2(1 - p)
}
ctrl_entropy <- binary_entropy(ctrl_beta)
ampu_entropy <- binary_entropy(ampu_beta)
delta_entropy <- ampu_entropy - ctrl_entropy

cat(sprintf("  CpGs: %s | Mean ctrl entropy: %.4f | Mean ampu entropy: %.4f\n",
    format(length(ctrl_entropy), big.mark = ","),
    mean(ctrl_entropy, na.rm = TRUE), mean(ampu_entropy, na.rm = TRUE)))

# GFF
gff <- load_gff()
genes <- gff[gff$type == "gene"]
genes <- genes[as.character(seqnames(genes)) %in% keep_chr]
exons <- gff[gff$type == "exon"]
exons <- exons[as.character(seqnames(exons)) %in% keep_chr]

# Expression
trans <- readRDS(CACHE$transcriptome)
vsd_mat <- trans$vst_clean

# Tail control expression (mean VST across control samples)
# VST matrix uses label columns like "Tail_Control_1", not raw filenames
tail_ctrl_cols <- grep("^Tail_Control", colnames(vsd_mat), value = TRUE)
if (length(tail_ctrl_cols) < 2) {
  # Fallback: use meta labels
  meta_clean <- trans$meta_clean
  tail_ctrl_labels <- meta_clean$label[meta_clean$tissue == "tail" & meta_clean$condition == "control"]
  tail_ctrl_cols <- intersect(colnames(vsd_mat), tail_ctrl_labels)
}
gene_expr <- data.table(gene_id = rownames(vsd_mat),
                         expression = rowMeans(vsd_mat[, tail_ctrl_cols, drop = FALSE], na.rm = TRUE))
gene_expr <- gene_expr[gene_id %in% genes$ID]
cat(sprintf("  Genes with expression: %s (using %d tail control columns)\n",
    format(nrow(gene_expr), big.mark = ","), length(tail_ctrl_cols)))

# =============================================================================
# 2. SETUP: gene structure classification + metagene machinery
# =============================================================================
cat("[2/6] Setting up metagene infrastructure...\n")

gene_strand <- as.character(strand(genes))
gene_starts <- start(genes); gene_ends <- end(genes)
gene_tss <- ifelse(gene_strand == "+", gene_starts, gene_ends)
gene_tts <- ifelse(gene_strand == "+", gene_ends, gene_starts)
gene_chr <- as.character(seqnames(genes))

exon_hits <- findOverlaps(exons, genes)
exon_to_gene <- data.table(exon_idx = queryHits(exon_hits), gene_idx = subjectHits(exon_hits),
  exon_start = start(exons)[queryHits(exon_hits)], exon_end = end(exons)[queryHits(exon_hits)])

N_BINS <- 20; FLANK <- 5000

exon_per_gene <- exon_to_gene[, .(n_exons = uniqueN(exon_idx)), by = gene_idx]
genes_with_expr <- which(genes$ID %in% gene_expr$gene_id)
idx_3plus <- intersect(exon_per_gene[n_exons >= 3, gene_idx], genes_with_expr)
idx_2exon <- intersect(exon_per_gene[n_exons == 2, gene_idx], genes_with_expr)
idx_1exon <- intersect(exon_per_gene[n_exons == 1, gene_idx], genes_with_expr)

cat(sprintf("  Gene counts — all: %d, 3+ exon: %d, 2-exon: %d, single-exon: %d\n",
    length(genes_with_expr), length(idx_3plus), length(idx_2exon), length(idx_1exon)))

# Segment orders
SEG_4 <- c("Distal upstream\n(5-2 kb)", "Promoter\n(2 kb)", "Gene body", "TTS downstream\n(5 kb)")
SEG_6 <- c("Distal upstream\n(5-2 kb)", "Promoter\n(2 kb)", "First exon", "Intron", "Last exon", "TTS downstream\n(5 kb)")
SEG_8 <- c("Distal upstream\n(5-2 kb)", "Promoter\n(2 kb)", "First exon", "First intron", "Body",
           "Last intron", "Last exon", "TTS downstream\n(5 kb)")

# Expression bin builders (reused from batch04)
make_expr_bins <- function(dt, n_bins) {
  d <- copy(dt)
  probs  <- seq(0, 1, length.out = n_bins + 1)
  raw_labels <- if (n_bins == 5) {
    c("Low (Qu1)", "Mid-low (Qu2)", "Mid (Qu3)", "Mid-high (Qu4)", "High (Qu5)")
  } else {
    paste0("D", sprintf("%02d", 1:10))
  }
  d[, expr_group := cut(expression,
                        breaks = quantile(expression, probs = probs, na.rm = TRUE),
                        labels = raw_labels,
                        include.lowest = TRUE)]
  d[, expr_group := factor(expr_group, levels = rev(raw_labels))]
  d
}

make_expr_colors <- function(n_bins) {
  if (n_bins == 5) {
    c("Low (Qu1)" = "#27AE60", "Mid-low (Qu2)" = "#2E86C1",
      "Mid (Qu3)" = "#F1C40F", "Mid-high (Qu4)" = "#E67E22",
      "High (Qu5)" = "#C0392B")
  } else {
    pal <- grDevices::colorRampPalette(c("#27AE60", "#2471A3", "#F1C40F",
                                          "#E67E22", "#C0392B"))(10)
    setNames(pal, paste0("D", sprintf("%02d", 1:10)))
  }
}

# Segment builder (from batch04)
build_segments <- function(gene_indices, group_lookup, mode = "3seg") {
  seg_list <- vector("list", length(gene_indices) * 9L)
  k <- 0L; n_used <- 0L
  for (gi in gene_indices) {
    gid <- genes$ID[gi]
    eg <- group_lookup[gene_id == gid, expr_group]
    if (length(eg) == 0 || is.na(eg[1])) next
    g_chr <- gene_chr[gi]; g_strand <- gene_strand[gi]
    g_tss <- gene_tss[gi]; g_tts <- gene_tts[gi]
    g_start <- gene_starts[gi]; g_end <- gene_ends[gi]

    if (mode == "3seg") {
      segs <- list(
        `Distal upstream\n(5-2 kb)` = c(max(1L, g_tss - FLANK), g_tss - 2001L),
        `Promoter\n(2 kb)` = c(max(1L, g_tss - 2000L), g_tss - 1L),
        `Gene body` = c(g_start, g_end),
        `TTS downstream\n(5 kb)` = c(g_tts + 1L, g_tts + FLANK))
    } else {
      g_exons <- exon_to_gene[gene_idx == gi][order(exon_start)]
      ne <- nrow(g_exons)
      if (mode == "5seg" && ne != 2) next
      if (mode == "7seg" && ne < 3) next

      if (mode == "5seg") {
        if (g_strand == "+") {
          fe <- c(g_exons$exon_start[1], g_exons$exon_end[1])
          le <- c(g_exons$exon_start[2], g_exons$exon_end[2])
          intr <- c(g_exons$exon_end[1] + 1L, g_exons$exon_start[2] - 1L)
        } else {
          fe <- c(g_exons$exon_start[2], g_exons$exon_end[2])
          le <- c(g_exons$exon_start[1], g_exons$exon_end[1])
          intr <- c(g_exons$exon_end[1] + 1L, g_exons$exon_start[2] - 1L)
        }
        segs <- list(
          `Distal upstream\n(5-2 kb)` = c(max(1L, g_tss - FLANK), g_tss - 2001L),
          `Promoter\n(2 kb)` = c(max(1L, g_tss - 2000L), g_tss - 1L),
          `First exon` = fe, Intron = intr, `Last exon` = le,
          `TTS downstream\n(5 kb)` = c(g_tts + 1L, g_tts + FLANK))
      } else {
        if (g_strand == "+") {
          fe <- c(g_exons$exon_start[1], g_exons$exon_end[1])
          le <- c(g_exons$exon_start[ne], g_exons$exon_end[ne])
          fi <- c(g_exons$exon_end[1] + 1L, g_exons$exon_start[2] - 1L)
          li <- c(g_exons$exon_end[ne - 1] + 1L, g_exons$exon_start[ne] - 1L)
          body_c <- c(g_exons$exon_start[2], g_exons$exon_end[ne - 1])
        } else {
          fe <- c(g_exons$exon_start[ne], g_exons$exon_end[ne])
          le <- c(g_exons$exon_start[1], g_exons$exon_end[1])
          fi <- c(g_exons$exon_end[ne - 1] + 1L, g_exons$exon_start[ne] - 1L)
          li <- c(g_exons$exon_end[1] + 1L, g_exons$exon_start[2] - 1L)
          body_c <- c(g_exons$exon_start[2], g_exons$exon_end[ne - 1])
        }
        segs <- list(
          `Distal upstream\n(5-2 kb)` = c(max(1L, g_tss - FLANK), g_tss - 2001L),
          `Promoter\n(2 kb)` = c(max(1L, g_tss - 2000L), g_tss - 1L),
          `First exon` = fe, `First intron` = fi, Body = body_c,
          `Last intron` = li, `Last exon` = le,
          `TTS downstream\n(5 kb)` = c(g_tts + 1L, g_tts + FLANK))
      }
    }

    n_used <- n_used + 1L
    for (sn in names(segs)) {
      s <- segs[[sn]]
      if (is.na(s[1]) || is.na(s[2]) || s[1] >= s[2] || s[1] < 1L) next
      k <- k + 1L
      seg_list[[k]] <- data.table(
        chr = g_chr, start = s[1], end = s[2], strand = g_strand,
        gene_id = gid, segment = sn, expr_group = as.character(eg[1]))
    }
  }
  list(segments = rbindlist(seg_list[1:k]), n_genes = n_used)
}

# Entropy metagene: same as run_metagene but uses entropy instead of beta*100
run_entropy_metagene <- function(entropy_vec, seg_result, seg_order) {
  seg_dt <- seg_result$segments
  seg_gr <- GRanges(seqnames = seg_dt$chr,
                     ranges = IRanges(start = seg_dt$start, end = seg_dt$end))
  cpg_gr <- GRanges(seqnames = bs_chr, ranges = IRanges(start = bs_pos, width = 1))
  ov <- findOverlaps(cpg_gr, seg_gr)

  ov_dt <- data.table(
    cpg_idx = queryHits(ov), seg_idx = subjectHits(ov),
    entropy = entropy_vec[queryHits(ov)], pos = bs_pos[queryHits(ov)])
  ov_dt <- ov_dt[!is.na(entropy)]
  ov_dt[, `:=`(seg_start = seg_dt$start[seg_idx], seg_end = seg_dt$end[seg_idx],
               seg_strand = seg_dt$strand[seg_idx], segment = seg_dt$segment[seg_idx],
               expr_group = seg_dt$expr_group[seg_idx],
               gene_id = seg_dt$gene_id[seg_idx])]

  ov_dt[, rel_pos := (pos - seg_start) / (seg_end - seg_start)]
  ov_dt[seg_strand == "-" & !segment %in% c("Body", "Gene body"), rel_pos := 1 - rel_pos]
  ov_dt[, bin := pmin(N_BINS, pmax(1L, ceiling(rel_pos * N_BINS)))]

  gene_bins <- ov_dt[, .(mean_entropy = mean(entropy)),
                     by = .(gene_id, segment, bin, expr_group)]
  agg <- gene_bins[, .(mean_entropy = mean(mean_entropy), se = sd(mean_entropy) / sqrt(.N)),
                    by = .(segment, bin, expr_group)]
  agg[, segment := factor(segment, levels = seg_order)]
  agg[, x_pos := (as.integer(segment) - 1) * N_BINS + bin]
  rm(ov_dt, gene_bins, ov); gc(verbose = FALSE)
  agg
}

# Plot entropy metagene
plot_entropy_metagene <- function(mg_agg, title, subtitle, colors, seg_order) {
  n_segs <- length(seg_order)
  seg_breaks <- (1:(n_segs - 1)) * N_BINS + 0.5
  seg_labels_pos <- (0:(n_segs - 1)) * N_BINS + N_BINS / 2

  ggplot(mg_agg, aes(x = x_pos, y = mean_entropy, color = expr_group)) +
    geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 1) +
    geom_vline(xintercept = seg_breaks, linetype = "dashed", color = "gray70", linewidth = 0.3) +
    scale_color_manual(values = colors) +
    scale_x_continuous(breaks = seg_labels_pos, labels = seg_order) +
    labs(title = title, subtitle = subtitle,
         x = NULL, y = "Binary Shannon Entropy (bits)", color = "Expression\ngroup") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank())
}

# =============================================================================
# 3. ENTROPY METAGENE PROFILES (5 + 10 deciles × 3 gene types × ctrl + ampu)
# =============================================================================
cat("[3/6] Generating entropy metagene profiles...\n")

metagene_specs <- list(
  list(name = "3plus_exon",   label = "3+ exon genes",            idx = idx_3plus, mode = "7seg", seg = SEG_8, ampu = TRUE),
  list(name = "2exon",        label = "2-exon genes",             idx = idx_2exon, mode = "5seg", seg = SEG_6, ampu = TRUE),
  list(name = "single_exon",  label = "Single-exon (intronless)", idx = idx_1exon, mode = "3seg", seg = SEG_4, ampu = TRUE)
)

binnings <- list(
  list(n = 5,  tag = "bin05", label = "expression quintiles (5 bins)"),
  list(n = 10, tag = "bin10", label = "expression deciles (10 bins)")
)

for (binning in binnings) {
  cat(sprintf("\n  === Entropy metagenes %s ===\n", binning$label))
  gene_bins_dt <- make_expr_bins(gene_expr, binning$n)
  colors_b <- make_expr_colors(binning$n)

  fig_num <- 0L
  for (spec in metagene_specs) {
    gi <- spec$idx
    if (length(gi) < 30) {
      cat(sprintf("    Skipping %s — too few genes (%d)\n", spec$label, length(gi))); next
    }

    cat(sprintf("    Entropy metagene: %s (n=%d)...\n", spec$label, length(gi)))
    seg_b <- build_segments(gi, gene_bins_dt, mode = spec$mode)
    n_genes <- seg_b$n_genes

    # Control entropy
    fig_num <- fig_num + 1L
    mg_ctrl <- run_entropy_metagene(ctrl_entropy, seg_b, spec$seg)
    p <- plot_entropy_metagene(mg_ctrl,
      sprintf("Entropy metagene — %s (Control, %s)", spec$label, binning$label),
      sprintf("Controls (C1+C2), n = %s genes | H = -p*log2(p) - (1-p)*log2(1-p)", format(n_genes, big.mark = ",")),
      colors_b, spec$seg)
    save_fig(p, BATCH_DIR,
             sprintf("fig13a_%s_%02d_ctrl_%s", binning$tag, fig_num, spec$name),
             w = 11, h = 6)
    rm(mg_ctrl); gc(verbose = FALSE)

    # Amputated entropy
    fig_num <- fig_num + 1L
    mg_ampu <- run_entropy_metagene(ampu_entropy, seg_b, spec$seg)
    p <- plot_entropy_metagene(mg_ampu,
      sprintf("Entropy metagene — %s (Amputated, %s)", spec$label, binning$label),
      sprintf("Amputated (A1+A2), n = %s genes | H = -p*log2(p) - (1-p)*log2(1-p)", format(n_genes, big.mark = ",")),
      colors_b, spec$seg)
    save_fig(p, BATCH_DIR,
             sprintf("fig13a_%s_%02d_ampu_%s", binning$tag, fig_num, spec$name),
             w = 11, h = 6)
    rm(mg_ampu); gc(verbose = FALSE)

    # Delta entropy metagene (ampu - ctrl)
    fig_num <- fig_num + 1L
    mg_delta <- run_entropy_metagene(delta_entropy, seg_b, spec$seg)
    # Override y-axis label for delta
    p <- plot_entropy_metagene(mg_delta,
      sprintf("\u0394 Entropy metagene — %s (%s)", spec$label, binning$label),
      sprintf("Amputated minus Control, n = %s genes | Positive = entropy gain", format(n_genes, big.mark = ",")),
      colors_b, spec$seg)
    p <- p + labs(y = "\u0394 Entropy (ampu - ctrl, bits)")
    save_fig(p, BATCH_DIR,
             sprintf("fig13a_%s_%02d_delta_%s", binning$tag, fig_num, spec$name),
             w = 11, h = 6)
    rm(mg_delta, seg_b); gc(verbose = FALSE)
  }
  rm(gene_bins_dt); gc(verbose = FALSE)
}

# =============================================================================
# 4. PER-GENE ENTROPY: identify genes with increased/decreased entropy
# =============================================================================
cat("[4/6] Computing per-gene entropy changes...\n")

# Assign CpGs to genes (gene body only, matching batch10 approach)
cpg_gr <- GRanges(seqnames = bs_chr, ranges = IRanges(start = bs_pos, width = 1))
gene_gr <- genes[genes$ID %in% gene_expr$gene_id]
ov <- findOverlaps(cpg_gr, gene_gr)

cpg_gene_dt <- data.table(
  cpg_idx = queryHits(ov),
  gene_id = gene_gr$ID[subjectHits(ov)],
  ctrl_entropy = ctrl_entropy[queryHits(ov)],
  ampu_entropy = ampu_entropy[queryHits(ov)],
  delta_entropy = delta_entropy[queryHits(ov)],
  ctrl_beta = ctrl_beta[queryHits(ov)],
  ampu_beta = ampu_beta[queryHits(ov)]
)
cpg_gene_dt <- cpg_gene_dt[!is.na(ctrl_entropy) & !is.na(ampu_entropy)]

# Per-gene summary: mean entropy, delta, Wilcoxon test
gene_entropy <- cpg_gene_dt[, .(
  n_cpgs = .N,
  mean_ctrl_entropy = mean(ctrl_entropy),
  mean_ampu_entropy = mean(ampu_entropy),
  mean_delta_entropy = mean(delta_entropy),
  mean_ctrl_beta = mean(ctrl_beta),
  mean_ampu_beta = mean(ampu_beta),
  p_wilcox = if (.N >= 5) suppressWarnings(wilcox.test(ctrl_entropy, ampu_entropy, paired = TRUE)$p.value) else NA_real_
), by = gene_id]

gene_entropy <- merge(gene_entropy, gene_expr, by = "gene_id")
gene_entropy[, padj := p.adjust(p_wilcox, method = "BH")]

# Add gene names from EviAnn annotations
if (file.exists(OG$annot)) {
  annot <- fread(OG$annot)
  name_col <- intersect(c("gene_name", "Name", "name"), names(annot))
  id_col   <- intersect(c("gene_id", "ID", "#gene_id"), names(annot))
  if (length(name_col) > 0 && length(id_col) > 0) {
    annot_sub <- annot[, c(id_col[1], name_col[1]), with = FALSE]
    setnames(annot_sub, c("gene_id", "gene_name"))
    gene_entropy <- merge(gene_entropy, annot_sub, by = "gene_id", all.x = TRUE)
  }
}

# Classify entropy direction
gene_entropy[, entropy_class := fifelse(
  padj < 0.05 & mean_delta_entropy > 0, "Entropy_Up",
  fifelse(padj < 0.05 & mean_delta_entropy < 0, "Entropy_Down", "NS"))]

n_up   <- sum(gene_entropy$entropy_class == "Entropy_Up", na.rm = TRUE)
n_down <- sum(gene_entropy$entropy_class == "Entropy_Down", na.rm = TRUE)
n_ns   <- sum(gene_entropy$entropy_class == "NS", na.rm = TRUE)
cat(sprintf("  Entropy genes: %d up, %d down, %d NS (total %d)\n", n_up, n_down, n_ns, nrow(gene_entropy)))

save_data(gene_entropy, BATCH_DIR, "gene_entropy_summary")

# Add exon count for gene structure classification
gene_entropy <- merge(gene_entropy,
  exon_per_gene[, .(gene_idx, n_exons)][, gene_id := genes$ID[gene_idx]][, .(gene_id, n_exons)],
  by = "gene_id", all.x = TRUE)
gene_entropy[, gene_class := fifelse(n_exons >= 3, "3+ exon",
                              fifelse(n_exons == 2, "2-exon", "Intronless"))]

# Save separate lists
entropy_up   <- gene_entropy[entropy_class == "Entropy_Up"][order(mean_delta_entropy, decreasing = TRUE)]
entropy_down <- gene_entropy[entropy_class == "Entropy_Down"][order(mean_delta_entropy)]
save_data(entropy_up, BATCH_DIR, "genes_entropy_increased")
save_data(entropy_down, BATCH_DIR, "genes_entropy_decreased")

# =============================================================================
# 5. ENTROPY GENE PLOTS
# =============================================================================
cat("[5/6] Plotting entropy gene-level results...\n")

# Fig 13b: Volcano-style: delta entropy vs -log10(padj)
gene_entropy[, neg_log10p := -log10(pmin(padj, 1))]
gene_entropy[is.na(neg_log10p), neg_log10p := 0]
gene_entropy[, plot_color := fifelse(entropy_class == "Entropy_Up", "#C0392B",
                              fifelse(entropy_class == "Entropy_Down", "#2471A3", "grey70"))]

p13b <- ggplot(gene_entropy[!is.na(padj)],
               aes(x = mean_delta_entropy, y = neg_log10p, color = plot_color)) +
  geom_point(alpha = 0.3, size = 0.8) +
  scale_color_identity() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Per-Gene Entropy Change During Regeneration",
       subtitle = sprintf("Entropy Up: %d | Entropy Down: %d | NS: %d (Wilcoxon paired, BH-adjusted)",
                          n_up, n_down, n_ns),
       x = "\u0394 Entropy (amputated - control, bits)",
       y = "-log10(adjusted p-value)") +
  theme_minimal(base_size = 12)
save_fig(p13b, BATCH_DIR, "fig13b_entropy_volcano", w = 9, h = 7)

# Fig 13c: Delta entropy by gene structure class
p13c <- ggplot(gene_entropy[!is.na(gene_class)],
               aes(x = gene_class, y = mean_delta_entropy, fill = gene_class)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("3+ exon" = "#2471A3", "2-exon" = "#E67E22", "Intronless" = "#27AE60"),
                    guide = "none") +
  labs(title = "Entropy Change by Gene Structure",
       subtitle = "Paired Wilcoxon per gene (ctrl vs ampu CpGs)",
       x = "Gene class", y = "\u0394 Entropy (bits)") +
  theme_minimal(base_size = 12)
save_fig(p13c, BATCH_DIR, "fig13c_entropy_by_gene_class", w = 7, h = 6)

# Fig 13d: Delta entropy vs expression (scatter)
p13d <- ggplot(gene_entropy, aes(x = expression, y = mean_delta_entropy)) +
  geom_point(alpha = 0.05, size = 0.5, color = "#2471A3") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = "#C0392B", linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  labs(title = "Entropy Change vs Expression Level",
       subtitle = sprintf("n = %s genes | GAM fit", format(nrow(gene_entropy), big.mark = ",")),
       x = "Expression (VST, tail controls)", y = "\u0394 Entropy (bits)") +
  theme_minimal(base_size = 12)
save_fig(p13d, BATCH_DIR, "fig13d_entropy_vs_expression", w = 8, h = 7)

# Fig 13e: Entropy vs methylation level — shows the concavity relationship
# (sites near 0% or 100% have low entropy by definition)
ctrl_sample_idx <- sample(length(ctrl_beta), min(500000, length(ctrl_beta)))
sample_dt <- data.table(beta = ctrl_beta[ctrl_sample_idx], entropy = ctrl_entropy[ctrl_sample_idx])
p13e <- ggplot(sample_dt, aes(x = beta * 100, y = entropy)) +
  geom_hex(bins = 100) +
  scale_fill_viridis_c(trans = "log10", name = "Count") +
  stat_function(fun = function(x) binary_entropy(x / 100), color = "red",
                linewidth = 1, linetype = "dashed") +
  labs(title = "Methylation Level vs Entropy (Controls)",
       subtitle = "Red = theoretical curve | Sites near 0%/100% have low entropy by definition",
       x = "Methylation (%)", y = "Binary Shannon Entropy (bits)") +
  theme_minimal(base_size = 12)
save_fig(p13e, BATCH_DIR, "fig13e_beta_vs_entropy_hex", w = 9, h = 7)

# Fig 13f: Entropy change by expression decile (10 bins, boxplot)
gene_entropy_dec <- make_expr_bins(gene_entropy, 10)
p13f <- ggplot(gene_entropy_dec, aes(x = expr_group, y = mean_delta_entropy, fill = expr_group)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_fill_manual(values = make_expr_colors(10), guide = "none") +
  labs(title = "Entropy Change by Expression Decile",
       subtitle = "D01=lowest expression, D10=highest",
       x = "Expression decile", y = "\u0394 Entropy (ampu - ctrl, bits)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p13f, BATCH_DIR, "fig13f_entropy_by_decile", w = 9, h = 6)

# =============================================================================
# 6. MOTIF ENRICHMENT IN DIFFERENTIAL ENTROPY GENES
# =============================================================================
cat("[6/6] Motif enrichment at differential entropy genes...\n")

# Load DMP annotations if available
dmp_file <- file.path(PIPE_DIR, "batch06/data/dmps_annotated.tsv")
has_dmps <- file.exists(dmp_file)
if (has_dmps) {
  dmps <- fread(dmp_file)
  cat(sprintf("  Loaded %s DMPs\n", format(nrow(dmps), big.mark = ",")))

  # Count DMPs per gene
  if ("gene_id" %in% names(dmps)) {
    dmp_counts <- dmps[!is.na(gene_id), .N, by = gene_id]
    setnames(dmp_counts, "N", "n_dmps")
    gene_entropy <- merge(gene_entropy, dmp_counts, by = "gene_id", all.x = TRUE)
    gene_entropy[is.na(n_dmps), n_dmps := 0L]
  }
}

# Load DESeq2 results for tail if available
has_de <- FALSE
if (!is.null(trans$res_tail)) {
  res_tail <- as.data.table(trans$res_tail, keep.rownames = "gene_id")
  if (!"gene_id" %in% names(res_tail) && "rn" %in% names(res_tail)) setnames(res_tail, "rn", "gene_id")
  gene_entropy <- merge(gene_entropy, res_tail[, .(gene_id, log2FoldChange, padj_de = padj)],
                        by = "gene_id", all.x = TRUE)
  has_de <- TRUE
  cat("  Merged DESeq2 tail results\n")
}

# Load motif hits from batch 1.5 cache
motif_file <- CACHE$extended
has_motifs <- file.exists(motif_file)
if (has_motifs) {
  cat("  Loading genome-wide motif hits...\n")
  ext_regions <- readRDS(motif_file)

  # Check if this is a motif hits table or extended regions
  if (is.data.frame(ext_regions) || is.data.table(ext_regions)) {
    motif_hits <- as.data.table(ext_regions)
    cat(sprintf("  Motif hits loaded: %s rows\n", format(nrow(motif_hits), big.mark = ",")))
  } else {
    has_motifs <- FALSE
    cat("  Extended regions cache is not motif hits — skipping motif analysis\n")
  }
}

# Check batch 1.5 motif data directly
motif_data_file <- file.path(PIPE_DIR, "batch1.5/data/genome_wide_motif_hits.tsv")
if (!has_motifs && file.exists(motif_data_file)) {
  cat("  Loading motif hits from batch 1.5 TSV...\n")
  motif_hits <- fread(motif_data_file)
  has_motifs <- TRUE
  cat(sprintf("  Motif hits: %s rows\n", format(nrow(motif_hits), big.mark = ",")))
}

if (has_motifs && "gene_id" %in% names(motif_hits)) {
  # Count unique motifs per gene
  motif_per_gene <- motif_hits[, .(n_motifs = uniqueN(motif_id)), by = gene_id]
  gene_entropy <- merge(gene_entropy, motif_per_gene, by = "gene_id", all.x = TRUE)
  gene_entropy[is.na(n_motifs), n_motifs := 0L]

  # Fisher's exact: are entropy-up/down genes enriched for motif-containing genes?
  gene_entropy[, has_motif := n_motifs > 0]
  for (cls in c("Entropy_Up", "Entropy_Down")) {
    genes_cls <- gene_entropy[entropy_class == cls]
    genes_bg  <- gene_entropy[entropy_class == "NS"]
    if (nrow(genes_cls) < 5 || nrow(genes_bg) < 5) next

    ct <- matrix(c(
      sum(genes_cls$has_motif), sum(!genes_cls$has_motif),
      sum(genes_bg$has_motif),  sum(!genes_bg$has_motif)
    ), nrow = 2, byrow = TRUE)
    ft <- fisher.test(ct)
    cat(sprintf("  %s motif enrichment: OR=%.2f, p=%.2e\n", cls, ft$estimate, ft$p.value))
  }

  # Per-motif enrichment in entropy-up genes (top enriched TF motifs)
  if ("motif_id" %in% names(motif_hits)) {
    entropy_up_ids <- gene_entropy[entropy_class == "Entropy_Up", gene_id]
    entropy_ns_ids <- gene_entropy[entropy_class == "NS", gene_id]

    motif_enrichment <- motif_hits[gene_id %in% c(entropy_up_ids, entropy_ns_ids)][,
      .(in_up = sum(gene_id %in% entropy_up_ids),
        in_ns = sum(gene_id %in% entropy_ns_ids)), by = motif_id]
    motif_enrichment <- motif_enrichment[in_up >= 3 & in_ns >= 3]

    if (nrow(motif_enrichment) > 0) {
      n_up_tot <- length(entropy_up_ids)
      n_ns_tot <- length(entropy_ns_ids)
      motif_enrichment[, p_fisher := mapply(function(a, b) {
        fisher.test(matrix(c(a, n_up_tot - a, b, n_ns_tot - b), nrow = 2))$p.value
      }, in_up, in_ns)]
      motif_enrichment[, OR := mapply(function(a, b) {
        fisher.test(matrix(c(a, n_up_tot - a, b, n_ns_tot - b), nrow = 2))$estimate
      }, in_up, in_ns)]
      motif_enrichment[, padj_motif := p.adjust(p_fisher, method = "BH")]
      motif_enrichment <- motif_enrichment[order(p_fisher)]
      save_data(motif_enrichment, BATCH_DIR, "motif_enrichment_entropy_up")

      # Fig 13g: Top motifs enriched in entropy-up genes
      top_motifs <- head(motif_enrichment[padj_motif < 0.1], 20)
      if (nrow(top_motifs) > 0) {
        p13g <- ggplot(top_motifs, aes(x = reorder(motif_id, -log10(p_fisher)),
                                        y = -log10(p_fisher), fill = log2(OR))) +
          geom_col() +
          scale_fill_gradient2(low = "#2471A3", mid = "white", high = "#C0392B", midpoint = 0,
                               name = "log2(OR)") +
          geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
          labs(title = "TFBS Motifs Enriched in Entropy-Up Genes",
               subtitle = sprintf("%d entropy-up vs %d NS genes", n_up_tot, n_ns_tot),
               x = "Motif", y = "-log10(Fisher p)") +
          theme_minimal(base_size = 11) +
          theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
        save_fig(p13g, BATCH_DIR, "fig13g_motif_enrichment_entropy_up", w = 11, h = 7)
      }
    }
  }
} else {
  cat("  No motif data available — skipping motif enrichment\n")
}

# Fig 13h: Entropy-up/down genes overlaid on DMP count (if available)
if (has_dmps && "n_dmps" %in% names(gene_entropy)) {
  p13h <- ggplot(gene_entropy[entropy_class != "NS"],
                 aes(x = n_dmps, y = mean_delta_entropy, color = entropy_class)) +
    geom_point(alpha = 0.3, size = 1) +
    scale_color_manual(values = c("Entropy_Up" = "#C0392B", "Entropy_Down" = "#2471A3"),
                       name = "Entropy\ndirection") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(title = "Entropy Change vs DMP Count per Gene",
         subtitle = "Do genes with more DMPs show more entropy change?",
         x = "Number of DMPs in gene", y = "\u0394 Entropy (bits)") +
    theme_minimal(base_size = 12)
  save_fig(p13h, BATCH_DIR, "fig13h_entropy_vs_dmp_count", w = 9, h = 7)
}

# Fig 13i: Four-way classification if DE data available
if (has_de) {
  gene_entropy[, de_class := fifelse(padj_de < 0.05 & log2FoldChange > 0, "DE_Up",
                              fifelse(padj_de < 0.05 & log2FoldChange < 0, "DE_Down", "NS_DE"))]
  gene_entropy[, combined := paste(entropy_class, de_class, sep = "_")]

  quad_counts <- gene_entropy[entropy_class != "NS" & de_class != "NS_DE",
    .N, by = .(entropy_class, de_class)]
  cat("  Entropy × DE quadrant counts:\n")
  print(quad_counts)
  save_data(quad_counts, BATCH_DIR, "entropy_de_quadrant_counts")

  # Concordant gene lists (entropy + DE)
  concordant <- gene_entropy[entropy_class != "NS" & de_class != "NS_DE"]
  if (nrow(concordant) > 0) {
    save_data(concordant[order(padj)], BATCH_DIR, "concordant_entropy_de_genes")
    cat(sprintf("  Concordant entropy-DE genes: %d\n", nrow(concordant)))
  }

  # Plot: entropy vs LFC colored by significance
  gene_entropy[, quad_color := fifelse(
    entropy_class != "NS" & de_class != "NS_DE", "#C0392B",
    fifelse(entropy_class != "NS", "#E67E22",
    fifelse(de_class != "NS_DE", "#2471A3", "grey80")))]

  p13i <- ggplot(gene_entropy[!is.na(log2FoldChange)],
                 aes(x = log2FoldChange, y = mean_delta_entropy, color = quad_color)) +
    geom_point(alpha = 0.2, size = 0.6) +
    scale_color_identity() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "Entropy Change vs Expression Change (Tail Regeneration)",
         subtitle = sprintf("Red: both sig | Orange: entropy only | Blue: DE only | Grey: NS | Concordant: %d genes",
                            nrow(concordant)),
         x = "log2(Fold Change), amputated vs control",
         y = "\u0394 Entropy (bits)") +
    theme_minimal(base_size = 12)
  save_fig(p13i, BATCH_DIR, "fig13i_entropy_vs_lfc_quadrant", w = 9, h = 8)
}

# =============================================================================
# DONE
# =============================================================================
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\n=== Batch 13 complete. Elapsed: %.1f min ===\n", elapsed / 60))
