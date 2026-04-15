#!/usr/bin/env Rscript
# =============================================================================
# Batch 12: Promoter Methylation Spike Characterization
# Question: What causes the methylation spike in the promoter of highly expressed
#           3+ exon genes? Is it a sequence feature, nucleosome, or TFBS?
# Requires: CACHE$bsseq, CACHE$transcriptome, CACHE$genome, batch1.5 motif hits
# =============================================================================

source("methylation_pipeline/_config.R")
t0 <- proc.time()

suppressPackageStartupMessages({
  library(bsseq)
  library(GenomicRanges)
  library(IRanges)
  library(data.table)
  library(ggplot2)
  library(Biostrings)
})

BATCH_DIR <- file.path(PIPE_DIR, "batch12")

# Clean old output
unlink(list.files(file.path(BATCH_DIR, "figures"), full.names = TRUE))
unlink(list.files(file.path(BATCH_DIR, "data"), full.names = TRUE))
dir.create(file.path(BATCH_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BATCH_DIR, "data"), showWarnings = FALSE, recursive = TRUE)

cat("=== Batch 12: Promoter Spike Characterization ===\n\n")

# Constants
TSS_UP   <- 5000L   # upstream of TSS
TSS_DOWN <- 2000L   # downstream of TSS
BIN_W    <- 10L      # aggregation bin width (bp)

# Decile color palette (same as batch04)
DCOLS <- grDevices::colorRampPalette(
  c("#27AE60", "#2471A3", "#F1C40F", "#E67E22", "#C0392B"))(10)
names(DCOLS) <- paste0("D", sprintf("%02d", 1:10))

# =============================================================================
# [1/8] LOAD DATA
# =============================================================================
cat("[1/8] Loading data...\n")

# BSseq
bs_obj <- readRDS(CACHE$bsseq)
ctrl_idx <- which(colnames(bs_obj) %in% c("C1", "C2"))
ampu_idx <- which(colnames(bs_obj) %in% c("A1", "A2"))
ctrl_beta <- rowMeans(getMeth(bs_obj[, ctrl_idx], type = "raw"), na.rm = TRUE)
ampu_beta <- rowMeans(getMeth(bs_obj[, ampu_idx], type = "raw"), na.rm = TRUE)
cpg_gr <- GRanges(seqnames = seqnames(bs_obj), ranges = IRanges(start = start(bs_obj), width = 1))
cat(sprintf("  CpGs: %s\n", format(length(cpg_gr), big.mark = ",")))

# GFF
gff <- load_gff()
genes <- gff[gff$type == "gene"]
exons <- gff[gff$type == "exon"]

# Expression + deciles
if (file.exists(CACHE$transcriptome)) {
  tc <- readRDS(CACHE$transcriptome)
  vst_mat <- tc$vst_clean
  ctrl_cols <- grep("Tail_Control|^C\\d+S\\d+", colnames(vst_mat), value = TRUE)
  if (length(ctrl_cols) == 0) ctrl_cols <- colnames(vst_mat)
  expr_vec <- rowMeans(vst_mat[, ctrl_cols, drop = FALSE])
  rm(tc, vst_mat); gc(verbose = FALSE)
} else {
  stop("CACHE$transcriptome not found — run batch 0.5 first")
}

# Per-gene methylation + expression
hits_gb <- findOverlaps(cpg_gr, genes)
gene_dt <- data.table(gene_id = genes$ID[subjectHits(hits_gb)],
                       beta = ctrl_beta[queryHits(hits_gb)])
gene_dt <- gene_dt[!is.na(beta)]
gene_means <- gene_dt[, .(mean_beta = mean(beta), n_cpgs = .N), by = gene_id]
gene_means[, expression := expr_vec[match(gene_id, names(expr_vec))]]
gene_means <- gene_means[!is.na(expression)]
rm(gene_dt); gc(verbose = FALSE)

# Assign deciles
probs <- seq(0, 1, length.out = 11)
d_labels <- paste0("D", sprintf("%02d", 1:10))
gene_means[, expr_decile := cut(expression,
  breaks = quantile(expression, probs = probs, na.rm = TRUE),
  labels = d_labels, include.lowest = TRUE)]

# 3+ exon gene indices
exon_hits <- findOverlaps(exons, genes)
exon_to_gene <- data.table(exon_idx = queryHits(exon_hits), gene_idx = subjectHits(exon_hits))
exon_per_gene <- exon_to_gene[, .(n_exons = uniqueN(exon_idx)), by = gene_idx]
genes_with_expr <- which(genes$ID %in% gene_means$gene_id)
idx_3plus <- intersect(exon_per_gene[n_exons >= 3, gene_idx], genes_with_expr)
cat(sprintf("  3+ exon genes with expression: %d\n", length(idx_3plus)))

# TSS positions
gene_strand <- as.character(strand(genes))
gene_starts <- start(genes); gene_ends <- end(genes)
gene_tss <- ifelse(gene_strand == "+", gene_starts, gene_ends)
gene_chr <- as.character(seqnames(genes))

# Build gene info table for 3+ exon genes
gi <- data.table(
  gene_idx = idx_3plus,
  gene_id  = genes$ID[idx_3plus],
  chr      = gene_chr[idx_3plus],
  tss      = gene_tss[idx_3plus],
  str      = gene_strand[idx_3plus]
)
gi <- merge(gi, gene_means[, .(gene_id, expr_decile)], by = "gene_id")
cat(sprintf("  Genes per decile: %s\n",
    paste(gi[, .N, by = expr_decile][order(expr_decile), sprintf("%s=%d", expr_decile, N)],
          collapse = ", ")))

# =============================================================================
# [2/8] SINGLE-BP RESOLUTION: METHYLATION vs TSS DISTANCE
# =============================================================================
cat("\n[2/8] Computing methylation vs TSS distance...\n")

# Build TSS windows for all 3+ exon genes
tss_windows <- GRanges(seqnames = gi$chr,
  ranges = IRanges(
    start = ifelse(gi$str == "+", gi$tss - TSS_UP, gi$tss - TSS_DOWN),
    end   = ifelse(gi$str == "+", gi$tss + TSS_DOWN, gi$tss + TSS_UP)
  ))
tss_windows <- trim(tss_windows)

# Bulk overlap
cat("  Finding CpG-TSS overlaps...\n")
ov <- findOverlaps(cpg_gr, tss_windows)
cat(sprintf("  Overlaps: %s\n", format(length(ov), big.mark = ",")))

# Compute strand-aware distance
ov_dt <- data.table(
  cpg_idx  = queryHits(ov),
  gene_row = subjectHits(ov)
)
ov_dt[, cpg_pos := start(cpg_gr)[cpg_idx]]
ov_dt[, tss := gi$tss[gene_row]]
ov_dt[, str := gi$str[gene_row]]
ov_dt[, dist := ifelse(str == "+", cpg_pos - tss, tss - cpg_pos)]
ov_dt <- ov_dt[dist >= -TSS_UP & dist <= TSS_DOWN]
ov_dt[, beta := ctrl_beta[cpg_idx]]
ov_dt[, beta_a := ampu_beta[cpg_idx]]
ov_dt[, decile := gi$expr_decile[gene_row]]
ov_dt <- ov_dt[!is.na(beta)]
cat(sprintf("  CpG-gene pairs in window: %s\n", format(nrow(ov_dt), big.mark = ",")))

# Bin into BIN_W-bp windows
ov_dt[, bin := floor(dist / BIN_W) * BIN_W]
prof <- ov_dt[, .(mean_beta = mean(beta) * 100, mean_beta_a = mean(beta_a) * 100,
                   n = .N), by = .(decile, bin)]
save_data(prof, BATCH_DIR, "tss_distance_methylation_by_decile")

# Fig 12A: Full window
cat("  Fig 12A: Methylation vs TSS distance (full window)...\n")
prof[, decile := factor(decile, levels = rev(d_labels))]
p12a <- ggplot(prof, aes(x = bin, y = mean_beta, color = decile)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = DCOLS, name = "Expression\ngroup") +
  geom_vline(xintercept = 0, color = "brown", linewidth = 0.6) +
  labs(title = "CpG Methylation Around TSS by Expression Decile (3+ exon genes)",
       subtitle = sprintf("10 bp bins | n = %s genes", format(nrow(gi), big.mark = ",")),
       x = "Distance to TSS (bp)", y = "Mean methylation (%)") +
  theme_minimal(base_size = 12)
save_fig(p12a, BATCH_DIR, "fig12a_meth_vs_tss_full", w = 12, h = 7)

# Fig 12B: Zoomed to promoter region
cat("  Fig 12B: Zoomed promoter region...\n")
p12b <- p12a + coord_cartesian(xlim = c(-1500, 500)) +
  labs(title = "Promoter Methylation Spike (zoomed)",
       subtitle = "10 bp bins | Vertical line = TSS")
save_fig(p12b, BATCH_DIR, "fig12b_meth_vs_tss_zoom", w = 12, h = 7)

# Locate spike: local max in D09+D10 between -1000 and -50
cat("  Locating spike position...\n")
top2 <- prof[decile %in% c("D09", "D10") & bin >= -1000 & bin <= -50,
             .(mean_beta = mean(mean_beta)), by = bin]
# Smooth to find robust peak
top2[, smooth := stats::filter(mean_beta, rep(1/5, 5), sides = 2)]
top2 <- top2[!is.na(smooth)]
spike_bin <- top2$bin[which.max(top2$smooth)]

# Flanking minima
left_of_spike <- top2[bin < spike_bin]
right_of_spike <- top2[bin > spike_bin & bin < 0]
left_min  <- if (nrow(left_of_spike) > 0)  left_of_spike$bin[which.min(left_of_spike$smooth)]  else spike_bin - 200
right_min <- if (nrow(right_of_spike) > 0) right_of_spike$bin[which.min(right_of_spike$smooth)] else spike_bin + 200

spike_info <- data.table(
  spike_center = spike_bin,
  spike_left_min = left_min,
  spike_right_min = right_min,
  spike_height_pct = top2$smooth[top2$bin == spike_bin],
  left_min_pct = top2$smooth[top2$bin == left_min],
  right_min_pct = top2$smooth[top2$bin == right_min]
)
spike_info[, amplitude := spike_height_pct - (left_min_pct + right_min_pct) / 2]
save_data(spike_info, BATCH_DIR, "spike_parameters")
cat(sprintf("  Spike center: %d bp | Left min: %d bp | Right min: %d bp | Amplitude: %.1f%%\n",
    spike_bin, left_min, right_min, spike_info$amplitude))

SPIKE_L <- left_min
SPIKE_R <- right_min

# =============================================================================
# [3/8] CpG DENSITY AROUND TSS
# =============================================================================
cat("\n[3/8] CpG density around TSS...\n")

# Count CpGs per 100bp sliding window (step=10bp) around each gene's TSS
# Use the overlap data we already computed
cpg_density <- ov_dt[, .(n_cpgs = .N), by = .(gene_row, bin)]
cpg_density[, decile := gi$expr_decile[gene_row]]
cpg_dens_agg <- cpg_density[, .(mean_cpg_per_bin = mean(n_cpgs)), by = .(decile, bin)]
# Convert to CpGs per kb
cpg_dens_agg[, cpg_per_kb := mean_cpg_per_bin / BIN_W * 1000]
cpg_dens_agg[, decile := factor(decile, levels = rev(d_labels))]

cat("  Fig 12C: CpG density around TSS...\n")
p12c <- ggplot(cpg_dens_agg, aes(x = bin, y = cpg_per_kb, color = decile)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = DCOLS, name = "Expression\ngroup") +
  geom_vline(xintercept = 0, color = "brown", linewidth = 0.6) +
  annotate("rect", xmin = SPIKE_L, xmax = SPIKE_R, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "red") +
  coord_cartesian(xlim = c(-1500, 500)) +
  labs(title = "CpG Density Around TSS by Expression Decile",
       subtitle = "CpGs per kb in 10 bp bins | Red shading = spike region",
       x = "Distance to TSS (bp)", y = "CpG density (CpGs / kb)") +
  theme_minimal(base_size = 12)
save_fig(p12c, BATCH_DIR, "fig12c_cpg_density_tss", w = 12, h = 7)
save_data(cpg_dens_agg, BATCH_DIR, "cpg_density_by_tss_distance")

# =============================================================================
# [4/8] SEQUENCE COMPOSITION AT THE SPIKE
# =============================================================================
cat("\n[4/8] Sequence composition around TSS...\n")

genome <- load_genome()

# Extract TSS-flanking sequences for 3+ exon genes (strand-aware)
# Use Views on genome for efficiency
SEQWIN <- 2000L  # +/- 2kb for sequence analysis
seq_starts <- ifelse(gi$str == "+", gi$tss - SEQWIN, gi$tss - SEQWIN)
seq_ends   <- ifelse(gi$str == "+", gi$tss + SEQWIN, gi$tss + SEQWIN)
# Trim to chromosome bounds
chr_lens <- seqlengths(bs_obj)[gi$chr]
if (is.null(chr_lens) || all(is.na(chr_lens))) {
  chr_lens <- sapply(gi$chr, function(ch) nchar(genome[[ch]]))
}
seq_starts <- pmax(seq_starts, 1L)
seq_ends   <- pmin(seq_ends, chr_lens)

cat("  Extracting promoter sequences...\n")
seqs <- mapply(function(ch, s, e, strand) {
  sq <- subseq(genome[[ch]], s, e)
  if (strand == "-") sq <- reverseComplement(sq)
  sq
}, gi$chr, seq_starts, seq_ends, gi$str, SIMPLIFY = FALSE)
seqs <- DNAStringSet(seqs)
names(seqs) <- gi$gene_id
cat(sprintf("  Extracted %d sequences (median length: %d bp)\n",
    length(seqs), median(width(seqs))))

# Sliding window GC% and CpG O/E
SW <- 100L; STEP <- 10L
max_len <- min(width(seqs))
positions <- seq(1, max_len - SW + 1, by = STEP)
# Center positions relative to TSS: SEQWIN is center
pos_relative <- positions + SW/2 - SEQWIN - 1

cat("  Computing GC% and CpG O/E in sliding windows...\n")
gc_mat <- matrix(NA_real_, nrow = length(seqs), ncol = length(positions))
oe_mat <- matrix(NA_real_, nrow = length(seqs), ncol = length(positions))

# Vectorize: extract all windows at once using Views
for (j in seq_along(positions)) {
  win_seqs <- subseq(seqs, start = positions[j], width = SW)
  af <- alphabetFrequency(win_seqs, as.prob = TRUE)
  gc_mat[, j] <- af[, "G"] + af[, "C"]
  df <- dinucleotideFrequency(win_seqs, as.prob = TRUE)
  p_c <- af[, "C"]; p_g <- af[, "G"]
  expected_cg <- p_c * p_g
  observed_cg <- df[, "CG"]
  oe_mat[, j] <- ifelse(expected_cg > 0, observed_cg / expected_cg, NA_real_)
}

# Aggregate by decile
gi_decile <- gi$expr_decile
seq_comp <- rbindlist(lapply(seq_along(positions), function(j) {
  data.table(
    pos = pos_relative[j],
    decile = gi_decile,
    gc = gc_mat[, j],
    oe = oe_mat[, j]
  )
}))
seq_agg <- seq_comp[, .(gc_pct = mean(gc, na.rm = TRUE) * 100,
                         cpg_oe = mean(oe, na.rm = TRUE)), by = .(decile, pos)]
seq_agg[, decile := factor(decile, levels = rev(d_labels))]
save_data(seq_agg, BATCH_DIR, "sequence_composition_around_tss")

cat("  Fig 12D: GC% around TSS...\n")
p12d <- ggplot(seq_agg, aes(x = pos, y = gc_pct, color = decile)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = DCOLS, name = "Expression\ngroup") +
  geom_vline(xintercept = 0, color = "brown", linewidth = 0.6) +
  annotate("rect", xmin = SPIKE_L, xmax = SPIKE_R, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "red") +
  coord_cartesian(xlim = c(-1500, 500)) +
  labs(title = "GC Content Around TSS by Expression Decile",
       subtitle = "100 bp sliding windows, step 10 bp | Red = spike region",
       x = "Distance to TSS (bp)", y = "GC (%)") +
  theme_minimal(base_size = 12)
save_fig(p12d, BATCH_DIR, "fig12d_gc_around_tss", w = 12, h = 7)

cat("  Fig 12E: CpG O/E around TSS...\n")
p12e <- ggplot(seq_agg, aes(x = pos, y = cpg_oe, color = decile)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = DCOLS, name = "Expression\ngroup") +
  geom_vline(xintercept = 0, color = "brown", linewidth = 0.6) +
  annotate("rect", xmin = SPIKE_L, xmax = SPIKE_R, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "red") +
  coord_cartesian(xlim = c(-1500, 500)) +
  labs(title = "CpG O/E Around TSS by Expression Decile",
       subtitle = "100 bp sliding windows | Red = spike region",
       x = "Distance to TSS (bp)", y = "CpG O/E") +
  theme_minimal(base_size = 12)
save_fig(p12e, BATCH_DIR, "fig12e_cpg_oe_around_tss", w = 12, h = 7)

rm(gc_mat, oe_mat, seq_comp); gc(verbose = FALSE)

# =============================================================================
# [5/8] CORE PROMOTER ELEMENT SEARCH
# =============================================================================
cat("\n[5/8] Core promoter element search...\n")

# Extract narrow TSS region (-100 to +50) for motif scanning
narrow_start <- pmax(1L, SEQWIN + 1 - 100)  # -100 relative to TSS in the extracted sequence
narrow_end   <- pmin(width(seqs), SEQWIN + 1 + 50)  # +50
narrow_seqs  <- subseq(seqs, start = narrow_start, end = narrow_end)

# Define core promoter elements using IUPAC ambiguity codes
# Biostrings vmatchPattern with fixed=FALSE supports IUPAC natively
cpe_patterns <- list(
  TATA = list(pattern = DNAString("TATAAA"), canonical_pos = -30),
  Inr  = list(pattern = DNAString("YYANWYY"), canonical_pos = 0),
  DPE  = list(pattern = DNAString("RGWYV"), canonical_pos = 30),
  BRE  = list(pattern = DNAString("SSRCGCC"), canonical_pos = -35)
)

# Scan each element
cpe_results <- rbindlist(lapply(names(cpe_patterns), function(motif_name) {
  pat <- cpe_patterns[[motif_name]]$pattern
  canon <- cpe_patterns[[motif_name]]$canonical_pos
  hits_list <- vmatchPattern(pat, narrow_seqs, fixed = FALSE)
  has_hit <- lengths(hits_list) > 0
  data.table(gene_id = gi$gene_id, motif = motif_name, present = has_hit,
             canonical_pos = canon)
}))
cpe_results <- merge(cpe_results, gi[, .(gene_id, expr_decile)], by = "gene_id")
save_data(cpe_results, BATCH_DIR, "core_promoter_elements")

# Enrichment per decile
cpe_enrich <- cpe_results[, .(
  n_with = sum(present), n_total = .N,
  pct = round(100 * mean(present), 1)
), by = .(motif, expr_decile)]
save_data(cpe_enrich, BATCH_DIR, "core_promoter_element_by_decile")

# Fisher: top (D09+D10) vs bottom (D01+D02)
cpe_fisher <- rbindlist(lapply(unique(cpe_results$motif), function(m) {
  top <- cpe_results[motif == m & expr_decile %in% c("D09", "D10")]
  bot <- cpe_results[motif == m & expr_decile %in% c("D01", "D02")]
  mat <- matrix(c(sum(top$present), sum(!top$present),
                   sum(bot$present), sum(!bot$present)), nrow = 2)
  ft <- fisher.test(mat)
  data.table(motif = m, top_pct = 100 * mean(top$present),
             bot_pct = 100 * mean(bot$present),
             OR = ft$estimate, p = ft$p.value)
}))
cpe_fisher[, padj := p.adjust(p, method = "BH")]
cat("  Core promoter element enrichment (top vs bottom decile):\n")
print(cpe_fisher)
save_data(cpe_fisher, BATCH_DIR, "core_promoter_element_enrichment")

cat("  Fig 12F: Core promoter elements by decile...\n")
cpe_enrich[, expr_decile := factor(expr_decile, levels = d_labels)]
p12f <- ggplot(cpe_enrich, aes(x = expr_decile, y = pct, fill = motif)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = "Element") +
  labs(title = "Core Promoter Element Frequency by Expression Decile",
       subtitle = "Scanned -100 to +50 relative to TSS | 3+ exon genes",
       x = "Expression decile", y = "Genes with element (%)") +
  theme_minimal(base_size = 12)
save_fig(p12f, BATCH_DIR, "fig12f_core_promoter_elements", w = 10, h = 6)

rm(narrow_seqs); gc(verbose = FALSE)

# =============================================================================
# [6/8] TFBS MOTIF ENRICHMENT AT SPIKE POSITION
# =============================================================================
cat("\n[6/8] TFBS motif enrichment at spike...\n")

motif_file <- file.path(PIPE_DIR, "batch1.5/data/motif_hits_extended.tsv.gz")
if (file.exists(motif_file)) {
  mh <- fread(motif_file)
  cat(sprintf("  Motif hits loaded: %s\n", format(nrow(mh), big.mark = ",")))

  # Filter to promoter hits for our genes
  mh <- mh[gene_id %in% gi$gene_id & region_type == "promoter"]
  cat(sprintf("  Promoter hits for 3+ exon genes: %s\n", format(nrow(mh), big.mark = ",")))

  if (nrow(mh) > 0 && "dist_to_tss" %in% names(mh)) {
    # Classify: in spike vs flanking
    mh[, in_spike := dist_to_tss >= SPIKE_L & dist_to_tss <= SPIKE_R]
    flank_size <- SPIKE_R - SPIKE_L
    mh[, in_flank := (dist_to_tss >= (SPIKE_L - flank_size) & dist_to_tss < SPIKE_L) |
                     (dist_to_tss > SPIKE_R & dist_to_tss <= (SPIKE_R + flank_size))]

    # Per-motif enrichment
    motif_counts <- mh[in_spike | in_flank, .(
      n_spike = sum(in_spike), n_flank = sum(in_flank)
    ), by = motif_name]
    motif_counts <- motif_counts[n_spike + n_flank >= 10]

    if (nrow(motif_counts) > 0) {
      motif_counts[, total_spike := sum(n_spike)]
      motif_counts[, total_flank := sum(n_flank)]
      motif_fisher <- motif_counts[, {
        mat <- matrix(c(n_spike, total_spike - n_spike,
                         n_flank, total_flank - n_flank), nrow = 2)
        ft <- fisher.test(mat)
        .(OR = ft$estimate, p = ft$p.value, n_spike = n_spike, n_flank = n_flank)
      }, by = motif_name]
      motif_fisher[, padj := p.adjust(p, method = "BH")]
      motif_fisher[, log2OR := log2(pmax(OR, 0.01))]
      motif_fisher[, neglog10p := -log10(pmax(padj, 1e-300))]
      motif_fisher <- motif_fisher[order(p)]
      save_data(motif_fisher, BATCH_DIR, "tfbs_spike_enrichment")

      cat(sprintf("  Significant motifs (padj<0.05): %d / %d\n",
          sum(motif_fisher$padj < 0.05), nrow(motif_fisher)))

      cat("  Fig 12G: TFBS enrichment at spike...\n")
      top_label <- motif_fisher[1:min(15, nrow(motif_fisher))]
      p12g <- ggplot(motif_fisher, aes(x = log2OR, y = neglog10p)) +
        geom_point(aes(size = n_spike), alpha = 0.5, color = "#2471A3") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
        ggrepel::geom_text_repel(data = top_label, aes(label = motif_name),
                                  size = 3, max.overlaps = 20) +
        scale_size_continuous(range = c(1, 6), name = "Hits in\nspike") +
        labs(title = "TFBS Motif Enrichment at Promoter Spike",
             subtitle = sprintf("Spike region: %d to %d bp from TSS | Fisher's exact vs flanking",
                                SPIKE_L, SPIKE_R),
             x = "log2(Odds Ratio)", y = "-log10(adjusted p-value)") +
        theme_minimal(base_size = 12)
      save_fig(p12g, BATCH_DIR, "fig12g_tfbs_spike_enrichment", w = 10, h = 7)

      # Fig 12H: Positional density of top enriched motifs
      cat("  Fig 12H: Top motif positional density...\n")
      top_motifs <- motif_fisher[padj < 0.05][1:min(6, sum(padj < 0.05)), motif_name]
      if (length(top_motifs) > 0) {
        mh_top <- mh[motif_name %in% top_motifs & gene_id %in%
                      gi[expr_decile %in% c("D09", "D10"), gene_id]]
        if (nrow(mh_top) > 0) {
          p12h <- ggplot(mh_top, aes(x = dist_to_tss, fill = motif_name)) +
            geom_histogram(binwidth = 20, alpha = 0.7, position = "identity") +
            facet_wrap(~ motif_name, scales = "free_y", ncol = 2) +
            annotate("rect", xmin = SPIKE_L, xmax = SPIKE_R, ymin = -Inf, ymax = Inf,
                     alpha = 0.15, fill = "red") +
            geom_vline(xintercept = 0, color = "brown") +
            coord_cartesian(xlim = c(-1500, 500)) +
            labs(title = "Top Enriched Motifs: Position Around TSS (D09+D10 genes)",
                 x = "Distance to TSS (bp)", y = "Motif hit count") +
            theme_minimal(base_size = 11) + theme(legend.position = "none")
          save_fig(p12h, BATCH_DIR, "fig12h_top_motif_positions", w = 10, h = 8)
        }
      }
    } else {
      cat("  Not enough motif hits in spike region for enrichment test.\n")
    }
  } else {
    cat("  No dist_to_tss column or no promoter hits — skipping TFBS analysis.\n")
  }
  rm(mh); gc(verbose = FALSE)
} else {
  cat("  batch1.5 motif hits not found — skipping TFBS analysis.\n")
}

# =============================================================================
# [7/8] NUCLEOSOME POSITIONING PREDICTION
# =============================================================================
cat("\n[7/8] Nucleosome positioning prediction...\n")

# WW dinucleotide (AA/TT/AT/TA) frequency in 147bp sliding windows
# High WW = nucleosome-favoring sequence
NUC_W <- 147L; NUC_STEP <- 10L
nuc_positions <- seq(1, min(width(seqs)) - NUC_W + 1, by = NUC_STEP)
nuc_pos_rel <- nuc_positions + NUC_W / 2 - SEQWIN - 1

cat("  Computing WW dinucleotide frequency...\n")
ww_mat <- matrix(NA_real_, nrow = length(seqs), ncol = length(nuc_positions))
for (j in seq_along(nuc_positions)) {
  win_seqs <- subseq(seqs, start = nuc_positions[j], width = NUC_W)
  df <- dinucleotideFrequency(win_seqs, as.prob = TRUE)
  ww_mat[, j] <- df[, "AA"] + df[, "TT"] + df[, "AT"] + df[, "TA"]
}

nuc_dt <- rbindlist(lapply(seq_along(nuc_positions), function(j) {
  data.table(pos = nuc_pos_rel[j], decile = gi$expr_decile, ww = ww_mat[, j])
}))
nuc_agg <- nuc_dt[, .(ww_score = mean(ww, na.rm = TRUE)), by = .(decile, pos)]
nuc_agg[, decile := factor(decile, levels = rev(d_labels))]
save_data(nuc_agg, BATCH_DIR, "nucleosome_prediction_by_tss_distance")

# Overlay with methylation for top decile
cat("  Fig 12I: Nucleosome prediction + methylation overlay...\n")
# Merge nuc and meth for D09+D10
nuc_top <- nuc_agg[decile %in% c("D09", "D10"), .(ww_score = mean(ww_score)), by = pos]
meth_top <- prof[decile %in% c("D09", "D10"), .(mean_beta = mean(mean_beta)), by = bin]

# Two-panel plot
p_nuc <- ggplot(nuc_top, aes(x = pos, y = ww_score)) +
  geom_line(color = "#8E44AD", linewidth = 0.8) +
  geom_vline(xintercept = 0, color = "brown", linewidth = 0.6) +
  annotate("rect", xmin = SPIKE_L, xmax = SPIKE_R, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "red") +
  coord_cartesian(xlim = c(-1500, 500)) +
  labs(y = "WW dinuc. freq.\n(nucleosome score)") +
  theme_minimal(base_size = 11) + theme(axis.title.x = element_blank())

p_meth <- ggplot(meth_top, aes(x = bin, y = mean_beta)) +
  geom_line(color = "#C0392B", linewidth = 0.8) +
  geom_vline(xintercept = 0, color = "brown", linewidth = 0.6) +
  annotate("rect", xmin = SPIKE_L, xmax = SPIKE_R, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "red") +
  coord_cartesian(xlim = c(-1500, 500)) +
  labs(x = "Distance to TSS (bp)", y = "Methylation (%)") +
  theme_minimal(base_size = 11)

p12i <- patchwork::wrap_plots(p_nuc, p_meth, ncol = 1) +
  patchwork::plot_annotation(
    title = "Nucleosome Positioning and Methylation Around TSS (D09+D10)",
    subtitle = "Top: WW dinucleotide score (higher = nucleosome-favoring) | Bottom: CpG methylation | Red = spike",
    theme = theme(plot.title = element_text(size = 14, face = "bold")))
save_fig(p12i, BATCH_DIR, "fig12i_nucleosome_vs_methylation", w = 12, h = 8)

rm(ww_mat, nuc_dt, seqs); gc(verbose = FALSE)

# =============================================================================
# [8/8] CONTROL vs AMPUTATED AT SPIKE + SUMMARY
# =============================================================================
cat("\n[8/8] Control vs amputated comparison + summary...\n")

# Fig 12J: ctrl vs ampu for D09+D10
top_prof <- ov_dt[decile %in% c("D09", "D10"),
  .(ctrl = mean(beta) * 100, ampu = mean(beta_a) * 100), by = bin]
top_long <- melt(top_prof, id.vars = "bin", variable.name = "condition", value.name = "meth")
top_long[, condition := ifelse(condition == "ctrl", "Control", "Amputated")]

p12j <- ggplot(top_long, aes(x = bin, y = meth, color = condition)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = c("Control" = "#2471A3", "Amputated" = "#C0392B")) +
  geom_vline(xintercept = 0, color = "brown", linewidth = 0.6) +
  annotate("rect", xmin = SPIKE_L, xmax = SPIKE_R, ymin = -Inf, ymax = Inf,
           alpha = 0.15, fill = "red") +
  coord_cartesian(xlim = c(-1500, 500)) +
  labs(title = "Promoter Spike: Control vs Amputated (D09+D10 genes)",
       subtitle = "Does the spike change during regeneration?",
       x = "Distance to TSS (bp)", y = "Methylation (%)", color = "Condition") +
  theme_minimal(base_size = 12)
save_fig(p12j, BATCH_DIR, "fig12j_spike_ctrl_vs_ampu", w = 10, h = 6)

# Summary
elapsed <- (proc.time() - t0)[3]
n_figs <- length(list.files(file.path(BATCH_DIR, "figures"), pattern = "\\.png$"))
n_data <- length(list.files(file.path(BATCH_DIR, "data"), pattern = "\\.tsv$"))
cat(sprintf("\n=== Batch 12 complete (%.0f s) ===\n", elapsed))
cat(sprintf("  Figures: %d | Data files: %d\n", n_figs, n_data))
cat(sprintf("  Spike at %d bp from TSS (amplitude: %.1f%%)\n",
    spike_info$spike_center, spike_info$amplitude))
