#!/usr/bin/env Rscript
# =============================================================================
# Batch 04: Baseline Methylation Landscape + Gene Body Meth vs Expression
# Question: Where is methylation? Does gene body meth correlate with expression?
# Output: data/ (6 tables) + figures/ (7 plots)
# Requires: BSseq cache (or builds from CpG reports), DESeq2 counts
# =============================================================================

source("methylation_pipeline/_config.R")
t0 <- proc.time()

suppressPackageStartupMessages({
  library(bsseq)
  library(GenomicRanges)
  library(rtracklayer)
  library(DESeq2)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(scales)
})

BATCH_DIR <- file.path(PIPE_DIR, "batch04")

# Clean old output
unlink(list.files(file.path(BATCH_DIR, "figures"), full.names = TRUE))
unlink(list.files(file.path(BATCH_DIR, "data"), full.names = TRUE))
dir.create(file.path(BATCH_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BATCH_DIR, "data"), showWarnings = FALSE, recursive = TRUE)

cat("=== Batch 04: Baseline Methylation Landscape ===\n\n")

# =============================================================================
# 1. LOAD OR BUILD BSseq
# =============================================================================
cat("[1/12] Loading BSseq...\n")
if (file.exists(CACHE$bsseq)) {
  bs_obj <- readRDS(CACHE$bsseq)
} else {
  cat("  Building BSseq from CpG reports (strand-collapsed, cov>=5)...\n")
  sample_info_bs <- data.frame(
    sample_id = c("C1", "C2", "A1", "A2"),
    condition = c("Control", "Control", "Amputated", "Amputated"),
    file_path = c(OG$cpg_C1, OG$cpg_C2, OG$cpg_A1, OG$cpg_A2))
  min_cov <- 5; tmp_dir <- tempdir(); common_dt <- NULL
  for (i in seq_len(nrow(sample_info_bs))) {
    sn <- sample_info_bs$sample_id[i]
    cat(sprintf("  Reading %s...\n", sn))
    dt <- fread(sample_info_bs$file_path[i], header = FALSE, sep = "\t",
                select = c(1,2,3,4,5,6),
                col.names = c("chr","pos","strand","meth","unmeth","context"))
    dt <- dt[context == "CG" & chr %in% keep_chr]
    dt[, context := NULL]
    dt[strand == "-", pos := pos - 1L]
    dt[, strand := NULL]
    dt <- dt[, .(meth = sum(meth), unmeth = sum(unmeth)), by = .(chr, pos)]
    dt <- dt[meth + unmeth >= min_cov]
    cat(sprintf("    After collapse + cov>=%d: %s\n", min_cov, format(nrow(dt), big.mark=",")))
    saveRDS(dt, file.path(tmp_dir, paste0("cpg_", sn, ".rds")))
    site_dt <- dt[, .(chr, pos)]
    if (is.null(common_dt)) { common_dt <- site_dt } else {
      setkeyv(common_dt, c("chr","pos")); setkeyv(site_dt, c("chr","pos"))
      common_dt <- fintersect(common_dt, site_dt)
    }
    cat(sprintf("    Running common: %s\n", format(nrow(common_dt), big.mark=",")))
    rm(dt, site_dt); gc(verbose = FALSE)
  }
  cat(sprintf("  Common sites (cov>=%d, all 4): %s\n", min_cov, format(nrow(common_dt), big.mark=",")))
  setkeyv(common_dt, c("chr","pos"))
  n_sites <- nrow(common_dt)
  M_mat <- matrix(0L, nrow = n_sites, ncol = 4)
  Cov_mat <- matrix(0L, nrow = n_sites, ncol = 4)
  colnames(M_mat) <- colnames(Cov_mat) <- sample_info_bs$sample_id
  gr <- NULL
  for (i in 1:4) {
    sn <- sample_info_bs$sample_id[i]
    dt <- readRDS(file.path(tmp_dir, paste0("cpg_", sn, ".rds")))
    setkeyv(dt, c("chr","pos")); dt <- dt[common_dt, nomatch = 0L]
    if (is.null(gr)) gr <- GRanges(seqnames = dt$chr, ranges = IRanges(start = dt$pos, width = 1))
    M_mat[, i] <- dt$meth; Cov_mat[, i] <- dt$meth + dt$unmeth
    rm(dt); gc(verbose = FALSE)
    file.remove(file.path(tmp_dir, paste0("cpg_", sn, ".rds")))
  }
  rm(common_dt); gc(verbose = FALSE)
  bs_obj <- BSseq(gr = gr, M = M_mat, Cov = Cov_mat)
  pData(bs_obj)$condition <- sample_info_bs$condition
  pData(bs_obj)$sample_id <- sample_info_bs$sample_id
  rm(gr, M_mat, Cov_mat); gc(verbose = FALSE)
  saveRDS(bs_obj, CACHE$bsseq)
  cat(sprintf("  BSseq cached: %s\n", CACHE$bsseq))
}
cat(sprintf("  BSseq: %s sites x %d samples\n\n", format(nrow(bs_obj), big.mark = ","), ncol(bs_obj)))

# =============================================================================
# 2. ANNOTATIONS
# =============================================================================
cat("[2/12] Loading annotations...\n")
gff <- load_gff()
genes <- gff[gff$type == "gene"]
exons <- gff[gff$type == "exon"]
promoters <- if (file.exists(CACHE$promoters)) readRDS(CACHE$promoters) else {
  p <- GenomicRanges::promoters(genes, upstream = 2000, downstream = 0)
  p <- GenomicRanges::trim(p)
  saveRDS(p, CACHE$promoters); p
}

# =============================================================================
# 3. PER-SAMPLE + PER-REGION METHYLATION
# =============================================================================
cat("[3/12] Computing per-sample and per-region methylation...\n")
ctrl_idx <- which(colnames(bs_obj) %in% c("C1", "C2"))
ampu_idx <- which(colnames(bs_obj) %in% c("A1", "A2"))
ctrl_beta <- rowMeans(getMeth(bs_obj[, ctrl_idx], type = "raw"), na.rm = TRUE)
ampu_beta <- rowMeans(getMeth(bs_obj[, ampu_idx], type = "raw"), na.rm = TRUE)

bs_chr <- as.character(seqnames(bs_obj))
bs_pos <- start(bs_obj)
gr_all <- GRanges(seqnames = bs_chr, ranges = IRanges(start = bs_pos, width = 1))

# Per-sample stats
cat("\n  Per-sample methylation:\n")
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

# Region annotation
regions <- annotate_regions(bs_chr, bs_pos, promoters, exons, genes)
dt_meth <- data.table(region = regions, ctrl = ctrl_beta, ampu = ampu_beta)
dt_meth <- dt_meth[!is.na(ctrl)]

reg_stats <- dt_meth[, .(n_cpgs = .N, mean_ctrl = round(mean(ctrl), 4),
                          median_ctrl = round(median(ctrl), 4),
                          mean_ampu = round(mean(ampu), 4),
                          median_ampu = round(median(ampu), 4)),
                      by = region]
reg_stats$region <- factor(reg_stats$region, levels = c("Promoter", "Exon", "Intron", "Intergenic"))
reg_stats <- reg_stats[order(region)]
cat("\n  Region methylation:\n"); print(reg_stats)
save_data(reg_stats, BATCH_DIR, "region_methylation_stats")

# =============================================================================
# 4. FIGURES — METHYLATION LANDSCAPE
# =============================================================================
cat("\n[4/12] Generating methylation landscape figures...\n")

# Fig 4A: Global methylation per sample bar
p4a <- ggplot(sample_stats, aes(x = sample, y = mean_beta * 100,
                                 fill = ifelse(grepl("^C", sample), "Control", "Amputated"))) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f%%", mean_beta * 100)), vjust = -0.3, size = 3.5) +
  scale_fill_manual(values = c(Control = "#2471A3", Amputated = "#C0392B")) +
  labs(x = "Sample", y = "Mean CpG methylation (%)", fill = "Condition",
       title = "Global CpG methylation per sample") +
  theme_minimal(base_size = 12)
save_fig(p4a, BATCH_DIR, "fig4a_global_methylation_per_sample", w = 7, h = 5)

# Fig 4C: Region methylation bar (control vs amputated)
cat("  Fig 4C: region methylation bar...\n")
reg_plot <- data.table::melt(reg_stats[, .(region, Control = mean_ctrl, Amputated = mean_ampu)],
                 id.vars = "region", variable.name = "condition", value.name = "mean_beta")
p4c <- ggplot(reg_plot, aes(x = region, y = mean_beta * 100, fill = condition)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c(Control = "#2471A3", Amputated = "#C0392B")) +
  labs(title = "Mean CpG methylation by genomic region",
       x = NULL, y = "Mean methylation (%)", fill = "Condition") +
  theme_minimal(base_size = 12)
save_fig(p4c, BATCH_DIR, "fig4c_region_methylation_bar", w = 8, h = 5)

# =============================================================================
# 5. GENE BODY METH vs EXPRESSION
# =============================================================================
cat("[5/12] Computing per-gene methylation + expression...\n")

if (file.exists(CACHE$transcriptome)) {
  cat("  Loading expression from CACHE$transcriptome...\n")
  tc <- readRDS(CACHE$transcriptome)
  vst_mat <- tc$vst_clean
  # Use tail control columns for baseline expression
  ctrl_cols <- grep("Tail_Control|^C\\d+S\\d+", colnames(vst_mat), value = TRUE)
  if (length(ctrl_cols) == 0) ctrl_cols <- colnames(vst_mat)[grepl("control", tc$meta_clean$condition, ignore.case = TRUE) &
                                                              grepl("tail", tc$meta_clean$tissue, ignore.case = TRUE)]
  if (length(ctrl_cols) > 0) {
    expr_vec <- rowMeans(vst_mat[, ctrl_cols, drop = FALSE])
  } else {
    expr_vec <- rowMeans(vst_mat)
  }
  cat(sprintf("  Genes with expression: %d (from %d tail control samples)\n",
              length(expr_vec), length(ctrl_cols)))
  rm(tc, vst_mat); gc(verbose = FALSE)
} else {
  cat("  Cache not found — running DESeq2 for tail controls...\n")
  tail_ctrl_files <- list.files(OG$counts_dir, pattern = "^C\\d+S\\d+_", ignore.case = TRUE)
  cat(sprintf("  Tail control files: %d\n", length(tail_ctrl_files)))
  sample_table <- data.frame(sampleName = tail_ctrl_files, fileName = tail_ctrl_files,
                              condition = "control", stringsAsFactors = FALSE)
  dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table, directory = OG$counts_dir, design = ~ 1)
  dds <- dds[rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds)
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  expr_vec <- rowMeans(assay(vsd))
  cat(sprintf("  Genes with expression: %d\n", length(expr_vec)))
  rm(dds, vsd); gc(verbose = FALSE)
}

# Gene body methylation via findOverlaps
hits <- findOverlaps(gr_all, genes)
gene_dt <- data.table(gene_id = genes$ID[subjectHits(hits)], beta = ctrl_beta[queryHits(hits)])
gene_dt <- gene_dt[!is.na(beta) & !is.nan(beta)]
gene_means <- gene_dt[, .(mean_beta = mean(beta), n_cpgs = .N), by = gene_id]
rm(gene_dt); gc(verbose = FALSE)

gene_means[, expression := expr_vec[match(gene_id, names(expr_vec))]]
gene_means <- gene_means[!is.na(expression)]

cor_test <- cor.test(gene_means$mean_beta, gene_means$expression, method = "spearman")
cat(sprintf("  Gene body meth vs expression: rho = %.3f, p = %s, n = %d\n",
            cor_test$estimate, format(cor_test$p.value, digits = 3), nrow(gene_means)))
save_data(gene_means, BATCH_DIR, "gene_body_meth_vs_expression")

# Per-region: promoter, exon methylation
cat("  Computing per-region methylation per gene...\n")
hits_prom <- findOverlaps(gr_all, promoters)
prom_dt <- data.table(gene_id = promoters$gene_id[subjectHits(hits_prom)],
                       beta = ctrl_beta[queryHits(hits_prom)])
prom_dt <- prom_dt[!is.na(beta)][, .(prom_beta = mean(beta)), by = gene_id]

hits_exon <- findOverlaps(gr_all, exons)
exon_gene <- data.table(exon_idx = seq_along(exons),
  gene_id = genes$ID[findOverlaps(exons, genes, select = "first")])
exon_dt <- data.table(exon_idx = subjectHits(hits_exon), beta = ctrl_beta[queryHits(hits_exon)])
exon_dt <- merge(exon_dt, exon_gene, by = "exon_idx")
exon_dt <- exon_dt[!is.na(beta)][, .(exon_beta = mean(beta)), by = gene_id]

# Downstream: 2kb downstream of TTS
downstream <- GenomicRanges::flank(genes, width = 2000, start = FALSE)
downstream <- GenomicRanges::trim(downstream)
hits_down <- findOverlaps(gr_all, downstream)
down_dt <- data.table(gene_id = genes$ID[subjectHits(hits_down)],
                       beta = ctrl_beta[queryHits(hits_down)])
down_dt <- down_dt[!is.na(beta)][, .(down_beta = mean(beta)), by = gene_id]

# Intron: gene body minus exons
introns <- GenomicRanges::setdiff(genes, exons)
hits_intron <- findOverlaps(gr_all, introns)
# Map introns back to genes
intron_gene_ov <- findOverlaps(introns, genes, select = "first")
intron_dt <- data.table(intron_idx = subjectHits(hits_intron), beta = ctrl_beta[queryHits(hits_intron)])
intron_dt[, gene_id := genes$ID[intron_gene_ov[intron_idx]]]
intron_dt <- intron_dt[!is.na(beta) & !is.na(gene_id)][, .(intron_beta = mean(beta)), by = gene_id]

intron_cor <- cor.test(intron_dt$intron_beta,
  gene_means$expression[match(intron_dt$gene_id, gene_means$gene_id)],
  method = "spearman", use = "complete.obs")
cat(sprintf("  Intron meth vs expression:    rho = %.3f (p = %s)\n",
            intron_cor$estimate, format(intron_cor$p.value, digits = 3)))

# Merge all regions
gene_full <- merge(gene_means, prom_dt, by = "gene_id", all.x = TRUE)
gene_full <- merge(gene_full, exon_dt, by = "gene_id", all.x = TRUE)
gene_full <- merge(gene_full, down_dt, by = "gene_id", all.x = TRUE)
gene_full <- merge(gene_full, intron_dt, by = "gene_id", all.x = TRUE)

prom_cor <- cor.test(gene_full$prom_beta, gene_full$expression, method = "spearman", use = "complete.obs")
exon_cor <- cor.test(gene_full$exon_beta, gene_full$expression, method = "spearman", use = "complete.obs")
down_cor <- cor.test(gene_full$down_beta, gene_full$expression, method = "spearman", use = "complete.obs")
cat(sprintf("  Promoter meth vs expression: rho = %.3f (p = %s)\n",
            prom_cor$estimate, format(prom_cor$p.value, digits = 3)))
cat(sprintf("  Exon meth vs expression:     rho = %.3f (p = %s)\n",
            exon_cor$estimate, format(exon_cor$p.value, digits = 3)))
cat(sprintf("  Downstream meth vs expression: rho = %.3f (p = %s)\n",
            down_cor$estimate, format(down_cor$p.value, digits = 3)))

save_data(gene_full, BATCH_DIR, "gene_region_meth_vs_expression")

# =============================================================================
# 6. METAGENE PROFILES BY GENE TYPE (9 plots)
# =============================================================================
cat("[6/12] Metagene profiles by gene type...\n")

gene_strand <- as.character(strand(genes))
gene_starts <- start(genes); gene_ends <- end(genes)
gene_tss <- ifelse(gene_strand == "+", gene_starts, gene_ends)
gene_tts <- ifelse(gene_strand == "+", gene_ends, gene_starts)
gene_chr <- as.character(seqnames(genes))

exon_hits <- findOverlaps(exons, genes)
exon_to_gene <- data.table(exon_idx = queryHits(exon_hits), gene_idx = subjectHits(exon_hits),
  exon_start = start(exons)[queryHits(exon_hits)], exon_end = end(exons)[queryHits(exon_hits)])

N_BINS <- 20; FLANK <- 5000

# -- Exon counts per gene for subsetting --
exon_per_gene <- exon_to_gene[, .(n_exons = uniqueN(exon_idx)), by = gene_idx]
genes_with_expr <- which(genes$ID %in% gene_full$gene_id)
idx_3plus <- intersect(exon_per_gene[n_exons >= 3, gene_idx], genes_with_expr)
idx_2exon <- intersect(exon_per_gene[n_exons == 2, gene_idx], genes_with_expr)
idx_1exon <- intersect(exon_per_gene[n_exons == 1, gene_idx], genes_with_expr)
idx_xloc  <- intersect(which(grepl("^XLOC", genes$ID)), genes_with_expr)

cat(sprintf("  Gene counts — all: %d, 3+ exon: %d, 2-exon: %d, single-exon: %d, XLOC/ncRNA: %d\n",
    length(genes_with_expr), length(idx_3plus), length(idx_2exon),
    length(idx_1exon), length(idx_xloc)))

# -- Expression bin builder (4=quartile, 5=quintile, 10=decile) --
# Legend always sorted High -> Low (top-of-legend = highest expression group),
# matching the visual top-to-bottom order of the lines in the plot.
make_expr_bins <- function(dt, n_bins) {
  d <- copy(dt)
  probs  <- seq(0, 1, length.out = n_bins + 1)
  prefix <- if (n_bins == 4) "Q" else if (n_bins == 5) "Qu" else "D"
  raw_labels <- if (n_bins == 4) {
    c("Low (Q1)", "Mid-low (Q2)", "Mid-high (Q3)", "High (Q4)")
  } else if (n_bins == 5) {
    c("Low (Qu1)", "Mid-low (Qu2)", "Mid (Qu3)", "Mid-high (Qu4)", "High (Qu5)")
  } else {
    paste0("D", sprintf("%02d", 1:10))   # D01 = lowest, D10 = highest
  }
  d[, expr_group := cut(expression,
                        breaks = quantile(expression, probs = probs, na.rm = TRUE),
                        labels = raw_labels,
                        include.lowest = TRUE)]
  # Reverse the factor levels so the legend prints High first
  d[, expr_group := factor(expr_group, levels = rev(raw_labels))]
  d
}

# Color palette: red (high) -> green (low), descending. Legend top = High.
make_expr_colors <- function(n_bins) {
  raw <- if (n_bins == 4) {
    c("Low (Q1)" = "#27AE60", "Mid-low (Q2)" = "#2471A3",
      "Mid-high (Q3)" = "#F39C12", "High (Q4)" = "#C0392B")
  } else if (n_bins == 5) {
    c("Low (Qu1)" = "#27AE60", "Mid-low (Qu2)" = "#2E86C1",
      "Mid (Qu3)" = "#F1C40F", "Mid-high (Qu4)" = "#E67E22",
      "High (Qu5)" = "#C0392B")
  } else {
    pal <- grDevices::colorRampPalette(c("#27AE60", "#2471A3", "#F1C40F",
                                          "#E67E22", "#C0392B"))(10)
    setNames(pal, paste0("D", sprintf("%02d", 1:10)))
  }
  raw
}

set.seed(42)
# NEVER subsample — use all genes for every metagene. Memory: feedback_no_subsampling.
cap_sample <- function(idx, n = NA) idx

# -- Segment orders (promoter 2kb + distal 3kb split) --
SEG_4 <- c("Distal upstream\n(5-2 kb)", "Promoter\n(2 kb)", "Gene body", "TTS downstream\n(5 kb)")
SEG_6 <- c("Distal upstream\n(5-2 kb)", "Promoter\n(2 kb)", "First exon", "Intron", "Last exon", "TTS downstream\n(5 kb)")
SEG_8 <- c("Distal upstream\n(5-2 kb)", "Promoter\n(2 kb)", "First exon", "First intron", "Body",
           "Last intron", "Last exon", "TTS downstream\n(5 kb)")

# -- Generalized segment builder (mode: "3seg", "5seg", "7seg") --
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

# -- Compute metagene from segments --
run_metagene <- function(beta_vec, seg_result, seg_order) {
  seg_dt <- seg_result$segments
  seg_gr <- GRanges(seqnames = seg_dt$chr,
                     ranges = IRanges(start = seg_dt$start, end = seg_dt$end))
  cpg_gr <- GRanges(seqnames = bs_chr, ranges = IRanges(start = bs_pos, width = 1))
  ov <- findOverlaps(cpg_gr, seg_gr)

  ov_dt <- data.table(
    cpg_idx = queryHits(ov), seg_idx = subjectHits(ov),
    beta = beta_vec[queryHits(ov)], pos = bs_pos[queryHits(ov)])
  ov_dt <- ov_dt[!is.na(beta)]
  ov_dt[, `:=`(seg_start = seg_dt$start[seg_idx], seg_end = seg_dt$end[seg_idx],
               seg_strand = seg_dt$strand[seg_idx], segment = seg_dt$segment[seg_idx],
               expr_group = seg_dt$expr_group[seg_idx],
               gene_id = seg_dt$gene_id[seg_idx])]

  ov_dt[, rel_pos := (pos - seg_start) / (seg_end - seg_start)]
  ov_dt[seg_strand == "-" & !segment %in% c("Body", "Gene body"), rel_pos := 1 - rel_pos]
  ov_dt[, bin := pmin(N_BINS, pmax(1L, ceiling(rel_pos * N_BINS)))]

  gene_bins <- ov_dt[, .(mean_beta = mean(beta) * 100),
                     by = .(gene_id, segment, bin, expr_group)]
  agg <- gene_bins[, .(mean_meth = mean(mean_beta), se = sd(mean_beta) / sqrt(.N)),
                    by = .(segment, bin, expr_group)]
  agg[, segment := factor(segment, levels = seg_order)]
  agg[, x_pos := (as.integer(segment) - 1) * N_BINS + bin]
  rm(ov_dt, gene_bins, ov); gc(verbose = FALSE)
  agg
}

# -- Plot metagene (adapts to any segment count) --
plot_metagene <- function(mg_agg, title, subtitle, colors, seg_order) {
  n_segs <- length(seg_order)
  seg_breaks <- (1:(n_segs - 1)) * N_BINS + 0.5
  seg_labels_pos <- (0:(n_segs - 1)) * N_BINS + N_BINS / 2

  ggplot(mg_agg, aes(x = x_pos, y = mean_meth, color = expr_group)) +
    geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 1) +
    geom_vline(xintercept = seg_breaks, linetype = "dashed", color = "gray70", linewidth = 0.3) +
    scale_color_manual(values = colors) +
    scale_x_continuous(breaks = seg_labels_pos, labels = seg_order) +
    labs(title = title, subtitle = subtitle,
         x = NULL, y = "Methylation level (%)", color = "Expression\ngroup") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(size = 8), panel.grid.minor = element_blank())
}

# -- Generate metagene plots: 5 gene classes x (ampu, ctrl) x 3 binnings (4/5/10)
# Amputated comes BEFORE control in the numbering — show experimental result first.
# xloc_ncrna only has control (no amputated paired data).
metagene_specs <- list(
  list(name = "all_genes",    label = "All genes",                idx = genes_with_expr, mode = "3seg", seg = SEG_4, ampu = TRUE),
  list(name = "3plus_exon",   label = "3+ exon genes",            idx = idx_3plus,       mode = "7seg", seg = SEG_8, ampu = TRUE),
  list(name = "2exon",        label = "2-exon genes",             idx = idx_2exon,       mode = "5seg", seg = SEG_6, ampu = TRUE),
  list(name = "single_exon",  label = "Single-exon (intronless)", idx = idx_1exon,       mode = "3seg", seg = SEG_4, ampu = TRUE),
  list(name = "xloc_ncrna",   label = "XLOC (non-coding RNA)",    idx = idx_xloc,        mode = "3seg", seg = SEG_4, ampu = FALSE)
)

binnings <- list(
  list(n = 4,  tag = "bin04", label = "expression quartiles (4 bins)"),
  list(n = 5,  tag = "bin05", label = "expression quintiles (5 bins)"),
  list(n = 10, tag = "bin10", label = "expression deciles (10 bins)")
)

for (binning in binnings) {
  cat(sprintf("\n  === Metagenes %s ===\n", binning$label))
  gene_bins <- make_expr_bins(gene_full, binning$n)
  colors_b  <- make_expr_colors(binning$n)

  fig_num <- 0L
  for (spec in metagene_specs) {
    gi <- cap_sample(spec$idx)
    if (length(gi) < 30) {
      cat(sprintf("    Skipping %s — too few genes (%d)\n", spec$label, length(gi))); next
    }

    cat(sprintf("    Metagene: %s (n=%d)...\n", spec$label, length(gi)))
    seg_b <- build_segments(gi, gene_bins, mode = spec$mode)
    n_genes <- seg_b$n_genes

    # Amputated FIRST (if requested) — leads the comparison
    if (isTRUE(spec$ampu)) {
      fig_num <- fig_num + 1L
      mg_ampu <- run_metagene(ampu_beta, seg_b, spec$seg)
      p <- plot_metagene(mg_ampu,
        sprintf("Metagene methylation — %s (Amputated, %s)", spec$label, binning$label),
        sprintf("Amputated (A1+A2), n = %s genes", format(n_genes, big.mark = ",")),
        colors_b, spec$seg)
      save_fig(p, BATCH_DIR,
               sprintf("fig4d_%s_%02d_ampu_%s", binning$tag, fig_num, spec$name),
               w = 11, h = 6)
      rm(mg_ampu); gc(verbose = FALSE)
    }

    # Control SECOND (or only, for xloc_ncrna)
    fig_num <- fig_num + 1L
    mg_ctrl <- run_metagene(ctrl_beta, seg_b, spec$seg)
    p <- plot_metagene(mg_ctrl,
      sprintf("Metagene methylation — %s (Control, %s)", spec$label, binning$label),
      sprintf("Controls (C1+C2), n = %s genes", format(n_genes, big.mark = ",")),
      colors_b, spec$seg)
    save_fig(p, BATCH_DIR,
             sprintf("fig4d_%s_%02d_ctrl_%s", binning$tag, fig_num, spec$name),
             w = 11, h = 6)
    rm(mg_ctrl); gc(verbose = FALSE)
    rm(seg_b)
  }
  rm(gene_bins); gc(verbose = FALSE)
}
gc(verbose = FALSE)

# =============================================================================
# 7. METH vs EXPRESSION FIGURES
# =============================================================================
cat("[7/12] Meth vs expression figures...\n")

# Fig 4F: Per-region meth vs expression (promoter, exon, downstream, gene body)
region_cor_df <- data.frame(
  Region = c("Promoter", "Exon", "Intron", "Downstream", "Gene body"),
  rho = c(prom_cor$estimate, exon_cor$estimate, intron_cor$estimate, down_cor$estimate, cor_test$estimate),
  p = c(prom_cor$p.value, exon_cor$p.value, intron_cor$p.value, down_cor$p.value, cor_test$p.value),
  stringsAsFactors = FALSE)
save_data(region_cor_df, BATCH_DIR, "region_meth_expression_correlations")

p4f <- ggplot(region_cor_df, aes(x = reorder(Region, rho), y = rho, fill = rho)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("rho = %.3f", rho)), hjust = -0.1, size = 3.5) +
  scale_fill_gradient2(low = "#2471A3", mid = "white", high = "#C0392B", midpoint = 0, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  coord_flip() +
  labs(title = "Methylation-expression correlation by region",
       x = "Region", y = "Spearman rho") +
  theme_minimal(base_size = 12)
save_fig(p4f, BATCH_DIR, "fig4f_region_meth_expression_correlation", w = 8, h = 5)

# Fig 4G: Methylation by expression decile
cat("  Fig 4G: methylation by expression decile...\n")
gene_means[, expr_decile := cut(expression,
  breaks = quantile(expression, probs = 0:10/10, na.rm = TRUE),
  labels = paste0("D", 1:10), include.lowest = TRUE)]

decile_summary <- gene_means[!is.na(expr_decile), .(
  mean_meth = mean(mean_beta) * 100,
  se = sd(mean_beta) / sqrt(.N) * 100, n = .N
), by = expr_decile]

p4g <- ggplot(decile_summary, aes(x = expr_decile, y = mean_meth, fill = expr_decile)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = mean_meth - se, ymax = mean_meth + se), width = 0.2) +
  scale_fill_manual(values = colorRampPalette(c("#2471A3", "#1ABC9C", "#F39C12", "#C0392B"))(10),
                    guide = "none") +
  labs(x = "Expression decile (D1=lowest, D10=highest)", y = "Mean gene body methylation (%)",
       title = "Gene body methylation by expression decile",
       subtitle = sprintf("n = %s genes total",
                          format(sum(decile_summary$n), big.mark = ","))) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p4g, BATCH_DIR, "fig4g_meth_by_expression_decile", w = 8, h = 6)
save_data(decile_summary, BATCH_DIR, "meth_by_expression_decile")

# Fig 4H: Region-specific methylation by expression decile (boxplots)
cat("  Fig 4H: region-specific methylation by expression decile...\n")
gene_full[, expr_decile := cut(expression,
  breaks = quantile(expression, probs = 0:10/10, na.rm = TRUE),
  labels = paste0("D", 1:10), include.lowest = TRUE)]

region_decile_long <- rbindlist(list(
  gene_full[!is.na(prom_beta) & !is.na(expr_decile),
    .(gene_id, meth = prom_beta * 100, expr_decile, region = "Promoter")],
  gene_full[!is.na(exon_beta) & !is.na(expr_decile),
    .(gene_id, meth = exon_beta * 100, expr_decile, region = "Exon")],
  gene_full[!is.na(intron_beta) & !is.na(expr_decile),
    .(gene_id, meth = intron_beta * 100, expr_decile, region = "Intron")],
  gene_full[!is.na(mean_beta) & !is.na(expr_decile),
    .(gene_id, meth = mean_beta * 100, expr_decile, region = "Gene body")],
  gene_full[!is.na(down_beta) & !is.na(expr_decile),
    .(gene_id, meth = down_beta * 100, expr_decile, region = "Downstream")]
))
region_decile_long[, region := factor(region,
  levels = c("Promoter", "Exon", "Intron", "Gene body", "Downstream"))]

p4h <- ggplot(region_decile_long, aes(x = expr_decile, y = meth, fill = expr_decile)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2) +
  facet_wrap(~ region, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = colorRampPalette(c("#2471A3", "#1ABC9C", "#F39C12", "#C0392B"))(10),
                    guide = "none") +
  labs(x = "Expression decile (D1=lowest, D10=highest)",
       y = "Mean methylation (%)",
       title = "Region-specific methylation by expression decile",
       subtitle = sprintf("n = %s gene-region observations",
                          format(nrow(region_decile_long), big.mark = ","))) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
save_fig(p4h, BATCH_DIR, "fig4h_region_meth_by_decile", w = 14, h = 10)
save_data(region_decile_long[, .(mean_meth = mean(meth), median_meth = median(meth),
  sd_meth = sd(meth), n = .N), by = .(region, expr_decile)],
  BATCH_DIR, "region_meth_by_expression_decile")

# =============================================================================
# 8. GENOME-WIDE METHYLATION BY CHROMOSOME (1 Mb windows)
# =============================================================================
cat("[8/12] Genome-wide methylation by chromosome (1 Mb windows)...\n")
chr_dt <- data.table(chr = bs_chr, pos = bs_pos)
for (si in 1:ncol(bs_obj)) {
  bv <- getMeth(bs_obj[, si], type = "raw")[, 1]
  chr_dt[, (colnames(bs_obj)[si]) := bv]
}
chr_dt[, window := floor(pos / 1e6) * 1e6 + 5e5]  # 1 Mb window midpoints
chr_dt <- chr_dt[chr %in% keep_chr]

# Melt to long format for all 4 samples
win_long <- data.table::melt(chr_dt, id.vars = c("chr", "pos", "window"),
                              measure.vars = colnames(bs_obj),
                              variable.name = "sample", value.name = "beta")
win_long <- win_long[!is.na(beta)]
win_summary <- win_long[, .(mean_beta = mean(beta) * 100, n_cpgs = .N),
                         by = .(chr, window, sample)]
win_summary[, condition := ifelse(sample %in% c("C1", "C2"), "Control", "Amputated")]

# Plot all chromosomes — scatter with all 4 samples
p4i <- ggplot(win_summary, aes(x = window / 1e6, y = mean_beta, color = sample)) +
  geom_point(size = 0.3, alpha = 0.4) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, linewidth = 0.5) +
  facet_wrap(~ chr, scales = "free_x", ncol = 4) +
  scale_color_manual(values = c(C1 = "#2471A3", C2 = "#5DADE2",
                                A1 = "#C0392B", A2 = "#EC7063")) +
  labs(x = "Position (Mb)", y = "Mean methylation (%)",
       title = "Genome-wide CpG methylation (1 Mb windows)",
       subtitle = "All 4 samples, all 31 chromosomes", color = "Sample") +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(size = 7), panel.grid.minor = element_blank())
save_fig(p4i, BATCH_DIR, "fig4i_genomewide_methylation_1mb", w = 16, h = 20)
save_data(win_summary, BATCH_DIR, "genomewide_methylation_1mb_windows")
rm(chr_dt, win_long, win_summary); gc(verbose = FALSE)

# =============================================================================
# 9. METHYLATION PIE CHART BY REGION (control and amputated)
# =============================================================================
cat("[9/12] Methylation pie chart by region...\n")
region_counts <- dt_meth[, .N, by = region]
region_counts[, pct := round(100 * N / sum(N), 1)]
setorder(region_counts, -N)
region_counts[, region := factor(region, levels = region)]
region_counts[, label := sprintf("%s: %s (%.1f%%)", region,
                                  format(N, big.mark = ","), pct)]

# Replaced overlapping pie with horizontal bar (labels readable, sorted by count).
p4j <- ggplot(region_counts, aes(x = N, y = region, fill = region)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = label), hjust = -0.05, size = 3.4) +
  scale_fill_manual(values = COLORS$region, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(title = "CpG distribution by genomic region",
       subtitle = "Sorted by count; pie replaced by bar to avoid label overlap",
       x = "CpG count", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank())
save_fig(p4j, BATCH_DIR, "fig4j_cpg_region_pie", w = 9, h = 6)

# Mean methylation by region per condition (bar version already exists as fig4c, add pie)
ctrl_pie <- dt_meth[, .(mean_meth = round(mean(ctrl) * 100, 1)), by = region]
ampu_pie <- dt_meth[, .(mean_meth = round(mean(ampu) * 100, 1)), by = region]
ctrl_pie[, condition := "Control"]; ampu_pie[, condition := "Amputated"]
pie_both <- rbind(ctrl_pie, ampu_pie)
save_data(pie_both, BATCH_DIR, "region_methylation_by_condition")

# =============================================================================
# 10. (Single-exon metagene now produced in step 6 above)
# =============================================================================
cat("[10/12] Single-exon metagene — already generated in step 6.\n")
expanded <- build_expanded_regions(gff)

# =============================================================================
# 11. FIRST EXON vs LAST EXON BASELINE CORRELATION WITH EXPRESSION
# =============================================================================
cat("[11/12] First exon vs last exon meth-expression correlation...\n")

# First exon methylation per gene
hits_fe <- findOverlaps(gr_all, expanded$first_exons)
fe_dt <- data.table(gene_id = expanded$first_exons$gene_id[subjectHits(hits_fe)],
                     beta = ctrl_beta[queryHits(hits_fe)])
fe_dt <- fe_dt[!is.na(beta)][, .(first_exon_beta = mean(beta)), by = gene_id]

# Last exon methylation per gene
hits_le <- findOverlaps(gr_all, expanded$last_exons)
le_dt <- data.table(gene_id = expanded$last_exons$gene_id[subjectHits(hits_le)],
                     beta = ctrl_beta[queryHits(hits_le)])
le_dt <- le_dt[!is.na(beta)][, .(last_exon_beta = mean(beta)), by = gene_id]

# Single exon methylation per gene
hits_se <- findOverlaps(gr_all, expanded$single_exons)
se_dt <- data.table(gene_id = expanded$single_exons$gene_id[subjectHits(hits_se)],
                     beta = ctrl_beta[queryHits(hits_se)])
se_dt <- se_dt[!is.na(beta)][, .(single_exon_beta = mean(beta)), by = gene_id]

# First intron methylation per gene
hits_fi <- findOverlaps(gr_all, expanded$first_introns)
fi_dt <- data.table(gene_id = expanded$first_introns$gene_id[subjectHits(hits_fi)],
                     beta = ctrl_beta[queryHits(hits_fi)])
fi_dt <- fi_dt[!is.na(beta)][, .(first_intron_beta = mean(beta)), by = gene_id]

# Merge all expanded regions with expression
gene_expanded <- copy(gene_means[, .(gene_id, expression)])
gene_expanded <- merge(gene_expanded, fe_dt, by = "gene_id", all.x = TRUE)
gene_expanded <- merge(gene_expanded, le_dt, by = "gene_id", all.x = TRUE)
gene_expanded <- merge(gene_expanded, se_dt, by = "gene_id", all.x = TRUE)
gene_expanded <- merge(gene_expanded, fi_dt, by = "gene_id", all.x = TRUE)

# Compute correlations
expanded_cors <- list()
for (col_name in c("first_exon_beta", "last_exon_beta", "single_exon_beta", "first_intron_beta")) {
  vals <- gene_expanded[[col_name]]
  valid <- !is.na(vals) & !is.na(gene_expanded$expression)
  if (sum(valid) >= 30) {
    ct <- cor.test(vals[valid], gene_expanded$expression[valid], method = "spearman")
    expanded_cors[[col_name]] <- data.frame(
      region = gsub("_beta$", "", col_name),
      rho = ct$estimate, p = ct$p.value, n = sum(valid))
    cat(sprintf("  %s vs expression: rho = %.3f (p = %s, n = %d)\n",
                col_name, ct$estimate, format(ct$p.value, digits = 3), sum(valid)))
  }
}

# Combine with existing region correlations
expanded_cor_df <- do.call(rbind, expanded_cors)
if (nrow(expanded_cor_df) > 0) {
  all_region_cors <- rbind(
    region_cor_df[, c("Region", "rho", "p")],
    setNames(expanded_cor_df[, c("region", "rho", "p")], c("Region", "rho", "p"))
  )
  all_region_cors$Region <- gsub("_", " ", tools::toTitleCase(all_region_cors$Region))
  save_data(all_region_cors, BATCH_DIR, "expanded_region_meth_expression_correlations")

  # Place rho label INSIDE the bar for positive bars, OUTSIDE for negative,
  # so labels never get clipped by the plot edge regardless of rho sign/magnitude.
  all_region_cors$lab     <- sprintf("rho=%.3f", all_region_cors$rho)
  all_region_cors$lab_hj  <- ifelse(all_region_cors$rho >= 0, 1.1, -0.1)
  all_region_cors$lab_col <- ifelse(abs(all_region_cors$rho) >= 0.15, "white", "black")
  rho_max <- max(abs(all_region_cors$rho), na.rm = TRUE)
  p4l <- ggplot(all_region_cors, aes(x = reorder(Region, rho), y = rho, fill = rho)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = lab, hjust = lab_hj, color = lab_col), size = 3.4) +
    scale_color_identity() +
    scale_fill_gradient2(low = "#2471A3", mid = "white", high = "#C0392B",
                         midpoint = 0, guide = "none") +
    scale_y_continuous(limits = c(-rho_max * 1.25, rho_max * 1.25),
                       expand = expansion(mult = 0.02)) +
    coord_flip() +
    labs(title = "Methylation-expression correlation: expanded regions",
         subtitle = "First/last exon, single exon, first intron vs basic regions",
         x = NULL, y = "Spearman rho") +
    theme_minimal(base_size = 12)
  save_fig(p4l, BATCH_DIR, "fig4l_expanded_region_correlation", w = 11, h = 7)
}
rm(fe_dt, le_dt, se_dt, fi_dt, gene_expanded); gc(verbose = FALSE)

# =============================================================================
# 12. SUMMARY
# =============================================================================
cat(sprintf("\n[12/12] Summary:\n"))
cat(sprintf("  Gene body meth vs expression: rho = %.3f\n", cor_test$estimate))
cat(sprintf("  Promoter meth vs expression:  rho = %.3f\n", prom_cor$estimate))
cat(sprintf("  Exon meth vs expression:      rho = %.3f\n", exon_cor$estimate))
cat(sprintf("  Intron meth vs expression:     rho = %.3f\n", intron_cor$estimate))
cat(sprintf("  Downstream meth vs expression: rho = %.3f\n", down_cor$estimate))

elapsed <- (proc.time() - t0)[3]
cat(sprintf("\n=== Batch 04 complete (%.1f min) ===\n", elapsed / 60))
n_figs <- length(list.files(file.path(BATCH_DIR, "figures"), pattern = "\\.png$"))
n_data <- length(list.files(file.path(BATCH_DIR, "data"), pattern = "\\.tsv$"))
cat(sprintf("Figures: %d | Data files: %d\n", n_figs, n_data))
