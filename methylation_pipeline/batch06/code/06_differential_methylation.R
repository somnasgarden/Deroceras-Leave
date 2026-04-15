#!/usr/bin/env Rscript
# =============================================================================
# Batch 06: Differential Methylation (DMP/DMR)
# Question: What changes during regeneration?
# Output: data/ (DMPs, DMRs, enrichment, GO/KEGG) + figures/ (12+ plots)
# Requires: BSseq + DMLtest cache, STRING enrichment terms
# =============================================================================

source("methylation_pipeline/_config.R")
t0 <- proc.time()

suppressPackageStartupMessages({
  library(bsseq)
  library(DSS)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(pheatmap)
  library(clusterProfiler)
})

BATCH_DIR <- file.path(PIPE_DIR, "batch06")

# Clean old output
unlink(list.files(file.path(BATCH_DIR, "figures"), full.names = TRUE))
unlink(list.files(file.path(BATCH_DIR, "data"), full.names = TRUE))
dir.create(file.path(BATCH_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BATCH_DIR, "data"), showWarnings = FALSE, recursive = TRUE)

cat("=== Batch 06: Differential Methylation ===\n\n")

# =============================================================================
# 1. LOAD DATA
# =============================================================================
cat("[1/10] Loading DMLtest + BSseq...\n")

if (file.exists(CACHE$dmltest)) {
  dml_test <- readRDS(CACHE$dmltest)
} else {
  cat("  No DMLtest cache. Loading BSseq and running DMLtest (takes hours)...\n")
  if (!file.exists(CACHE$bsseq)) stop("BSseq cache not found. Run batch04 first.")
  bs_obj <- readRDS(CACHE$bsseq)
  dml_test <- DMLtest(bs_obj, group1 = c("C1","C2"), group2 = c("A1","A2"), smoothing = TRUE)
  saveRDS(dml_test, CACHE$dmltest)
  rm(bs_obj); gc(verbose = FALSE)
}
cat(sprintf("  DMLtest: %s sites\n", format(nrow(dml_test), big.mark = ",")))

bs_obj <- readRDS(CACHE$bsseq)

# =============================================================================
# 2. CALL DMPs AND DMRs
# =============================================================================
cat("[2/10] Calling DMPs and DMRs...\n")

dmp <- callDML(dml_test, p.threshold = 0.05, delta = 0.1)
dmr <- callDMR(dml_test, p.threshold = 0.05, delta = 0.1, minlen = 50, minCG = 3)
if (is.null(dmr) || nrow(dmr) == 0) {
  dmr <- data.frame(chr = character(), start = integer(), end = integer(),
                     length = integer(), nCG = integer(), areaStat = numeric(),
                     diff.Methy = numeric())
}

n_dmp <- nrow(dmp); n_dmr <- nrow(dmr)
n_hyper <- sum(dmp$diff > 0); n_hypo <- sum(dmp$diff < 0)
cat(sprintf("  DMPs: %s (hyper: %s, hypo: %s)\n", format(n_dmp, big.mark = ","),
            format(n_hyper, big.mark = ","), format(n_hypo, big.mark = ",")))
cat(sprintf("  DMRs: %s\n", format(n_dmr, big.mark = ",")))

# =============================================================================
# 3. ANNOTATE DMPs
# =============================================================================
cat("[3/10] Annotating DMPs...\n")

gff <- load_gff()
genes <- gff[gff$type == "gene"]
exons <- gff[gff$type == "exon"]
promoters <- if (file.exists(CACHE$promoters)) readRDS(CACHE$promoters) else {
  p <- GenomicRanges::trim(GenomicRanges::promoters(genes, upstream = 2000, downstream = 0))
  saveRDS(p, CACHE$promoters); p
}

dmp$annotation <- annotate_regions(dmp$chr, dmp$pos, promoters, exons, genes)
dmp$direction <- ifelse(dmp$diff > 0, "Hyper", "Hypo")

# TE overlap
te_data <- load_te()
dmp_gr <- GRanges(seqnames = dmp$chr, ranges = IRanges(start = dmp$pos, width = 1))
dmp$in_te <- overlapsAny(dmp_gr, te_data$te_gr)

# Nearest gene
nearest_idx <- nearest(dmp_gr, genes)
dmp$nearest_gene <- ifelse(!is.na(nearest_idx), genes$ID[nearest_idx], NA)
dmp$dist_to_gene <- mcols(distanceToNearest(dmp_gr, genes))$distance

# Load gene names from EviAnn annotation (deduplicate: keep first name per gene_id)
gene_annot <- if (file.exists(OG$annot)) {
  ga <- fread(OG$annot, header = FALSE, col.names = c("gene_id", "gene_name", "description"))
  ga[, .SD[1], by = gene_id]
} else { data.table(gene_id = character(), gene_name = character(), description = character()) }

# Add gene names to DMPs
dmp <- merge(as.data.table(dmp), gene_annot[, .(gene_id, gene_name)],
             by.x = "nearest_gene", by.y = "gene_id", all.x = TRUE, sort = FALSE)
cat(sprintf("  DMPs with gene names: %d / %d\n",
            sum(!is.na(dmp$gene_name) & dmp$gene_name != ""), nrow(dmp)))

# Annotate DMRs
if (n_dmr > 0) {
  dmr_gr <- GRanges(seqnames = dmr$chr, ranges = IRanges(start = dmr$start, end = dmr$end))
  dmr_mid <- (dmr$start + dmr$end) %/% 2L
  dmr$annotation <- annotate_regions(dmr$chr, dmr_mid, promoters, exons, genes)
  dmr$direction <- ifelse(dmr$diff.Methy > 0, "Hyper", "Hypo")
  dmr_nearest <- nearest(dmr_gr, genes)
  dmr$nearest_gene <- ifelse(!is.na(dmr_nearest), genes$ID[dmr_nearest], NA)
  dmr <- merge(as.data.table(dmr), gene_annot[, .(gene_id, gene_name)],
               by.x = "nearest_gene", by.y = "gene_id", all.x = TRUE, sort = FALSE)
}

save_data(dmp, BATCH_DIR, "dmps_annotated")
save_data(as.data.frame(dmr), BATCH_DIR, "dmrs_annotated")

# =============================================================================
# 4. REGION ENRICHMENT (Fisher's exact test)
# =============================================================================
cat("[4/10] Computing region enrichment...\n")

# Use dml_test positions for background (chunked for memory on login node)
cat("  Annotating 25M background CpGs (chunked)...\n")
chunk_size <- 5000000L
n_total <- nrow(dml_test)
all_anno <- character(n_total)
for (ci in seq(1, n_total, by = chunk_size)) {
  ce <- min(ci + chunk_size - 1L, n_total)
  all_anno[ci:ce] <- annotate_regions(dml_test$chr[ci:ce], dml_test$pos[ci:ce],
                                       promoters, exons, genes)
  cat(sprintf("    %s / %s\n", format(ce, big.mark = ","), format(n_total, big.mark = ",")))
}
rm(dml_test); gc(verbose = FALSE)  # Free 1.4GB
bg_counts <- as.data.frame(table(all_anno)); colnames(bg_counts) <- c("Region", "BG_Count")
rm(all_anno); gc(verbose = FALSE)
dmp_counts <- as.data.frame(table(dmp$annotation)); colnames(dmp_counts) <- c("Region", "DMP_Count")
enrichment <- merge(dmp_counts, bg_counts, by = "Region")
enrichment$DMP_Pct <- round(100 * enrichment$DMP_Count / sum(enrichment$DMP_Count), 1)
enrichment$BG_Pct  <- round(100 * enrichment$BG_Count / sum(enrichment$BG_Count), 1)
enrichment$Fold <- round(enrichment$DMP_Pct / enrichment$BG_Pct, 2)
enrichment$pvalue <- sapply(1:nrow(enrichment), function(i) {
  fisher.test(matrix(c(enrichment$DMP_Count[i], n_dmp - enrichment$DMP_Count[i],
                        enrichment$BG_Count[i], sum(bg_counts$BG_Count) - enrichment$BG_Count[i]),
                      nrow = 2))$p.value
})
enrichment$sig <- ifelse(enrichment$pvalue < 0.001, "***",
                  ifelse(enrichment$pvalue < 0.01, "**",
                  ifelse(enrichment$pvalue < 0.05, "*", "ns")))
cat("  Region enrichment:\n"); print(enrichment[, c("Region", "DMP_Pct", "BG_Pct", "Fold", "sig")])
save_data(enrichment, BATCH_DIR, "dmp_region_enrichment")

# =============================================================================
# 5. GO / KEGG / REACTOME ENRICHMENT (STRING-based)
# =============================================================================
cat("[5/10] GO/KEGG/Reactome enrichment (STRING-based)...\n")

if (file.exists(OG$string_enrich)) {
  string_raw <- fread(OG$string_enrich, sep = "\t", header = TRUE, quote = "")
  # Parse gene IDs from STRING protein IDs
  string_raw[, gene_id := sub("^.*\\.(LOC_\\d+)$", "\\1", `#string_protein_id`)]
  string_raw[!grepl("^LOC_", gene_id), gene_id := sub("^.*\\.(XLOC_\\S+)$", "\\1", `#string_protein_id`)]
  string_raw <- string_raw[grepl("^LOC_|^XLOC_", gene_id)]
  cat(sprintf("  STRING annotations: %s entries, %s genes\n",
              format(nrow(string_raw), big.mark = ","),
              format(uniqueN(string_raw$gene_id), big.mark = ",")))

  # Build term tables (as data.frame — clusterProfiler incompatible with data.table)
  string_df <- as.data.frame(string_raw)
  go_df <- string_df[grepl("Gene Ontology", string_df$category), ]
  go_df$go_cat <- ifelse(grepl("Biological Process", go_df$category), "BP",
                  ifelse(grepl("Molecular Function", go_df$category), "MF",
                  ifelse(grepl("Cellular Component", go_df$category), "CC", NA)))
  T2G_go <- unique(go_df[, c("term", "gene_id")])
  T2N_go <- unique(go_df[, c("term", "description")])

  reactome_df <- string_df[grepl("Reactome", string_df$category), ]
  T2G_react <- unique(reactome_df[, c("term", "gene_id")])
  T2N_react <- unique(reactome_df[, c("term", "description")])

  kegg_df <- string_df[grepl("KEGG", string_df$category), ]
  T2G_kegg <- unique(kegg_df[, c("term", "gene_id")])
  T2N_kegg <- unique(kegg_df[, c("term", "description")])
  rm(string_raw, string_df); gc(verbose = FALSE)

  cat(sprintf("  GO: %d terms | Reactome: %d | KEGG: %d\n",
              nrow(T2N_go), nrow(T2N_react), nrow(T2N_kegg)))

  # Enrichment helper
  run_enrichment <- function(gene_list, T2G, T2N, label, universe) {
    if (length(gene_list) < 10 || nrow(T2G) == 0) return(NULL)
    res <- tryCatch(enricher(gene = gene_list,
                              TERM2GENE = T2G,
                              TERM2NAME = T2N,
                              pvalueCutoff = 0.05, pAdjustMethod = "BH",
                              minGSSize = 5, maxGSSize = 500, universe = universe),
                    error = function(e) { cat(sprintf("    ERROR: %s\n", e$message)); NULL })
    if (is.null(res) || nrow(res@result[res@result$p.adjust < 0.05, ]) == 0) return(NULL)
    sig <- res@result[res@result$p.adjust < 0.05, ]
    cat(sprintf("    %s: %d significant terms\n", label, nrow(sig)))
    save_data(res@result, BATCH_DIR, paste0("enrichment_", tolower(gsub(" ", "_", label))))

    # Dotplot
    top_n <- min(20, nrow(sig))
    top <- sig %>% arrange(p.adjust) %>% head(top_n)
    p <- ggplot(top, aes(x = Count, y = reorder(Description, Count),
                          color = p.adjust, size = Count)) +
      geom_point() +
      scale_color_gradient(low = "#E74C3C", high = "#3498DB", name = "Adj. p") +
      scale_size_continuous(range = c(3, 8)) +
      labs(title = label, x = "Gene Count", y = NULL) +
      theme_minimal(base_size = 11) +
      theme(axis.text.y = element_text(size = 8))
    save_fig(p, BATCH_DIR, paste0("fig6_", tolower(gsub("[/ ]", "_", label)), "_dotplot"),
             w = 12, h = max(5, top_n * 0.35))
    res
  }

  dmp_genes <- unique(dmp$nearest_gene[!is.na(dmp$nearest_gene)])
  dmr_genes <- if (n_dmr > 0) unique(dmr$nearest_gene[!is.na(dmr$nearest_gene)]) else character(0)
  all_genes <- genes$ID
  hyper_genes <- unique(dmp$nearest_gene[!is.na(dmp$nearest_gene) & dmp$diff > 0])
  hypo_genes  <- unique(dmp$nearest_gene[!is.na(dmp$nearest_gene) & dmp$diff < 0])

  cat(sprintf("  DMP genes: %d | DMR genes: %d | Hyper: %d | Hypo: %d\n",
              length(dmp_genes), length(dmr_genes), length(hyper_genes), length(hypo_genes)))

  # GO by category
  for (gc in c("BP", "MF", "CC")) {
    cat_terms <- unique(go_df$term[go_df$go_cat == gc])
    t2g_cat <- T2G_go[T2G_go$term %in% cat_terms, ]
    t2n_cat <- T2N_go[T2N_go$term %in% cat_terms, ]
    run_enrichment(dmp_genes, t2g_cat, t2n_cat, paste0("DMP GO ", gc), all_genes)
    run_enrichment(dmr_genes, t2g_cat, t2n_cat, paste0("DMR GO ", gc), all_genes)
  }

  # GO combined
  run_enrichment(dmp_genes, T2G_go, T2N_go, "DMP GO All", all_genes)
  run_enrichment(dmr_genes, T2G_go, T2N_go, "DMR GO All", all_genes)

  # Reactome + KEGG
  run_enrichment(dmp_genes, T2G_react, T2N_react, "DMP Reactome", all_genes)
  run_enrichment(dmr_genes, T2G_react, T2N_react, "DMR Reactome", all_genes)
  run_enrichment(dmp_genes, T2G_kegg, T2N_kegg, "DMP KEGG", all_genes)
  run_enrichment(dmr_genes, T2G_kegg, T2N_kegg, "DMR KEGG", all_genes)

  # Hyper vs Hypo (GO BP only)
  bp_terms <- unique(go_df$term[go_df$go_cat == "BP"])
  t2g_bp <- T2G_go[T2G_go$term %in% bp_terms, ]
  t2n_bp <- T2N_go[T2N_go$term %in% bp_terms, ]
  run_enrichment(hyper_genes, t2g_bp, t2n_bp, "Hyper GO BP", all_genes)
  run_enrichment(hypo_genes,  t2g_bp, t2n_bp, "Hypo GO BP",  all_genes)

  rm(go_df, reactome_df, kegg_df); gc(verbose = FALSE)
} else {
  cat("  WARNING: STRING enrichment file not found at", OG$string_enrich, "\n")
  cat("  Skipping GO/KEGG/Reactome enrichment.\n")
}

# =============================================================================
# 6. FIGURES: PCA + SAMPLE CORRELATION
# =============================================================================
cat("[6/10] PCA + sample correlation...\n")

beta_pca <- getMeth(bs_obj, type = "raw")
beta_pca <- beta_pca[complete.cases(beta_pca), ]

# PCA
pca <- prcomp(t(beta_pca), center = TRUE, scale. = FALSE)
ve <- summary(pca)$importance[2, 1:2] * 100
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                      Sample = rownames(pca$x),
                      Condition = c("Control", "Control", "Amputated", "Amputated"))
p6a <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Condition, label = Sample)) +
  geom_point(size = 5) +
  geom_text(vjust = -1.2, size = 4, show.legend = FALSE) +
  scale_color_manual(values = COLORS$condition) +
  labs(title = "PCA of DNA methylation profiles",
       x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2])) +
  theme_minimal(base_size = 13)
save_fig(p6a, BATCH_DIR, "fig6a_pca", w = 8, h = 6)

# Sample correlation heatmap
cat("  Sample correlation heatmap...\n")
sample_cors <- cor(beta_pca, method = "spearman")
ann_df <- data.frame(Condition = c("Control", "Control", "Amputated", "Amputated"),
                      row.names = colnames(sample_cors))
ann_colors <- list(Condition = c(Control = "#2471A3", Amputated = "#C0392B"))

png(file.path(BATCH_DIR, "figures", "fig6b_correlation_heatmap.png"),
    width = 8, height = 6, units = "in", res = 300)
pheatmap(sample_cors, annotation_col = ann_df, annotation_row = ann_df,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main = "Sample-to-sample methylation correlations (Spearman)",
         fontsize = 12, display_numbers = TRUE, number_format = "%.3f")
dev.off()
tryCatch({
  cairo_pdf(file.path(BATCH_DIR, "figures", "fig6b_correlation_heatmap.pdf"), width = 8, height = 6)
  pheatmap(sample_cors, annotation_col = ann_df, annotation_row = ann_df,
           annotation_colors = ann_colors,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           main = "Sample-to-sample methylation correlations (Spearman)",
           fontsize = 12, display_numbers = TRUE, number_format = "%.3f")
  dev.off()
}, error = function(e) NULL)
cat("  Saved: fig6b_correlation_heatmap\n")

rm(beta_pca, pca); gc(verbose = FALSE)

# =============================================================================
# 7. GENOME-WIDE METHYLATION + MANHATTAN PLOT
# =============================================================================
cat("[7/10] Genome-wide visualization + Manhattan...\n")

# Manhattan plot using DMP data (dml_test already freed for memory)
cat("  Manhattan plot (DMPs only)...\n")
dmp_dt <- as.data.table(dmp)
dmp_dt[, chr_num := as.integer(sub("chr", "", chr))]

chr_info <- dmp_dt[, .(max_pos = max(pos)), by = .(chr, chr_num)][order(chr_num)]
chr_info[, cumulative := cumsum(c(0, head(max_pos, -1)))]
chr_info[, chr_center := cumulative + max_pos / 2]

dmp_dt <- merge(dmp_dt, chr_info[, .(chr, cumulative)], by = "chr")
dmp_dt[, genome_pos := pos + cumulative]

p6c <- ggplot(dmp_dt, aes(x = genome_pos, y = diff, color = direction)) +
  geom_point(alpha = 0.5, size = 0.5) +
  scale_color_manual(values = COLORS$direction, name = "Direction") +
  scale_x_continuous(breaks = chr_info$chr_center,
                      labels = sub("chr", "", chr_info$chr), expand = c(0.01, 0)) +
  geom_hline(yintercept = 0, color = "gray40", linewidth = 0.3) +
  geom_vline(xintercept = chr_info$cumulative[-1], color = "gray80", linewidth = 0.2) +
  labs(title = "Manhattan plot — DMPs across genome",
       subtitle = sprintf("%s DMPs (FDR<0.05, |diff|>10%%): %s hyper + %s hypo",
                          format(n_dmp, big.mark = ","),
                          format(n_hyper, big.mark = ","),
                          format(n_hypo, big.mark = ",")),
       x = "Chromosome", y = "Methylation difference (Amputated - Control)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())
save_fig(p6c, BATCH_DIR, "fig6c_manhattan", w = 14, h = 6)
rm(dmp_dt); gc(verbose = FALSE)

# =============================================================================
# 8. VOLCANO PLOT
# =============================================================================
cat("[8/10] Volcano plot...\n")

p6d <- ggplot(data.frame(diff = dmp$diff, pval = dmp$pval, dir = dmp$direction),
              aes(x = diff, y = -log10(pval), color = dir)) +
  geom_point(alpha = 0.4, size = 0.5) +
  scale_color_manual(values = COLORS$direction, name = "Direction") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray40") +
  labs(title = "DMP Volcano Plot",
       subtitle = sprintf("Hyper: %s | Hypo: %s", format(n_hyper, big.mark = ","),
                          format(n_hypo, big.mark = ",")),
       x = "Methylation Difference (Amputated - Control)",
       y = expression(-log[10](p-value))) +
  theme_minimal(base_size = 12)
save_fig(p6d, BATCH_DIR, "fig6d_dmp_volcano", w = 10, h = 7)

# DMP volcano with gene name labels (top DMPs)
suppressPackageStartupMessages(library(ggrepel))
top_label_dmp <- as.data.frame(dmp[order(-abs(dmp$diff)), ][1:min(30, nrow(dmp)), ])
top_label_dmp$label <- ifelse(!is.na(top_label_dmp$gene_name) & top_label_dmp$gene_name != "",
                               top_label_dmp$gene_name, top_label_dmp$nearest_gene)
p6d_lab <- ggplot(data.frame(diff = dmp$diff, pval = dmp$pval, dir = dmp$direction),
              aes(x = diff, y = -log10(pval), color = dir)) +
  geom_point(alpha = 0.3, size = 0.4) +
  scale_color_manual(values = COLORS$direction, name = "Direction") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray40") +
  geom_text_repel(data = top_label_dmp,
                  aes(x = diff, y = -log10(pval), label = label),
                  size = 2.5, max.overlaps = 25, inherit.aes = FALSE, color = "black") +
  labs(title = "DMP Volcano Plot — top genes labeled",
       x = "Methylation Difference (Amputated - Control)",
       y = expression(-log[10](p-value))) +
  theme_minimal(base_size = 12)
save_fig(p6d_lab, BATCH_DIR, "fig6d2_dmp_volcano_labeled", w = 12, h = 8)

# DMR volcano: log10(|areaStat|+1) as significance proxy — readable y range
if (n_dmr > 0) {
  cat("  DMR volcano (areaStat-based)...\n")
  dmr_df <- as.data.frame(dmr)
  dmr_df$sig_proxy <- log10(abs(dmr_df$areaStat) + 1)
  dmr_df$sig_cat <- ifelse(dmr_df$direction == "Hyper", "Hypermethylated",
                    ifelse(dmr_df$direction == "Hypo", "Hypomethylated", "Not Significant"))
  n_hyper <- sum(dmr_df$direction == "Hyper")
  n_hypo  <- sum(dmr_df$direction == "Hypo")

  # Top 10 DMRs by |areaStat| for labels
  dmr_df$label <- paste0(dmr_df$seqnames, ":", format(dmr_df$start, big.mark = ","))
  top_dmrs <- dmr_df[order(-abs(dmr_df$areaStat)), ][1:min(10, nrow(dmr_df)), ]

  # Base volcano (no labels)
  p6d_dmr <- ggplot(dmr_df, aes(x = diff.Methy, y = sig_proxy, color = sig_cat)) +
    geom_point(alpha = 0.6, size = 1.4) +
    scale_color_manual(values = c("Hypermethylated" = COLORS$direction[["Hyper"]],
                                  "Hypomethylated"  = COLORS$direction[["Hypo"]],
                                  "Not Significant" = "gray70"),
                       name = "Methylation\nChange") +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", color = "gray40") +
    labs(title = "Volcano Plot \u2014 Differentially Methylated Regions (DMRs)",
         subtitle = sprintf("Hypermethylated: %s | Hypomethylated: %s | Total: %s",
                            format(n_hyper, big.mark = ","),
                            format(n_hypo, big.mark = ","),
                            format(n_hyper + n_hypo, big.mark = ",")),
         x = "Methylation Difference (Amputated \u2212 Control)",
         y = expression(log[10]("|areaStat| + 1") ~ "[significance proxy]"),
         caption = "areaStat from DSS::callDMR \u2014 larger |value| = stronger signal") +
    theme_minimal(base_size = 12)
  save_fig(p6d_dmr, BATCH_DIR, "fig6d3_dmr_volcano", w = 10, h = 7)

  # Labeled version (top 10 DMRs)
  p6d_dmr_lab <- p6d_dmr +
    ggrepel::geom_text_repel(data = top_dmrs, aes(label = label),
                             size = 3, max.overlaps = 15, show.legend = FALSE)
  save_fig(p6d_dmr_lab, BATCH_DIR, "fig6d4_dmr_volcano_labeled", w = 10, h = 7)

}

# =============================================================================
# 9. ANNOTATION PLOTS + ENRICHMENT BAR
# =============================================================================
cat("[9/10] Annotation and enrichment plots...\n")

# Annotation pie
anno_tab <- as.data.frame(table(dmp$annotation))
anno_tab$pct <- round(100 * anno_tab$Freq / sum(anno_tab$Freq), 1)
anno_tab$label <- paste0(anno_tab$Var1, "\n(", anno_tab$pct, "%)")

p6e <- ggplot(anno_tab, aes(x = "", y = Freq, fill = Var1)) +
  geom_col(width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = COLORS$region) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3.5) +
  labs(title = "DMP genomic distribution", fill = "Region") +
  theme_void(base_size = 12)
save_fig(p6e, BATCH_DIR, "fig6e_dmp_annotation_pie", w = 8, h = 7)

# Direction by region stacked bar
dir_tab <- as.data.frame(table(dmp$annotation, dmp$direction))
colnames(dir_tab) <- c("Region", "Direction", "Count")
p6f <- ggplot(dir_tab, aes(x = Region, y = Count, fill = Direction)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = COLORS$direction) +
  labs(title = "DMP direction by genomic region", x = NULL, y = "Count") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p6f, BATCH_DIR, "fig6f_dmp_direction_by_region", w = 8, h = 6)

# DMR annotation pie + direction by region
if (n_dmr > 0) {
  cat("  DMR annotation plots...\n")
  dmr_anno_tab <- as.data.frame(table(dmr$annotation))
  dmr_anno_tab$pct <- round(100 * dmr_anno_tab$Freq / sum(dmr_anno_tab$Freq), 1)
  dmr_anno_tab$label <- paste0(dmr_anno_tab$Var1, "\n(", dmr_anno_tab$pct, "%)")

  p6f2 <- ggplot(dmr_anno_tab, aes(x = "", y = Freq, fill = Var1)) +
    geom_col(width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = COLORS$region) +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3.5) +
    labs(title = "DMR genomic distribution", fill = "Region") +
    theme_void(base_size = 12)
  save_fig(p6f2, BATCH_DIR, "fig6f2_dmr_annotation_pie", w = 8, h = 7)

  dmr_dir_tab <- as.data.frame(table(dmr$annotation, dmr$direction))
  colnames(dmr_dir_tab) <- c("Region", "Direction", "Count")
  p6f3 <- ggplot(dmr_dir_tab, aes(x = Region, y = Count, fill = Direction)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = COLORS$direction) +
    labs(title = "DMR direction by genomic region", x = NULL, y = "Count") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p6f3, BATCH_DIR, "fig6f3_dmr_direction_by_region", w = 8, h = 6)
}

# Enrichment bar
enrich_long <- rbind(
  data.frame(Region = enrichment$Region, Pct = enrichment$DMP_Pct, Source = "DMPs"),
  data.frame(Region = enrichment$Region, Pct = enrichment$BG_Pct, Source = "All CpGs"))
p6g <- ggplot(enrich_long, aes(x = Region, y = Pct, fill = Source)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("DMPs" = "#E74C3C", "All CpGs" = "#95A5A6")) +
  geom_text(data = enrichment, aes(x = Region, y = pmax(DMP_Pct, BG_Pct) + 1.5,
            label = paste0(Fold, "x ", sig)), inherit.aes = FALSE, size = 3.5) +
  labs(title = "DMP enrichment by genomic region (Fisher's exact test)",
       x = NULL, y = "% of positions") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p6g, BATCH_DIR, "fig6g_dmp_enrichment", w = 8, h = 6)

# DMR enrichment by region (Fisher's exact test)
if (n_dmr > 0) {
  cat("  DMR region enrichment (Fisher's exact)...\n")
  dmr_counts <- as.data.frame(table(dmr$annotation)); colnames(dmr_counts) <- c("Region", "DMR_Count")
  dmr_enrichment <- merge(dmr_counts, bg_counts, by = "Region")
  dmr_enrichment$DMR_Pct <- round(100 * dmr_enrichment$DMR_Count / sum(dmr_enrichment$DMR_Count), 1)
  dmr_enrichment$BG_Pct  <- round(100 * dmr_enrichment$BG_Count / sum(dmr_enrichment$BG_Count), 1)
  dmr_enrichment$Fold <- round(dmr_enrichment$DMR_Pct / dmr_enrichment$BG_Pct, 2)
  dmr_enrichment$pvalue <- sapply(1:nrow(dmr_enrichment), function(i) {
    fisher.test(matrix(c(dmr_enrichment$DMR_Count[i], n_dmr - dmr_enrichment$DMR_Count[i],
                          dmr_enrichment$BG_Count[i], sum(bg_counts$BG_Count) - dmr_enrichment$BG_Count[i]),
                        nrow = 2))$p.value
  })
  dmr_enrichment$sig <- ifelse(dmr_enrichment$pvalue < 0.001, "***",
                        ifelse(dmr_enrichment$pvalue < 0.01, "**",
                        ifelse(dmr_enrichment$pvalue < 0.05, "*", "ns")))
  cat("  DMR region enrichment:\n")
  print(dmr_enrichment[, c("Region", "DMR_Pct", "BG_Pct", "Fold", "sig")])
  save_data(dmr_enrichment, BATCH_DIR, "dmr_region_enrichment")

  dmr_enrich_long <- rbind(
    data.frame(Region = dmr_enrichment$Region, Pct = dmr_enrichment$DMR_Pct, Source = "DMRs"),
    data.frame(Region = dmr_enrichment$Region, Pct = dmr_enrichment$BG_Pct, Source = "All CpGs"))
  p6g2 <- ggplot(dmr_enrich_long, aes(x = Region, y = Pct, fill = Source)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(values = c("DMRs" = "#E74C3C", "All CpGs" = "#95A5A6")) +
    geom_text(data = dmr_enrichment, aes(x = Region, y = pmax(DMR_Pct, BG_Pct) + 1.5,
              label = paste0(Fold, "x ", sig)), inherit.aes = FALSE, size = 3.5) +
    labs(title = "DMR enrichment by genomic region (Fisher's exact test)",
         x = NULL, y = "% of regions") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p6g2, BATCH_DIR, "fig6g2_dmr_enrichment", w = 8, h = 6)
}

# =============================================================================
# 10. DMP + DMR HEATMAPS
# =============================================================================
cat("[10/10] DMP and DMR heatmaps...\n")

# Top 100 DMPs heatmap
n_top <- min(100, n_dmp)
top_dmps <- as.data.frame(dmp[order(-abs(dmp$diff)), ][1:n_top, ])
dmp_ids <- paste(top_dmps$chr, top_dmps$pos, sep = ":")
bs_ids  <- paste(as.character(seqnames(bs_obj)), start(bs_obj), sep = ":")
match_idx <- match(dmp_ids, bs_ids)
# Build labels with gene names
dmp_labels <- ifelse(!is.na(top_dmps$gene_name) & top_dmps$gene_name != "",
                     paste0(top_dmps$gene_name, " (", top_dmps$chr, ":", top_dmps$pos, ")"),
                     paste0(top_dmps$nearest_gene, " (", top_dmps$chr, ":", top_dmps$pos, ")"))
match_valid <- !is.na(match_idx)
match_idx <- match_idx[match_valid]
dmp_labels <- dmp_labels[match_valid]

if (length(match_idx) >= 2) {
  beta_mat <- getMeth(bs_obj[match_idx, ], type = "raw")
  rownames(beta_mat) <- dmp_labels
  valid_rows <- complete.cases(beta_mat) & apply(beta_mat, 1, var, na.rm = TRUE) > 0
  beta_mat <- beta_mat[valid_rows, ]

  if (nrow(beta_mat) >= 2) {
    col_ann <- data.frame(Condition = c("Control", "Control", "Amputated", "Amputated"),
                           row.names = colnames(beta_mat))
    ann_colors <- list(Condition = c(Control = "#2471A3", Amputated = "#C0392B"))

    png(file.path(BATCH_DIR, "figures", "fig6h_dmp_heatmap.png"),
        width = 8, height = 12, units = "in", res = 300)
    pheatmap(beta_mat, scale = "row", annotation_col = col_ann,
             annotation_colors = ann_colors,
             color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
             main = sprintf("Top %d DMPs — site-level methylation", nrow(beta_mat)),
             fontsize = 8, show_rownames = TRUE, fontsize_row = 5, show_colnames = TRUE,
             clustering_method = "complete")
    dev.off()
    tryCatch({
      cairo_pdf(file.path(BATCH_DIR, "figures", "fig6h_dmp_heatmap.pdf"), width = 8, height = 12)
      pheatmap(beta_mat, scale = "row", annotation_col = col_ann,
               annotation_colors = ann_colors,
               color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
               main = sprintf("Top %d DMPs — site-level methylation", nrow(beta_mat)),
               fontsize = 8, show_rownames = TRUE, fontsize_row = 5, show_colnames = TRUE,
               clustering_method = "complete")
      dev.off()
    }, error = function(e) NULL)
    cat("  Saved: fig6h_dmp_heatmap\n")
  }
}

# Top 50 DMRs heatmap
if (n_dmr > 0) {
  n_top_dmr <- min(50, n_dmr)
  top_dmrs <- dmr[order(-abs(dmr$diff.Methy)), ][1:n_top_dmr, ]
  bs_ranges <- GRanges(seqnames = seqnames(bs_obj), ranges = IRanges(start = start(bs_obj), width = 1))

  dmr_matrix <- matrix(NA_real_, nrow = n_top_dmr, ncol = ncol(bs_obj))
  colnames(dmr_matrix) <- colnames(bs_obj)
  top_dmrs_df <- as.data.frame(top_dmrs)
  rownames(dmr_matrix) <- ifelse(!is.na(top_dmrs_df$gene_name) & top_dmrs_df$gene_name != "",
    paste0(top_dmrs_df$gene_name, " (", top_dmrs_df$chr, ":", top_dmrs_df$start, "-", top_dmrs_df$end, ")"),
    paste0(top_dmrs_df$nearest_gene, " (", top_dmrs_df$chr, ":", top_dmrs_df$start, "-", top_dmrs_df$end, ")"))

  for (i in 1:n_top_dmr) {
    dmr_range <- GRanges(seqnames = top_dmrs$chr[i],
                          ranges = IRanges(start = top_dmrs$start[i], end = top_dmrs$end[i]))
    ov <- findOverlaps(bs_ranges, dmr_range)
    cpg_idx <- queryHits(ov)
    if (length(cpg_idx) > 0) {
      dmr_beta <- getMeth(bs_obj[cpg_idx, ], type = "raw")
      dmr_matrix[i, ] <- if (length(cpg_idx) == 1) as.numeric(dmr_beta) else colMeans(dmr_beta, na.rm = TRUE)
    }
  }

  valid_dmrs <- !apply(dmr_matrix, 1, function(x) any(is.na(x)))
  if (sum(valid_dmrs) >= 2) {
    dmr_matrix <- dmr_matrix[valid_dmrs, ]
    col_ann <- data.frame(Condition = c("Control", "Control", "Amputated", "Amputated"),
                           row.names = colnames(dmr_matrix))

    png(file.path(BATCH_DIR, "figures", "fig6i_dmr_heatmap.png"),
        width = 10, height = 12, units = "in", res = 300)
    pheatmap(dmr_matrix, scale = "row", annotation_col = col_ann,
             annotation_colors = ann_colors,
             color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
             main = sprintf("Top %d DMRs — region-level methylation", nrow(dmr_matrix)),
             fontsize = 10, show_rownames = TRUE, show_colnames = TRUE,
             clustering_method = "complete")
    dev.off()
    tryCatch({
      cairo_pdf(file.path(BATCH_DIR, "figures", "fig6i_dmr_heatmap.pdf"), width = 10, height = 12)
      pheatmap(dmr_matrix, scale = "row", annotation_col = col_ann,
               annotation_colors = ann_colors,
               color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
               main = sprintf("Top %d DMRs — region-level methylation", nrow(dmr_matrix)),
               fontsize = 10, show_rownames = TRUE, show_colnames = TRUE,
               clustering_method = "complete")
      dev.off()
    }, error = function(e) NULL)
    cat("  Saved: fig6i_dmr_heatmap\n")
  }
}

# Top DMP genes barplot (gene names already merged above)
cat("  Top genes by DMP count...\n")
top_genes <- as.data.table(dmp)[, .(n = .N, gene_name = gene_name[1]),
                                 by = nearest_gene][order(-n)][!is.na(nearest_gene)][1:min(20, .N)]
top_genes[, label := ifelse(!is.na(gene_name) & gene_name != "",
                            paste0(gene_name, " (", nearest_gene, ")"),
                            nearest_gene)]

p6j <- ggplot(top_genes, aes(x = reorder(label, n), y = n)) +
  geom_col(fill = "#C0392B") +
  coord_flip() +
  labs(title = "Top 20 genes by DMP count", x = NULL, y = "Number of DMPs") +
  theme_minimal(base_size = 11)
save_fig(p6j, BATCH_DIR, "fig6j_top_genes_dmp", w = 10, h = 7)

# Top DMR genes barplot
if (n_dmr > 0) {
  cat("  Top genes by DMR count...\n")
  top_dmr_genes <- as.data.table(dmr)[, .(n = .N, gene_name = gene_name[1]),
                                       by = nearest_gene][order(-n)][!is.na(nearest_gene)][1:min(20, .N)]
  top_dmr_genes[, label := ifelse(!is.na(gene_name) & gene_name != "",
                              paste0(gene_name, " (", nearest_gene, ")"),
                              nearest_gene)]
  p6k <- ggplot(top_dmr_genes, aes(x = reorder(label, n), y = n)) +
    geom_col(fill = "#2471A3") +
    coord_flip() +
    labs(title = "Top 20 genes by DMR count", x = NULL, y = "Number of DMRs") +
    theme_minimal(base_size = 11)
  save_fig(p6k, BATCH_DIR, "fig6k_top_genes_dmr", w = 10, h = 7)
}

# =============================================================================
# SUMMARY
# =============================================================================
rm(bs_obj); gc(verbose = FALSE)

elapsed <- (proc.time() - t0)[3]
n_figs <- length(list.files(file.path(BATCH_DIR, "figures"), pattern = "\\.png$"))
n_data <- length(list.files(file.path(BATCH_DIR, "data"), pattern = "\\.tsv$"))
cat(sprintf("\n=== Batch 06 complete (%.1f min) ===\n", elapsed / 60))
cat(sprintf("Figures: %d\nData files: %d\n", n_figs, n_data))
cat(sprintf("DMPs: %s | DMRs: %s\n", format(n_dmp, big.mark = ","), format(n_dmr, big.mark = ",")))
