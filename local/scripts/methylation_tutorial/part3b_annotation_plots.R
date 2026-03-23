#!/usr/bin/env Rscript
# =============================================================================
# Part 3B: Annotation, Enrichment & Plots (uses cached DMLtest results)
# Skips the 3.5-hour DMLtest — reads from cache, redoes annotation + plots
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(bsseq)
  library(DSS)
  library(GenomicRanges)
  library(rtracklayer)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(scales)
  library(gridExtra)
})

cat("=== Part 3B: Annotation + Plots (cached DMLtest) ===\n\n")

# --- PATHS ---
project_dir <- "C:/Users/rafae/Projects/STANDBY"
cache_dir   <- file.path(project_dir, "genome/cache")
out_dir     <- file.path(project_dir, "results/methylation_tutorial")
fig_pdf_dir <- file.path(out_dir, "pdf")
fig_png_dir <- file.path(out_dir, "png")
data_dir_out <- file.path(out_dir, "data")

dir.create(fig_pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_png_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir_out, recursive = TRUE, showWarnings = FALSE)

keep_chr <- paste0("chr", 1:31)

save_plot <- function(p, name, w = 8, h = 5) {
  png(file.path(fig_png_dir, paste0(name, ".png")),
      width = w, height = h, units = "in", res = 300)
  print(p); dev.off()
  tryCatch({
    cairo_pdf(file.path(fig_pdf_dir, paste0(name, ".pdf")), width = w, height = h)
    print(p); dev.off()
  }, error = function(e) {
    pdf(file.path(fig_pdf_dir, paste0(name, ".pdf")), width = w, height = h)
    print(p); dev.off()
  })
  cat("  Saved:", name, "\n")
}

# --- LOAD CACHED DATA ---
cat("Loading cached DMLtest...\n")
dml_cache <- file.path(cache_dir, "dmltest_full.rds")
if (!file.exists(dml_cache)) stop("No cached DMLtest found. Run part3 first.")
dml_test <- readRDS(dml_cache)
cat(sprintf("DMLtest: %s sites\n", format(nrow(dml_test), big.mark = ",")))

cat("Loading BSseq...\n")
bs_obj <- readRDS(file.path(cache_dir, "bsseq_tutorial.rds"))
cat(sprintf("BSseq: %s sites\n\n", format(nrow(bs_obj), big.mark = ",")))

# --- CALL DMPs/DMRs ---
cat("=== Calling DMPs ===\n")
dmp_results <- callDML(dml_test, p.threshold = 0.05, delta = 0.1)
n_dmp <- if (!is.null(dmp_results)) nrow(dmp_results) else 0
cat(sprintf("DMPs (FDR<0.05, |diff|>10%%): %s\n", format(n_dmp, big.mark = ",")))

dmr_results <- callDMR(dml_test, p.threshold = 0.05, delta = 0.1,
                       minlen = 50, minCG = 3)
if (is.null(dmr_results)) {
  dmr_results <- data.frame(chr = character(), start = integer(), end = integer(),
                            length = integer(), nCG = integer(), areaStat = numeric())
}
cat(sprintf("DMRs: %s\n", format(nrow(dmr_results), big.mark = ",")))

n_hyper <- sum(dmp_results$diff > 0)
n_hypo  <- sum(dmp_results$diff < 0)
cat(sprintf("Hyper: %s | Hypo: %s\n\n", format(n_hyper, big.mark = ","),
            format(n_hypo, big.mark = ",")))

rm(dml_test); gc(verbose = FALSE)

# --- LOAD ANNOTATIONS ---
cat("Loading annotations...\n")
gff <- readRDS(file.path(cache_dir, "gff_chr1_31.rds"))
genes <- gff[gff$type == "gene"]
exons <- gff[gff$type == "exon"]

te_rds <- file.path(cache_dir, "te_chr1_31.rds")
has_te <- file.exists(te_rds)
if (has_te) {
  te_data <- readRDS(te_rds)
  if (is(te_data, "list") && "te_gr" %in% names(te_data)) {
    te_gr <- te_data$te_gr
  } else if (is(te_data, "GRanges")) {
    te_gr <- te_data
  } else {
    te_gr <- te_data[[1]]
  }
  rm(te_data); gc(verbose = FALSE)
}

promoters_rds <- file.path(cache_dir, "promoters_2kb.rds")
if (file.exists(promoters_rds)) {
  promoters <- readRDS(promoters_rds)
} else {
  gene_strand <- as.character(strand(genes))
  gene_chr <- as.character(seqnames(genes))
  prom_start <- ifelse(gene_strand == "+", pmax(1L, start(genes) - 2000L), end(genes) + 1L)
  prom_end <- ifelse(gene_strand == "+", start(genes) - 1L, end(genes) + 2000L)
  valid <- prom_start < prom_end
  promoters <- GRanges(seqnames = gene_chr[valid],
                       ranges = IRanges(start = prom_start[valid], end = prom_end[valid]))
}

# Functional region annotation (NO TE — TEs are repeat annotation, not functional region)
annotate_positions <- function(chr, pos) {
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1))
  in_promoter <- overlapsAny(gr, promoters)
  in_exon     <- overlapsAny(gr, exons)
  in_gene     <- overlapsAny(gr, genes)
  annotation <- rep("Intergenic", length(gr))
  annotation[in_gene & !in_exon] <- "Intron"
  annotation[in_exon] <- "Exon"
  annotation[in_promoter] <- "Promoter"
  return(annotation)
}

annotate_te <- function(chr, pos) {
  if (!has_te) return(rep(FALSE, length(chr)))
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1))
  overlapsAny(gr, te_gr)
}

# --- ANNOTATE DMPs ---
cat("\n=== Annotating DMPs ===\n")
dmp_results$annotation <- annotate_positions(dmp_results$chr, dmp_results$pos)
dmp_results$in_te <- annotate_te(dmp_results$chr, dmp_results$pos)
dmp_results$direction <- ifelse(dmp_results$diff > 0, "Hyper", "Hypo")

dmp_anno_summary <- as.data.frame(table(dmp_results$annotation))
colnames(dmp_anno_summary) <- c("Region", "Count")
dmp_anno_summary$Percent <- round(100 * dmp_anno_summary$Count / sum(dmp_anno_summary$Count), 1)
dmp_anno_summary <- dmp_anno_summary[order(-dmp_anno_summary$Count), ]
cat("DMP annotation (functional regions):\n")
print(dmp_anno_summary)
cat(sprintf("DMPs overlapping TEs: %s / %s (%.1f%%)\n\n",
            format(sum(dmp_results$in_te), big.mark = ","),
            format(nrow(dmp_results), big.mark = ","),
            100 * mean(dmp_results$in_te)))


# =============================================================================
# ENRICHMENT: Are DMPs enriched/depleted in each region vs genome background?
# =============================================================================
cat("=== DMP Region Enrichment ===\n")

# Count CpGs per region from the BSseq object (all 25M tested CpGs)
all_chr <- as.character(seqnames(bs_obj))
all_pos <- start(bs_obj)
all_anno <- annotate_positions(all_chr, all_pos)
all_te <- annotate_te(all_chr, all_pos)

bg_counts <- as.data.frame(table(all_anno))
colnames(bg_counts) <- c("Region", "BG_Count")
bg_counts$BG_Percent <- round(100 * bg_counts$BG_Count / sum(bg_counts$BG_Count), 1)

# Merge DMP counts with background
enrichment <- merge(dmp_anno_summary, bg_counts, by = "Region")
enrichment$DMP_Pct <- enrichment$Percent
enrichment$BG_Pct <- enrichment$BG_Percent
enrichment$Fold <- round(enrichment$DMP_Pct / enrichment$BG_Pct, 2)

# Fisher's exact test per region
enrichment$pvalue <- sapply(1:nrow(enrichment), function(i) {
  reg <- enrichment$Region[i]
  dmp_in <- enrichment$Count[i]
  dmp_out <- n_dmp - dmp_in
  bg_in <- enrichment$BG_Count[i]
  bg_out <- sum(bg_counts$BG_Count) - bg_in
  fisher.test(matrix(c(dmp_in, dmp_out, bg_in, bg_out), nrow = 2))$p.value
})
enrichment$sig <- ifelse(enrichment$pvalue < 0.001, "***",
                  ifelse(enrichment$pvalue < 0.01, "**",
                  ifelse(enrichment$pvalue < 0.05, "*", "ns")))

cat("Region enrichment (DMPs vs all CpGs):\n")
print(enrichment[, c("Region", "DMP_Pct", "BG_Pct", "Fold", "sig")])

# TE enrichment
te_dmp_pct <- 100 * mean(dmp_results$in_te)
te_bg_pct <- 100 * mean(all_te)
te_fold <- round(te_dmp_pct / te_bg_pct, 2)
te_fisher <- fisher.test(matrix(c(sum(dmp_results$in_te), n_dmp - sum(dmp_results$in_te),
                                   sum(all_te), length(all_te) - sum(all_te)), nrow = 2))
cat(sprintf("\nTE overlap: DMPs=%.1f%% vs background=%.1f%% (fold=%.2f, p=%s)\n",
            te_dmp_pct, te_bg_pct, te_fold,
            ifelse(te_fisher$p.value < 2.2e-16, "< 2.2e-16",
                   format(te_fisher$p.value, digits = 3))))

rm(all_chr, all_pos, all_anno, all_te, bs_obj); gc(verbose = FALSE)

# --- ENRICHMENT PLOT ---
enrich_plot_df <- enrichment[, c("Region", "DMP_Pct", "BG_Pct", "Fold", "sig")]
enrich_long <- rbind(
  data.frame(Region = enrichment$Region, Pct = enrichment$DMP_Pct, Source = "DMPs"),
  data.frame(Region = enrichment$Region, Pct = enrichment$BG_Pct, Source = "All CpGs")
)
enrich_long$Region <- factor(enrich_long$Region,
                             levels = enrichment$Region[order(enrichment$Fold, decreasing = TRUE)])

region_colors_func <- c("Promoter" = "#8E44AD", "Exon" = "#2471A3",
                        "Intron" = "#1ABC9C", "Intergenic" = "#C0392B")

p_enrich <- ggplot(enrich_long, aes(x = Region, y = Pct, fill = Source)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("DMPs" = "#E74C3C", "All CpGs" = "#95A5A6")) +
  geom_text(data = enrichment, aes(x = Region, y = pmax(DMP_Pct, BG_Pct) + 1.5,
            label = paste0(Fold, "x ", sig)), inherit.aes = FALSE, size = 3.5) +
  labs(title = "DMP Enrichment by Genomic Region",
       subtitle = sprintf("D. laeve | %s DMPs vs %s background CpGs | Fold change + Fisher's test",
                          format(n_dmp, big.mark = ","),
                          format(sum(bg_counts$BG_Count), big.mark = ",")),
       x = NULL, y = "% of positions") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_enrich, "dmp_enrichment_by_region", w = 8, h = 6)


# --- PIE CHART (functional regions only, no TE) ---
p_anno <- ggplot(dmp_anno_summary, aes(x = "", y = Count, fill = Region)) +
  geom_col(width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = region_colors_func) +
  labs(title = "DMP Genomic Distribution",
       subtitle = sprintf("D. laeve | %s DMPs (TE overlap: %.1f%%)",
                          format(n_dmp, big.mark = ","),
                          100 * mean(dmp_results$in_te))) +
  theme_void(base_size = 12) +
  theme(legend.position = "right")

save_plot(p_anno, "dmp_annotation_pie", w = 8, h = 6)

# --- DIRECTION BAR ---
dir_anno <- as.data.frame(table(dmp_results$annotation, dmp_results$direction))
colnames(dir_anno) <- c("Region", "Direction", "Count")

p_dir <- ggplot(dir_anno, aes(x = Region, y = Count, fill = Direction)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Hyper" = "#C0392B", "Hypo" = "#2471A3")) +
  labs(title = "DMP Direction by Genomic Region",
       subtitle = "D. laeve | Amputated vs Control",
       x = NULL, y = "Number of DMPs") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p_dir, "dmp_direction_by_region", w = 8, h = 6)

# --- VOLCANO ---
cat("\nCreating volcano plot...\n")
volcano_data <- data.frame(
  diff = dmp_results$diff, pvalue = dmp_results$pval, fdr = dmp_results$fdr
)
volcano_data$category <- "Not Significant"
volcano_data$category[volcano_data$fdr < 0.05 & volcano_data$diff > 0.1] <- "Hypermethylated"
volcano_data$category[volcano_data$fdr < 0.05 & volcano_data$diff < -0.1] <- "Hypomethylated"
volcano_data$category <- factor(volcano_data$category,
                                levels = c("Not Significant", "Hypermethylated", "Hypomethylated"))

p_volcano <- ggplot(volcano_data, aes(x = diff, y = -log10(pvalue), color = category)) +
  geom_point(alpha = 0.5, size = 0.6) +
  scale_color_manual(values = c("Not Significant" = "gray70",
                                "Hypermethylated" = "#C0392B",
                                "Hypomethylated" = "#2471A3"),
                     name = "Methylation Change") +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", alpha = 0.7, color = "gray40") +
  labs(title = "Volcano Plot -- Differentially Methylated Positions",
       subtitle = sprintf("D. laeve | Hyper: %s | Hypo: %s",
                          format(n_hyper, big.mark = ","), format(n_hypo, big.mark = ",")),
       x = "Methylation Difference (Amputated - Control)",
       y = expression(-log[10](p-value))) +
  theme_minimal(base_size = 12)

save_plot(p_volcano, "volcano_plot_dmps", w = 10, h = 7)
rm(volcano_data); gc(verbose = FALSE)

# --- NEAREST GENE ---
cat("Finding nearest genes...\n")
dmp_gr <- GRanges(seqnames = dmp_results$chr,
                  ranges = IRanges(start = dmp_results$pos, width = 1))
nearest_idx <- nearest(dmp_gr, genes)
dmp_results$nearest_gene <- ifelse(!is.na(nearest_idx), genes$ID[nearest_idx], NA)
dmp_results$dist_to_gene <- ifelse(!is.na(nearest_idx),
                                   mcols(distanceToNearest(dmp_gr, genes))$distance, NA)

# --- DMR ANNOTATION ---
if (nrow(dmr_results) > 0) {
  dmr_mid <- (dmr_results$start + dmr_results$end) %/% 2
  dmr_results$annotation <- annotate_positions(dmr_results$chr, dmr_mid)
  dmr_results$in_te <- annotate_te(dmr_results$chr, dmr_mid)
  dmr_gr <- GRanges(seqnames = dmr_results$chr,
                    ranges = IRanges(start = dmr_results$start, end = dmr_results$end))
  nearest_idx <- nearest(dmr_gr, genes)
  dmr_results$nearest_gene <- ifelse(!is.na(nearest_idx), genes$ID[nearest_idx], NA)
}

# --- SAVE DATA ---
cat("Saving data...\n")
fwrite(dmp_results, file.path(data_dir_out, "dmps_annotated.tsv"), sep = "\t")
fwrite(as.data.frame(dmr_results), file.path(data_dir_out, "dmrs_annotated.tsv"), sep = "\t")
fwrite(enrichment, file.path(data_dir_out, "dmp_region_enrichment.tsv"), sep = "\t")

# --- HTML REPORT ---
cat("\n=== Generating HTML Report ===\n")
report_path <- file.path(out_dir, "part3_report.html")
sink(report_path)

cat("<!DOCTYPE html>\n<html><head><meta charset='utf-8'>\n")
cat("<title>Methylation Analysis: D. laeve Full-Genome</title>\n")
cat("<style>\n")
cat("body { font-family: 'Segoe UI', Arial, sans-serif; max-width: 1200px; margin: 40px auto; padding: 0 20px; color: #333; }\n")
cat("h1 { color: #2c3e50; border-bottom: 2px solid #2471A3; padding-bottom: 10px; }\n")
cat("h2 { color: #2471A3; margin-top: 30px; }\n")
cat("table { border-collapse: collapse; width: auto; margin: 15px 0; }\n")
cat("th { background: #2471A3; color: white; padding: 8px 12px; text-align: left; }\n")
cat("td { padding: 6px 12px; border-bottom: 1px solid #ddd; }\n")
cat("tr:nth-child(even) { background: #f8f9fa; }\n")
cat("img { max-width: 100%; border: 1px solid #ddd; margin: 10px 0; }\n")
cat(".metric { font-size: 1.2em; font-weight: bold; color: #2471A3; }\n")
cat("</style>\n</head><body>\n")

cat("<h1>Full-Genome Methylation Analysis -- D. laeve</h1>\n")
cat("<p>25,439,068 CpGs | DSS with smoothing | Generated:", format(Sys.Date()), "</p>\n\n")

cat("<h2>1. Summary</h2>\n")
cat("<table>\n")
cat(sprintf("<tr><td>Total CpGs tested</td><td class='metric'>%s</td></tr>\n", format(n_dmp + (25439068 - n_dmp), big.mark = ",")))
cat(sprintf("<tr><td>Significant DMPs</td><td class='metric'>%s</td></tr>\n", format(n_dmp, big.mark = ",")))
cat(sprintf("<tr><td>Hypermethylated</td><td class='metric'>%s</td></tr>\n", format(n_hyper, big.mark = ",")))
cat(sprintf("<tr><td>Hypomethylated</td><td class='metric'>%s</td></tr>\n", format(n_hypo, big.mark = ",")))
cat(sprintf("<tr><td>DMRs</td><td class='metric'>%s</td></tr>\n", format(nrow(dmr_results), big.mark = ",")))
cat("</table>\n\n")

cat("<h2>2. Volcano Plot</h2>\n")
cat("<img src='png/volcano_plot_dmps.png'>\n\n")

cat("<h2>3. DMP Region Enrichment</h2>\n")
cat("<img src='png/dmp_enrichment_by_region.png'>\n")
cat("<table><tr><th>Region</th><th>DMP %</th><th>Background %</th><th>Fold</th><th>Sig</th></tr>\n")
for (i in seq_len(nrow(enrichment))) {
  cat(sprintf("<tr><td>%s</td><td>%.1f</td><td>%.1f</td><td>%.2f</td><td>%s</td></tr>\n",
              enrichment$Region[i], enrichment$DMP_Pct[i], enrichment$BG_Pct[i],
              enrichment$Fold[i], enrichment$sig[i]))
}
cat("</table>\n")
cat(sprintf("<p>TE overlap: DMPs=%.1f%% vs background=%.1f%% (fold=%.2f)</p>\n\n",
            te_dmp_pct, te_bg_pct, te_fold))

cat("<h2>4. DMP Genomic Distribution</h2>\n")
cat("<img src='png/dmp_annotation_pie.png'>\n\n")

cat("<h2>5. DMP Direction by Region</h2>\n")
cat("<img src='png/dmp_direction_by_region.png'>\n\n")

cat("\n</body></html>\n")
sink()
cat("Report saved:", report_path, "\n")

cat("\n=== Part 3B Complete ===\n")
