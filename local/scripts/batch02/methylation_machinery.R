#!/usr/bin/env Rscript
# =============================================================================
# Batch 02: Methylation Machinery Expression Across Tissues
# Paper: Act 1, Section 1.1 — Does D. laeve have a DNA methylation toolkit?
# Question: Which methylation writers/erasers/readers are expressed, and where?
# Plot 1: Genome presence/absence (green = present, black = absent)
# Plot 2: Expression by tissue (boxplot, outliers filtered)
# =============================================================================

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(rtracklayer)
library(RColorBrewer)
library(scales)
library(grid)
library(pheatmap)

# -- Paths --
data_dir    <- "C:/Users/rafae/Projects/DATA"
project_dir <- "C:/Users/rafae/Projects/STANDBY"
dir_counts  <- file.path(data_dir, "counts_HTseq_EviAnn")
gff_cache   <- file.path(project_dir, "genome/cache/gff_chr1_31.rds")
annot_file  <- file.path(data_dir, "derLaeGenome_eviann_annotations.tsv")
out_dir     <- file.path(project_dir, "results/batch02")

dir.create(file.path(out_dir, "pdf"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "png"), recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. CURATED METHYLATION MACHINERY GENE LIST
# =============================================================================

cat("=== Loading GFF + annotations ===\n")
gff <- readRDS(gff_cache)
genes <- gff[gff$type == "gene"]
gff_gene_ids <- as.character(mcols(genes)$ID)

annot <- read.delim(annot_file, header = FALSE,
                    col.names = c("gene_id", "gene_name", "description"))

# User-specified gene list — verified from GFF + annotation search
gene_list <- data.frame(
  gene_id = c(
    "LOC_00009210",   # DNMT1
    "LOC_00008803",   # TET2
    "LOC_00017436",   # TET3
    "LOC_00005931",   # MBD2
    "LOC_00023929",   # MBD4
    "LOC_00011698",   # UHRF1_1
    "LOC_00017302",   # UHRF1_2
    "LOC_00018573",   # TDG_1
    "LOC_00024530",   # TDG_2
    "LOC_00012796",   # GADD45B
    "LOC_00002628",   # GADD45G_1
    "LOC_00012923",   # GADD45G_2
    "LOC_00002413",   # Gadd45gip1
    "LOC_00000327",   # apex1
    "LOC_00010497",   # HELLS_1
    "LOC_00013941",   # HELLS_2
    "LOC_00003192",   # HDAC2
    "LOC_00024130"    # CFP1 (cpf1)
  ),
  gene_name = c(
    "DNMT1", "TET2", "TET3", "MBD2", "MBD4",
    "UHRF1_1", "UHRF1_2",
    "TDG_1", "TDG_2",
    "GADD45B", "GADD45G_1", "GADD45G_2", "GADD45GIP1",
    "APEX1", "HELLS_1", "HELLS_2",
    "HDAC2", "CFP1"
  ),
  category = c(
    "Writer", "Eraser", "Eraser", "Reader", "Reader",
    "Recruiter", "Recruiter",
    "BER", "BER",
    "Demethylation", "Demethylation", "Demethylation", "Demethylation",
    "BER", "Remodeler", "Remodeler",
    "Deacetylase", "Reader"
  ),
  stringsAsFactors = FALSE
)

# Notable ABSENT genes for presence/absence plot
absent_genes <- data.frame(
  gene_name = c("DNMT3A", "DNMT3B", "DNMT3L", "TET1", "MeCP2", "MBD1", "MBD3"),
  category  = c("Writer", "Writer", "Writer", "Eraser", "Reader", "Reader", "Reader"),
  stringsAsFactors = FALSE
)

# Verify present genes exist in GFF
in_gff <- gene_list$gene_id %in% gff_gene_ids
cat("Genes found in GFF (chr1-31):", sum(in_gff), "/", nrow(gene_list), "\n")
if (any(!in_gff)) cat("Not on chr1-31:", gene_list$gene_name[!in_gff], "\n")

# =============================================================================
# 2. BUILD SAMPLE TABLE + DESeq2
# =============================================================================

cat("\n=== Building sample table ===\n")

files <- list.files(dir_counts, pattern = "htseq_gene_counts\\.txt$", full.names = TRUE)
files <- files[!grepl(" - Copy", files)]
stopifnot(length(files) > 0)
bn <- basename(files)

df <- data.frame(
  file      = files,
  fileName  = bn,
  sample_id = tools::file_path_sans_ext(bn),
  stringsAsFactors = FALSE
)

infer_row <- function(b) {
  experiment <- NA_character_; tissue <- NA_character_; condition <- NA_character_
  if (str_starts(b, "T\\d+S\\d+")) {
    experiment <- "tail_amputation"; tissue <- "tail"; condition <- "amputated"
  } else if (str_starts(b, "C\\d+S\\d+")) {
    experiment <- "tail_amputation"; tissue <- "tail"; condition <- "control"
  } else if (str_starts(b, "R\\d+_")) {
    experiment <- "eye_amputation"; tissue <- "eye"; condition <- "amputated"
  } else if (str_starts(b, "C\\d+_")) {
    experiment <- "eye_amputation"; tissue <- "eye"; condition <- "control"
  } else if (str_starts(b, "irrep")) {
    experiment <- "bodywall_irradiation"; tissue <- "bodywall"; condition <- "irradiated"
  } else if (str_starts(b, "dcrep")) {
    experiment <- "bodywall_irradiation"; tissue <- "bodywall"; condition <- "control"
  } else if (str_starts(b, "fungicide_l30")) {
    experiment <- "bodywall_fungicide"; tissue <- "bodywall"; condition <- "fungicide"
  } else if (str_starts(b, "fungicide_l0")) {
    experiment <- "bodywall_fungicide"; tissue <- "bodywall"; condition <- "control"
  } else if (str_detect(b, "rmoverrep_fpDLHead")) {
    experiment <- "baseline_tissues"; tissue <- "head"; condition <- "control"
  } else if (str_detect(b, "rmoverrep_fpDLJuv")) {
    experiment <- "baseline_tissues"; tissue <- "juvenile"; condition <- "control"
  } else if (str_detect(b, "rmoverrep_fpDLOvo")) {
    experiment <- "baseline_tissues"; tissue <- "ovotestis"; condition <- "control"
  }
  list(experiment = experiment, tissue = tissue, condition = condition)
}

parsed <- do.call(rbind, lapply(bn, function(x) as.data.frame(infer_row(x), stringsAsFactors = FALSE)))
df$experiment <- parsed$experiment
df$tissue     <- parsed$tissue
df$condition  <- parsed$condition
df$experiment <- factor(df$experiment,
  levels = c("baseline_tissues", "tail_amputation", "eye_amputation",
             "bodywall_irradiation", "bodywall_fungicide"))
df$condition <- factor(df$condition,
  levels = c("control", "irradiated", "amputated", "fungicide"))
df <- df %>% mutate(
  sample_id = str_remove(sample_id, "_htseq_gene_counts"),
  sample_id = str_remove(sample_id, ".Aligned.out.bam"),
  sample_id = str_remove(sample_id, "rmoverrep_fpDL")
)
rownames(df) <- df$sample_id

cat("Samples per tissue:\n")
print(table(df$tissue, df$condition))

# --- Run DESeq2 ---
cat("\n=== Running DESeq2 ===\n")
sampleTable <- df[, c("sample_id", "fileName", "condition", "tissue", "experiment")]
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                  directory = dir_counts,
                                  design = ~ tissue + condition)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
vst_mat <- assay(vsd)

cat("DESeq2 done. Genes retained:", nrow(vst_mat), "\n")

# =============================================================================
# 3. EXTRACT EXPRESSION DATA
# =============================================================================

present_in_counts <- gene_list$gene_id %in% rownames(vst_mat)
cat("\nGenes in count matrix:", sum(present_in_counts), "/", nrow(gene_list), "\n")
if (any(!present_in_counts)) cat("Missing:", gene_list$gene_name[!present_in_counts], "\n")
gene_list <- gene_list[present_in_counts, ]

# Control samples only, averaged by tissue
ctrl_samples <- df$sample_id[df$condition == "control"]
vst_ctrl <- vst_mat[gene_list$gene_id, ctrl_samples, drop = FALSE]

tissue_map <- df[ctrl_samples, "tissue"]
tissue_order <- c("head", "eye", "bodywall", "tail", "juvenile", "ovotestis")

avg_by_tissue <- do.call(cbind, lapply(tissue_order, function(t) {
  samps <- ctrl_samples[tissue_map == t]
  if (length(samps) > 0) rowMeans(vst_ctrl[, samps, drop = FALSE]) else rep(NA, nrow(vst_ctrl))
}))
colnames(avg_by_tissue) <- tissue_order
rownames(avg_by_tissue) <- gene_list$gene_name

# =============================================================================
# 4. PLOT 1: Genome Presence/Absence
# =============================================================================

cat("\n=== Generating Plot 1: Genome Presence/Absence ===\n")

# Combine present + absent genes
all_names <- c(gene_list$gene_name, absent_genes$gene_name)
all_cats  <- c(gene_list$category, absent_genes$category)
all_status <- c(rep("Present", nrow(gene_list)), rep("Absent", nrow(absent_genes)))

pres_df <- data.frame(
  gene     = factor(all_names, levels = rev(all_names)),
  category = all_cats,
  status   = factor(all_status, levels = c("Present", "Absent"))
)

# Category order for faceting
cat_order <- c("Writer", "Eraser", "Reader", "Recruiter",
               "BER", "Demethylation", "Remodeler", "Deacetylase")
pres_df$category <- factor(pres_df$category, levels = cat_order)

p1 <- ggplot(pres_df, aes(x = 1, y = gene, fill = status)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = ifelse(status == "Present", "\u2713", "\u2717")),
            size = 5, color = "white", fontface = "bold") +
  facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = c(Present = "#27AE60", Absent = "#2C3E50"),
                    name = "Genome Status") +
  labs(
    title = expression(paste("DNA Methylation Toolkit in ", italic("D. laeve"))),
    subtitle = "Gene presence in genome assembly (chr1\u201331)",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.text.y      = element_text(face = "italic", size = 11),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
    strip.placement  = "outside",
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 15),
    plot.subtitle    = element_text(size = 10, color = "grey40"),
    panel.spacing    = unit(0.3, "lines")
  )

ggsave(file.path(out_dir, "png/methylation_toolkit_presence.png"),
       p1, width = 7, height = 10, dpi = 300)
tmp1 <- file.path(tempdir(), "methylation_toolkit_presence.pdf")
ggsave(tmp1, p1, width = 7, height = 10, device = cairo_pdf)
if (!file.copy(tmp1, file.path(out_dir, "pdf/methylation_toolkit_presence.pdf"), overwrite = TRUE))
  cat("  WARNING: PDF locked, saved to:", tmp1, "\n")

cat("Plot 1 saved.\n")

# =============================================================================
# 5. PLOT 2: Expression by Tissue (boxplot, outliers filtered)
# =============================================================================

cat("=== Generating Plot 2: Expression by Tissue ===\n")

# Build long-format data from per-sample VST (control only)
expr_long <- as.data.frame(vst_ctrl)
expr_long$gene_id <- rownames(expr_long)
expr_long <- merge(expr_long, gene_list[, c("gene_id", "gene_name", "category")], by = "gene_id")
expr_long <- pivot_longer(expr_long,
  cols = all_of(ctrl_samples),
  names_to = "sample_id",
  values_to = "vst_expression"
)
expr_long$tissue <- df[expr_long$sample_id, "tissue"]
expr_long$tissue <- factor(expr_long$tissue, levels = tissue_order)

# Filter outliers per gene × tissue (IQR method)
expr_long <- expr_long %>%
  group_by(gene_name, tissue) %>%
  mutate(
    q1  = quantile(vst_expression, 0.25),
    q3  = quantile(vst_expression, 0.75),
    iqr = q3 - q1,
    is_outlier = vst_expression < (q1 - 1.5 * iqr) | vst_expression > (q3 + 1.5 * iqr)
  ) %>%
  filter(!is_outlier) %>%
  ungroup() %>%
  select(-q1, -q3, -iqr, -is_outlier)

cat("After outlier filter:", nrow(expr_long), "data points\n")

# Pretty labels
tissue_labels <- c(head = "Head", eye = "Eye", bodywall = "Body Wall",
                   tail = "Tail", juvenile = "Juvenile", ovotestis = "Ovotestis")
expr_long$tissue_label <- tissue_labels[as.character(expr_long$tissue)]
expr_long$tissue_label <- factor(expr_long$tissue_label, levels = tissue_labels[tissue_order])

# Gene order matching plot 1
expr_long$gene_name <- factor(expr_long$gene_name, levels = gene_list$gene_name)
expr_long$category <- factor(expr_long$category, levels = cat_order)

tissue_colors <- c(
  "Head" = "#8E44AD", "Eye" = "#2471A3", "Body Wall" = "#1ABC9C",
  "Tail" = "#C0392B", "Juvenile" = "#F39C12", "Ovotestis" = "#27AE60"
)

p2 <- ggplot(expr_long, aes(x = gene_name, y = vst_expression, fill = tissue_label)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.8),
               linewidth = 0.3) +
  facet_grid(~ category, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = tissue_colors, name = "Tissue") +
  labs(
    title = "Methylation Machinery Expression Across Tissues",
    subtitle = "VST-normalized counts (control samples, outliers removed)",
    x = NULL,
    y = "VST Expression"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, face = "italic", size = 9),
    strip.text          = element_text(face = "bold", size = 10),
    panel.grid.major.x  = element_blank(),
    legend.position     = "bottom",
    plot.title          = element_text(face = "bold", size = 14),
    plot.subtitle       = element_text(size = 10, color = "grey40")
  )

ggsave(file.path(out_dir, "png/methylation_machinery_expression.png"),
       p2, width = 16, height = 6, dpi = 300)
tmp2 <- file.path(tempdir(), "methylation_machinery_expression.pdf")
ggsave(tmp2, p2, width = 16, height = 6, device = cairo_pdf)
if (!file.copy(tmp2, file.path(out_dir, "pdf/methylation_machinery_expression.pdf"), overwrite = TRUE))
  cat("  WARNING: PDF locked, saved to:", tmp2, "\n")

cat("Plot 2 saved.\n")

# =============================================================================
# 5b. PLOT 3: Heatmap (z-scored VST, pheatmap)
# =============================================================================

cat("=== Generating Plot 3: Expression Heatmap ===\n")

mat_z <- t(scale(t(avg_by_tissue)))
colnames(mat_z) <- c("Head", "Eye", "Body Wall", "Tail", "Juvenile", "Ovotestis")

anno_row <- data.frame(
  Category = gene_list$category[match(rownames(mat_z), gene_list$gene_name)],
  row.names = rownames(mat_z)
)

cat_colors <- list(Category = c(
  Writer = "#C0392B", Eraser = "#2471A3", Reader = "#27AE60", Recruiter = "#F39C12",
  BER = "#8E44AD", Demethylation = "#1ABC9C", Remodeler = "#E67E22", Deacetylase = "#95A5A6"
))

p3 <- pheatmap(mat_z,
  annotation_row = anno_row,
  annotation_colors = cat_colors,
  cluster_cols = FALSE,
  color = colorRampPalette(c("#2471A3", "white", "#C0392B"))(100),
  main = "DNA Methylation Toolkit \u2014 Expression Across Tissues",
  fontsize = 11,
  fontsize_row = 10,
  silent = TRUE
)

png(file.path(out_dir, "png/methylation_toolkit_heatmap.png"), width = 8, height = 8, units = "in", res = 300)
grid::grid.draw(p3$gtable); dev.off()

tmp_pdf <- file.path(tempdir(), "methylation_toolkit_heatmap.pdf")
cairo_pdf(tmp_pdf, width = 8, height = 8)
grid::grid.draw(p3$gtable); dev.off()
ok <- file.copy(tmp_pdf, file.path(out_dir, "pdf/methylation_toolkit_heatmap.pdf"), overwrite = TRUE)
if (!ok) cat("  WARNING: could not overwrite PDF (file locked?) — saved to:", tmp_pdf, "\n")

cat("Plot 3 saved.\n")

# =============================================================================
# 6. SUMMARY TABLE + HTML REPORT
# =============================================================================

cat("\n=== Building summary table ===\n")

col_labels <- c(head = "Head", eye = "Eye", bodywall = "Body Wall",
                tail = "Tail", juvenile = "Juvenile", ovotestis = "Ovotestis")
summary_table <- as.data.frame(avg_by_tissue)
colnames(summary_table) <- col_labels[tissue_order]
summary_table$Gene <- rownames(summary_table)
summary_table$Category <- gene_list$category[match(rownames(summary_table), gene_list$gene_name)]
summary_table$Gene_ID <- gene_list$gene_id[match(rownames(summary_table), gene_list$gene_name)]
summary_table <- summary_table[, c("Gene", "Gene_ID", "Category", col_labels[tissue_order])]
summary_table <- summary_table[order(match(summary_table$Category, cat_order), summary_table$Gene), ]

cat("\nMethylation machinery expression (mean VST, control):\n")
print(summary_table, row.names = FALSE)

# --- HTML report ---
cat("\n=== Writing HTML report ===\n")

html <- paste0('<!DOCTYPE html>
<html><head>
<meta charset="UTF-8">
<title>Batch 02: Methylation Machinery</title>
<style>
  body { font-family: "Segoe UI", Arial, sans-serif; max-width: 1100px; margin: 40px auto; padding: 0 20px; color: #333; }
  h1 { color: #2c3e50; border-bottom: 2px solid #2471A3; padding-bottom: 10px; }
  h2 { color: #2471A3; margin-top: 30px; }
  table { border-collapse: collapse; width: 100%; margin: 15px 0; font-size: 13px; }
  th { background: #2471A3; color: white; padding: 8px 10px; text-align: left; }
  td { padding: 5px 10px; border-bottom: 1px solid #ddd; }
  tr:nth-child(even) { background: #f8f9fa; }
  img { max-width: 100%; border: 1px solid #ddd; margin: 10px 0; }
  .key-finding { background: #eaf2f8; padding: 12px; border-left: 4px solid #2471A3; margin: 15px 0; }
</style>
</head><body>
<h1>Batch 02: Methylation Machinery in <em>D. laeve</em></h1>
<p><strong>Question:</strong> Does <em>D. laeve</em> possess a functional DNA methylation toolkit, and how is it expressed across tissues?</p>

<h2>1. Genome Presence/Absence</h2>
<p>Presence (green) or absence (black) of methylation machinery genes in the <em>D. laeve</em> genome (chr1&ndash;31). Notable absences: DNMT3A/B/L (no de novo methyltransferase), TET1, MeCP2, MBD1, MBD3.</p>
<img src="png/methylation_toolkit_presence.png" alt="Presence/absence">

<h2>2. Expression Across Tissues</h2>
<p>VST-normalized expression of 18 methylation machinery genes across 6 tissues (control samples only, IQR outliers removed).</p>
<img src="png/methylation_machinery_expression.png" alt="Expression by tissue">

<h2>3. Expression Heatmap</h2>
<p>Z-scored VST expression of all 18 methylation machinery genes. Rows clustered by expression similarity, columns ordered anatomically.</p>
<img src="png/methylation_toolkit_heatmap.png" alt="Heatmap">

<h2>4. Expression Summary</h2>
<table>
<tr><th>Gene</th><th>ID</th><th>Category</th>',
paste0('<th>', col_labels[tissue_order], '</th>', collapse = ''), '</tr>\n')

for (i in seq_len(nrow(summary_table))) {
  row <- summary_table[i, ]
  html <- paste0(html, '<tr>',
    '<td><em>', row$Gene, '</em></td>',
    '<td>', row$Gene_ID, '</td>',
    '<td>', row$Category, '</td>',
    paste0('<td>', round(as.numeric(row[4:ncol(summary_table)]), 2), '</td>', collapse = ''),
    '</tr>\n')
}

html <- paste0(html, '</table>

<h2>5. Key Findings</h2>
<div class="key-finding">
<ul>
<li><strong>DNMT1</strong> present, <strong>no DNMT3</strong> &mdash; maintenance-only methylation, no de novo writer.</li>
<li><strong>TET2 + TET3</strong> &mdash; active demethylation pathway complete.</li>
<li><strong>TDG + APEX1 + GADD45 family</strong> &mdash; base excision repair (BER) demethylation pathway present.</li>
<li><strong>HELLS (×2)</strong> &mdash; chromatin remodeler that facilitates DNA methylation present.</li>
<li><strong>UHRF1 (×2) + CFP1</strong> &mdash; methylation readers/recruiters present.</li>
<li><strong>Absent: DNMT3A/B/L, TET1, MeCP2, MBD1, MBD3</strong> &mdash; streamlined toolkit compared to vertebrates.</li>
</ul>
</div>
</body></html>')

writeLines(html, file.path(out_dir, "methylation_machinery_report.html"))

cat("\n=== Batch 02 complete ===\n")
cat("Output:", out_dir, "\n")
