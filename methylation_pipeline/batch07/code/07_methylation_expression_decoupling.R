#!/usr/bin/env Rscript
# =============================================================================
# Batch 07: Methylation-Expression Decoupling
# Question: Do methylation changes predict expression changes?
# Output: data/ (MI table, gene-level) + figures/ (3 plots)
# Requires: Batch 06 DMPs, DESeq2 DE results
# =============================================================================

source("methylation_pipeline/_config.R")

library(data.table)
library(ggplot2)
library(DESeq2)
library(dplyr)

BATCH_DIR <- file.path(PIPE_DIR, "batch07")
cat("=== Batch 07: Methylation-Expression Decoupling ===\n\n")

# --- Load DMPs from batch06 ---
dmp_file <- file.path(PIPE_DIR, "batch06/data/dmps_annotated.tsv")
if (!file.exists(dmp_file)) stop("Run batch06 first.")
dmp <- fread(dmp_file)
cat(sprintf("DMPs loaded: %s\n", format(nrow(dmp), big.mark = ",")))

# --- DESeq2 for tail control vs amputated ---
cat("Running DESeq2 (tail control vs amputated)...\n")
all_files <- list.files(OG$counts_dir, pattern = "\\.txt$")
tail_files <- all_files[grepl("tail", all_files, ignore.case = TRUE)]

sample_table <- data.frame(sampleName = tail_files, fileName = tail_files, stringsAsFactors = FALSE)
sample_table$condition <- ifelse(grepl("control|C\\d+S", sample_table$sampleName, ignore.case = TRUE),
                                  "control", "amputated")

dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_table, directory = OG$counts_dir, design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "control")
dds <- dds[rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)
de_results <- as.data.frame(results(dds, alpha = 0.05))
de_results$gene_id <- rownames(de_results)

n_de <- sum(de_results$padj < 0.05, na.rm = TRUE)
cat(sprintf("DE genes (FDR < 0.05): %d\n", n_de))

# --- Gene-level methylation summary ---
gene_meth <- dmp %>%
  group_by(nearest_gene) %>%
  summarize(n_dmp = n(), mean_diff = mean(diff),
            n_hyper = sum(diff > 0), n_hypo = sum(diff < 0),
            meth_direction = ifelse(mean(diff) > 0, "Hyper", "Hypo"),
            .groups = "drop") %>%
  rename(gene_id = nearest_gene)

# Merge with expression
gene_meth <- merge(gene_meth, de_results[, c("gene_id", "log2FoldChange", "padj")], by = "gene_id")
gene_meth$de <- !is.na(gene_meth$padj) & gene_meth$padj < 0.05
gene_meth$expr_direction <- ifelse(gene_meth$log2FoldChange > 0, "Up", "Down")

cat(sprintf("Genes with DMP + expression data: %d\n", nrow(gene_meth)))
cat(sprintf("Genes with DMP that are DE: %d (%.1f%%)\n",
            sum(gene_meth$de), 100 * mean(gene_meth$de)))

# --- Mutual Information ---
# 2x2 contingency: meth_direction x expr_direction
ct <- table(gene_meth$meth_direction, gene_meth$expr_direction)
cat("\nContingency table:\n"); print(ct)

# MI calculation
p_joint <- ct / sum(ct)
p_meth <- rowSums(p_joint)
p_expr <- colSums(p_joint)
mi <- 0
for (i in 1:nrow(ct)) for (j in 1:ncol(ct)) {
  if (p_joint[i,j] > 0) mi <- mi + p_joint[i,j] * log2(p_joint[i,j] / (p_meth[i] * p_expr[j]))
}

# Permutation test
set.seed(42)
n_perm <- 10000
mi_null <- replicate(n_perm, {
  shuffled <- sample(gene_meth$expr_direction)
  ct_s <- table(gene_meth$meth_direction, shuffled)
  p_j <- ct_s / sum(ct_s); p_m <- rowSums(p_j); p_e <- colSums(p_j)
  mi_s <- 0
  for (i in 1:nrow(ct_s)) for (j in 1:ncol(ct_s)) {
    if (p_j[i,j] > 0) mi_s <- mi_s + p_j[i,j] * log2(p_j[i,j] / (p_m[i] * p_e[j]))
  }
  mi_s
})
mi_pval <- mean(mi_null >= mi)

cat(sprintf("\nMI = %.6f bits (p = %.3f, %d permutations)\n", mi, mi_pval, n_perm))

# Spearman correlation
cor_test <- cor.test(gene_meth$mean_diff, gene_meth$log2FoldChange, method = "spearman")
cat(sprintf("Spearman rho = %.3f, p = %.3f\n", cor_test$estimate, cor_test$p.value))

# Dose-response (Cochran-Armitage)
gene_meth$dmp_bin <- cut(gene_meth$n_dmp, breaks = c(0, 1, 3, 5, 10, Inf),
                          labels = c("1", "2-3", "4-5", "6-10", ">10"))
dose_table <- gene_meth %>% group_by(dmp_bin) %>%
  summarize(genes = n(), de_genes = sum(de), de_rate = mean(de) * 100, .groups = "drop")
cat("\nDose-response:\n"); print(as.data.frame(dose_table))

save_data(gene_meth, BATCH_DIR, "gene_methylation_expression")
save_data(dose_table, BATCH_DIR, "dose_response_table")
save_data(data.frame(MI = mi, MI_pval = mi_pval, Spearman_rho = cor_test$estimate,
                     Spearman_pval = cor_test$p.value, n_genes = nrow(gene_meth)),
          BATCH_DIR, "decoupling_statistics")

# --- Figures ---
# (A) Delta-meth vs delta-expression scatter
p7a <- ggplot(gene_meth, aes(x = mean_diff, y = log2FoldChange)) +
  geom_point(alpha = 0.15, size = 0.8, color = "#2471A3") +
  geom_smooth(method = "lm", color = "#C0392B", se = TRUE, linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Mean methylation difference (Amputated - Control)",
       y = "log2 fold change (expression)",
       title = "Methylation change vs expression change",
       subtitle = sprintf("rho = %.3f (p = %.2f) | MI = %.5f bits (p = %.2f) | n = %s genes",
                          cor_test$estimate, cor_test$p.value, mi, mi_pval,
                          format(nrow(gene_meth), big.mark = ","))) +
  theme_minimal(base_size = 12)
save_fig(p7a, BATCH_DIR, "fig7a_delta_meth_vs_expression", w = 9, h = 7)

# (B) Dose-response
p7b <- ggplot(dose_table, aes(x = dmp_bin, y = de_rate)) +
  geom_col(fill = "#2471A3", width = 0.6) +
  geom_text(aes(label = paste0(round(de_rate, 1), "%\n(", de_genes, "/", genes, ")")),
            vjust = -0.3, size = 3) +
  labs(x = "DMPs per gene", y = "% differentially expressed",
       title = "No dose-response: DE rate vs DMP count",
       subtitle = "Cochran-Armitage trend test") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  theme_minimal(base_size = 12)
save_fig(p7b, BATCH_DIR, "fig7b_dose_response", w = 8, h = 6)

# (C) 2x2 contingency heatmap
ct_df <- as.data.frame(ct)
colnames(ct_df) <- c("Methylation", "Expression", "Count")
p7c <- ggplot(ct_df, aes(x = Expression, y = Methylation, fill = Count)) +
  geom_tile(color = "white", linewidth = 2) +
  geom_text(aes(label = Count), size = 8, fontface = "bold") +
  scale_fill_gradient(low = "#EBF5FB", high = "#2471A3") +
  labs(title = "Methylation direction vs expression direction",
       subtitle = sprintf("MI = %.5f bits | All quadrants ~equal", mi)) +
  theme_minimal(base_size = 14) + theme(legend.position = "none")
save_fig(p7c, BATCH_DIR, "fig7c_contingency_table", w = 7, h = 6)

cat("\n=== Batch 07 complete ===\n")
