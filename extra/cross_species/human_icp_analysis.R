#!/usr/bin/env Rscript
# =============================================================================
# Cross-species ICP validation: Does the ICP methylation-expression pattern
# seen in D. laeve also exist in humans?
#
# Question: In D. laeve, 85% of genes are ICP. ICP gene body methylation
# correlates with expression (rho=0.407) but promoter methylation does NOT
# (rho=0.013). Is this a general ICP property or slug-specific?
#
# Approach:
# 1. Classify human genes into HCP/ICP/LCP using Weber criteria (same as batch03)
# 2. Download Roadmap Epigenomics WGBS + RNA-seq for ~33 tissues
# 3. Compute gene body + promoter methylation per gene per tissue
# 4. Correlate with expression, stratified by promoter class
#
# Data: Roadmap Epigenomics (Kundaje et al. 2015, Nature)
# Output: extra/cross_species/figures/ + data/
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(rtracklayer)
  library(ggplot2)
  library(Biostrings)
})

t0 <- proc.time()[3]
BASE <- "/mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave/extra/cross_species"
dir.create(file.path(BASE, "figures"), showWarnings = FALSE)
dir.create(file.path(BASE, "data"), showWarnings = FALSE)
dir.create(file.path(BASE, "bigwig"), showWarnings = FALSE)

save_fig <- function(p, name, w = 10, h = 7) {
  png(file.path(BASE, "figures", paste0(name, ".png")), width = w, height = h,
      units = "in", res = 300)
  print(p)
  dev.off()
  cairo_pdf(file.path(BASE, "figures", paste0(name, ".pdf")), width = w, height = h)
  print(p)
  dev.off()
  cat(sprintf("  Saved: %s\n", name))
}

# =============================================================================
# 1. CLASSIFY HUMAN PROMOTERS (Weber criteria, same as D. laeve batch03)
# =============================================================================
cat("=== Step 1: Classify human promoters ===\n")

genome <- BSgenome.Hsapiens.UCSC.hg19
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Get gene models — one per gene (longest transcript)
genes <- genes(txdb)
# Keep standard chromosomes only
keep_chr <- paste0("chr", c(1:22, "X"))
genes <- genes[seqnames(genes) %in% keep_chr]
cat(sprintf("  Genes on standard chromosomes: %d\n", length(genes)))

# Promoter regions: TSS ± 500bp (Weber uses ±500 for CpG O/E calculation)
promoters <- promoters(genes, upstream = 500, downstream = 500)

# Get promoter sequences
prom_seqs <- getSeq(genome, promoters)
cat(sprintf("  Promoter sequences extracted: %d\n", length(prom_seqs)))

# Compute CpG O/E and GC% for each promoter
cat("  Computing CpG O/E and GC%...\n")
df <- dinucleotideFrequency(prom_seqs)
nf <- oligonucleotideFrequency(prom_seqs, width = 1)

n_C <- nf[, "C"]
n_G <- nf[, "G"]
n_CG <- df[, "CG"]
seq_len <- width(prom_seqs)

gc_pct <- (n_C + n_G) / seq_len
cpg_obs <- n_CG / seq_len
cpg_exp <- (n_C / seq_len) * (n_G / seq_len)
cpg_oe <- ifelse(cpg_exp > 0, cpg_obs / cpg_exp, 0)

# Weber classification (same thresholds as D. laeve batch03)
prom_class <- ifelse(cpg_oe >= 0.75 & gc_pct >= 0.55, "HCP",
             ifelse(cpg_oe < 0.48, "LCP", "ICP"))

prom_dt <- data.table(
  gene_id = names(promoters),
  chr = as.character(seqnames(promoters)),
  cpg_oe = round(cpg_oe, 4),
  gc_pct = round(gc_pct, 4),
  prom_class = prom_class
)

cat("  Promoter classification:\n")
print(prom_dt[, .(.N, pct = round(100 * .N / nrow(prom_dt), 1)), by = prom_class][order(-N)])

# =============================================================================
# 2. DEFINE GENE BODY AND PROMOTER REGIONS FOR METHYLATION
# =============================================================================
cat("\n=== Step 2: Define gene body + promoter regions ===\n")

# Gene body = entire gene minus promoter
gene_body <- genes
# Promoter for methylation = TSS - 2kb to TSS
prom_meth <- promoters(genes, upstream = 2000, downstream = 0)

cat(sprintf("  Gene bodies: %d | Promoter regions: %d\n",
            length(gene_body), length(prom_meth)))

# =============================================================================
# 3. LOAD RNA-SEQ DATA
# =============================================================================
cat("\n=== Step 3: Load Roadmap RNA-seq ===\n")

rna <- fread(cmd = paste0("zcat ", file.path(BASE, "roadmap_rpkm.pc.gz")))
setnames(rna, 1, "gene_id")
cat(sprintf("  RNA matrix: %d genes x %d samples\n", nrow(rna), ncol(rna) - 1))

# Load epigenome names
eg_names <- fread(file.path(BASE, "epigenome_names.txt"), header = FALSE,
                  col.names = c("eid", "tissue"))

# Epigenomes with both WGBS and RNA
both <- fread(file.path(BASE, "both_epigenomes.txt"), header = FALSE)$V1
cat(sprintf("  Epigenomes with WGBS + RNA: %d\n", length(both)))

# Map Ensembl IDs to Entrez (gene_id in TxDb is Entrez)
# RNA data uses Ensembl IDs; TxDb uses Entrez
# We need biomaRt or org.Hs.eg.db for mapping
# Try a simpler approach: use the Ensembl ID directly and match later
rna[, ensembl := sub("\\..*$", "", gene_id)]  # strip version

# =============================================================================
# 4. DOWNLOAD WGBS BIGWIGS AND COMPUTE GENE BODY METHYLATION
# =============================================================================
cat("\n=== Step 4: Download WGBS and compute methylation ===\n")

# Pick a representative subset of diverse tissues (not all 33 — too slow)
# Select: stem cell, brain, liver, heart, lung, intestine, muscle, ovary, spleen
target_eids <- c("E003", "E066", "E071", "E095", "E096",
                 "E084", "E100", "E097", "E113", "E058")
target_eids <- intersect(target_eids, both)
cat(sprintf("  Target tissues: %d\n", length(target_eids)))
for (e in target_eids) {
  cat(sprintf("    %s: %s\n", e, eg_names[eid == e, tissue]))
}

# Download bigwigs (if not already present)
bw_base <- "https://egg2.wustl.edu/roadmap/data/byDataType/dnamethylation/WGBS/FractionalMethylation_bigwig"
for (eid in target_eids) {
  bw_file <- file.path(BASE, "bigwig", paste0(eid, "_WGBS_FractionalMethylation.bigwig"))
  if (!file.exists(bw_file)) {
    url <- paste0(bw_base, "/", eid, "_WGBS_FractionalMethylation.bigwig")
    cat(sprintf("  Downloading %s...\n", eid))
    tryCatch(
      download.file(url, bw_file, mode = "wb", quiet = TRUE),
      error = function(e) cat(sprintf("    FAILED: %s\n", e$message))
    )
  } else {
    cat(sprintf("  %s already downloaded (%.0f MB)\n", eid, file.size(bw_file) / 1e6))
  }
}

# Compute mean methylation per gene body and promoter for each tissue
cat("\n  Computing gene body + promoter methylation...\n")

# Reduce gene set to those in prom_dt
gene_ids_use <- prom_dt$gene_id

meth_results <- list()
for (eid in target_eids) {
  bw_file <- file.path(BASE, "bigwig", paste0(eid, "_WGBS_FractionalMethylation.bigwig"))
  if (!file.exists(bw_file) || file.size(bw_file) < 1000) {
    cat(sprintf("  Skipping %s (file missing)\n", eid))
    next
  }
  cat(sprintf("  Processing %s (%s)...\n", eid, eg_names[eid == eid, tissue][1]))

  bw <- tryCatch(BigWigFile(bw_file), error = function(e) { cat("    BigWig error\n"); NULL })
  if (is.null(bw)) next

  # Gene body methylation (mean over gene body)
  gb_meth <- tryCatch({
    import(bw, which = gene_body, as = "NumericList")
  }, error = function(e) NULL)

  if (is.null(gb_meth)) {
    # Fallback: use summary
    gb_summary <- tryCatch(
      summary(bw, gene_body, type = "mean"),
      error = function(e) NULL
    )
    if (!is.null(gb_summary)) {
      gb_means <- unlist(lapply(gb_summary, function(x) x$score[1]))
    } else {
      cat("    Failed to read gene body methylation\n")
      next
    }
  } else {
    gb_means <- vapply(gb_meth, mean, numeric(1), na.rm = TRUE)
  }

  # Promoter methylation
  pm_meth <- tryCatch({
    import(bw, which = prom_meth, as = "NumericList")
  }, error = function(e) NULL)

  if (is.null(pm_meth)) {
    pm_summary <- tryCatch(
      summary(bw, prom_meth, type = "mean"),
      error = function(e) NULL
    )
    if (!is.null(pm_summary)) {
      pm_means <- unlist(lapply(pm_summary, function(x) x$score[1]))
    } else {
      pm_means <- rep(NA_real_, length(gene_body))
    }
  } else {
    pm_means <- vapply(pm_meth, mean, numeric(1), na.rm = TRUE)
  }

  meth_results[[eid]] <- data.table(
    gene_id = names(gene_body),
    eid = eid,
    tissue = eg_names[eid == eid, tissue][1],
    gb_meth = gb_means,
    prom_meth = pm_means
  )
  cat(sprintf("    Done: %d genes with methylation data\n", sum(!is.na(gb_means))))
}

if (length(meth_results) == 0) {
  cat("ERROR: No methylation data extracted. Check bigwig files.\n")
  quit(status = 1)
}

meth_dt <- rbindlist(meth_results)
cat(sprintf("\n  Total methylation measurements: %s\n", format(nrow(meth_dt), big.mark = ",")))

# =============================================================================
# 5. MERGE METHYLATION + EXPRESSION + PROMOTER CLASS
# =============================================================================
cat("\n=== Step 5: Merge all data ===\n")

# We need to map between Entrez (TxDb gene_id) and Ensembl (RNA gene_id)
# Use biomaRt
library(biomaRt)
cat("  Querying biomaRt for Entrez-Ensembl mapping...\n")
mart <- tryCatch(
  useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://grch37.ensembl.org"),
  error = function(e) {
    cat("  biomaRt failed, trying alternative...\n")
    useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }
)

id_map <- tryCatch(
  getBM(attributes = c("entrezgene_id", "ensembl_gene_id"),
        filters = "entrezgene_id",
        values = unique(prom_dt$gene_id),
        mart = mart),
  error = function(e) { cat("  biomaRt query failed\n"); NULL }
)

if (!is.null(id_map) && nrow(id_map) > 0) {
  id_map <- as.data.table(id_map)
  setnames(id_map, c("gene_id", "ensembl"))
  id_map[, gene_id := as.character(gene_id)]
  cat(sprintf("  ID mapping: %d Entrez → Ensembl pairs\n", nrow(id_map)))
} else {
  cat("  WARNING: No ID mapping available. Trying numeric matching.\n")
  id_map <- data.table(gene_id = character(0), ensembl = character(0))
}

# Merge methylation with promoter class
meth_class <- merge(meth_dt, prom_dt[, .(gene_id, prom_class, cpg_oe, gc_pct)],
                    by = "gene_id", all.x = TRUE)

# Merge with expression via ID mapping
meth_class <- merge(meth_class, id_map, by = "gene_id", all.x = TRUE, allow.cartesian = TRUE)

# For each tissue-gene pair, get expression
rna_long <- melt(rna, id.vars = c("gene_id", "ensembl"),
                 variable.name = "eid", value.name = "rpkm")
rna_long[, eid := as.character(eid)]

final <- merge(meth_class, rna_long[, .(ensembl, eid, rpkm)],
               by = c("ensembl", "eid"), all.x = TRUE)
final <- final[!is.na(rpkm) & !is.na(gb_meth) & !is.na(prom_class)]
cat(sprintf("  Final merged dataset: %s gene-tissue pairs\n",
            format(nrow(final), big.mark = ",")))
cat(sprintf("  Unique genes: %d | Tissues: %d\n",
            uniqueN(final$gene_id), uniqueN(final$tissue)))

# Add log expression
final[, log_rpkm := log2(rpkm + 1)]

# Save merged data
fwrite(final, file.path(BASE, "data/human_icp_meth_expression.tsv"), sep = "\t")

# =============================================================================
# 6. COMPUTE CORRELATIONS BY PROMOTER CLASS
# =============================================================================
cat("\n=== Step 6: Correlations by promoter class ===\n")

# Per-tissue, per-class correlations
cor_results <- final[, {
  gb_cor <- tryCatch(cor.test(gb_meth, log_rpkm, method = "spearman"),
                     error = function(e) list(estimate = NA, p.value = NA))
  pm_cor <- tryCatch(cor.test(prom_meth, log_rpkm, method = "spearman"),
                     error = function(e) list(estimate = NA, p.value = NA))
  .(gb_rho = as.numeric(gb_cor$estimate),
    gb_p = gb_cor$p.value,
    pm_rho = as.numeric(pm_cor$estimate),
    pm_p = pm_cor$p.value,
    n_genes = .N)
}, by = .(tissue, prom_class)]

cat("\nGene body methylation vs expression (Spearman rho):\n")
gb_wide <- dcast(cor_results, tissue ~ prom_class, value.var = "gb_rho")
print(gb_wide)

cat("\nPromoter methylation vs expression (Spearman rho):\n")
pm_wide <- dcast(cor_results, tissue ~ prom_class, value.var = "pm_rho")
print(pm_wide)

fwrite(cor_results, file.path(BASE, "data/human_icp_correlations.tsv"), sep = "\t")

# =============================================================================
# 7. FIGURES
# =============================================================================
cat("\n=== Step 7: Figures ===\n")

# Fig 1: Gene body meth vs expression, faceted by promoter class
# Use one representative tissue (Adult_Liver as it's well-characterized)
liver <- final[grepl("Liver", tissue)]
if (nrow(liver) > 0) {
  p1 <- ggplot(liver[prom_class %in% c("HCP", "ICP", "LCP")],
               aes(x = gb_meth, y = log_rpkm)) +
    geom_point(alpha = 0.05, size = 0.3) +
    geom_smooth(method = "gam", color = "red", se = FALSE) +
    facet_wrap(~prom_class, scales = "free") +
    labs(x = "Gene body methylation (fractional)",
         y = "Expression (log2 RPKM + 1)",
         title = "Human Adult Liver: Gene body methylation vs expression",
         subtitle = "Stratified by Weber promoter class (HCP/ICP/LCP)") +
    theme_minimal(base_size = 12)
  save_fig(p1, "fig_human_gb_meth_vs_expr_by_class", w = 14, h = 5)
}

# Fig 2: Promoter meth vs expression, faceted by class
if (nrow(liver) > 0) {
  p2 <- ggplot(liver[prom_class %in% c("HCP", "ICP", "LCP")],
               aes(x = prom_meth, y = log_rpkm)) +
    geom_point(alpha = 0.05, size = 0.3) +
    geom_smooth(method = "gam", color = "blue", se = FALSE) +
    facet_wrap(~prom_class, scales = "free") +
    labs(x = "Promoter methylation (fractional, 2kb upstream)",
         y = "Expression (log2 RPKM + 1)",
         title = "Human Adult Liver: Promoter methylation vs expression",
         subtitle = "Stratified by Weber promoter class") +
    theme_minimal(base_size = 12)
  save_fig(p2, "fig_human_prom_meth_vs_expr_by_class", w = 14, h = 5)
}

# Fig 3: Correlation barplot across tissues — ICP vs HCP vs LCP
# Gene body
p3 <- ggplot(cor_results[!is.na(gb_rho)],
             aes(x = reorder(tissue, gb_rho), y = gb_rho, fill = prom_class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("HCP" = "#E74C3C", "ICP" = "#3498DB", "LCP" = "#2ECC71")) +
  coord_flip() +
  labs(x = NULL, y = "Spearman rho (gene body meth vs expression)",
       title = "Gene body methylation-expression correlation across human tissues",
       subtitle = "Roadmap Epigenomics WGBS + RNA-seq | Stratified by promoter class",
       fill = "Promoter\nclass") +
  theme_minimal(base_size = 11)
save_fig(p3, "fig_human_gb_correlation_by_tissue_and_class", w = 12, h = 8)

# Fig 4: Promoter meth correlation barplot
p4 <- ggplot(cor_results[!is.na(pm_rho)],
             aes(x = reorder(tissue, pm_rho), y = pm_rho, fill = prom_class)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("HCP" = "#E74C3C", "ICP" = "#3498DB", "LCP" = "#2ECC71")) +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "Spearman rho (promoter meth vs expression)",
       title = "Promoter methylation-expression correlation across human tissues",
       subtitle = "Roadmap Epigenomics | 2kb upstream of TSS",
       fill = "Promoter\nclass") +
  theme_minimal(base_size = 11)
save_fig(p4, "fig_human_prom_correlation_by_tissue_and_class", w = 12, h = 8)

# Fig 5: Side-by-side comparison with D. laeve
# D. laeve values (from batch04):
dlaeve <- data.table(
  species = "D. laeve",
  region = rep(c("Gene body", "Promoter"), each = 3),
  prom_class = rep(c("HCP", "ICP", "LCP"), 2),
  rho = c(NA, 0.407, NA,    # gene body (only ICP computed as dominant class)
          NA, 0.013, NA)     # promoter
)
# Actually use the overall values since D. laeve is ~all ICP
dlaeve <- data.table(
  species = "D. laeve (85% ICP)",
  region = c("Gene body", "Promoter"),
  rho = c(0.407, 0.013)
)

# Human mean across tissues for ICP genes
human_means <- cor_results[prom_class == "ICP",
  .(rho = mean(gb_rho, na.rm = TRUE)), by = .(region = "Gene body")]
human_means <- rbind(human_means,
  cor_results[prom_class == "ICP",
    .(rho = mean(pm_rho, na.rm = TRUE)), by = .(region = "Promoter")])
human_means[, species := "Human (ICP genes)"]

comparison <- rbindlist(list(dlaeve, human_means[, .(species, region, rho)]))

p5 <- ggplot(comparison, aes(x = region, y = rho, fill = species)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c("D. laeve (85% ICP)" = "#E67E22",
                                "Human (ICP genes)" = "#2980B9")) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "Spearman rho (methylation vs expression)",
       title = "ICP genes: methylation-expression correlation",
       subtitle = "D. laeve whole genome vs Human ICP-classified genes (Roadmap Epigenomics)",
       fill = "Species") +
  theme_minimal(base_size = 14)
save_fig(p5, "fig_cross_species_icp_comparison", w = 8, h = 6)

# Fig 6: Human promoter class distribution (pie chart like batch03)
class_counts <- prom_dt[, .N, by = prom_class]
class_counts[, pct := round(100 * N / sum(N), 1)]
p6 <- ggplot(class_counts, aes(x = "", y = N, fill = prom_class)) +
  geom_col(width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c("HCP" = "#E74C3C", "ICP" = "#3498DB", "LCP" = "#2ECC71"),
                    labels = paste0(class_counts$prom_class, " (", class_counts$pct, "%)")) +
  labs(title = "Human promoter classification (Weber criteria, hg19)",
       fill = "Class") +
  theme_void(base_size = 14)
save_fig(p6, "fig_human_promoter_class_distribution", w = 7, h = 6)

# =============================================================================
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\n=== Done (%.1f minutes) ===\n", elapsed / 60))
cat(sprintf("  %d tissues analyzed\n", length(meth_results)))
cat(sprintf("  %s gene-tissue pairs in final dataset\n", format(nrow(final), big.mark = ",")))
cat(sprintf("  Figures in: %s/figures/\n", BASE))
