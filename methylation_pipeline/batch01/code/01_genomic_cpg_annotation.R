#!/usr/bin/env Rscript
# =============================================================================
# Batch 01: Genomic CpG Annotation
# Question: What does the D. laeve CpG landscape look like?
# Output: data/ (region stats, TE stats) + figures/ (9 plots)
# =============================================================================

source("methylation_pipeline/_config.R")

library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
library(ggplot2)
library(data.table)

BATCH_DIR <- file.path(PIPE_DIR, "batch01")

cat("=== Batch 01: Genomic CpG Annotation ===\n\n")

# --- Load data ---
genome <- load_genome()
cat("Genome:", length(genome), "chromosomes\n")

gff <- load_gff()
genes <- gff[gff$type == "gene"]
exons <- gff[gff$type == "exon"]

te_data <- load_te()
te_gr <- te_data$te_gr; te_regions <- te_data$te_regions

# --- Genome stats ---
chr_sizes <- width(genome)
names(chr_sizes) <- names(genome)
total_size <- sum(as.numeric(chr_sizes))
gc_per_chr <- letterFrequency(genome, letters = "GC", as.prob = TRUE)[, 1]
overall_gc <- sum(letterFrequency(genome, letters = "GC")) / total_size

cat(sprintf("Genome: %s bp, GC: %.2f%%\n", format(total_size, big.mark = ","), overall_gc * 100))

# --- CpG counts genome-wide ---
dinuc <- dinucleotideFrequency(genome)
total_cpg <- sum(dinuc[, "CG"])
cpg_density <- total_cpg / (total_size / 1000)

# O/E from base composition
base_freq <- colSums(letterFrequency(genome, letters = c("A", "C", "G", "T")))
base_freq <- base_freq / sum(base_freq)
expected_cpg <- base_freq["C"] * base_freq["G"] * total_size
cpg_oe <- total_cpg / expected_cpg

cat(sprintf("CpGs: %s (%.2f/kb, O/E: %.4f)\n\n",
            format(total_cpg, big.mark = ","), cpg_density, cpg_oe))

# --- Region definitions ---
gene_cov <- 100 * sum(as.numeric(width(reduce(genes)))) / total_size
te_cov <- 100 * sum(as.numeric(width(te_regions))) / total_size

chr_lengths <- setNames(as.integer(chr_sizes), names(genome))
seqlengths(gff) <- chr_lengths[seqlevels(gff)]
gene_promoters <- trim(promoters(genes, upstream = 2000, downstream = 0))
gene_bodies <- reduce(genes)
exon_regions <- reduce(exons)
intron_regions <- GenomicRanges::setdiff(gene_bodies, exon_regions)
genic_plus_prom <- reduce(c(gene_bodies, gene_promoters))
strand(genic_plus_prom) <- "*"
genic_plus_prom <- reduce(genic_plus_prom)
seqlengths(genic_plus_prom) <- chr_lengths[seqlevels(genic_plus_prom)]
all_gaps <- gaps(genic_plus_prom)
intergenic_regions <- all_gaps[strand(all_gaps) == "*"]

regions_list <- list(Promoter = gene_promoters, `Gene body` = gene_bodies,
                     Exon = exon_regions, Intron = intron_regions,
                     TE = te_regions, Intergenic = intergenic_regions)

# --- CpG sites genome-wide ---
cat("Locating all CpG sites...\n")
cpg_all <- do.call(c, lapply(names(genome), function(chr_name) {
  hits <- matchPattern("CG", genome[[chr_name]])
  GRanges(seqnames = chr_name, ranges = IRanges(start = start(hits), width = 1))
}))
cat(sprintf("Total CpG sites: %s\n\n", format(length(cpg_all), big.mark = ",")))

# --- Per-region CpG stats ---
cat("Computing per-region CpG stats...\n")
cpg_by_chr <- split(cpg_all, seqnames(cpg_all))

region_stats <- do.call(rbind, lapply(names(regions_list), function(nm) {
  reg <- trim(regions_list[[nm]]); reg <- reg[width(reg) > 0]
  len <- sum(as.numeric(width(reg)))
  n_cpg <- sum(countOverlaps(reg, cpg_all) > 0)  # approximate
  # Exact count using Views
  total_cpg_r <- 0; total_c <- 0; total_g <- 0
  for (chr_name in unique(as.character(seqnames(reg)))) {
    chr_reg <- reg[seqnames(reg) == chr_name]
    if (length(chr_reg) == 0) next
    chr_cpg <- cpg_by_chr[[chr_name]]
    if (!is.null(chr_cpg)) total_cpg_r <- total_cpg_r + sum(overlapsAny(chr_cpg, chr_reg))
    v <- Views(genome[[chr_name]], start = start(chr_reg), end = end(chr_reg))
    seqs <- as(v, "DNAStringSet")
    total_c <- total_c + sum(letterFrequency(seqs, letters = "C"))
    total_g <- total_g + sum(letterFrequency(seqs, letters = "G"))
    rm(seqs); gc(verbose = FALSE)
  }
  oe <- ifelse(total_c > 0 & total_g > 0,
               (as.numeric(total_cpg_r) * len) / (as.numeric(total_c) * as.numeric(total_g)), NA)
  data.frame(region = nm, total_bp = len, cpg_count = total_cpg_r,
             cpg_per_kb = total_cpg_r / (len / 1000), cpg_oe = oe,
             enrichment = ifelse(oe > cpg_oe, "OVER", "UNDER"))
}))

cat("Region stats:\n"); print(region_stats)
save_data(region_stats, BATCH_DIR, "region_cpg_stats")

# --- Plots ---
region_colors <- COLORS$region
region_stats$region_name <- gsub("_", " ", region_stats$region)

# Fig 1A: Dinucleotide frequency
obs_counts <- colSums(dinuc); obs_freq <- obs_counts / sum(obs_counts)
combos <- expand.grid(first = names(base_freq), second = names(base_freq), stringsAsFactors = FALSE)
combos$dinuc <- paste0(combos$first, combos$second)
combos$expected <- base_freq[combos$first] * base_freq[combos$second]
exp_freq <- combos$expected[match(names(obs_freq), combos$dinuc)]

plot_df <- data.frame(dinucleotide = names(obs_freq), observed = as.numeric(obs_freq) * 100,
                      expected = as.numeric(exp_freq) * 100)
plot_df <- plot_df[order(-plot_df$observed), ]
plot_df$dinucleotide <- factor(plot_df$dinucleotide, levels = plot_df$dinucleotide)
plot_df$is_cpg <- ifelse(plot_df$dinucleotide == "CG", "CpG", "other")

p1a <- ggplot(plot_df, aes(x = dinucleotide, y = observed, fill = is_cpg)) +
  geom_col(width = 0.7) +
  geom_point(aes(y = expected), shape = 18, size = 3, color = "black") +
  scale_fill_manual(values = c("CpG" = "#C0392B", "other" = "#2471A3"), guide = "none") +
  labs(x = "Dinucleotide", y = "Frequency (%)",
       title = "Dinucleotide frequencies", subtitle = "Bars = observed, diamonds = expected") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p1a, BATCH_DIR, "fig1a_dinucleotide_freq")

# Fig 1B: CpG O/E by region
region_stats$region_f <- factor(region_stats$region_name,
                                levels = region_stats$region_name[order(region_stats$cpg_oe)])
p1b <- ggplot(region_stats, aes(x = region_f, y = cpg_oe, fill = region_name)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = cpg_oe, linetype = "dotted", color = "black") +
  annotate("text", x = 0.7, y = cpg_oe + 0.02, label = paste0("genome: ", round(cpg_oe, 3)), size = 3, hjust = 0) +
  scale_fill_manual(values = c(region_colors, "Gene body" = "#27AE60", TE = "#F39C12"), guide = "none") +
  labs(x = NULL, y = "CpG O/E ratio", title = "CpG O/E by genomic region") +
  theme_minimal()
save_fig(p1b, BATCH_DIR, "fig1b_cpg_oe_by_region", w = 7, h = 5)

# Fig 1C: CpG density by region
region_stats$region_f2 <- factor(region_stats$region_name,
                                 levels = region_stats$region_name[order(region_stats$cpg_per_kb)])
p1c <- ggplot(region_stats, aes(x = region_f2, y = cpg_per_kb, fill = region_name)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = cpg_density, linetype = "dotted", color = "black") +
  scale_fill_manual(values = c(region_colors, "Gene body" = "#27AE60", TE = "#F39C12"), guide = "none") +
  labs(x = NULL, y = "CpG per kb", title = "CpG density by genomic region") +
  theme_minimal()
save_fig(p1c, BATCH_DIR, "fig1c_cpg_density_by_region", w = 7, h = 5)

# Fig 1D: Per-chromosome O/E distributions (outlier-filtered: >= 10kb per region)
cat("Computing per-chromosome distributions...\n")
chr_region_stats <- do.call(rbind, lapply(names(regions_list), function(nm) {
  reg <- trim(regions_list[[nm]]); reg <- reg[width(reg) > 0]
  do.call(rbind, lapply(names(genome), function(chr_name) {
    chr_reg <- reg[seqnames(reg) == chr_name]
    if (length(chr_reg) == 0) return(NULL)
    len <- sum(as.numeric(width(chr_reg)))
    if (len < 10000) return(NULL)  # outlier filter
    chr_cpg <- cpg_by_chr[[chr_name]]
    n_cpg <- if (!is.null(chr_cpg)) sum(overlapsAny(chr_cpg, chr_reg)) else 0L
    v <- Views(genome[[chr_name]], start = start(chr_reg), end = end(chr_reg))
    seqs <- as(v, "DNAStringSet")
    cc <- sum(letterFrequency(seqs, letters = "C"))
    gg <- sum(letterFrequency(seqs, letters = "G"))
    oe <- ifelse(cc > 0 & gg > 0, (as.numeric(n_cpg) * len) / (as.numeric(cc) * as.numeric(gg)), NA)
    rm(seqs); gc(verbose = FALSE)
    data.frame(region = gsub("_", " ", nm), chr = chr_name, cpg_per_kb = n_cpg / (len / 1000), cpg_oe = oe)
  }))
}))
chr_region_stats <- chr_region_stats[!is.na(chr_region_stats$cpg_oe), ]
chr_region_stats$region <- factor(chr_region_stats$region,
                                  levels = c("Promoter", "Gene body", "Exon", "Intron", "TE", "Intergenic"))

p1d <- ggplot(chr_region_stats, aes(x = cpg_oe, fill = region, color = region)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = cpg_oe, linetype = "dotted", color = "black") +
  scale_fill_manual(values = c(region_colors, "Gene body" = "#27AE60", TE = "#F39C12")) +
  scale_color_manual(values = c(region_colors, "Gene body" = "#27AE60", TE = "#F39C12")) +
  labs(x = "CpG O/E ratio", y = "Density", title = "CpG O/E distribution by region",
       subtitle = "Per-chromosome (>= 10kb). Dashed = expected, dotted = genome avg") +
  theme_minimal() + theme(legend.title = element_blank())
save_fig(p1d, BATCH_DIR, "fig1d_cpg_oe_distribution")

rm(cpg_all, cpg_by_chr, chr_region_stats); gc(verbose = FALSE)

# --- TE class analysis ---
cat("Computing TE class CpG stats...\n")
te_gr$te_class <- sub("/.*", "", te_gr$class_family)
te_class_counts <- table(te_gr$te_class)
keep_classes <- names(te_class_counts[te_class_counts > 100])

te_family_stats <- do.call(rbind, lapply(keep_classes, function(cls) {
  cls_gr <- reduce(te_gr[te_gr$te_class == cls]); cls_gr <- trim(cls_gr); cls_gr <- cls_gr[width(cls_gr) > 0]
  len <- sum(as.numeric(width(cls_gr)))
  total_cpg_cls <- 0; total_c <- 0; total_g <- 0
  for (chr_name in unique(as.character(seqnames(cls_gr)))) {
    chr_reg <- cls_gr[seqnames(cls_gr) == chr_name]
    v <- Views(genome[[chr_name]], start = start(chr_reg), end = end(chr_reg))
    seqs <- as(v, "DNAStringSet")
    total_cpg_cls <- total_cpg_cls + sum(dinucleotideFrequency(seqs)[, "CG"])
    total_c <- total_c + sum(letterFrequency(seqs, letters = "C"))
    total_g <- total_g + sum(letterFrequency(seqs, letters = "G"))
    rm(seqs); gc(verbose = FALSE)
  }
  oe <- ifelse(total_c > 0 & total_g > 0,
               (as.numeric(total_cpg_cls) * len) / (as.numeric(total_c) * as.numeric(total_g)), NA)
  data.frame(te_class = cls, total_bp = len, cpg_count = total_cpg_cls,
             cpg_per_kb = total_cpg_cls / (len / 1000), cpg_oe = oe,
             n_elements = as.integer(te_class_counts[cls]))
}))
te_family_stats <- te_family_stats[order(te_family_stats$cpg_oe), ]
te_family_stats$te_class <- factor(te_family_stats$te_class, levels = te_family_stats$te_class)
save_data(te_family_stats, BATCH_DIR, "te_class_cpg_stats")

# TE O/E bar
p_te <- ggplot(te_family_stats, aes(x = te_class, y = cpg_oe, fill = cpg_oe)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = cpg_oe, linetype = "dotted") +
  scale_fill_gradient(low = "#2471A3", high = "#C0392B", guide = "none") +
  labs(x = NULL, y = "CpG O/E", title = "CpG O/E by TE class") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_te, BATCH_DIR, "fig1e_te_cpg_oe_by_class")

# TE per-chr distributions (outlier-filtered)
te_chr_stats <- do.call(rbind, lapply(keep_classes, function(cls) {
  cls_gr <- reduce(te_gr[te_gr$te_class == cls]); cls_gr <- trim(cls_gr); cls_gr <- cls_gr[width(cls_gr) > 0]
  do.call(rbind, lapply(names(genome), function(chr_name) {
    chr_reg <- cls_gr[seqnames(cls_gr) == chr_name]
    if (length(chr_reg) == 0) return(NULL)
    len <- sum(as.numeric(width(chr_reg)))
    if (len < 10000) return(NULL)
    v <- Views(genome[[chr_name]], start = start(chr_reg), end = end(chr_reg))
    seqs <- as(v, "DNAStringSet")
    n_cpg <- sum(dinucleotideFrequency(seqs)[, "CG"])
    cc <- sum(letterFrequency(seqs, letters = "C"))
    gg <- sum(letterFrequency(seqs, letters = "G"))
    oe <- ifelse(cc > 0 & gg > 0, (as.numeric(n_cpg) * len) / (as.numeric(cc) * as.numeric(gg)), NA)
    rm(seqs); gc(verbose = FALSE)
    data.frame(te_class = cls, chr = chr_name, cpg_per_kb = n_cpg / (len / 1000), cpg_oe = oe)
  }))
}))
te_chr_stats <- te_chr_stats[!is.na(te_chr_stats$cpg_oe), ]
te_chr_stats$te_class <- factor(te_chr_stats$te_class, levels = levels(te_family_stats$te_class))

p_te_dist <- ggplot(te_chr_stats, aes(x = cpg_oe, fill = te_class, color = te_class)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_vline(xintercept = cpg_oe, linetype = "dotted") +
  labs(x = "CpG O/E", y = "Density", title = "CpG O/E distribution by TE class",
       subtitle = "Per-chromosome (>= 10kb)") +
  theme_minimal() + theme(legend.title = element_blank())
save_fig(p_te_dist, BATCH_DIR, "fig1f_te_cpg_oe_distribution")

# Save genome summary
genome_summary <- data.frame(
  Metric = c("Total size (bp)", "Chromosomes", "GC content (%)", "Total CpGs",
             "CpG density (/kb)", "CpG O/E", "Genes", "Gene coverage (%)", "TE coverage (%)"),
  Value = c(format(total_size, big.mark = ","), length(genome), round(overall_gc * 100, 2),
            format(total_cpg, big.mark = ","), round(cpg_density, 2), round(cpg_oe, 4),
            length(genes), round(gene_cov, 1), round(te_cov, 1))
)
save_data(genome_summary, BATCH_DIR, "genome_summary")

cat("\n=== Batch 01 complete ===\n")
cat("Figures:", length(list.files(file.path(BATCH_DIR, "figures"), pattern = "\\.png$")), "\n")
cat("Data files:", length(list.files(file.path(BATCH_DIR, "data"), pattern = "\\.tsv$")), "\n")
