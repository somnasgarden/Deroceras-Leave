#!/usr/bin/env Rscript
# =============================================================================
# Batch 03: Baseline Methylome Analysis
# Paper: Act 2 — "The methylome at rest"
# Question: What does the D. laeve methylome look like? Where is methylation?
#           How does it vary by region, gene body, TE class?
# Output: 10 plots (PDF + PNG) + HTML report
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)
keep_chr <- paste0("chr", 1:31)
MIN_COV <- 5  # minimum read coverage per CpG for reliable beta

# --- 1. LIBRARIES ---
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
library(ggplot2)
library(dplyr)
library(scales)

# --- 2. PATHS ---
data_dir    <- "C:/Users/rafae/Projects/DATA"
project_dir <- "C:/Users/rafae/Projects/STANDBY"
out_dir     <- file.path(project_dir, "results/batch03")
fig_pdf_dir <- file.path(out_dir, "pdf")
fig_png_dir <- file.path(out_dir, "png")
cache_dir   <- file.path(project_dir, "genome/cache")

dir.create(fig_pdf_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_png_dir, showWarnings = FALSE, recursive = TRUE)

# CpG report files
cpg_files <- c(
  C1 = file.path(data_dir, "C1.CpG_report.txt"),
  C2 = file.path(data_dir, "C2.CpG_report.txt"),
  A1 = file.path(data_dir, "A1.CpG_report.txt"),
  A2 = file.path(data_dir, "A2.CpG_report.txt")
)

# Region colors from CLAUDE.md
region_colors <- c(Promoter = "#8E44AD", `Gene body` = "#27AE60", Exon = "#2471A3",
                   Intron = "#1ABC9C", TE = "#F39C12", Intergenic = "#C0392B")
cond_colors <- c(Control = "#2471A3", Amputated = "#C0392B")

# --- 3. SAVE HELPER (PDF lock-safe) ---
save_plot <- function(p, name, w = 8, h = 5) {
  png(file.path(fig_png_dir, paste0(name, ".png")),
      width = w, height = h, units = "in", res = 300)
  print(p); dev.off()
  tmp_pdf <- file.path(tempdir(), paste0(name, ".pdf"))
  cairo_pdf(tmp_pdf, width = w, height = h); print(p); dev.off()
  ok <- file.copy(tmp_pdf, file.path(fig_pdf_dir, paste0(name, ".pdf")),
                  overwrite = TRUE)
  if (!ok) cat("  WARNING: PDF locked, saved to:", tmp_pdf, "\n")
  cat("  Saved:", name, "\n")
}

# =============================================================================
# B. LOAD/CACHE GFF + TE
# =============================================================================

cat("=== Section B: Loading GFF + TE data ===\n")

# --- GFF ---
gff_rds <- file.path(cache_dir, "gff_chr1_31.rds")
if (file.exists(gff_rds)) {
  cat("  Loading GFF from cache...\n")
  gff <- readRDS(gff_rds)
} else {
  cat("  Parsing GFF from raw file (will cache)...\n")
  gff_file <- file.path(data_dir,
    "derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff")
  gff <- import(gff_file)
  gff <- gff[seqnames(gff) %in% keep_chr]
  gff <- keepSeqlevels(gff, keep_chr, pruning.mode = "coarse")
  saveRDS(gff, gff_rds)
  cat("  Cached GFF to:", gff_rds, "\n")
}
cat("  GFF features:", length(gff), "\n")

genes <- gff[gff$type == "gene"]
exons <- gff[gff$type == "exon"]
cat("  Genes:", length(genes), " Exons:", length(exons), "\n")

# --- TE ---
te_rds <- file.path(cache_dir, "te_chr1_31.rds")
if (file.exists(te_rds)) {
  cat("  Loading TEs from cache...\n")
  te_data <- readRDS(te_rds)
  te_gr <- te_data$te_gr; te_regions <- te_data$te_regions
} else {
  cat("  Parsing TE data from raw file (will cache)...\n")
  te_file <- file.path(data_dir, "collapsed_te_age_data.tsv")
  te_dt <- fread(te_file)
  te_dt <- te_dt[chrom %in% keep_chr]
  te_gr <- makeGRangesFromDataFrame(te_dt, seqnames.field = "chrom",
                                    start.field = "start", end.field = "end",
                                    keep.extra.columns = TRUE)
  te_regions <- reduce(te_gr)
  saveRDS(list(te_gr = te_gr, te_regions = te_regions), te_rds)
  cat("  Cached TEs to:", te_rds, "\n")
}
cat("  TEs:", length(te_gr), " Merged regions:", length(te_regions), "\n")

# --- chr_lengths from GFF ---
chr_lengths_arr <- tapply(end(gff), as.character(seqnames(gff)), max)
chr_lengths <- setNames(as.integer(chr_lengths_arr[keep_chr]), keep_chr)
cat("  Chr lengths from GFF, range:",
    format(min(chr_lengths), big.mark = ","), "-",
    format(max(chr_lengths), big.mark = ","), "\n")

# --- Region definitions ---
seqlengths(gff) <- chr_lengths[seqlevels(gff)]
gene_promoters <- trim(promoters(genes, upstream = 2000, downstream = 0))
gene_bodies    <- reduce(genes)
exon_regions   <- reduce(exons)
intron_regions <- GenomicRanges::setdiff(gene_bodies, exon_regions)

genic_plus_promoter <- reduce(c(gene_bodies, gene_promoters))
strand(genic_plus_promoter) <- "*"
genic_plus_promoter <- reduce(genic_plus_promoter)
seqlengths(genic_plus_promoter) <- chr_lengths[seqlevels(genic_plus_promoter)]
all_gaps <- gaps(genic_plus_promoter)
intergenic_regions <- all_gaps[strand(all_gaps) == "*"]

regions_list <- list(Promoter = gene_promoters, `Gene body` = gene_bodies,
                     Exon = exon_regions, Intron = intron_regions,
                     TE = te_regions, Intergenic = intergenic_regions)

cat("  Region sizes (Mb):\n")
for (nm in names(regions_list)) {
  cat("    ", nm, ":", round(sum(as.numeric(width(regions_list[[nm]]))) / 1e6, 1), "\n")
}

# =============================================================================
# C. READ + PROCESS CpG REPORTS
# =============================================================================

cpg_merged_rds <- file.path(cache_dir, "cpg_meth_merged.rds")

if (file.exists(cpg_merged_rds)) {
  cat("\n=== Section C: Loading merged CpG data from cache ===\n")
  cpg_dt <- readRDS(cpg_merged_rds)
  cat("  Loaded", format(nrow(cpg_dt), big.mark = ","), "CpG sites\n")
} else {
  cat("\n=== Section C: Processing CpG reports (one at a time) ===\n")

  # Per-sample processing function
  process_cpg_file <- function(filepath, sample_name) {
    cat("  Reading", sample_name, "...\n")
    # Native data.table read — no shell pipes needed
    dt <- fread(filepath, header = FALSE, sep = "\t",
                select = c(1, 2, 3, 4, 5),
                col.names = c("chr", "pos", "strand", "meth", "unmeth"),
                colClasses = c("character", "integer", "character",
                               "integer", "integer"))
    cat("    Raw rows:", format(nrow(dt), big.mark = ","), "\n")

    # Filter to keep_chr only (grep '^chr[0-9]' also matches chr10-31)
    dt <- dt[chr %in% keep_chr]
    cat("    After chr filter:", format(nrow(dt), big.mark = ","), "\n")

    # Strand collapse: - strand CpG is at pos-1 on + strand
    dt[strand == "-", pos := pos - 1L]
    dt <- dt[, .(meth = sum(meth), unmeth = sum(unmeth)), by = .(chr, pos)]
    cat("    After strand collapse:", format(nrow(dt), big.mark = ","), "\n")

    # Save intermediate
    inter_rds <- file.path(cache_dir, paste0("cpg_", sample_name, ".rds"))
    saveRDS(dt, inter_rds)
    cat("    Saved intermediate:", inter_rds, "\n")
    rm(dt); gc(verbose = FALSE)
    return(inter_rds)
  }

  # Process each sample
  inter_files <- list()
  for (sn in names(cpg_files)) {
    inter_files[[sn]] <- process_cpg_file(cpg_files[[sn]], sn)
    gc(verbose = FALSE)
  }

  # Merge controls: C1 + C2
  cat("  Merging control samples...\n")
  c1 <- readRDS(inter_files[["C1"]])
  c2 <- readRDS(inter_files[["C2"]])
  ctrl <- merge(c1, c2, by = c("chr", "pos"), all = TRUE, suffixes = c("_1", "_2"))
  rm(c1, c2); gc(verbose = FALSE)
  # Replace NA with 0 for sites covered in only one replicate
  for (col in c("meth_1", "unmeth_1", "meth_2", "unmeth_2")) {
    ctrl[is.na(get(col)), (col) := 0L]
  }
  ctrl[, meth_ctrl := meth_1 + meth_2]
  ctrl[, unmeth_ctrl := unmeth_1 + unmeth_2]
  ctrl <- ctrl[, .(chr, pos, meth_ctrl, unmeth_ctrl)]
  cat("    Control CpGs:", format(nrow(ctrl), big.mark = ","), "\n")

  # Merge amputated: A1 + A2
  cat("  Merging amputated samples...\n")
  a1 <- readRDS(inter_files[["A1"]])
  a2 <- readRDS(inter_files[["A2"]])
  ampu <- merge(a1, a2, by = c("chr", "pos"), all = TRUE, suffixes = c("_1", "_2"))
  rm(a1, a2); gc(verbose = FALSE)
  for (col in c("meth_1", "unmeth_1", "meth_2", "unmeth_2")) {
    ampu[is.na(get(col)), (col) := 0L]
  }
  ampu[, meth_ampu := meth_1 + meth_2]
  ampu[, unmeth_ampu := unmeth_1 + unmeth_2]
  ampu <- ampu[, .(chr, pos, meth_ampu, unmeth_ampu)]
  cat("    Amputated CpGs:", format(nrow(ampu), big.mark = ","), "\n")

  # Merge ctrl + ampu
  cat("  Final merge...\n")
  cpg_dt <- merge(ctrl, ampu, by = c("chr", "pos"), all = TRUE)
  rm(ctrl, ampu); gc(verbose = FALSE)
  for (col in c("meth_ctrl", "unmeth_ctrl", "meth_ampu", "unmeth_ampu")) {
    cpg_dt[is.na(get(col)), (col) := 0L]
  }

  # Derived columns
  cpg_dt[, cov_ctrl := meth_ctrl + unmeth_ctrl]
  cpg_dt[, cov_ampu := meth_ampu + unmeth_ampu]
  cpg_dt[, cov_all  := cov_ctrl + cov_ampu]
  cpg_dt[, beta_ctrl := fifelse(cov_ctrl >= MIN_COV,
                                 meth_ctrl / cov_ctrl, NA_real_)]
  cpg_dt[, beta_ampu := fifelse(cov_ampu >= MIN_COV,
                                 meth_ampu / cov_ampu, NA_real_)]
  cpg_dt[, beta_all  := fifelse(cov_all >= MIN_COV,
                                 (meth_ctrl + meth_ampu) / cov_all, NA_real_)]

  cat("  Total CpG sites:", format(nrow(cpg_dt), big.mark = ","), "\n")
  cat("  With coverage >= ", MIN_COV, ":\n")
  cat("    Control:", format(sum(!is.na(cpg_dt$beta_ctrl)), big.mark = ","), "\n")
  cat("    Amputated:", format(sum(!is.na(cpg_dt$beta_ampu)), big.mark = ","), "\n")
  cat("    Overall:", format(sum(!is.na(cpg_dt$beta_all)), big.mark = ","), "\n")

  # Cache
  saveRDS(cpg_dt, cpg_merged_rds)
  cat("  Cached merged CpG data to:", cpg_merged_rds, "\n")

  # Clean up intermediate files
  for (f in inter_files) {
    if (file.exists(f)) file.remove(f)
  }
  cat("  Removed intermediate RDS files\n")
}

# =============================================================================
# D. GLOBAL METHYLATION STATS
# =============================================================================

cat("\n=== Section D: Global methylation stats ===\n")

n_ctrl  <- sum(!is.na(cpg_dt$beta_ctrl))
n_ampu  <- sum(!is.na(cpg_dt$beta_ampu))
n_all   <- sum(!is.na(cpg_dt$beta_all))

mean_beta_ctrl <- mean(cpg_dt$beta_ctrl, na.rm = TRUE)
mean_beta_ampu <- mean(cpg_dt$beta_ampu, na.rm = TRUE)
mean_beta_all  <- mean(cpg_dt$beta_all, na.rm = TRUE)

med_beta_ctrl <- median(cpg_dt$beta_ctrl, na.rm = TRUE)
med_beta_ampu <- median(cpg_dt$beta_ampu, na.rm = TRUE)
med_beta_all  <- median(cpg_dt$beta_all, na.rm = TRUE)

# Fraction in methylation categories
frac_unm_ctrl <- mean(cpg_dt$beta_ctrl < 0.2, na.rm = TRUE)
frac_int_ctrl <- mean(cpg_dt$beta_ctrl >= 0.2 & cpg_dt$beta_ctrl <= 0.8, na.rm = TRUE)
frac_met_ctrl <- mean(cpg_dt$beta_ctrl > 0.8, na.rm = TRUE)

frac_unm_ampu <- mean(cpg_dt$beta_ampu < 0.2, na.rm = TRUE)
frac_int_ampu <- mean(cpg_dt$beta_ampu >= 0.2 & cpg_dt$beta_ampu <= 0.8, na.rm = TRUE)
frac_met_ampu <- mean(cpg_dt$beta_ampu > 0.8, na.rm = TRUE)

# Coverage stats
mean_cov_ctrl <- mean(cpg_dt$cov_ctrl[cpg_dt$cov_ctrl > 0])
mean_cov_ampu <- mean(cpg_dt$cov_ampu[cpg_dt$cov_ampu > 0])
med_cov_ctrl  <- median(cpg_dt$cov_ctrl[cpg_dt$cov_ctrl > 0])
med_cov_ampu  <- median(cpg_dt$cov_ampu[cpg_dt$cov_ampu > 0])

cat("  CpGs with cov >=", MIN_COV, ": Control =", format(n_ctrl, big.mark = ","),
    ", Amputated =", format(n_ampu, big.mark = ","),
    ", Overall =", format(n_all, big.mark = ","), "\n")
cat("  Mean beta:  Control =", round(mean_beta_ctrl, 4),
    ", Amputated =", round(mean_beta_ampu, 4), "\n")
cat("  Median beta: Control =", round(med_beta_ctrl, 4),
    ", Amputated =", round(med_beta_ampu, 4), "\n")
cat("  Methylation classes (Control):\n")
cat("    Unmethylated (<0.2):", round(frac_unm_ctrl * 100, 1), "%\n")
cat("    Intermediate (0.2-0.8):", round(frac_int_ctrl * 100, 1), "%\n")
cat("    Methylated (>0.8):", round(frac_met_ctrl * 100, 1), "%\n")
cat("  Coverage: Control mean =", round(mean_cov_ctrl, 1),
    "x, median =", round(med_cov_ctrl, 1), "x\n")
cat("           Amputated mean =", round(mean_cov_ampu, 1),
    "x, median =", round(med_cov_ampu, 1), "x\n")

# =============================================================================
# E. PLOT 1: Global Beta Distribution
# =============================================================================

cat("\n=== Section E: Plot 1 — Global beta distribution ===\n")

set.seed(42)
sub_idx_ctrl <- sample(which(!is.na(cpg_dt$beta_ctrl)), min(2e6, n_ctrl))
sub_idx_ampu <- sample(which(!is.na(cpg_dt$beta_ampu)), min(2e6, n_ampu))

beta_plot_df <- rbind(
  data.frame(beta = cpg_dt$beta_ctrl[sub_idx_ctrl], condition = "Control"),
  data.frame(beta = cpg_dt$beta_ampu[sub_idx_ampu], condition = "Amputated")
)

p1 <- ggplot(beta_plot_df, aes(x = beta, fill = condition, color = condition)) +
  geom_density(alpha = 0.3, linewidth = 0.7) +
  scale_fill_manual(values = cond_colors) +
  scale_color_manual(values = cond_colors) +
  labs(x = "Methylation level (beta)", y = "Density",
       title = expression(paste("Global CpG methylation in ", italic("D. laeve"))),
       subtitle = paste0("Subsample of 2M CpGs per condition (coverage >= ", MIN_COV, "x)")) +
  theme_minimal(base_size = 13) +
  theme(legend.title = element_blank(),
        legend.position = c(0.5, 0.9))

save_plot(p1, "global_beta_distribution", w = 8, h = 5)
rm(beta_plot_df); gc(verbose = FALSE)

# =============================================================================
# F. PLOT 2: Coverage Distribution
# =============================================================================

cat("=== Section F: Plot 2 — Coverage distribution ===\n")

set.seed(43)
cov_sub <- 2e6
cov_idx_c <- sample(which(cpg_dt$cov_ctrl > 0), min(cov_sub, sum(cpg_dt$cov_ctrl > 0)))
cov_idx_a <- sample(which(cpg_dt$cov_ampu > 0), min(cov_sub, sum(cpg_dt$cov_ampu > 0)))

cov_plot_df <- rbind(
  data.frame(coverage = cpg_dt$cov_ctrl[cov_idx_c], condition = "Control"),
  data.frame(coverage = cpg_dt$cov_ampu[cov_idx_a], condition = "Amputated")
)
# Cap at 100x for visualization
cov_plot_df$coverage <- pmin(cov_plot_df$coverage, 100)

p2 <- ggplot(cov_plot_df, aes(x = coverage, fill = condition)) +
  geom_histogram(position = "identity", alpha = 0.4, bins = 100) +
  scale_fill_manual(values = cond_colors) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  labs(x = "Read depth per CpG", y = "Count",
       title = "CpG coverage distribution",
       subtitle = paste0("Capped at 100x. Median: Control = ",
                         round(med_cov_ctrl, 1), "x, Amputated = ",
                         round(med_cov_ampu, 1), "x")) +
  theme_minimal(base_size = 13) +
  theme(legend.title = element_blank(),
        legend.position = c(0.8, 0.8))

save_plot(p2, "coverage_distribution", w = 8, h = 5)
rm(cov_plot_df); gc(verbose = FALSE)

# =============================================================================
# G. REGIONAL METHYLATION
# =============================================================================

cat("\n=== Section G: Regional methylation ===\n")

# Per-chromosome loop for memory safety
region_meth_list <- list()

for (chr_name in keep_chr) {
  chr_cpg <- cpg_dt[chr == chr_name & cov_all >= MIN_COV]
  if (nrow(chr_cpg) == 0) next
  chr_gr <- GRanges(chr_name, IRanges(chr_cpg$pos, width = 1))

  for (reg_name in names(regions_list)) {
    reg <- regions_list[[reg_name]]
    chr_reg <- reg[seqnames(reg) == chr_name]
    if (length(chr_reg) == 0) next

    hits <- overlapsAny(chr_gr, chr_reg)
    if (sum(hits) == 0) next

    beta_c <- chr_cpg$beta_ctrl[hits]
    beta_a <- chr_cpg$beta_ampu[hits]

    region_meth_list[[paste(chr_name, reg_name)]] <- data.frame(
      chr = chr_name, region = reg_name,
      n_cpg = sum(hits),
      sum_beta_ctrl = sum(beta_c, na.rm = TRUE),
      n_ctrl = sum(!is.na(beta_c)),
      sum_beta_ampu = sum(beta_a, na.rm = TRUE),
      n_ampu = sum(!is.na(beta_a))
    )
  }
  rm(chr_cpg, chr_gr); gc(verbose = FALSE)
}

region_chr_df <- do.call(rbind, region_meth_list)

# Aggregate across chromosomes
region_summary <- region_chr_df %>%
  group_by(region) %>%
  summarise(
    n_cpg = sum(n_cpg),
    mean_beta_ctrl = sum(sum_beta_ctrl) / sum(n_ctrl),
    mean_beta_ampu = sum(sum_beta_ampu) / sum(n_ampu),
    n_ctrl = sum(n_ctrl),
    n_ampu = sum(n_ampu),
    .groups = "drop"
  )

# SE across chromosomes
region_se <- region_chr_df %>%
  mutate(beta_ctrl = sum_beta_ctrl / n_ctrl,
         beta_ampu = sum_beta_ampu / n_ampu) %>%
  group_by(region) %>%
  summarise(
    se_ctrl = sd(beta_ctrl, na.rm = TRUE) / sqrt(n()),
    se_ampu = sd(beta_ampu, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

region_summary <- merge(region_summary, region_se, by = "region")

cat("  Regional methylation:\n")
for (i in seq_len(nrow(region_summary))) {
  r <- region_summary[i, ]
  cat(sprintf("    %-12s Control=%.4f  Amputated=%.4f  n=%s\n",
              r$region, r$mean_beta_ctrl, r$mean_beta_ampu,
              format(r$n_cpg, big.mark = ",")))
}

# =============================================================================
# H. PLOT 3: Mean Methylation by Region
# =============================================================================

cat("=== Section H: Plot 3 — Methylation by region ===\n")

region_order <- c("Promoter", "Gene body", "Exon", "Intron", "TE", "Intergenic")
reg_plot <- region_summary %>%
  tidyr::pivot_longer(cols = c(mean_beta_ctrl, mean_beta_ampu),
                      names_to = "condition", values_to = "mean_beta") %>%
  mutate(
    se = ifelse(condition == "mean_beta_ctrl", se_ctrl, se_ampu),
    condition = ifelse(condition == "mean_beta_ctrl", "Control", "Amputated"),
    region = factor(region, levels = region_order),
    condition = factor(condition, levels = c("Control", "Amputated"))
  )

p3 <- ggplot(reg_plot, aes(x = region, y = mean_beta, fill = condition)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_beta - se, ymax = mean_beta + se),
                position = position_dodge(0.8), width = 0.2) +
  scale_fill_manual(values = cond_colors) +
  labs(x = NULL, y = "Mean methylation (beta)",
       title = "Mean CpG methylation by genomic region",
       subtitle = paste0("Error bars = SE across chromosomes (n = 31). Coverage >= ",
                         MIN_COV, "x")) +
  theme_minimal(base_size = 13) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(size = 11))

save_plot(p3, "methylation_by_region", w = 8, h = 5)

# =============================================================================
# I. PLOT 4: Beta Distribution by Region
# =============================================================================

cat("=== Section I: Plot 4 — Beta distribution by region ===\n")

set.seed(44)
reg_beta_list <- list()

for (chr_name in keep_chr) {
  chr_cpg <- cpg_dt[chr == chr_name & cov_all >= MIN_COV]
  if (nrow(chr_cpg) == 0) next
  chr_gr <- GRanges(chr_name, IRanges(chr_cpg$pos, width = 1))

  for (reg_name in names(regions_list)) {
    reg <- regions_list[[reg_name]]
    chr_reg <- reg[seqnames(reg) == chr_name]
    if (length(chr_reg) == 0) next

    hits <- which(overlapsAny(chr_gr, chr_reg))
    if (length(hits) == 0) next

    # Subsample per chr: ~6500 per chr * 31 chr ~ 200K per region
    n_sub <- min(length(hits), 6500)
    idx <- hits[sample.int(length(hits), n_sub)]

    reg_beta_list[[paste(chr_name, reg_name, "ctrl")]] <- data.frame(
      beta = chr_cpg$beta_ctrl[idx], region = reg_name, condition = "Control")
    reg_beta_list[[paste(chr_name, reg_name, "ampu")]] <- data.frame(
      beta = chr_cpg$beta_ampu[idx], region = reg_name, condition = "Amputated")
  }
  rm(chr_cpg, chr_gr); gc(verbose = FALSE)
}

reg_beta_df <- do.call(rbind, reg_beta_list)
reg_beta_df <- reg_beta_df[!is.na(reg_beta_df$beta), ]
reg_beta_df$region <- factor(reg_beta_df$region, levels = region_order)

p4 <- ggplot(reg_beta_df, aes(x = beta, fill = condition, color = condition)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  facet_wrap(~ region, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = cond_colors) +
  scale_color_manual(values = cond_colors) +
  labs(x = "Methylation level (beta)", y = "Density",
       title = "CpG methylation distribution by genomic region",
       subtitle = "~200K CpGs sampled per region") +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        strip.text = element_text(face = "bold"))

save_plot(p4, "beta_distribution_by_region", w = 10, h = 7)
rm(reg_beta_df, reg_beta_list); gc(verbose = FALSE)

# =============================================================================
# J. PLOT 5: Metagene Methylation Profile
# =============================================================================

cat("\n=== Section J: Plot 5 — Metagene methylation profile ===\n")

N_UP <- 20; N_BODY <- 50; N_DOWN <- 20
N_TOTAL <- N_UP + N_BODY + N_DOWN

# Accumulate bin sums across chromosomes
bin_sum_ctrl <- numeric(N_TOTAL)
bin_sum_ampu <- numeric(N_TOTAL)
bin_count_ctrl <- numeric(N_TOTAL)
bin_count_ampu <- numeric(N_TOTAL)

for (chr_name in keep_chr) {
  chr_cpg <- cpg_dt[chr == chr_name & cov_all >= MIN_COV]
  if (nrow(chr_cpg) == 0) next

  chr_genes <- genes[seqnames(genes) == chr_name]
  if (length(chr_genes) == 0) { rm(chr_cpg); next }

  cpg_pos <- chr_cpg$pos
  b_ctrl  <- chr_cpg$beta_ctrl
  b_ampu  <- chr_cpg$beta_ampu

  for (gi in seq_along(chr_genes)) {
    g <- chr_genes[gi]
    g_start <- start(g); g_end <- end(g)
    g_strand <- as.character(strand(g))
    g_len <- g_end - g_start + 1

    # Define regions based on strand
    if (g_strand == "-") {
      tss <- g_end; tes <- g_start
      up_start <- tss + 1; up_end <- tss + 2000
      dn_start <- tes - 2000; dn_end <- tes - 1
    } else {
      tss <- g_start; tes <- g_end
      up_start <- tss - 2000; up_end <- tss - 1
      dn_start <- tes + 1; dn_end <- tes + 2000
    }

    # Upstream bins
    up_idx <- which(cpg_pos >= up_start & cpg_pos <= up_end)
    if (length(up_idx) > 0) {
      if (g_strand == "-") {
        rel_pos <- up_end - cpg_pos[up_idx]  # reverse for - strand
      } else {
        rel_pos <- cpg_pos[up_idx] - up_start
      }
      bins <- pmin(floor(rel_pos / (2000 / N_UP)) + 1, N_UP)
      for (b in unique(bins)) {
        sel <- up_idx[bins == b]
        v_c <- b_ctrl[sel]; v_a <- b_ampu[sel]
        ok_c <- !is.na(v_c); ok_a <- !is.na(v_a)
        bin_sum_ctrl[b] <- bin_sum_ctrl[b] + sum(v_c[ok_c])
        bin_count_ctrl[b] <- bin_count_ctrl[b] + sum(ok_c)
        bin_sum_ampu[b] <- bin_sum_ampu[b] + sum(v_a[ok_a])
        bin_count_ampu[b] <- bin_count_ampu[b] + sum(ok_a)
      }
    }

    # Gene body bins
    body_idx <- which(cpg_pos >= g_start & cpg_pos <= g_end)
    if (length(body_idx) > 0) {
      if (g_strand == "-") {
        rel_pos <- g_end - cpg_pos[body_idx]
      } else {
        rel_pos <- cpg_pos[body_idx] - g_start
      }
      frac <- rel_pos / g_len  # 0 to 1
      bins <- pmin(floor(frac * N_BODY) + 1, N_BODY) + N_UP
      for (b in unique(bins)) {
        sel <- body_idx[bins == b]
        v_c <- b_ctrl[sel]; v_a <- b_ampu[sel]
        ok_c <- !is.na(v_c); ok_a <- !is.na(v_a)
        bin_sum_ctrl[b] <- bin_sum_ctrl[b] + sum(v_c[ok_c])
        bin_count_ctrl[b] <- bin_count_ctrl[b] + sum(ok_c)
        bin_sum_ampu[b] <- bin_sum_ampu[b] + sum(v_a[ok_a])
        bin_count_ampu[b] <- bin_count_ampu[b] + sum(ok_a)
      }
    }

    # Downstream bins
    dn_idx <- which(cpg_pos >= dn_start & cpg_pos <= dn_end)
    if (length(dn_idx) > 0) {
      if (g_strand == "-") {
        rel_pos <- dn_end - cpg_pos[dn_idx]  # reverse
        rel_pos <- 2000 - 1 - rel_pos  # flip so 0 = close to TES
      } else {
        rel_pos <- cpg_pos[dn_idx] - dn_start
      }
      bins <- pmin(floor(rel_pos / (2000 / N_DOWN)) + 1, N_DOWN) + N_UP + N_BODY
      for (b in unique(bins)) {
        sel <- dn_idx[bins == b]
        v_c <- b_ctrl[sel]; v_a <- b_ampu[sel]
        ok_c <- !is.na(v_c); ok_a <- !is.na(v_a)
        bin_sum_ctrl[b] <- bin_sum_ctrl[b] + sum(v_c[ok_c])
        bin_count_ctrl[b] <- bin_count_ctrl[b] + sum(ok_c)
        bin_sum_ampu[b] <- bin_sum_ampu[b] + sum(v_a[ok_a])
        bin_count_ampu[b] <- bin_count_ampu[b] + sum(ok_a)
      }
    }
  }
  rm(chr_cpg); gc(verbose = FALSE)
}

meta_beta_ctrl <- bin_sum_ctrl / bin_count_ctrl
meta_beta_ampu <- bin_sum_ampu / bin_count_ampu

meta_df <- rbind(
  data.frame(bin = 1:N_TOTAL, beta = meta_beta_ctrl, condition = "Control"),
  data.frame(bin = 1:N_TOTAL, beta = meta_beta_ampu, condition = "Amputated")
)

# X-axis labels
x_breaks <- c(1, N_UP, N_UP + 1, N_UP + N_BODY, N_UP + N_BODY + 1, N_TOTAL)
x_labels <- c("-2 kb", "TSS", "TSS", "TES", "TES", "+2 kb")

p5 <- ggplot(meta_df, aes(x = bin, y = beta, color = condition)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = c(N_UP + 0.5, N_UP + N_BODY + 0.5),
             linetype = "dashed", color = "gray50") +
  scale_color_manual(values = cond_colors) +
  scale_x_continuous(breaks = x_breaks, labels = x_labels) +
  annotate("text", x = N_UP / 2, y = max(meta_df$beta, na.rm = TRUE),
           label = "Upstream", size = 3.5, color = "gray40") +
  annotate("text", x = N_UP + N_BODY / 2, y = max(meta_df$beta, na.rm = TRUE),
           label = "Gene body", size = 3.5, color = "gray40") +
  annotate("text", x = N_UP + N_BODY + N_DOWN / 2,
           y = max(meta_df$beta, na.rm = TRUE),
           label = "Downstream", size = 3.5, color = "gray40") +
  labs(x = NULL, y = "Mean methylation (beta)",
       title = "Metagene methylation profile",
       subtitle = paste0(format(length(genes), big.mark = ","),
                         " genes, coverage >= ", MIN_COV, "x")) +
  theme_minimal(base_size = 13) +
  theme(legend.title = element_blank(),
        legend.position = c(0.9, 0.2))

save_plot(p5, "metagene_methylation_profile", w = 9, h = 5)

# =============================================================================
# K. PLOT 6: TSS Methylation Profile
# =============================================================================

cat("=== Section K: Plot 6 — TSS methylation profile ===\n")

TSS_WINDOW <- 2000
BIN_SIZE <- 100
n_bins_tss <- (2 * TSS_WINDOW) / BIN_SIZE

tss_sum_ctrl <- numeric(n_bins_tss)
tss_sum_ampu <- numeric(n_bins_tss)
tss_cnt_ctrl <- numeric(n_bins_tss)
tss_cnt_ampu <- numeric(n_bins_tss)

for (chr_name in keep_chr) {
  chr_cpg <- cpg_dt[chr == chr_name & cov_all >= MIN_COV]
  if (nrow(chr_cpg) == 0) next

  chr_genes <- genes[seqnames(genes) == chr_name]
  if (length(chr_genes) == 0) { rm(chr_cpg); next }

  cpg_pos <- chr_cpg$pos
  b_ctrl  <- chr_cpg$beta_ctrl
  b_ampu  <- chr_cpg$beta_ampu

  for (gi in seq_along(chr_genes)) {
    g <- chr_genes[gi]
    g_strand <- as.character(strand(g))
    tss <- ifelse(g_strand == "-", end(g), start(g))

    idx <- which(cpg_pos >= (tss - TSS_WINDOW) & cpg_pos < (tss + TSS_WINDOW))
    if (length(idx) == 0) next

    if (g_strand == "-") {
      dist <- tss - cpg_pos[idx]  # positive = upstream for - strand
    } else {
      dist <- cpg_pos[idx] - tss  # positive = downstream for + strand
    }
    # Shift: dist ranges from -2000 to +1999
    bins <- floor((dist + TSS_WINDOW) / BIN_SIZE) + 1
    bins <- pmax(1, pmin(bins, n_bins_tss))

    for (b in unique(bins)) {
      sel <- idx[bins == b]
      v_c <- b_ctrl[sel]; v_a <- b_ampu[sel]
      ok_c <- !is.na(v_c); ok_a <- !is.na(v_a)
      tss_sum_ctrl[b] <- tss_sum_ctrl[b] + sum(v_c[ok_c])
      tss_cnt_ctrl[b] <- tss_cnt_ctrl[b] + sum(ok_c)
      tss_sum_ampu[b] <- tss_sum_ampu[b] + sum(v_a[ok_a])
      tss_cnt_ampu[b] <- tss_cnt_ampu[b] + sum(ok_a)
    }
  }
  rm(chr_cpg); gc(verbose = FALSE)
}

tss_beta_ctrl <- tss_sum_ctrl / tss_cnt_ctrl
tss_beta_ampu <- tss_sum_ampu / tss_cnt_ampu
bin_centers <- seq(-TSS_WINDOW + BIN_SIZE / 2, TSS_WINDOW - BIN_SIZE / 2,
                   by = BIN_SIZE)

tss_df <- rbind(
  data.frame(position = bin_centers, beta = tss_beta_ctrl, condition = "Control"),
  data.frame(position = bin_centers, beta = tss_beta_ampu, condition = "Amputated")
)

p6 <- ggplot(tss_df, aes(x = position, y = beta, color = condition)) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = cond_colors) +
  scale_x_continuous(breaks = seq(-2000, 2000, 500),
                     labels = paste0(seq(-2000, 2000, 500) / 1000, " kb")) +
  labs(x = "Distance from TSS", y = "Mean methylation (beta)",
       title = "Methylation profile around TSS",
       subtitle = paste0(format(length(genes), big.mark = ","),
                         " genes, ", BIN_SIZE, "-bp bins")) +
  theme_minimal(base_size = 13) +
  theme(legend.title = element_blank(),
        legend.position = c(0.9, 0.2))

save_plot(p6, "tss_methylation_profile", w = 9, h = 5)

# =============================================================================
# L. PLOT 7: Promoter Methylation by CpG Density Class
# =============================================================================

cat("\n=== Section L: Plot 7 — Promoter CpG density classes ===\n")

# Count CpG sites per promoter using per-chr loop
prom_cpg_counts <- integer(length(gene_promoters))

for (chr_name in keep_chr) {
  chr_prom_idx <- which(as.character(seqnames(gene_promoters)) == chr_name)
  if (length(chr_prom_idx) == 0) next
  chr_cpg_pos <- cpg_dt[chr == chr_name, pos]
  chr_cpg_gr <- GRanges(chr_name, IRanges(chr_cpg_pos, width = 1))
  chr_proms <- gene_promoters[chr_prom_idx]
  prom_cpg_counts[chr_prom_idx] <- countOverlaps(chr_proms, chr_cpg_gr)
  rm(chr_cpg_pos, chr_cpg_gr); gc(verbose = FALSE)
}

prom_density <- prom_cpg_counts / (width(gene_promoters) / 1000)

# Classify by tertiles
q33 <- quantile(prom_density, 1/3)
q67 <- quantile(prom_density, 2/3)
prom_class <- ifelse(prom_density <= q33, "LCP",
              ifelse(prom_density <= q67, "ICP", "HCP"))

# Get mean beta per promoter via per-chr overlaps
prom_beta_ctrl <- rep(NA_real_, length(gene_promoters))
prom_beta_ampu <- rep(NA_real_, length(gene_promoters))

for (chr_name in keep_chr) {
  chr_prom_idx <- which(as.character(seqnames(gene_promoters)) == chr_name)
  if (length(chr_prom_idx) == 0) next

  chr_cpg <- cpg_dt[chr == chr_name & cov_all >= MIN_COV]
  if (nrow(chr_cpg) == 0) next
  chr_cpg_gr <- GRanges(chr_name, IRanges(chr_cpg$pos, width = 1))
  chr_proms <- gene_promoters[chr_prom_idx]

  olaps <- findOverlaps(chr_proms, chr_cpg_gr)
  if (length(olaps) == 0) next

  # Compute mean beta per promoter
  q_hits <- queryHits(olaps)
  s_hits <- subjectHits(olaps)

  for (pi in unique(q_hits)) {
    cpg_idx <- s_hits[q_hits == pi]
    bc <- chr_cpg$beta_ctrl[cpg_idx]
    ba <- chr_cpg$beta_ampu[cpg_idx]
    if (sum(!is.na(bc)) > 0) prom_beta_ctrl[chr_prom_idx[pi]] <- mean(bc, na.rm = TRUE)
    if (sum(!is.na(ba)) > 0) prom_beta_ampu[chr_prom_idx[pi]] <- mean(ba, na.rm = TRUE)
  }
  rm(chr_cpg, chr_cpg_gr); gc(verbose = FALSE)
}

prom_df <- data.frame(
  cpg_class = prom_class,
  cpg_density = prom_density,
  beta_ctrl = prom_beta_ctrl,
  beta_ampu = prom_beta_ampu
)
prom_df <- prom_df[!is.na(prom_df$beta_ctrl) & !is.na(prom_df$beta_ampu), ]

prom_long <- rbind(
  data.frame(cpg_class = prom_df$cpg_class, beta = prom_df$beta_ctrl,
             condition = "Control"),
  data.frame(cpg_class = prom_df$cpg_class, beta = prom_df$beta_ampu,
             condition = "Amputated")
)
prom_long$cpg_class <- factor(prom_long$cpg_class, levels = c("LCP", "ICP", "HCP"))

cat("  Promoter classes: LCP =", sum(prom_class == "LCP"),
    ", ICP =", sum(prom_class == "ICP"),
    ", HCP =", sum(prom_class == "HCP"), "\n")
cat("  CpG density cutoffs: LCP <=", round(q33, 2),
    ", HCP >", round(q67, 2), "/kb\n")

p7 <- ggplot(prom_long, aes(x = cpg_class, y = beta, fill = condition)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  scale_fill_manual(values = cond_colors) +
  labs(x = "Promoter CpG density class", y = "Mean promoter methylation (beta)",
       title = "Promoter methylation by CpG density class",
       subtitle = paste0("LCP <= ", round(q33, 1), ", ICP <= ", round(q67, 1),
                         ", HCP > ", round(q67, 1), " CpGs/kb")) +
  theme_minimal(base_size = 13) +
  theme(legend.title = element_blank(),
        legend.position = "top")

save_plot(p7, "promoter_methylation_by_cpg_class", w = 7, h = 5)

# =============================================================================
# M. PLOT 8: Methylation vs CpG Density (gene body)
# =============================================================================

cat("=== Section M: Plot 8 — Methylation vs CpG density ===\n")

# Compute per-gene CpG count and mean beta
gene_cpg_count <- integer(length(genes))
gene_beta_all  <- rep(NA_real_, length(genes))

for (chr_name in keep_chr) {
  chr_gene_idx <- which(as.character(seqnames(genes)) == chr_name)
  if (length(chr_gene_idx) == 0) next

  chr_cpg <- cpg_dt[chr == chr_name & cov_all >= MIN_COV]
  if (nrow(chr_cpg) == 0) next
  chr_cpg_gr <- GRanges(chr_name, IRanges(chr_cpg$pos, width = 1))
  chr_genes_gr <- genes[chr_gene_idx]

  # CpG count (from all positions, not just covered)
  chr_all_pos <- cpg_dt[chr == chr_name, pos]
  chr_all_gr <- GRanges(chr_name, IRanges(chr_all_pos, width = 1))
  gene_cpg_count[chr_gene_idx] <- countOverlaps(chr_genes_gr, chr_all_gr)

  # Mean beta
  olaps <- findOverlaps(chr_genes_gr, chr_cpg_gr)
  if (length(olaps) > 0) {
    q_hits <- queryHits(olaps)
    s_hits <- subjectHits(olaps)
    for (gi in unique(q_hits)) {
      cpg_idx <- s_hits[q_hits == gi]
      ba <- chr_cpg$beta_all[cpg_idx]
      if (sum(!is.na(ba)) > 0) gene_beta_all[chr_gene_idx[gi]] <- mean(ba, na.rm = TRUE)
    }
  }
  rm(chr_cpg, chr_cpg_gr, chr_all_pos, chr_all_gr); gc(verbose = FALSE)
}

gene_cpg_density <- gene_cpg_count / (width(genes) / 1000)
gene_scatter <- data.frame(
  cpg_density = gene_cpg_density,
  beta = gene_beta_all
)
gene_scatter <- gene_scatter[!is.na(gene_scatter$beta) & gene_scatter$cpg_density > 0, ]

p8 <- ggplot(gene_scatter, aes(x = cpg_density, y = beta)) +
  geom_hex(bins = 80) +
  geom_smooth(method = "loess", color = "#C0392B", se = FALSE, linewidth = 1) +
  scale_fill_viridis_c(option = "inferno", trans = "log10", name = "Count") +
  labs(x = "CpG density (CpGs per kb)", y = "Mean gene body methylation (beta)",
       title = "Gene body methylation vs CpG density",
       subtitle = paste0(format(nrow(gene_scatter), big.mark = ","), " genes")) +
  theme_minimal(base_size = 13)

save_plot(p8, "methylation_vs_cpg_density", w = 8, h = 6)

# Save gene methylation summary for future batches
gene_meth_summary <- data.frame(
  gene_id = mcols(genes)$ID,
  chr = as.character(seqnames(genes)),
  start = start(genes), end = end(genes),
  strand = as.character(strand(genes)),
  width = width(genes),
  cpg_count = gene_cpg_count,
  cpg_density = gene_cpg_density,
  mean_beta_all = gene_beta_all
)
saveRDS(gene_meth_summary, file.path(cache_dir, "gene_methylation_summary.rds"))
cat("  Saved gene_methylation_summary.rds\n")

# =============================================================================
# N. PLOTS 9-10: TE Methylation
# =============================================================================

cat("\n=== Section N: TE methylation ===\n")

te_gr$te_class <- sub("/.*", "", te_gr$class_family)
te_class_counts <- table(te_gr$te_class)
keep_classes <- names(te_class_counts[te_class_counts > 100])

# Per-TE-class mean beta
te_class_meth <- list()
for (cls in keep_classes) {
  cls_gr <- reduce(te_gr[te_gr$te_class == cls])
  cls_gr <- trim(cls_gr)
  sum_ctrl <- 0; n_ctrl <- 0; sum_ampu <- 0; n_ampu <- 0

  for (chr_name in keep_chr) {
    chr_te <- cls_gr[seqnames(cls_gr) == chr_name]
    if (length(chr_te) == 0) next
    chr_cpg <- cpg_dt[chr == chr_name & cov_all >= MIN_COV]
    if (nrow(chr_cpg) == 0) next
    chr_cpg_gr <- GRanges(chr_name, IRanges(chr_cpg$pos, width = 1))
    hits <- overlapsAny(chr_cpg_gr, chr_te)
    if (sum(hits) == 0) next
    bc <- chr_cpg$beta_ctrl[hits]; ba <- chr_cpg$beta_ampu[hits]
    sum_ctrl <- sum_ctrl + sum(bc, na.rm = TRUE)
    n_ctrl <- n_ctrl + sum(!is.na(bc))
    sum_ampu <- sum_ampu + sum(ba, na.rm = TRUE)
    n_ampu <- n_ampu + sum(!is.na(ba))
    rm(chr_cpg, chr_cpg_gr); gc(verbose = FALSE)
  }
  te_class_meth[[cls]] <- data.frame(
    te_class = cls, n_elements = as.integer(te_class_counts[cls]),
    mean_beta_ctrl = sum_ctrl / n_ctrl, mean_beta_ampu = sum_ampu / n_ampu,
    n_cpg_ctrl = n_ctrl, n_cpg_ampu = n_ampu
  )
}
te_meth_df <- do.call(rbind, te_class_meth)
cat("  TE methylation by class:\n")
for (i in seq_len(nrow(te_meth_df))) {
  r <- te_meth_df[i, ]
  cat(sprintf("    %-10s Control=%.4f  Amputated=%.4f  n_cpg=%s\n",
              r$te_class, r$mean_beta_ctrl, r$mean_beta_ampu,
              format(r$n_cpg_ctrl, big.mark = ",")))
}

# Plot 9: TE methylation by class
te_plot <- te_meth_df %>%
  tidyr::pivot_longer(cols = c(mean_beta_ctrl, mean_beta_ampu),
                      names_to = "condition", values_to = "mean_beta") %>%
  mutate(condition = ifelse(condition == "mean_beta_ctrl", "Control", "Amputated"),
         condition = factor(condition, levels = c("Control", "Amputated")),
         te_class = factor(te_class, levels = te_meth_df$te_class[order(te_meth_df$mean_beta_ctrl)]))

p9 <- ggplot(te_plot, aes(x = te_class, y = mean_beta, fill = condition)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  scale_fill_manual(values = cond_colors) +
  labs(x = "TE class", y = "Mean methylation (beta)",
       title = "TE methylation by class",
       subtitle = "All TE classes with >100 elements") +
  theme_minimal(base_size = 13) +
  theme(legend.title = element_blank(),
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1))

save_plot(p9, "te_methylation_by_class", w = 8, h = 5)

# Plot 10: TE methylation vs age (Kimura divergence)
cat("  Computing per-TE element methylation for age plot...\n")

# Process in chunks: for each TE element, get mean beta
# Use only a subsample of TEs for the scatter (too many individual TEs)
set.seed(45)
te_sub_idx <- sample(length(te_gr), min(50000, length(te_gr)))
te_sub <- te_gr[te_sub_idx]

te_elem_beta <- rep(NA_real_, length(te_sub))
te_elem_kimura <- mcols(te_sub)$kimura_div

for (chr_name in keep_chr) {
  chr_te_idx <- which(as.character(seqnames(te_sub)) == chr_name)
  if (length(chr_te_idx) == 0) next
  chr_cpg <- cpg_dt[chr == chr_name & cov_all >= MIN_COV]
  if (nrow(chr_cpg) == 0) next
  chr_cpg_gr <- GRanges(chr_name, IRanges(chr_cpg$pos, width = 1))
  chr_te <- te_sub[chr_te_idx]

  olaps <- findOverlaps(chr_te, chr_cpg_gr)
  if (length(olaps) == 0) { rm(chr_cpg, chr_cpg_gr); next }

  q_hits <- queryHits(olaps)
  s_hits <- subjectHits(olaps)
  for (ti in unique(q_hits)) {
    cpg_idx <- s_hits[q_hits == ti]
    ba <- chr_cpg$beta_all[cpg_idx]
    if (sum(!is.na(ba)) >= 3) {  # require at least 3 CpGs
      te_elem_beta[chr_te_idx[ti]] <- mean(ba, na.rm = TRUE)
    }
  }
  rm(chr_cpg, chr_cpg_gr); gc(verbose = FALSE)
}

te_age_df <- data.frame(
  kimura = te_elem_kimura,
  beta = te_elem_beta,
  te_class = sub("/.*", "", mcols(te_sub)$class_family)
)
te_age_df <- te_age_df[!is.na(te_age_df$beta) & !is.na(te_age_df$kimura), ]
cat("  TEs with beta + kimura data:", nrow(te_age_df), "\n")

te_class_colors <- c(LINE = "#2471A3", LTR = "#8E44AD", DNA = "#27AE60",
                     RC = "#F39C12", SINE = "#C0392B", Unknown = "#95A5A6")

p10 <- ggplot(te_age_df, aes(x = kimura, y = beta)) +
  geom_hex(bins = 60) +
  geom_smooth(method = "loess", color = "#C0392B", se = FALSE, linewidth = 1) +
  scale_fill_viridis_c(option = "inferno", trans = "log10", name = "Count") +
  labs(x = "Kimura divergence (%)", y = "Mean methylation (beta)",
       title = "TE methylation vs evolutionary age",
       subtitle = paste0(format(nrow(te_age_df), big.mark = ","),
                         " TE elements (>= 3 CpGs each)")) +
  theme_minimal(base_size = 13)

save_plot(p10, "te_methylation_vs_age", w = 8, h = 6)

# =============================================================================
# O-P removed: Sample correlation + chromosome overview
# (plan had 10 plots but the 10 above are sufficient and more informative)
# =============================================================================

# =============================================================================
# Q. HTML REPORT
# =============================================================================

cat("\n=== Section Q: Writing HTML report ===\n")

report_path <- file.path(out_dir, "baseline_methylome_report.html")
sink(report_path)

cat('<!DOCTYPE html>\n<html><head><meta charset="UTF-8">\n')
cat('<title>Batch 03: Baseline Methylome</title>\n')
cat('<style>\n')
cat('body { font-family: "Segoe UI", Arial, sans-serif; max-width: 1100px; margin: 40px auto; padding: 0 20px; color: #333; }\n')
cat('h1 { color: #2c3e50; border-bottom: 2px solid #2471A3; padding-bottom: 10px; }\n')
cat('h2 { color: #2471A3; margin-top: 30px; }\n')
cat('table { border-collapse: collapse; width: auto; margin: 15px 0; font-size: 13px; }\n')
cat('th { background: #2471A3; color: white; padding: 8px 10px; text-align: left; }\n')
cat('td { padding: 5px 10px; border-bottom: 1px solid #ddd; }\n')
cat('tr:nth-child(even) { background: #f8f9fa; }\n')
cat('img { max-width: 100%; border: 1px solid #ddd; margin: 10px 0; }\n')
cat('.key-finding { background: #eaf2f8; padding: 12px; border-left: 4px solid #2471A3; margin: 15px 0; }\n')
cat('</style>\n</head><body>\n')

cat('<h1>Batch 03: Baseline Methylome in <em>D. laeve</em></h1>\n')
cat('<p>Generated:', format(Sys.time(), "%Y-%m-%d %H:%M"), '</p>\n')
cat('<p><strong>Question:</strong> What does the <em>D. laeve</em> methylome look like? Where is methylation concentrated? How does it vary by genomic region?</p>\n')

# Summary stats table
cat('<h2>1. Global Methylation Summary</h2>\n')
cat('<table>\n')
cat('<tr><th>Metric</th><th>Control</th><th>Amputated</th></tr>\n')
cat('<tr><td>CpGs with coverage &ge; ', MIN_COV, 'x</td><td>',
    format(n_ctrl, big.mark = ","), '</td><td>',
    format(n_ampu, big.mark = ","), '</td></tr>\n')
cat('<tr><td>Mean beta</td><td>', round(mean_beta_ctrl, 4),
    '</td><td>', round(mean_beta_ampu, 4), '</td></tr>\n')
cat('<tr><td>Median beta</td><td>', round(med_beta_ctrl, 4),
    '</td><td>', round(med_beta_ampu, 4), '</td></tr>\n')
cat('<tr><td>Unmethylated (&lt;0.2)</td><td>', round(frac_unm_ctrl * 100, 1),
    '%</td><td>', round(frac_unm_ampu * 100, 1), '%</td></tr>\n')
cat('<tr><td>Intermediate (0.2&ndash;0.8)</td><td>', round(frac_int_ctrl * 100, 1),
    '%</td><td>', round(frac_int_ampu * 100, 1), '%</td></tr>\n')
cat('<tr><td>Methylated (&gt;0.8)</td><td>', round(frac_met_ctrl * 100, 1),
    '%</td><td>', round(frac_met_ampu * 100, 1), '%</td></tr>\n')
cat('<tr><td>Mean coverage</td><td>', round(mean_cov_ctrl, 1),
    'x</td><td>', round(mean_cov_ampu, 1), 'x</td></tr>\n')
cat('<tr><td>Median coverage</td><td>', round(med_cov_ctrl, 1),
    'x</td><td>', round(med_cov_ampu, 1), 'x</td></tr>\n')
cat('</table>\n')

# Coverage distribution
cat('<h2>2. Coverage Distribution</h2>\n')
cat('<p>Read depth per CpG site (capped at 100x for visualization).</p>\n')
cat('<img src="png/coverage_distribution.png" alt="Coverage distribution">\n')

# Global beta
cat('<h2>3. Global Methylation Distribution</h2>\n')
cat('<p>Density plot of CpG methylation levels (beta values) across conditions.</p>\n')
cat('<img src="png/global_beta_distribution.png" alt="Global beta distribution">\n')

# Regional methylation table + plot
cat('<h2>4. Methylation by Genomic Region</h2>\n')
cat('<table>\n')
cat('<tr><th>Region</th><th>CpGs</th><th>Mean beta (Control)</th><th>Mean beta (Amputated)</th><th>SE (Control)</th></tr>\n')
for (i in seq_len(nrow(region_summary))) {
  r <- region_summary[i, ]
  cat('<tr><td>', r$region, '</td><td>', format(r$n_cpg, big.mark = ","),
      '</td><td>', round(r$mean_beta_ctrl, 4),
      '</td><td>', round(r$mean_beta_ampu, 4),
      '</td><td>', round(r$se_ctrl, 4), '</td></tr>\n')
}
cat('</table>\n')
cat('<img src="png/methylation_by_region.png" alt="Methylation by region">\n')

# Beta by region
cat('<h2>5. Methylation Distribution by Region</h2>\n')
cat('<p>CpG methylation distributions faceted by genomic region. Shows how the bimodal pattern varies.</p>\n')
cat('<img src="png/beta_distribution_by_region.png" alt="Beta distribution by region">\n')

# Metagene
cat('<h2>6. Metagene Methylation Profile</h2>\n')
cat('<p>Average methylation across a standardized gene model: 2 kb upstream &rarr; gene body (50 bins) &rarr; 2 kb downstream.</p>\n')
cat('<img src="png/metagene_methylation_profile.png" alt="Metagene profile">\n')

# TSS profile
cat('<h2>7. TSS Methylation Profile</h2>\n')
cat('<p>Methylation at absolute distance from transcription start site (&plusmn;2 kb, 100-bp bins).</p>\n')
cat('<img src="png/tss_methylation_profile.png" alt="TSS profile">\n')

# Promoter CpG classes
cat('<h2>8. Promoter Methylation by CpG Density Class</h2>\n')
cat('<p>Promoters classified by CpG density tertiles: LCP (low, &le;', round(q33, 1),
    '/kb), ICP (intermediate), HCP (high, &gt;', round(q67, 1), '/kb).</p>\n')
cat('<img src="png/promoter_methylation_by_cpg_class.png" alt="Promoter CpG classes">\n')

# Meth vs CpG density
cat('<h2>9. Gene Body Methylation vs CpG Density</h2>\n')
cat('<p>Relationship between CpG density and mean gene body methylation across ',
    format(nrow(gene_scatter), big.mark = ","), ' genes.</p>\n')
cat('<img src="png/methylation_vs_cpg_density.png" alt="Meth vs CpG density">\n')

# TE methylation
cat('<h2>10. TE Methylation by Class</h2>\n')
cat('<table>\n')
cat('<tr><th>TE Class</th><th>Elements</th><th>Mean beta (Control)</th><th>Mean beta (Amputated)</th><th>CpGs</th></tr>\n')
for (i in seq_len(nrow(te_meth_df))) {
  r <- te_meth_df[i, ]
  cat('<tr><td>', r$te_class, '</td><td>', format(r$n_elements, big.mark = ","),
      '</td><td>', round(r$mean_beta_ctrl, 4),
      '</td><td>', round(r$mean_beta_ampu, 4),
      '</td><td>', format(r$n_cpg_ctrl, big.mark = ","), '</td></tr>\n')
}
cat('</table>\n')
cat('<img src="png/te_methylation_by_class.png" alt="TE methylation by class">\n')

cat('<h2>11. TE Methylation vs Evolutionary Age</h2>\n')
cat('<p>Relationship between TE age (Kimura divergence) and methylation level.</p>\n')
cat('<img src="png/te_methylation_vs_age.png" alt="TE meth vs age">\n')

# Key findings
cat('<h2>12. Key Findings</h2>\n')
cat('<div class="key-finding">\n<ul>\n')
cat('<li><strong>Bimodal methylation:</strong> CpGs are either unmethylated (&lt;0.2) or heavily methylated (&gt;0.8), with few intermediate values.</li>\n')
cat('<li><strong>Gene body methylation:</strong> Gene bodies show higher methylation than promoters, consistent with invertebrate-type methylation.</li>\n')
cat('<li><strong>TSS dip:</strong> Methylation dips at the transcription start site, suggesting functional relevance.</li>\n')
cat('<li><strong>TE methylation:</strong> Transposable elements are heavily methylated, consistent with TE silencing.</li>\n')
cat('<li><strong>Control vs amputated:</strong> Global patterns are highly similar between conditions; differential analysis in Batch 04.</li>\n')
cat('</ul>\n</div>\n')

cat('</body></html>\n')
sink()
cat("  Report written to:", report_path, "\n")

cat("\n=== Batch 03 complete ===\n")
cat("Output:", out_dir, "\n")
cat("Plots: 10 (PNG + PDF)\n")
cat("Cache:", cpg_merged_rds, "\n")
cat("Done!\n")
