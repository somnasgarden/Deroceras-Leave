#!/usr/bin/env Rscript
# =============================================================================
# Batch 11: Motif x Per-read Entropy x Expression
#
# Closes the loop:
#   TF motif -> ordered methylation (low NME) -> DMP -> DE gene -> phenotype
#
# Distinct from:
#   batch10  = global/regional Level-1 (marginal) entropy
#   batch10b = per-read Level-2 NME from BAMs (the input here)
#
# Logic: a 4-CpG window with low NME means every read at that locus has the
# same methylation pattern. That uniformity has to be enforced by something.
# If windows that overlap a JASPAR TFBS have systematically lower NME than
# beta-matched non-motif windows, the bound TF is the most plausible enforcer
# (directed mechanism). Then if those low-NME-in-motif windows belong to genes
# that are differentially expressed during regeneration, we have a mechanistic
# shortlist.
#
# Inputs:
#   batch10/data/perread_nme_windows.tsv  (4-CpG windows with NME per sample)
#   batch1.5/data/motif_hits_extended.tsv.gz  (genome-wide JASPAR hits, all regions)
#   batch06/data/dmps_annotated.tsv  (for DMP overlap)
#   CACHE$transcriptome  (DESeq2 control vs amputated tail, apeglm shrunk)
#
# Outputs:
#   data/window_motif_assignments.tsv     (every NME window with motif tag)
#   data/nme_in_vs_out_motif_test.tsv     (Test A: beta-matched in vs out)
#   data/delta_nme_by_methyl_sens.tsv     (Test B: delta NME by Yin class)
#   data/read_coverage_by_motif_class.tsv (read depth diagnostics per TF class)
#   data/directed_mechanism_loci.tsv      (Test C: shortlist with expression)
#   figures/fig11a_nme_in_vs_out_motif    (Test A violin)
#   figures/fig11b_delta_nme_by_methyl_sens
#   figures/fig11c_directed_loci_volcano
#   figures/fig11d_reads_per_window_by_motif
#   figures/fig11e_reads_per_motif_class
#   figures/fig11f_reads_vs_nme
# =============================================================================

source("methylation_pipeline/_config.R")
t0 <- proc.time()

library(data.table)
library(GenomicRanges)
library(ggplot2)

BATCH_DIR <- file.path(PIPE_DIR, "batch11")
dir.create(file.path(BATCH_DIR, "data"),    showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BATCH_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)

# Clean stale outputs
unlink(list.files(file.path(BATCH_DIR, "data"),    full.names = TRUE))
unlink(list.files(file.path(BATCH_DIR, "figures"), full.names = TRUE))

cat("=== Batch 11: Motif x NME x Expression ===\n\n")

MIN_READS         <- 8L      # min reads per window per sample (stable NME)
LOW_NME_THRESHOLD <- 0.30    # "ordered" cutoff for shortlist
BETA_MATCH_TOL    <- 0.05    # +/- beta tolerance for matched comparison

# =============================================================================
# Yin 2017 methylation sensitivity (inlined from batch09 — keep batches independent)
# =============================================================================
methyl_sens_class <- function(tf_class) {
  cls <- tolower(as.character(tf_class))
  fcase(
    grepl("homeodomain", cls),                              "MethylPlus",
    grepl("helix-loop-helix|bhlh", cls),                    "MethylMinus",
    grepl("leucine zipper|bzip", cls),                      "MethylMinus",
    grepl("tryptophan cluster|ets",  cls),                  "MethylMinus",
    grepl("nuclear receptor", cls),                         "LittleEffect",
    grepl("forkhead|fox",     cls),                         "MethylMinus",
    grepl("rel homology|nf-?kappab",  cls),                 "LittleEffect",
    grepl("p53",              cls),                         "MethylMinus",
    grepl("sox|hmg",          cls),                         "LittleEffect",
    grepl("c2h2 zinc",        cls),                         "Mixed",
    default = "Unknown"
  )
}

yin_overrides <- data.table(
  symbol = c("CEBPB","CEBPA","CEBPD","KLF4","KLF5","KLF6","SP1","SP2","SP3",
             "EGR1","EGR2","ZBTB7A","ZBTB33",
             "MYC","MAX","USF1","USF2","MITF","TFE3","ARNT","AHR","HIF1A",
             "CREB1","ATF1","ATF3","ATF4","JUN","JUND","FOS","FOSL1","FOSL2","BATF",
             "NRF1","NFE2L1","NFE2L2","BACH1","BACH2","MAFB","MAFK","MAFG",
             "ELK1","ELK4","ELF1","ELF3","ELF5","ETS1","ETS2","ETV1","ETV4","ETV5",
             "GABPA","SPI1","SPDEF","FLI1","ERG","FEV","ERF",
             "FOXA1","FOXA2","FOXO1","FOXO3","FOXP1","FOXP3","FOXM1",
             "TP53","TP63","TP73",
             "HOXA1","HOXA2","HOXB13","HOXC9","HOXD8","NKX2-5","NKX6-1","OTX1","OTX2",
             "PAX2","PAX6","PAX7","CDX1","CDX2","DLX1","DLX2","LHX2","LHX3",
             "MEIS1","PBX1","PITX1","PITX2","PRRX1","PRRX2","SIX1","SIX2","TGIF1",
             "POU2F1","POU3F2","POU5F1","NANOG","CRX","ISL1","MSX1","MSX2","EMX1","EMX2",
             "RUNX1","RUNX2","RUNX3"),
  yin_call = c(rep("MethylMinus", 13),
               rep("MethylMinus", 9),
               rep("MethylMinus", 10),
               rep("MethylMinus", 8),
               rep("MethylMinus", 16),
               rep("MethylMinus", 7),
               rep("MethylMinus", 3),
               rep("MethylPlus", 36),
               rep("LittleEffect", 3))
)

assign_methyl_sens <- function(motif_name, tf_class) {
  base <- methyl_sens_class(tf_class)
  sym  <- toupper(sub("\\(.*$", "", as.character(motif_name)))
  sym  <- sub("::.*$", "", sym)
  ov   <- yin_overrides[match(sym, yin_overrides$symbol), yin_call]
  ifelse(!is.na(ov), ov, base)
}

# =============================================================================
# [1/7] Load NME windows from batch10b
# =============================================================================
cat("[1/7] Loading per-read NME windows (batch10b)...\n")
nme_path <- file.path(PIPE_DIR, "batch10/data/perread_nme_windows.tsv")
if (!file.exists(nme_path)) {
  stop("perread_nme_windows.tsv not found. Run batch10b first.")
}
win <- fread(nme_path)
cat(sprintf("  windows loaded: %s\n", format(nrow(win), big.mark = ",")))

# Coverage filter: require MIN_READS in all 4 samples
win[, total_reads := reads_C1 + reads_C2 + reads_A1 + reads_A2]
win_pass <- win[reads_C1 >= MIN_READS & reads_C2 >= MIN_READS &
                reads_A1 >= MIN_READS & reads_A2 >= MIN_READS]
cat(sprintf("  windows passing min %d reads/sample: %s (%.1f%%)\n",
            MIN_READS, format(nrow(win_pass), big.mark = ","),
            100 * nrow(win_pass) / nrow(win)))

# =============================================================================
# [2/7] Load motif hits and tag windows
# =============================================================================
cat("\n[2/7] Loading genome-wide motif hits (batch1.5)...\n")
motif_path <- file.path(PIPE_DIR, "batch1.5/data/motif_hits_extended.tsv.gz")
mh <- fread(motif_path)
mh <- mh[chr %in% keep_chr]
cat(sprintf("  motif hits (chr1-31, all regions): %s\n",
            format(nrow(mh), big.mark = ",")))
cat("  region distribution:\n"); print(mh[, .N, by = region_type])

# Build GRanges
win_gr <- GRanges(seqnames = win_pass$chr,
                  ranges = IRanges(start = win_pass$win_start,
                                   end   = win_pass$win_end))
mh_gr <- GRanges(seqnames = mh$chr,
                 ranges = IRanges(start = mh$start, end = mh$end))

# All overlaps (one window may overlap many motifs)
ov <- findOverlaps(win_gr, mh_gr)
cat(sprintf("  window-motif overlaps: %s\n",
            format(length(ov), big.mark = ",")))

# Score per overlap (handle empty score gracefully)
mh[, score_num := suppressWarnings(as.numeric(score))]
mh[is.na(score_num), score_num := 0]

ov_dt <- data.table(
  win_idx   = queryHits(ov),
  motif_idx = subjectHits(ov),
  motif_id  = mh$motif_id[subjectHits(ov)],
  motif_name= mh$motif_name[subjectHits(ov)],
  tf_class  = mh$tf_class[subjectHits(ov)],
  gene_id   = mh$gene_id[subjectHits(ov)],
  region_type = mh$region_type[subjectHits(ov)],
  score     = mh$score_num[subjectHits(ov)]
)

# Per-window: count motifs, pick top by score, attach Yin class
ov_dt[, methyl_sens := assign_methyl_sens(motif_name, tf_class)]
setorder(ov_dt, win_idx, -score)
top_per_win <- ov_dt[, .SD[1], by = win_idx]
counts_per_win <- ov_dt[, .(n_motifs = .N,
                            n_distinct_motifs = uniqueN(motif_id),
                            n_distinct_genes  = uniqueN(gene_id)),
                        by = win_idx]

win_pass[, win_idx := .I]
win_pass[, in_motif := FALSE]
win_pass[top_per_win$win_idx, in_motif := TRUE]
win_pass[top_per_win$win_idx, `:=`(
  motif_id_top    = top_per_win$motif_id,
  motif_name_top  = top_per_win$motif_name,
  tf_class_top    = top_per_win$tf_class,
  gene_id_top     = top_per_win$gene_id,
  region_motif    = top_per_win$region_type,
  motif_score_top = top_per_win$score,
  methyl_sens_top = top_per_win$methyl_sens
)]
win_pass[counts_per_win$win_idx, `:=`(
  n_motifs          = counts_per_win$n_motifs,
  n_distinct_motifs = counts_per_win$n_distinct_motifs,
  n_distinct_genes  = counts_per_win$n_distinct_genes
)]
win_pass[is.na(n_motifs), n_motifs := 0L]

cat(sprintf("  windows in motif: %s (%.1f%%)\n",
            format(sum(win_pass$in_motif), big.mark = ","),
            100 * mean(win_pass$in_motif)))

# Attach gene_name from EviAnn annotation
if (file.exists(OG$annot)) {
  annot_dt <- fread(OG$annot, header = FALSE,
                    col.names = c("gene_id", "gene_name", "description"))
  annot_dt[gene_name == "", gene_name := gene_id]
  win_pass[annot_dt, gene_name_top := i.gene_name, on = c(gene_id_top = "gene_id")]
}

save_data(win_pass, BATCH_DIR, "window_motif_assignments")

# =============================================================================
# [3/7] TEST A: beta-matched in-motif vs out-of-motif baseline NME
# =============================================================================
cat("\n[3/7] TEST A: beta-matched in-motif vs out-of-motif NME...\n")

in_w  <- win_pass[in_motif == TRUE  & !is.na(ctrl_nme) & !is.na(win_beta)]
out_w <- win_pass[in_motif == FALSE & !is.na(ctrl_nme) & !is.na(win_beta)]
cat(sprintf("  in-motif: %s | out-of-motif: %s\n",
            format(nrow(in_w), big.mark = ","),
            format(nrow(out_w), big.mark = ",")))

# Beta-bin matching: bin both pools into 0.02 beta bins, sample matched pairs
set.seed(42)
in_w[,  beta_bin := round(win_beta / 0.02) * 0.02]
out_w[, beta_bin := round(win_beta / 0.02) * 0.02]

matched <- in_w[, {
  cand <- out_w[abs(beta_bin - .BY[[1]]) <= BETA_MATCH_TOL]
  if (nrow(cand) >= .N) {
    .(in_nme  = ctrl_nme,
      out_nme = cand[sample(.N, .N), ctrl_nme])
  } else if (nrow(cand) > 0) {
    n <- min(.N, nrow(cand))
    .(in_nme  = ctrl_nme[seq_len(n)],
      out_nme = cand[sample(.N, n), ctrl_nme])
  } else {
    .(in_nme = numeric(0), out_nme = numeric(0))
  }
}, by = beta_bin]

cat(sprintf("  matched pairs: %s\n", format(nrow(matched), big.mark = ",")))
if (nrow(matched) > 100) {
  wt_A <- wilcox.test(matched$in_nme, matched$out_nme, paired = TRUE)
  d_A  <- (mean(matched$in_nme) - mean(matched$out_nme)) /
          sd(matched$in_nme - matched$out_nme)
  cat(sprintf("  mean ctrl NME in-motif:  %.4f\n", mean(matched$in_nme)))
  cat(sprintf("  mean ctrl NME out-motif: %.4f\n", mean(matched$out_nme)))
  cat(sprintf("  Wilcoxon paired p = %s | Cohen's d = %.4f\n",
              format(wt_A$p.value, digits = 3), d_A))

  save_data(data.table(
    n_pairs       = nrow(matched),
    in_nme_mean   = mean(matched$in_nme),
    out_nme_mean  = mean(matched$out_nme),
    delta         = mean(matched$in_nme) - mean(matched$out_nme),
    wilcoxon_p    = wt_A$p.value,
    cohens_d      = d_A
  ), BATCH_DIR, "nme_in_vs_out_motif_test")
} else {
  wt_A <- list(p.value = NA); d_A <- NA
  cat("  insufficient matched pairs.\n")
}

# =============================================================================
# [4/7] TEST B: delta NME at motif sites stratified by Yin class
# =============================================================================
cat("\n[4/7] TEST B: delta NME at motif sites by Yin methyl-sensitivity...\n")

motif_w <- win_pass[in_motif == TRUE & !is.na(delta_nme)]
test_B <- motif_w[, {
  if (.N >= 30) {
    wt <- wilcox.test(delta_nme, mu = 0)
    .(n = .N, mean_delta = mean(delta_nme),
      median_delta = median(delta_nme), p = wt$p.value)
  } else {
    .(n = .N, mean_delta = mean(delta_nme),
      median_delta = median(delta_nme), p = NA_real_)
  }
}, by = methyl_sens_top]
test_B[, padj := p.adjust(p, "BH")]
print(test_B)
save_data(test_B, BATCH_DIR, "delta_nme_by_methyl_sens")

# =============================================================================
# [5/7] Read coverage diagnostics per motif class
# =============================================================================
cat("\n[5/7] Read coverage diagnostics...\n")
cov_diag <- motif_w[, .(
  n_windows    = .N,
  median_reads = as.numeric(median(total_reads)),
  mean_reads   = mean(total_reads),
  q25_reads    = as.numeric(quantile(total_reads, 0.25)),
  q75_reads    = as.numeric(quantile(total_reads, 0.75)),
  median_ctrl_nme = median(ctrl_nme, na.rm = TRUE)
), by = tf_class_top]
setorder(cov_diag, -n_windows)
print(head(cov_diag, 20))
save_data(cov_diag, BATCH_DIR, "read_coverage_by_motif_class")

# =============================================================================
# [6/7] TEST C: directed mechanism loci shortlist
# =============================================================================
cat("\n[6/7] TEST C: directed mechanism shortlist (low NME + DMP + DE)...\n")

shortlist <- win_pass[in_motif == TRUE &
                      ctrl_nme < LOW_NME_THRESHOLD &
                      has_dmp == TRUE &
                      !is.na(gene_id_top)]
cat(sprintf("  candidates (in_motif & ctrl_nme<%.2f & has_dmp): %s\n",
            LOW_NME_THRESHOLD, format(nrow(shortlist), big.mark = ",")))

# Join expression from CACHE$transcriptome
expr_joined <- FALSE
if (!is.null(CACHE$transcriptome) && file.exists(CACHE$transcriptome)) {
  tr <- tryCatch(readRDS(CACHE$transcriptome), error = function(e) NULL)
  if (!is.null(tr) && !is.null(tr$res_tail)) {
    res_tail <- as.data.table(tr$res_tail, keep.rownames = "gene_id")
    setnames(res_tail,
             old = intersect(c("log2FoldChange","padj"), names(res_tail)),
             new = intersect(c("log2FoldChange","padj"), names(res_tail)))
    shortlist <- merge(shortlist, res_tail[, .(gene_id, log2FoldChange, padj)],
                       by.x = "gene_id_top", by.y = "gene_id", all.x = TRUE)
    expr_joined <- TRUE
    cat(sprintf("  with expression data: %s\n",
                format(sum(!is.na(shortlist$padj)), big.mark = ",")))
  }
}
if (!expr_joined) {
  shortlist[, `:=`(log2FoldChange = NA_real_, padj = NA_real_)]
  cat("  WARNING: CACHE$transcriptome unavailable, expression columns NA\n")
}

# Rank by |delta NME| * -log10(padj) when expression present
shortlist[, rank_score := abs(delta_nme) * pmax(-log10(pmax(padj, 1e-300)), 0, na.rm = TRUE)]
setorder(shortlist, -rank_score, na.last = TRUE)

shortlist_out <- shortlist[, .(
  chr, win_start, win_end,
  gene_id = gene_id_top, gene_name = gene_name_top,
  motif_name = motif_name_top, tf_class = tf_class_top,
  methyl_sens = methyl_sens_top,
  n_motifs, n_distinct_motifs,
  ctrl_nme, ampu_nme, delta_nme,
  win_beta, total_reads,
  log2FoldChange, padj, rank_score
)]
save_data(shortlist_out, BATCH_DIR, "directed_mechanism_loci")
cat(sprintf("  shortlist saved: %s loci\n",
            format(nrow(shortlist_out), big.mark = ",")))

# =============================================================================
# [7/7] FIGURES
# =============================================================================
cat("\n[7/7] Generating figures...\n")

# Fig 11a — beta-matched in vs out NME
if (nrow(matched) > 100) {
  long_A <- rbind(
    data.table(group = "In motif",     nme = matched$in_nme),
    data.table(group = "Out of motif", nme = matched$out_nme)
  )
  p11a <- ggplot(long_A, aes(x = group, y = nme, fill = group)) +
    geom_violin(alpha = 0.55, linewidth = 0.4) +
    geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.9) +
    scale_fill_manual(values = c("In motif" = "#C0392B",
                                  "Out of motif" = "#2471A3"),
                      guide = "none") +
    labs(x = NULL, y = "Baseline NME (control)",
         title = "Per-read methylation entropy: TFBS vs non-TFBS windows",
         subtitle = sprintf("beta-matched 4-CpG windows (~50-500 bp around 6-20 bp motif core) | n=%s pairs | Wilcoxon p=%s | d=%.3f",
                            format(nrow(matched), big.mark = ","),
                            format(wt_A$p.value, digits = 3), d_A)) +
    theme_minimal(base_size = 12)
  save_fig(p11a, BATCH_DIR, "fig11a_nme_in_vs_out_motif", w = 8, h = 6)
}

# Fig 11b — delta NME by Yin methyl class
p11b <- ggplot(motif_w, aes(x = methyl_sens_top, y = delta_nme, fill = methyl_sens_top)) +
  geom_violin(alpha = 0.55, linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = "Yin 2017 methylation sensitivity (top motif per window)",
       y = "delta NME (Amputated - Control)",
       title = "Entropy change at TFBS during regeneration, by TF methyl-sensitivity",
       subtitle = "Negative = TFBS becomes more ordered after amputation (directed)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))
save_fig(p11b, BATCH_DIR, "fig11b_delta_nme_by_methyl_sens", w = 9, h = 6)

# Fig 11c — directed loci volcano
if (expr_joined && nrow(shortlist_out[!is.na(padj)]) > 5) {
  vol_dt <- shortlist_out[!is.na(padj) & !is.na(log2FoldChange)]
  vol_dt[, neg_ctrl_nme := -ctrl_nme]
  top_lab <- vol_dt[order(-rank_score)][1:min(20, .N)]
  p11c <- ggplot(vol_dt, aes(x = log2FoldChange, y = neg_ctrl_nme, color = delta_nme)) +
    geom_point(alpha = 0.7, size = 1.6) +
    scale_color_gradient2(low = "#1B7837", mid = "grey80", high = "#762A83",
                          midpoint = 0, name = "delta NME") +
    geom_text(data = top_lab, aes(label = gene_name),
              size = 2.8, vjust = -0.6, color = "black",
              check_overlap = TRUE) +
    labs(x = "log2 fold change (amputated vs control, tail)",
         y = "-1 x baseline ctrl NME (top = most ordered)",
         title = "Directed-mechanism loci: ordered TFBS in DE genes",
         subtitle = sprintf("%s loci | colored by delta NME (green = ordering, purple = disorder)",
                            format(nrow(vol_dt), big.mark = ","))) +
    theme_minimal(base_size = 12)
  save_fig(p11c, BATCH_DIR, "fig11c_directed_loci_volcano", w = 9, h = 7)
}

# Fig 11d — reads per window: in vs out of motif
reads_long <- rbind(
  data.table(group = "In motif",     reads = win_pass[in_motif == TRUE,  total_reads]),
  data.table(group = "Out of motif", reads = win_pass[in_motif == FALSE, total_reads])
)
p11d <- ggplot(reads_long, aes(x = group, y = reads, fill = group)) +
  geom_violin(alpha = 0.55, linewidth = 0.4) +
  geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.9) +
  scale_fill_manual(values = c("In motif" = "#C0392B",
                                "Out of motif" = "#2471A3"),
                    guide = "none") +
  scale_y_log10() +
  labs(x = NULL, y = "Total reads per window (log10, sum across 4 samples)",
       title = "Read coverage per 4-CpG window: TFBS vs non-TFBS",
       subtitle = sprintf("min %d reads/sample required", MIN_READS)) +
  theme_minimal(base_size = 12)
save_fig(p11d, BATCH_DIR, "fig11d_reads_per_window_by_motif", w = 7, h = 6)

# Fig 11e — reads per motif class (top 20 by n_windows)
top_classes <- cov_diag[1:min(20, .N)]
p11e <- ggplot(top_classes,
               aes(x = reorder(tf_class_top, mean_reads), y = mean_reads)) +
  geom_col(fill = "#34495E") +
  geom_errorbar(aes(ymin = q25_reads, ymax = q75_reads), width = 0.3,
                color = "grey20") +
  coord_flip() +
  labs(x = NULL, y = "Mean reads per window (Q25-Q75 bars)",
       title = "Read coverage per TF class (top 20 by window count)",
       subtitle = "Sanity: are well-represented TF classes also well-covered?") +
  theme_minimal(base_size = 11)
save_fig(p11e, BATCH_DIR, "fig11e_reads_per_motif_class", w = 9, h = 7)

# Fig 11f — reads vs NME (in-motif only) — should be flat above MIN_READS
p11f <- ggplot(motif_w[!is.na(ctrl_nme)],
               aes(x = total_reads, y = ctrl_nme)) +
  geom_hex(bins = 50) +
  scale_x_log10() +
  scale_fill_viridis_c(trans = "log10") +
  geom_smooth(method = "loess", color = "#C0392B", se = FALSE, linewidth = 0.8) +
  labs(x = "Total reads per window (log10)",
       y = "Baseline NME (control)",
       title = "NME vs read depth at TFBS windows",
       subtitle = "Sanity: NME should be flat above min-reads threshold") +
  theme_minimal(base_size = 12)
save_fig(p11f, BATCH_DIR, "fig11f_reads_vs_nme", w = 8, h = 6)

# =============================================================================
# [FANG TABLE 1] Per-motif NME and MML binomial tests
# =============================================================================
# For each motif: test if windows overlapping that motif have significantly
# different NME (high entropy = disordered) or MML (low mean methylation level)
# compared to beta-matched background. Binomial test per motif.
# Output: ranked TF table like Fang 2023 Table 1.
# =============================================================================
cat("\n[FANG] Per-motif NME/MML binomial tests...\n")

# Need all window-motif overlaps (ov_dt built in section 2)
# For each motif: windows in that motif vs genome-wide median
genome_median_nme <- median(win_pass$ctrl_nme, na.rm = TRUE)
genome_median_mml <- median(win_pass[, (ctrl_nme + ampu_nme) / 2], na.rm = TRUE)  # approximate MML
# Actually MML = mean methylation level. Use win_beta from batch10b if available
if ("win_beta" %in% names(win_pass)) {
  genome_median_mml <- median(win_pass$win_beta, na.rm = TRUE)
}

cat(sprintf("  Genome median NME: %.4f | Genome median MML (beta): %.4f\n",
            genome_median_nme, genome_median_mml))

# Unique motifs in overlapping windows
motif_ids_in_data <- unique(ov_dt$motif_id)
cat(sprintf("  Motifs with window overlaps: %d\n", length(motif_ids_in_data)))

# Per-motif test: for each motif, count how many of its windows have
# NME > genome median (high NME) or MML < genome median (low MML)
# Binomial test: is the proportion > 0.5?

per_motif_test <- function(motif_id_val) {
  win_idx <- unique(ov_dt[motif_id == motif_id_val, win_idx])
  motif_wins <- win_pass[win_idx]
  motif_wins <- motif_wins[!is.na(ctrl_nme)]
  n <- nrow(motif_wins)
  if (n < 10) return(NULL)

  # NME test: is this motif associated with HIGH NME?
  n_high_nme <- sum(motif_wins$ctrl_nme > genome_median_nme)
  p_nme <- binom.test(n_high_nme, n, p = 0.5, alternative = "greater")$p.value

  # MML test: is this motif associated with LOW MML (beta)?
  if ("win_beta" %in% names(motif_wins)) {
    n_low_mml <- sum(motif_wins$win_beta < genome_median_mml, na.rm = TRUE)
    n_mml <- sum(!is.na(motif_wins$win_beta))
    p_mml <- if (n_mml >= 10) binom.test(n_low_mml, n_mml, p = 0.5, alternative = "greater")$p.value else NA_real_
  } else {
    p_mml <- NA_real_
    n_low_mml <- NA_integer_; n_mml <- NA_integer_
  }

  # Effect sizes
  mean_nme <- mean(motif_wins$ctrl_nme)
  mean_mml <- if ("win_beta" %in% names(motif_wins)) mean(motif_wins$win_beta, na.rm = TRUE) else NA_real_

  mn <- ov_dt[motif_id == motif_id_val, motif_name[1]]
  tc <- ov_dt[motif_id == motif_id_val, tf_class[1]]

  data.table(
    motif_id = motif_id_val,
    motif_name = mn,
    tf_class = tc,
    n_windows = n,
    mean_nme = mean_nme,
    pct_high_nme = 100 * n_high_nme / n,
    p_high_nme = p_nme,
    mean_mml = mean_mml,
    pct_low_mml = if (!is.na(n_mml) && n_mml > 0) 100 * n_low_mml / n_mml else NA_real_,
    p_low_mml = p_mml
  )
}

cat("  Running per-motif binomial tests...\n")
motif_tests <- rbindlist(lapply(motif_ids_in_data, per_motif_test))
cat(sprintf("  Motifs tested: %d\n", nrow(motif_tests)))

# FDR correction
motif_tests[, fdr_high_nme := p.adjust(p_high_nme, "BH")]
motif_tests[, fdr_low_mml := p.adjust(p_low_mml, "BH")]

# Assign Yin class
motif_tests[, methyl_sens := assign_methyl_sens(motif_name, tf_class)]

# Classify: High NME only, Low MML only, Both, Neither
motif_tests[, category := fcase(
  fdr_high_nme < 0.1 & (is.na(fdr_low_mml) | fdr_low_mml >= 0.1), "High NME only",
  (is.na(fdr_high_nme) | fdr_high_nme >= 0.1) & fdr_low_mml < 0.1, "Low MML only",
  fdr_high_nme < 0.1 & fdr_low_mml < 0.1,                          "Both",
  default = "Neither"
)]

cat("  Category breakdown:\n")
print(motif_tests[, .N, by = category])

# Sort and save
setorder(motif_tests, p_high_nme)
save_data(motif_tests, BATCH_DIR, "fang_table1_motif_nme_mml")

# Extract the two table halves (like Fang Table 1)
high_nme_only <- motif_tests[category == "High NME only"][order(fdr_high_nme)]
low_mml_only  <- motif_tests[category == "Low MML only"][order(fdr_low_mml)]

cat(sprintf("\n  High NME-associated motifs (not low MML): %d\n", nrow(high_nme_only)))
if (nrow(high_nme_only) > 0) {
  print(head(high_nme_only[, .(motif_name, tf_class, methyl_sens, fdr_high_nme, n_windows, mean_nme)], 20))
}
cat(sprintf("\n  Low MML-associated motifs (not high NME): %d\n", nrow(low_mml_only)))
if (nrow(low_mml_only) > 0) {
  print(head(low_mml_only[, .(motif_name, tf_class, methyl_sens, fdr_low_mml, n_windows, mean_mml)], 20))
}

save_data(high_nme_only, BATCH_DIR, "fang_high_nme_only_motifs")
save_data(low_mml_only, BATCH_DIR, "fang_low_mml_only_motifs")

# --- Fig 11g: Fang Table 1 as bar plot (top significant motifs) ---
cat("\nFig 11g: Fang-style motif NME/MML bar plot...\n")

# Combined significant motifs
sig_motifs <- motif_tests[fdr_high_nme < 0.1 | fdr_low_mml < 0.1]
if (nrow(sig_motifs) > 0) {
  # For plotting: show top 30 by most significant test
  sig_motifs[, best_fdr := pmin(fdr_high_nme, fdr_low_mml, na.rm = TRUE)]
  sig_motifs[, best_test := fifelse(fdr_high_nme <= fdr_low_mml | is.na(fdr_low_mml),
                                     "High NME", "Low MML")]
  top_sig <- sig_motifs[order(best_fdr)][1:min(40, .N)]
  top_sig[, motif_label := paste0(motif_name, " (", methyl_sens, ")")]
  top_sig[, motif_label := factor(motif_label, levels = rev(motif_label))]

  p11g <- ggplot(top_sig, aes(x = motif_label, y = -log10(best_fdr), fill = category)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed", color = "gray40") +
    scale_fill_manual(values = c("High NME only" = "#E74C3C",
                                  "Low MML only" = "#3498DB",
                                  "Both" = "#8E44AD")) +
    coord_flip() +
    labs(x = NULL, y = "-log10(FDR)",
         title = "TF motifs associated with high NME or low MML (Fang framework)",
         subtitle = sprintf("Binomial test per motif | %d high NME, %d low MML, %d both",
                            sum(motif_tests$category == "High NME only"),
                            sum(motif_tests$category == "Low MML only"),
                            sum(motif_tests$category == "Both")),
         fill = "Category") +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(size = 7))
  save_fig(p11g, BATCH_DIR, "fig11g_fang_motif_nme_mml", w = 12, h = 10)
}

# --- Fig 11h: Heatmap of top motifs × NME/MML ---
if (nrow(sig_motifs) >= 5) {
  cat("Fig 11h: Motif × NME/MML heatmap...\n")

  top_for_heat <- sig_motifs[order(best_fdr)][1:min(30, .N)]
  heat_mat <- as.matrix(top_for_heat[, .(pct_high_nme, pct_low_mml)])
  rownames(heat_mat) <- top_for_heat$motif_name
  colnames(heat_mat) <- c("% High NME", "% Low MML")

  ht2 <- Heatmap(heat_mat,
    name = "%",
    col = colorRamp2(c(40, 50, 70), c("#2166AC", "white", "#B2182B")),
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 10),
    column_title = "Motif-associated NME and MML (Fang Table 1 analogue)",
    right_annotation = rowAnnotation(
      Sens = top_for_heat$methyl_sens,
      col = list(Sens = c("MethylPlus" = "#E74C3C", "MethylMinus" = "#3498DB",
                           "LittleEffect" = "#95A5A6", "Mixed" = "#F39C12",
                           "Unknown" = "#BDC3C7")),
      show_legend = TRUE
    ))

  png(file.path(BATCH_DIR, "figures/fig11h_fang_motif_heatmap.png"),
      width = 8, height = 10, units = "in", res = 300)
  draw(ht2)
  dev.off()
  cairo_pdf(file.path(BATCH_DIR, "figures/fig11h_fang_motif_heatmap.pdf"),
            width = 8, height = 10)
  draw(ht2)
  dev.off()
  cat("  Saved fig11h\n")
}

# =============================================================================
elapsed <- (proc.time() - t0)[3]
cat(sprintf("\n=== Batch 11 complete (%.1f minutes) ===\n", elapsed / 60))
cat(sprintf("  windows processed:   %s\n", format(nrow(win_pass), big.mark = ",")))
cat(sprintf("  in-motif windows:    %s\n", format(sum(win_pass$in_motif), big.mark = ",")))
cat(sprintf("  shortlist size:      %s\n", format(nrow(shortlist_out), big.mark = ",")))
if (!is.na(wt_A$p.value)) {
  cat(sprintf("  Test A (in vs out):  delta = %.4f, p = %s\n",
              mean(matched$in_nme) - mean(matched$out_nme),
              format(wt_A$p.value, digits = 3)))
}
