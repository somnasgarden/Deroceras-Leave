#!/usr/bin/env Rscript
# =============================================================================
# Full-genome DMLtest with smoothing — Cluster job
# Reads BSseq from CpG reports (or cache), runs DMLtest, saves DMPs/DMRs
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(bsseq)
  library(DSS)
  library(GenomicRanges)
  library(data.table)
})

cat("=== Full-Genome DMLtest (Cluster) ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

# --- Paths ---
CLUSTER_ROOT <- Sys.getenv("CLUSTER_ROOT", ".")
cache_dir    <- file.path(CLUSTER_ROOT, "genome/cache")
out_dir      <- file.path(CLUSTER_ROOT, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# CpG report files on cluster (compressed)
cpg_dir <- "/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/09_methylation_calls"

sample_info <- data.frame(
  sample_id = c("C1", "C2", "A1", "A2"),
  condition = c("Control", "Control", "Amputated", "Amputated"),
  file_path = file.path(cpg_dir, paste0(
    c("C1", "C2", "A1", "A2"), ".CpG_report.txt.gz"
  ))
)

keep_chr <- paste0("chr", 1:31)

# --- Load or build BSseq ---
bsseq_cache <- file.path(cache_dir, "bsseq_full.rds")

if (file.exists(bsseq_cache)) {
  cat("Loading cached BSseq...\n")
  bs_obj <- readRDS(bsseq_cache)
  pData(bs_obj)$condition <- sample_info$condition
  pData(bs_obj)$sample_id <- sample_info$sample_id
  cat(sprintf("BSseq: %s sites x %d samples\n\n",
              format(nrow(bs_obj), big.mark = ","), ncol(bs_obj)))
} else {
  cat("Building BSseq from CpG reports...\n\n")

  min_cov <- 5
  tmp_dir <- tempdir()
  common_dt <- NULL

  for (i in seq_len(nrow(sample_info))) {
    sn <- sample_info$sample_id[i]
    fp <- sample_info$file_path[i]
    cat(sprintf("  Reading %s from %s...\n", sn, fp))
    t0 <- proc.time()

    dt <- fread(fp, header = FALSE, sep = "\t",
                select = c(1, 2, 3, 4, 5, 6),
                col.names = c("chr", "pos", "strand", "meth", "unmeth", "context"))

    cat(sprintf("    Raw rows: %s\n", format(nrow(dt), big.mark = ",")))

    # Filter: CG context, chr1-31
    dt <- dt[context == "CG" & chr %in% keep_chr]
    dt[, context := NULL]

    # Strand collapse
    dt[strand == "-", pos := pos - 1L]
    dt[, strand := NULL]
    dt <- dt[, .(meth = sum(meth), unmeth = sum(unmeth)), by = .(chr, pos)]

    # Coverage filter
    dt <- dt[meth + unmeth >= min_cov]
    cat(sprintf("    After collapse + cov>=%d: %s\n",
                min_cov, format(nrow(dt), big.mark = ",")))

    # Save temp
    tmp_rds <- file.path(tmp_dir, paste0("cpg_", sn, ".rds"))
    saveRDS(dt, tmp_rds)

    # Common sites
    site_dt <- dt[, .(chr, pos)]
    if (is.null(common_dt)) {
      common_dt <- site_dt
    } else {
      setkeyv(common_dt, c("chr", "pos"))
      setkeyv(site_dt, c("chr", "pos"))
      common_dt <- fintersect(common_dt, site_dt)
    }
    cat(sprintf("    Running common: %s\n", format(nrow(common_dt), big.mark = ",")))

    elapsed <- (proc.time() - t0)[3]
    cat(sprintf("    Done in %.0f sec\n\n", elapsed))
    rm(dt, site_dt); gc(verbose = FALSE)
  }

  cat(sprintf("Common sites (cov >= %d, all 4 samples): %s\n\n",
              min_cov, format(nrow(common_dt), big.mark = ",")))

  # Build BSseq
  setkeyv(common_dt, c("chr", "pos"))
  n_sites <- nrow(common_dt)
  M_mat <- matrix(0L, nrow = n_sites, ncol = 4)
  Cov_mat <- matrix(0L, nrow = n_sites, ncol = 4)
  colnames(M_mat) <- colnames(Cov_mat) <- sample_info$sample_id
  gr <- NULL

  for (i in 1:4) {
    sn <- sample_info$sample_id[i]
    cat(sprintf("  Building matrix for %s...\n", sn))
    dt <- readRDS(file.path(tmp_dir, paste0("cpg_", sn, ".rds")))
    setkeyv(dt, c("chr", "pos"))
    dt <- dt[common_dt, nomatch = 0L]
    if (is.null(gr)) {
      gr <- GRanges(seqnames = dt$chr, ranges = IRanges(start = dt$pos, width = 1))
    }
    M_mat[, i] <- dt$meth
    Cov_mat[, i] <- dt$meth + dt$unmeth
    rm(dt); gc(verbose = FALSE)
    file.remove(file.path(tmp_dir, paste0("cpg_", sn, ".rds")))
  }
  rm(common_dt); gc(verbose = FALSE)

  bs_obj <- BSseq(gr = gr, M = M_mat, Cov = Cov_mat)
  pData(bs_obj)$condition <- sample_info$condition
  pData(bs_obj)$sample_id <- sample_info$sample_id
  rm(gr, M_mat, Cov_mat); gc(verbose = FALSE)

  cat(sprintf("BSseq: %s sites x %d samples\n",
              format(nrow(bs_obj), big.mark = ","), ncol(bs_obj)))

  cat("Caching BSseq...\n")
  saveRDS(bs_obj, bsseq_cache)
  cat(sprintf("Cached: %s (%.1f MB)\n\n",
              bsseq_cache, file.size(bsseq_cache) / 1e6))
}

# --- Per-sample stats ---
cat("=== Per-sample methylation stats ===\n")
for (i in 1:ncol(bs_obj)) {
  bv <- getMeth(bs_obj[, i], type = "raw")[, 1]
  bv <- bv[!is.nan(bv) & !is.na(bv)]
  cat(sprintf("  %s: %s sites, mean=%.4f, median=%.4f, high(>0.8)=%s, low(<0.2)=%s\n",
              colnames(bs_obj)[i], format(length(bv), big.mark = ","),
              mean(bv), median(bv),
              format(sum(bv > 0.8), big.mark = ","),
              format(sum(bv < 0.2), big.mark = ",")))
  rm(bv); gc(verbose = FALSE)
}

# --- DMLtest ---
cat("\n=== DMLtest with smoothing ===\n")
group1 <- c("C1", "C2")
group2 <- c("A1", "A2")
cat(sprintf("Testing: %s vs %s\n", paste(group1, collapse = "+"), paste(group2, collapse = "+")))
cat(sprintf("Sites: %s\n", format(nrow(bs_obj), big.mark = ",")))
cat("Start:", format(Sys.time()), "\n\n")

t0 <- proc.time()
dml_test <- DMLtest(bs_obj, group1 = group1, group2 = group2, smoothing = TRUE)
dml_time <- (proc.time() - t0)[3]
cat(sprintf("\nDMLtest done in %.1f minutes\n", dml_time / 60))
cat(sprintf("Sites tested: %s\n", format(nrow(dml_test), big.mark = ",")))

# Cache DMLtest
dml_cache <- file.path(out_dir, "dmltest_full.rds")
saveRDS(dml_test, dml_cache)
cat(sprintf("Cached: %s (%.1f MB)\n", dml_cache, file.size(dml_cache) / 1e6))

rm(bs_obj); gc(verbose = FALSE)

# --- DMPs ---
cat("\n=== Calling DMPs ===\n")
dmp_strict <- callDML(dml_test, p.threshold = 0.05, delta = 0.1)
dmp_relaxed <- callDML(dml_test, p.threshold = 0.01)

n_strict <- if (!is.null(dmp_strict)) nrow(dmp_strict) else 0
n_relaxed <- if (!is.null(dmp_relaxed)) nrow(dmp_relaxed) else 0

cat(sprintf("DMPs (FDR<0.05, |diff|>10%%): %s\n", format(n_strict, big.mark = ",")))
cat(sprintf("DMPs (p<0.01, no delta):      %s\n", format(n_relaxed, big.mark = ",")))

if (n_strict > 0) {
  n_hyper <- sum(dmp_strict$diff > 0)
  n_hypo <- sum(dmp_strict$diff < 0)
  cat(sprintf("  Hyper: %s  Hypo: %s\n", format(n_hyper, big.mark = ","),
              format(n_hypo, big.mark = ",")))
  fwrite(dmp_strict, file.path(out_dir, "dmps_strict.tsv"), sep = "\t")
}
if (n_relaxed > 0) {
  fwrite(dmp_relaxed, file.path(out_dir, "dmps_relaxed.tsv"), sep = "\t")
}

# --- DMRs ---
cat("\n=== Calling DMRs ===\n")
dmr_results <- callDMR(dml_test, p.threshold = 0.05, delta = 0.1,
                       minlen = 50, minCG = 3)

if (is.null(dmr_results) || nrow(dmr_results) == 0) {
  cat("No DMRs found.\n")
  dmr_results <- data.frame(chr = character(), start = integer(), end = integer(),
                            length = integer(), nCG = integer(), areaStat = numeric())
} else {
  n_dmr <- nrow(dmr_results)
  n_hyper_dmr <- sum(dmr_results$areaStat > 0)
  n_hypo_dmr <- sum(dmr_results$areaStat < 0)
  cat(sprintf("DMRs: %s (hyper: %s, hypo: %s)\n",
              format(n_dmr, big.mark = ","),
              format(n_hyper_dmr, big.mark = ","),
              format(n_hypo_dmr, big.mark = ",")))
  cat(sprintf("Mean length: %.0f bp, mean CpGs: %.1f\n",
              mean(dmr_results$length), mean(dmr_results$nCG)))
}
fwrite(as.data.frame(dmr_results), file.path(out_dir, "dmrs.tsv"), sep = "\t")

rm(dml_test); gc(verbose = FALSE)

# --- Summary ---
cat("\n=== Summary ===\n")
summary_df <- data.frame(
  Metric = c("Total CpGs tested", "DMPs (strict)", "DMPs (relaxed)",
             "DMRs", "DMLtest time (min)"),
  Value = c(format(n_relaxed + n_strict, big.mark = ","),
            format(n_strict, big.mark = ","),
            format(n_relaxed, big.mark = ","),
            format(nrow(dmr_results), big.mark = ","),
            round(dml_time / 60, 1))
)
fwrite(summary_df, file.path(out_dir, "dmltest_summary.tsv"), sep = "\t")
print(summary_df)

cat("\nOutput files:\n")
for (f in list.files(out_dir, pattern = "\\.(tsv|rds)$")) {
  cat(sprintf("  %s (%.1f MB)\n", f, file.size(file.path(out_dir, f)) / 1e6))
}
cat("\nDone:", format(Sys.time()), "\n")
