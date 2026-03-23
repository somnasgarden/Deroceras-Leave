# =============================================================================
# Shared configuration for all methylation pipeline batches
# Auto-detects local (Windows) vs cluster (HPC Linux) environment
# Source this at the top of every batch script:
#   source("methylation_pipeline/_config.R")
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)
keep_chr <- paste0("chr", 1:31)

# --- Auto-detect environment ---
ON_CLUSTER <- dir.exists("/mnt/data/alfredvar")

if (ON_CLUSTER) {
  # ---- CLUSTER PATHS (HPC) ----
  cat("=== CLUSTER environment detected ===\n")
  BASE_DATA   <- "/mnt/data/alfredvar"
  # Find repo root (works whether sourced from methylation_pipeline/ or repo root)
  PROJECT_DIR <- Sys.getenv("PROJECT_DIR",
    if (file.exists("methylation_pipeline/_config.R")) getwd()
    else if (file.exists("_config.R")) dirname(getwd())
    else stop("Cannot find project root. Set PROJECT_DIR env var or cd to repo root."))
  DATA_DIR    <- PROJECT_DIR  # OG data comes from cluster paths below
  CACHE_DIR   <- file.path(PROJECT_DIR, "genome/cache")
  PIPE_DIR    <- file.path(PROJECT_DIR, "methylation_pipeline")

  OG <- list(
    cpg_C1     = file.path(BASE_DATA, "jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/C1.CpG_report.txt.gz"),
    cpg_C2     = file.path(BASE_DATA, "jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/C2.CpG_report.txt.gz"),
    cpg_A1     = file.path(BASE_DATA, "jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/A1.CpG_report.txt.gz"),
    cpg_A2     = file.path(BASE_DATA, "jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/A2.CpG_report.txt.gz"),
    gff        = file.path(BASE_DATA, "30-Genoma/31-Alternative_Annotation_EviAnn/derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff"),
    annot      = file.path(BASE_DATA, "30-Genoma/31-Alternative_Annotation_EviAnn/derLaeGenome_eviann_annotations.tsv"),
    te_age     = file.path(BASE_DATA, "30-Genoma/32-Repeats/age_of_transposons/collapsed_te_age_data.tsv"),
    counts_dir = file.path(BASE_DATA, "jmiranda/20-Transcriptomic_Bulk/25-metaAnalysisTranscriptome/counts_HTseq_EviAnn"),
    deeptf     = file.path(BASE_DATA, "rlopezt/DeepFactor1/DeepFactorV1/deeptfactor/result/prediction_result.txt"),
    string     = file.path(BASE_DATA, "rlopezt/Preliminary/STRG0A31YWK.protein.orthology.v12.0.txt"),
    genie3     = file.path(BASE_DATA, "wgutierrez/genie3_2/genie3_all_links.tsv")
  )
} else {
  # ---- LOCAL PATHS (Windows) ----
  cat("=== LOCAL environment detected ===\n")
  PROJECT_DIR <- "C:/Users/rafae/Projects/STANDBY"
  DATA_DIR    <- "C:/Users/rafae/Projects/DATA"
  CACHE_DIR   <- file.path(PROJECT_DIR, "genome/cache")
  PIPE_DIR    <- file.path(PROJECT_DIR, "methylation_pipeline")

  OG <- list(
    cpg_C1     = file.path(DATA_DIR, "C1.CpG_report.txt"),
    cpg_C2     = file.path(DATA_DIR, "C2.CpG_report.txt"),
    cpg_A1     = file.path(DATA_DIR, "A1.CpG_report.txt"),
    cpg_A2     = file.path(DATA_DIR, "A2.CpG_report.txt"),
    gff        = file.path(DATA_DIR, "derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff"),
    annot      = file.path(DATA_DIR, "derLaeGenome_eviann_annotations.tsv"),
    te_age     = file.path(DATA_DIR, "collapsed_te_age_data.tsv"),
    counts_dir = file.path(DATA_DIR, "counts_HTseq_EviAnn"),
    deeptf     = file.path(DATA_DIR, "prediction_result.txt"),
    string     = file.path(DATA_DIR, "STRG0A31YWK.protein.orthology.v12.0.txt"),
    genie3     = file.path(DATA_DIR, "genie3_top500k.tsv")
  )
}

# --- Verify critical paths exist ---
for (nm in c("gff", "te_age", "counts_dir", "deeptf")) {
  if (!file.exists(OG[[nm]])) cat("  WARNING: missing", nm, "at", OG[[nm]], "\n")
}

# --- Cache dir ---
dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Cached objects ---
CACHE <- list(
  genome     = file.path(CACHE_DIR, "genome_chr1_31.rds"),
  gff        = file.path(CACHE_DIR, "gff_chr1_31.rds"),
  te         = file.path(CACHE_DIR, "te_chr1_31.rds"),
  promoters  = file.path(CACHE_DIR, "promoters_2kb.rds"),
  bsseq      = file.path(CACHE_DIR, "bsseq_tutorial.rds"),
  dmltest    = file.path(CACHE_DIR, "dmltest_full.rds"),
  extended   = file.path(CACHE_DIR, "extended_regions_10kb.rds")
)

# --- Plotting conventions ---
COLORS <- list(
  direction = c(Hyper = "#C0392B", Hypo = "#2471A3"),
  region    = c(Promoter = "#8E44AD", Exon = "#2471A3",
                Intron = "#1ABC9C", Intergenic = "#C0392B"),
  te        = c(LTR = "#2471A3", RC = "#F39C12", LINE = "#8E44AD",
                DNA = "#2471A3", Unknown = "#27AE60", SINE = "#C0392B"),
  condition = c(Control = "#2471A3", Amputated = "#C0392B")
)

# --- Save helper: figures only (PDF + PNG, no HTML) ---
save_fig <- function(p, batch_dir, name, w = 8, h = 5) {
  fig_dir <- file.path(batch_dir, "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  png(file.path(fig_dir, paste0(name, ".png")),
      width = w, height = h, units = "in", res = 300)
  print(p); dev.off()
  tryCatch({
    cairo_pdf(file.path(fig_dir, paste0(name, ".pdf")), width = w, height = h)
    print(p); dev.off()
  }, error = function(e) {
    pdf(file.path(fig_dir, paste0(name, ".pdf")), width = w, height = h)
    print(p); dev.off()
  })
  cat("  Saved:", name, "\n")
}

# --- Save data helper ---
save_data <- function(df, batch_dir, name) {
  data_dir <- file.path(batch_dir, "data")
  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  fpath <- file.path(data_dir, paste0(name, ".tsv"))
  data.table::fwrite(as.data.frame(df), fpath, sep = "\t")
  cat("  Data saved:", name, ".tsv\n")
}

# --- Load genome (BSgenome fallback) ---
load_genome <- function() {
  if (file.exists(CACHE$genome)) {
    cat("Loading genome from cache...\n")
    readRDS(CACHE$genome)
  } else {
    cat("Loading from BSgenome...\n")
    library(BSgenome.Dlaeve.NCBI.dlgm)
    bsg <- BSgenome.Dlaeve.NCBI.dlgm
    seqs <- getSeq(bsg, intersect(keep_chr, seqnames(bsg)))
    names(seqs) <- intersect(keep_chr, seqnames(bsg))
    cat("Caching genome RDS...\n")
    saveRDS(seqs, CACHE$genome)
    rm(bsg); gc(verbose = FALSE)
    seqs
  }
}

# --- Load GFF ---
load_gff <- function() {
  if (file.exists(CACHE$gff)) {
    readRDS(CACHE$gff)
  } else {
    gff <- rtracklayer::import(OG$gff)
    gff <- gff[seqnames(gff) %in% keep_chr]
    gff <- keepSeqlevels(gff, keep_chr, pruning.mode = "coarse")
    saveRDS(gff, CACHE$gff)
    gff
  }
}

# --- Load TE ---
load_te <- function() {
  if (file.exists(CACHE$te)) {
    te_data <- readRDS(CACHE$te)
    if (is(te_data, "list") && "te_gr" %in% names(te_data)) {
      list(te_gr = te_data$te_gr, te_regions = te_data$te_regions)
    } else {
      list(te_gr = te_data, te_regions = GenomicRanges::reduce(te_data))
    }
  } else {
    te_dt <- data.table::fread(OG$te_age)
    te_dt <- te_dt[chrom %in% keep_chr]
    te_gr <- GenomicRanges::makeGRangesFromDataFrame(te_dt,
      seqnames.field = "chrom", start.field = "start", end.field = "end",
      keep.extra.columns = TRUE)
    te_regions <- GenomicRanges::reduce(te_gr)
    saveRDS(list(te_gr = te_gr, te_regions = te_regions), CACHE$te)
    list(te_gr = te_gr, te_regions = te_regions)
  }
}

# --- Functional region annotation (mutually exclusive, no TE) ---
annotate_regions <- function(chr, pos, promoters, exons, genes) {
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = pos, width = 1))
  ann <- rep("Intergenic", length(gr))
  ann[overlapsAny(gr, genes) & !overlapsAny(gr, exons)] <- "Intron"
  ann[overlapsAny(gr, exons)] <- "Exon"
  ann[overlapsAny(gr, promoters)] <- "Promoter"
  ann
}

cat("Config loaded. Environment:", ifelse(ON_CLUSTER, "CLUSTER", "LOCAL"), "\n")
cat("  PROJECT_DIR:", PROJECT_DIR, "\n")
cat("  CACHE_DIR:", CACHE_DIR, "\n")
cat("  PIPE_DIR:", PIPE_DIR, "\n\n")
