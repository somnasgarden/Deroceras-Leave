# =============================================================================
# Shared configuration for all methylation pipeline batches
# Source this at the top of every batch script:
#   source("methylation_pipeline/_config.R")
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

# --- Paths ---
PROJECT_DIR <- "C:/Users/rafae/Projects/STANDBY"
DATA_DIR    <- "C:/Users/rafae/Projects/DATA"
CACHE_DIR   <- file.path(PROJECT_DIR, "genome/cache")
PIPE_DIR    <- file.path(PROJECT_DIR, "methylation_pipeline")

keep_chr <- paste0("chr", 1:31)

# --- OG data files ---
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
    cairo_pdf(file.path(fig_dir, paste0(name, ".pdf")),
              width = w, height = h)
    print(p); dev.off()
  }, error = function(e) {
    pdf(file.path(fig_dir, paste0(name, ".pdf")),
        width = w, height = h)
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

cat("Config loaded. PROJECT_DIR:", PROJECT_DIR, "\n")
