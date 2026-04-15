#!/usr/bin/env Rscript
# =============================================================================
# Cluster: Expanded Motif Annotation — JASPAR + HOMER
# =============================================================================
# Supersedes batch1.5_motif_positions.R (which only scanned 2kb promoters).
#
# What this does:
#   A. Load genome + GFF, define extended regions (10kb upstream + gene body +
#      10kb downstream per gene)
#   B. Export promoter FASTA + convert JASPAR to HOMER format (for HOMER steps)
#   C. Run HOMER de novo motif discovery on promoters (system call)
#   D. JASPAR motif scanning on extended regions (per-chromosome, memory-safe)
#   E. Annotate every hit with: gene_id, region_type, distance_to_tss
#   F. Save comprehensive output
#
# REQUIRES: ~64 GB RAM, HOMER in PATH (run setup_homer.sh first)
#
# Input:
#   - Genome FASTA or genome_chr1_31.rds in genome/cache/
#   - GFF (from cluster path or gff_chr1_31.rds in genome/cache/)
#   - JASPAR 2024 (downloaded via BiocFileCache)
#   - HOMER installed (for de novo — optional, script continues without it)
#
# Output (all in cluster/results/):
#   - motif_hits_extended.tsv.gz         Full annotated motif hit table
#   - motif_summary_by_region.tsv        Per-motif per-region counts
#   - motif_region_overview.tsv          Region-level summary
#   - promoters_for_homer.fa             Promoter sequences (FASTA)
#   - jaspar_motifs_homer.txt            JASPAR in HOMER format
#   - homer_denovo/                      HOMER de novo results
# =============================================================================

source("methylation_pipeline/_config.R")

# -- Cache dir for JASPAR/BiocFileCache --
cache_home <- Sys.getenv("XDG_CACHE_HOME", file.path(PROJECT_DIR, "cluster/.cache"))
Sys.setenv(XDG_CACHE_HOME = cache_home)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(Biostrings)
  library(GenomeInfoDb)
  library(TFBSTools)
  library(JASPAR2024)
  library(motifmatchr)
  library(universalmotif)
})

# -- Paths (from _config.R: PROJECT_DIR, CACHE_DIR, CACHE, OG, keep_chr) --
out_dir <- file.path(PROJECT_DIR, "cluster/results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Expanded Motif Annotation (Cluster) ===\n")
cat("Project dir:", PROJECT_DIR, "\n")
cat("Cache dir:", CACHE_DIR, "\n")
cat("Output dir:", out_dir, "\n")
cat("Available RAM:", system("free -h | grep Mem | awk '{print $2}'", intern = TRUE), "\n\n")


# =============================================================================
# A. LOAD GENOME + GFF, DEFINE EXTENDED REGIONS
# =============================================================================
cat("=== A. Loading data and defining regions ===\n")

# -- A1. Load genome (via _config.R loader — caches automatically) --
genome <- load_genome()
chr_lengths <- setNames(width(genome), names(genome))
cat(sprintf("Genome: %d chromosomes, %s bp total\n",
            length(genome), format(sum(chr_lengths), big.mark = ",")))

# -- A2. Load GFF (via _config.R loader — caches automatically) --
gff <- load_gff()

genes <- gff[gff$type == "gene"]
cat("Genes:", length(genes), "\n")

# -- A3. Define extended regions --
regions_rds <- CACHE$extended
if (file.exists(regions_rds)) {
  cat("Loading extended regions from cache...\n")
  regions <- readRDS(regions_rds)
} else {
  cat("Defining extended regions (10kb upstream + gene body + 10kb downstream)...\n")

  gene_strand <- as.character(strand(genes))
  gene_chr    <- as.character(seqnames(genes))
  gene_ids    <- genes$ID
  chrlen      <- as.integer(chr_lengths[gene_chr])

  tss <- ifelse(gene_strand == "+", start(genes), end(genes))
  tes <- ifelse(gene_strand == "+", end(genes), start(genes))

  build_gr <- function(s, e, type) {
    valid <- !is.na(s) & !is.na(e) & s <= e & s >= 1L
    GRanges(
      seqnames    = gene_chr[valid],
      ranges      = IRanges(start = s[valid], end = e[valid]),
      strand      = gene_strand[valid],
      gene_id     = gene_ids[valid],
      region_type = type
    )
  }

  # Promoter: 2kb upstream of TSS
  prom_s <- ifelse(gene_strand == "+", pmax(1L, tss - 2000L), tss + 1L)
  prom_e <- ifelse(gene_strand == "+", tss - 1L, pmin(chrlen, tss + 2000L))

  # Upstream distal: 2-10kb upstream
  up_s <- ifelse(gene_strand == "+", pmax(1L, tss - 10000L), tss + 2001L)
  up_e <- ifelse(gene_strand == "+", pmax(1L, tss - 2001L), pmin(chrlen, tss + 10000L))

  # Gene body
  gb_s <- start(genes)
  gb_e <- end(genes)

  # Downstream: 10kb past TES
  dn_s <- ifelse(gene_strand == "+", tes + 1L, pmax(1L, tes - 10000L))
  dn_e <- ifelse(gene_strand == "+", pmin(chrlen, tes + 10000L), tes - 1L)

  regions <- c(
    build_gr(prom_s, prom_e, "promoter"),
    build_gr(up_s,   up_e,   "upstream_distal"),
    build_gr(gb_s,   gb_e,   "gene_body"),
    build_gr(dn_s,   dn_e,   "downstream")
  )
  seqlevels(regions) <- keep_chr
  seqlengths(regions) <- chr_lengths[keep_chr]
  regions <- trim(regions)
  saveRDS(regions, regions_rds)
}

for (rt in c("promoter", "upstream_distal", "gene_body", "downstream")) {
  n <- sum(regions$region_type == rt)
  cat(sprintf("  %-18s %6d regions\n", rt, n))
}

# -- A4. TSS lookup table (for distance calculation) --
gene_strand <- as.character(strand(genes))
tss_vals    <- ifelse(gene_strand == "+", start(genes), end(genes))
tss_lookup  <- data.frame(
  gene_id = genes$ID,
  chr     = as.character(seqnames(genes)),
  tss     = tss_vals,
  strand  = gene_strand,
  stringsAsFactors = FALSE
)


# =============================================================================
# B. EXPORT FOR HOMER
# =============================================================================
cat("\n=== B. Exporting data for HOMER ===\n")

# -- B1. Export genome FASTA for HOMER (chr1-31 only) --
homer_genome_fa <- file.path(PROJECT_DIR, "genome", "derLaeGenome_chr1_31.fasta")
if (!file.exists(homer_genome_fa)) {
  cat("Exporting genome FASTA for HOMER (chr1-31)...\n")
  dir.create(dirname(homer_genome_fa), showWarnings = FALSE, recursive = TRUE)
  writeXStringSet(genome, homer_genome_fa)
  cat("Genome FASTA exported:", homer_genome_fa, "\n")
} else {
  cat("Genome FASTA already exists:", homer_genome_fa, "\n")
}

# -- B2. Export promoter sequences as FASTA --
prom_regions <- regions[regions$region_type == "promoter"]
prom_seqs <- getSeq(genome, prom_regions)
names(prom_seqs) <- prom_regions$gene_id

prom_fasta <- file.path(out_dir, "promoters_for_homer.fa")
writeXStringSet(prom_seqs, prom_fasta)
cat("Promoter FASTA:", length(prom_seqs), "sequences →", prom_fasta, "\n")

# -- B2. Convert JASPAR motifs to HOMER format --
cat("Loading JASPAR 2024 motifs...\n")
jaspar_sqlite <- list.files(cache_home, pattern = "JASPAR.*sqlite",
                            recursive = TRUE, full.names = TRUE)
if (length(jaspar_sqlite) == 0) {
  cat("Downloading JASPAR2024 database...\n")
  jaspar_db <- JASPAR2024()
  jaspar_sqlite <- list.files(cache_home, pattern = "JASPAR.*sqlite",
                              recursive = TRUE, full.names = TRUE)
}
jaspar_sqlite <- jaspar_sqlite[1]

all_pfms <- getMatrixSet(jaspar_sqlite, list(collection = "CORE", all_versions = FALSE))
metazoan_groups <- c("insects", "nematodes", "urochordates", "vertebrates")
is_metazoa <- sapply(all_pfms, function(x) {
  tg <- tags(x)$tax_group
  !is.null(tg) && tg %in% metazoan_groups
})
jaspar_pfms <- all_pfms[is_metazoa]
rm(all_pfms); gc(verbose = FALSE)
cat("JASPAR metazoan motifs:", length(jaspar_pfms), "\n")

# Build metadata
motif_meta <- data.frame(
  motif_id   = sapply(jaspar_pfms, ID),
  motif_name = sapply(jaspar_pfms, name),
  tf_class   = vapply(jaspar_pfms, function(x) {
    tgs <- TFBSTools::tags(x)
    if ("family" %in% names(tgs) && length(tgs$family) > 0) tgs$family[1]
    else NA_character_
  }, character(1)),
  stringsAsFactors = FALSE
)

# Convert to universalmotif and export as HOMER format
cat("Converting JASPAR to HOMER format...\n")
umotifs <- convert_motifs(jaspar_pfms)
homer_motif_file <- file.path(out_dir, "jaspar_motifs_homer.txt")
write_homer(umotifs, homer_motif_file)
cat("HOMER motif file:", homer_motif_file, "\n")


# =============================================================================
# C. HOMER DE NOVO MOTIF DISCOVERY
# =============================================================================
cat("\n=== C. HOMER de novo motif discovery ===\n")

homer_available <- Sys.which("findMotifs.pl") != ""
homer_denovo_dir <- file.path(out_dir, "homer_denovo")

if (homer_available) {
  dir.create(homer_denovo_dir, showWarnings = FALSE, recursive = TRUE)
  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))

  cat("Running HOMER de novo on promoter sequences...\n")
  cat(sprintf("  %d promoters, %d cores, motif lengths 8,10,12\n",
              length(prom_seqs), n_cores))

  homer_cmd <- sprintf(
    "findMotifs.pl %s fasta %s -len 8,10,12 -p %d -S 25 2>&1",
    prom_fasta, homer_denovo_dir, n_cores
  )
  cat("Command:", homer_cmd, "\n")
  homer_start <- proc.time()
  homer_out <- system(homer_cmd, intern = TRUE)
  homer_time <- (proc.time() - homer_start)[3]
  cat(sprintf("HOMER de novo complete in %.1f seconds\n", homer_time))

  # Check for de novo motifs
  denovo_file <- file.path(homer_denovo_dir, "homerMotifs.all.motifs")
  if (file.exists(denovo_file)) {
    denovo_motifs <- read_homer(denovo_file)
    cat(sprintf("De novo motifs discovered: %d\n", length(denovo_motifs)))
  } else {
    cat("WARNING: No de novo motif file found\n")
  }

  # Check for known motif enrichment
  known_file <- file.path(homer_denovo_dir, "knownResults.txt")
  if (file.exists(known_file)) {
    known_results <- fread(known_file, sep = "\t")
    cat(sprintf("Known motifs tested: %d\n", nrow(known_results)))
    # Save top enriched known motifs
    top_known <- known_results[1:min(50, nrow(known_results)), ]
    fwrite(top_known, file.path(out_dir, "homer_top_known_enriched.tsv"), sep = "\t")
  }
} else {
  cat("HOMER not in PATH — skipping de novo discovery.\n")
  cat("Run setup_homer.sh first to install HOMER.\n")
}

# Free promoter sequences (no longer needed)
rm(prom_seqs); gc(verbose = FALSE)


# =============================================================================
# D. JASPAR MOTIF SCANNING ON EXTENDED REGIONS
# =============================================================================
cat("\n=== D. JASPAR motif scanning on extended regions ===\n")
cat(sprintf("  %d motifs x %d regions\n", length(jaspar_pfms), length(regions)))

# -- D0. Diagnostic: verify matchMotifs works --
cat("\n  [DIAGNOSTIC] Testing matchMotifs on 10 promoters x 5 motifs...\n")
test_reg <- regions[regions$region_type == "promoter"][1:10]
test_seqs <- getSeq(genome, test_reg)
names(test_seqs) <- test_reg$gene_id
test_result <- matchMotifs(jaspar_pfms[1:5], test_seqs, out = "positions", p.cutoff = 5e-5)
cat(sprintf("  [DIAGNOSTIC] Result class: %s, length: %d\n", class(test_result)[1], length(test_result)))
cat(sprintf("  [DIAGNOSTIC] Element class: %s\n", class(test_result[[1]])[1]))
test_nhits <- sum(sapply(seq_along(test_result), function(i) {
  x <- test_result[[i]]
  if (is(x, "IRangesList")) sum(lengths(x))
  else if (is(x, "GRanges") || is(x, "GRangesList")) length(unlist(x))
  else { cat(sprintf("  [DIAGNOSTIC] Unexpected class: %s\n", class(x)[1])); 0 }
}))
cat(sprintf("  [DIAGNOSTIC] Test hits: %d (expect > 0)\n", test_nhits))
if (test_nhits == 0) {
  cat("  [DIAGNOSTIC] WARNING: 0 test hits! Trying p.cutoff=1e-4...\n")
  test_result2 <- matchMotifs(jaspar_pfms[1:5], test_seqs, out = "positions", p.cutoff = 1e-4)
  test_nhits2 <- sum(sapply(seq_along(test_result2), function(i) {
    x <- test_result2[[i]]
    if (is(x, "IRangesList")) sum(lengths(x)) else length(unlist(x))
  }))
  cat(sprintf("  [DIAGNOSTIC] p=1e-4 hits: %d\n", test_nhits2))
}
rm(test_reg, test_seqs, test_result); gc(verbose = FALSE)

# -- Helper: extract hits handling both IRangesList and GRanges --
extract_hits <- function(pos_result, motif_idx, reg_lookup, rtype, motif_meta) {
  x <- pos_result[[motif_idx]]
  if (is(x, "IRangesList")) {
    if (sum(lengths(x)) == 0) return(NULL)
    ir <- unlist(x)
    hit_genes <- names(ir)
    rel_s <- start(ir); rel_e <- end(ir)
  } else if (is(x, "GRanges")) {
    if (length(x) == 0) return(NULL)
    hit_genes <- as.character(seqnames(x))
    rel_s <- start(x); rel_e <- end(x)
  } else if (is(x, "GRangesList")) {
    x <- unlist(x)
    if (length(x) == 0) return(NULL)
    hit_genes <- names(x)
    rel_s <- start(x); rel_e <- end(x)
  } else { return(NULL) }

  ridx <- match(hit_genes, reg_lookup$gene_id)
  valid <- !is.na(ridx)
  if (!any(valid)) return(NULL)
  ridx <- ridx[valid]; rel_s <- rel_s[valid]; rel_e <- rel_e[valid]
  hit_genes <- hit_genes[valid]
  rs <- reg_lookup$reg_strand[ridx]
  ps <- reg_lookup$reg_start[ridx]; pe <- reg_lookup$reg_end[ridx]
  data.frame(
    chr = reg_lookup$reg_chr[ridx],
    start = ifelse(rs == "+", ps + rel_s - 1L, pe - rel_e + 1L),
    end   = ifelse(rs == "+", ps + rel_e - 1L, pe - rel_s + 1L),
    strand = ifelse(rs == "+", "+", "-"),
    motif_id = motif_meta$motif_id[motif_idx],
    motif_name = motif_meta$motif_name[motif_idx],
    tf_class = motif_meta$tf_class[motif_idx],
    gene_id = hit_genes, region_type = rtype,
    score = NA_real_, source = "JASPAR", stringsAsFactors = FALSE)
}

# -- D1. Scan per region type (all chromosomes batched in one call) --
region_types <- c("promoter", "upstream_distal", "gene_body", "downstream")
all_chunk_files <- c()
total_hits <- 0L
scan_t0 <- proc.time()[3]

for (rtype in region_types) {
  type_regions <- regions[regions$region_type == rtype]
  n_reg <- length(type_regions)
  cat(sprintf("\n  [%s] Extracting %d sequences...\n", rtype, n_reg))
  flush.console()

  t0 <- proc.time()[3]
  type_seqs <- getSeq(genome, type_regions)
  names(type_seqs) <- type_regions$gene_id
  cat(sprintf("  [%s] Sequences extracted (%.0fs), total bp: %s\n",
              rtype, proc.time()[3] - t0, format(sum(width(type_seqs)), big.mark = ",")))
  flush.console()

  reg_lookup <- data.frame(
    gene_id = type_regions$gene_id,
    reg_chr = as.character(seqnames(type_regions)),
    reg_start = start(type_regions), reg_end = end(type_regions),
    reg_strand = as.character(strand(type_regions)), stringsAsFactors = FALSE)

  cat(sprintf("  [%s] Running matchMotifs (%d motifs x %d seqs)...\n",
              rtype, length(jaspar_pfms), n_reg))
  flush.console()
  t0 <- proc.time()[3]

  pos_result <- tryCatch(
    matchMotifs(jaspar_pfms, type_seqs, out = "positions", p.cutoff = 5e-5),
    error = function(e) {
      cat(sprintf("  [%s] ERROR in matchMotifs: %s\n", rtype, e$message))
      NULL
    })

  scan_time <- proc.time()[3] - t0
  cat(sprintf("  [%s] matchMotifs done (%.0fs = %.1f min)\n", rtype, scan_time, scan_time/60))
  flush.console()

  if (is.null(pos_result)) {
    cat(sprintf("  [%s] SKIPPED due to error\n", rtype))
    rm(type_seqs); gc(verbose = FALSE)
    next
  }

  cat(sprintf("  [%s] Extracting hits from %d motifs...\n", rtype, length(pos_result)))
  flush.console()
  chunk_hits <- vector("list", length(pos_result))
  for (i in seq_along(pos_result)) {
    chunk_hits[[i]] <- extract_hits(pos_result, i, reg_lookup, rtype, motif_meta)
    if (i %% 200 == 0) {
      n_so_far <- sum(sapply(chunk_hits[1:i], function(x) if(is.null(x)) 0L else nrow(x)))
      cat(sprintf("    [%s] %d/%d motifs, %s hits so far\n",
                  rtype, i, length(pos_result), format(n_so_far, big.mark = ",")))
      flush.console()
    }
  }

  dt <- rbindlist(chunk_hits, fill = TRUE)
  rtype_hits <- nrow(dt)
  rm(chunk_hits, pos_result, type_seqs); gc(verbose = FALSE)

  if (rtype_hits > 0) {
    chunk_file <- file.path(out_dir, sprintf("_chunk_%s.tsv", rtype))
    fwrite(dt, chunk_file, sep = "\t")
    all_chunk_files <- c(all_chunk_files, chunk_file)
    total_hits <- total_hits + rtype_hits
  }
  rm(dt); gc(verbose = FALSE)

  elapsed <- proc.time()[3] - scan_t0
  cat(sprintf("  [%s] DONE: %s hits | Running total: %s | Elapsed: %.1f min\n",
              rtype, format(rtype_hits, big.mark = ","),
              format(total_hits, big.mark = ","), elapsed/60))
  flush.console()
}

cat(sprintf("\nTotal JASPAR scanning: %s hits in %.1f min\n",
            format(total_hits, big.mark = ","), (proc.time()[3] - scan_t0)/60))


# =============================================================================
# E. COMBINE AND ANNOTATE ALL HITS
# =============================================================================
cat("\n=== E. Combining and annotating all hits ===\n")

# -- E1. Read all JASPAR chunks --
cat("Reading", length(all_chunk_files), "chunk files...\n")
hits_list <- lapply(all_chunk_files, fread, sep = "\t")
hits_dt <- rbindlist(hits_list, fill = TRUE)
rm(hits_list); gc(verbose = FALSE)
cat(sprintf("Combined JASPAR hits: %s\n", format(nrow(hits_dt), big.mark = ",")))

# -- E2. Add distance to TSS --
cat("Computing distance to TSS...\n")
tss_idx <- match(hits_dt$gene_id, tss_lookup$gene_id)
hit_tss <- tss_lookup$tss[tss_idx]
hit_str <- tss_lookup$strand[tss_idx]

# Midpoint of motif hit
hit_mid <- (hits_dt$start + hits_dt$end) %/% 2L

# Signed distance: negative = upstream, positive = downstream
# For + strand: distance = hit_mid - TSS
# For - strand: distance = TSS - hit_mid (reversed)
hits_dt[, dist_to_tss := ifelse(hit_str == "+",
                                hit_mid - hit_tss,
                                hit_tss - hit_mid)]
rm(hit_tss, hit_str, hit_mid, tss_idx); gc(verbose = FALSE)

# -- E3. Clean up temp files --
cat("Cleaning up chunk files...\n")
file.remove(all_chunk_files)


# =============================================================================
# F. SAVE COMPREHENSIVE OUTPUT
# =============================================================================
cat("\n=== F. Saving output ===\n")

# -- F1. Full annotated hit table --
out_hits <- file.path(out_dir, "motif_hits_extended.tsv.gz")
fwrite(hits_dt, out_hits, sep = "\t", compress = "gzip")
cat(sprintf("Annotated hits: %s (%s rows, %.1f MB)\n",
            out_hits, format(nrow(hits_dt), big.mark = ","),
            file.size(out_hits) / 1e6))

# -- F2. Per-motif per-region summary --
cat("Computing per-motif summaries...\n")
motif_region_summary <- hits_dt[, .(
  n_hits      = .N,
  n_genes     = uniqueN(gene_id),
  mean_score  = mean(score, na.rm = TRUE),
  mean_dist   = mean(dist_to_tss, na.rm = TRUE)
), by = .(motif_id, motif_name, tf_class, region_type, source)]

out_summary <- file.path(out_dir, "motif_summary_by_region.tsv")
fwrite(motif_region_summary, out_summary, sep = "\t")
cat(sprintf("Summary by region: %d rows → %s\n",
            nrow(motif_region_summary), out_summary))

# -- F3. Region-level overview --
region_overview <- hits_dt[, .(
  total_hits     = .N,
  unique_motifs  = uniqueN(motif_id),
  unique_genes   = uniqueN(gene_id),
  mean_hits_per_gene = .N / uniqueN(gene_id)
), by = .(region_type)]
region_overview <- region_overview[match(
  c("promoter", "upstream_distal", "gene_body", "downstream"),
  region_overview$region_type
)]

out_overview <- file.path(out_dir, "motif_region_overview.tsv")
fwrite(region_overview, out_overview, sep = "\t")
cat("\nRegion overview:\n")
print(region_overview)

# -- F4. Motif metadata (for downstream use) --
out_meta <- file.path(out_dir, "motif_metadata.tsv")
fwrite(motif_meta, out_meta, sep = "\t")


# =============================================================================
# G. VERIFICATION
# =============================================================================
cat("\n=== G. Verification ===\n")
cat(sprintf("  Total hits:        %s\n", format(nrow(hits_dt), big.mark = ",")))
cat(sprintf("  Unique motifs:     %d / %d\n",
            uniqueN(hits_dt$motif_id), nrow(motif_meta)))
cat(sprintf("  Unique genes:      %d / %d\n",
            uniqueN(hits_dt$gene_id), length(genes)))
cat(sprintf("  Region types:      %s\n",
            paste(unique(hits_dt$region_type), collapse = ", ")))
cat(sprintf("  Output files:      %s\n", out_dir))

# Spot check: promoter hits should have negative dist_to_tss
prom_hits <- hits_dt[region_type == "promoter"]
if (nrow(prom_hits) > 0) {
  pct_neg <- 100 * sum(prom_hits$dist_to_tss < 0, na.rm = TRUE) / nrow(prom_hits)
  cat(sprintf("  Promoter sanity:   %.1f%% have negative dist_to_tss (expected ~100%%)\n",
              pct_neg))
}

# Memory report
mem_info <- gc()
cat(sprintf("  Peak memory:       %.1f MB\n", sum(mem_info[, 6])))

rm(hits_dt, motif_region_summary, region_overview); gc(verbose = FALSE)


# =============================================================================
# H. HOMER DE NOVO RESULTS SUMMARY
# =============================================================================
if (homer_available && dir.exists(homer_denovo_dir)) {
  cat("\n=== H. HOMER de novo summary ===\n")

  denovo_file <- file.path(homer_denovo_dir, "homerMotifs.all.motifs")
  if (file.exists(denovo_file)) {
    denovo_motifs <- read_homer(denovo_file)
    cat(sprintf("  De novo motifs: %d\n", length(denovo_motifs)))

    # Compare de novo motifs to JASPAR
    cat("  Comparing de novo motifs to JASPAR database...\n")
    tryCatch({
      # Compare each de novo motif to all JASPAR motifs
      for (i in seq_along(denovo_motifs)) {
        sim <- compare_motifs(c(denovo_motifs[i], umotifs),
                              method = "PCC", min.overlap = 6)
        # Best match (excluding self-comparison)
        best_sim <- sort(sim[1, -1], decreasing = TRUE)[1]
        best_name <- names(best_sim)
        cat(sprintf("    Motif %d: best JASPAR match = %s (PCC = %.3f)\n",
                    i, best_name, best_sim))
      }
    }, error = function(e) {
      cat("  WARNING: Motif comparison failed:", e$message, "\n")
    })
  }
}


cat("\n=== Cluster job complete ===\n")
cat("Output files in:", out_dir, "\n")
cat("Copy results back to local:\n")
cat("  cluster/results/motif_hits_extended.tsv.gz     → results/batch1.5/\n")
cat("  cluster/results/motif_summary_by_region.tsv    → results/batch1.5/\n")
cat("  cluster/results/motif_region_overview.tsv       → results/batch1.5/\n")
cat("  cluster/results/motif_metadata.tsv              → results/batch1.5/\n")
cat("  cluster/results/homer_denovo/                   → results/batch1.5/homer_denovo/\n")
cat("  cluster/results/homer_top_known_enriched.tsv    → results/batch1.5/\n")
