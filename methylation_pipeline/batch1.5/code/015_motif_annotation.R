#!/usr/bin/env Rscript
# =============================================================================
# Batch 1.5: TFBS Motif Annotation — JASPAR 2024 + HOMER de novo
# Scans 10kb upstream + gene body + 10kb downstream per gene for TFBS motifs.
# Output: data/ (motif hits, summaries) — used by batch09
# Requires: ~64 GB RAM, HOMER optional
# =============================================================================

source("methylation_pipeline/_config.R")

cache_home <- Sys.getenv("XDG_CACHE_HOME",
  if (ON_CLUSTER) file.path(PROJECT_DIR, "cluster/.cache") else "~/.cache")
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

BATCH_DIR <- file.path(PIPE_DIR, "batch1.5")

cat("=== Batch 1.5: TFBS Motif Annotation ===\n\n")

# --- A. Load genome + GFF, define extended regions ---
genome <- load_genome()
chr_lengths <- setNames(width(genome), names(genome))
cat(sprintf("Genome: %d chromosomes, %s bp total\n",
            length(genome), format(sum(chr_lengths), big.mark = ",")))

gff <- load_gff()
genes <- gff[gff$type == "gene"]
cat("Genes:", length(genes), "\n")

# Extended regions (10kb flanks + gene body)
regions_rds <- CACHE$extended
if (file.exists(regions_rds)) {
  regions <- readRDS(regions_rds)
} else {
  gene_strand <- as.character(strand(genes))
  gene_chr    <- as.character(seqnames(genes))
  gene_ids    <- genes$ID
  chrlen      <- as.integer(chr_lengths[gene_chr])

  tss <- ifelse(gene_strand == "+", start(genes), end(genes))
  tes <- ifelse(gene_strand == "+", end(genes), start(genes))

  build_gr <- function(s, e, type) {
    valid <- !is.na(s) & !is.na(e) & s <= e & s >= 1L
    GRanges(seqnames = gene_chr[valid], ranges = IRanges(start = s[valid], end = e[valid]),
            strand = gene_strand[valid], gene_id = gene_ids[valid], region_type = type)
  }

  prom_s <- ifelse(gene_strand == "+", pmax(1L, tss - 2000L), tss + 1L)
  prom_e <- ifelse(gene_strand == "+", tss - 1L, pmin(chrlen, tss + 2000L))
  up_s   <- ifelse(gene_strand == "+", pmax(1L, tss - 10000L), tss + 2001L)
  up_e   <- ifelse(gene_strand == "+", pmax(1L, tss - 2001L), pmin(chrlen, tss + 10000L))
  gb_s   <- start(genes); gb_e <- end(genes)
  dn_s   <- ifelse(gene_strand == "+", tes + 1L, pmax(1L, tes - 10000L))
  dn_e   <- ifelse(gene_strand == "+", pmin(chrlen, tes + 10000L), tes - 1L)

  regions <- c(build_gr(prom_s, prom_e, "promoter"),
               build_gr(up_s, up_e, "upstream_distal"),
               build_gr(gb_s, gb_e, "gene_body"),
               build_gr(dn_s, dn_e, "downstream"))
  seqlevels(regions) <- keep_chr
  seqlengths(regions) <- chr_lengths[keep_chr]
  regions <- trim(regions)
  saveRDS(regions, regions_rds)
}

for (rt in c("promoter", "upstream_distal", "gene_body", "downstream")) {
  cat(sprintf("  %-18s %6d regions\n", rt, sum(regions$region_type == rt)))
}

# TSS lookup
gene_strand <- as.character(strand(genes))
tss_vals    <- ifelse(gene_strand == "+", start(genes), end(genes))
tss_lookup  <- data.frame(gene_id = genes$ID, chr = as.character(seqnames(genes)),
                           tss = tss_vals, strand = gene_strand, stringsAsFactors = FALSE)

# --- B. JASPAR 2024 motifs ---
cat("\nLoading JASPAR 2024 motifs...\n")
jaspar_sqlite <- list.files(cache_home, pattern = "JASPAR.*sqlite", recursive = TRUE, full.names = TRUE)
if (length(jaspar_sqlite) == 0) {
  jaspar_db <- JASPAR2024()
  jaspar_sqlite <- list.files(cache_home, pattern = "JASPAR.*sqlite", recursive = TRUE, full.names = TRUE)
}

all_pfms <- getMatrixSet(jaspar_sqlite[1], list(collection = "CORE", all_versions = FALSE))
metazoan_groups <- c("insects", "nematodes", "urochordates", "vertebrates")
is_metazoa <- sapply(all_pfms, function(x) {
  tg <- tags(x)$tax_group; !is.null(tg) && tg %in% metazoan_groups
})
jaspar_pfms <- all_pfms[is_metazoa]
rm(all_pfms); gc(verbose = FALSE)
cat("JASPAR metazoan motifs:", length(jaspar_pfms), "\n")

motif_meta <- data.frame(
  motif_id   = sapply(jaspar_pfms, ID),
  motif_name = sapply(jaspar_pfms, name),
  tf_class   = vapply(jaspar_pfms, function(x) {
    tgs <- TFBSTools::tags(x)
    if ("family" %in% names(tgs) && length(tgs$family) > 0) tgs$family[1] else NA_character_
  }, character(1)),
  stringsAsFactors = FALSE)

# --- C. HOMER de novo (optional) ---
cat("\nHOMER de novo motif discovery...\n")
homer_available <- Sys.which("findMotifs.pl") != ""
data_dir <- file.path(BATCH_DIR, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

prom_regions <- regions[regions$region_type == "promoter"]
prom_seqs <- getSeq(genome, prom_regions)
names(prom_seqs) <- prom_regions$gene_id
prom_fasta <- file.path(data_dir, "promoters_for_homer.fa")
writeXStringSet(prom_seqs, prom_fasta)

homer_denovo_dir <- file.path(data_dir, "homer_denovo")
denovo_file <- file.path(homer_denovo_dir, "homerMotifs.all.motifs")
if (file.exists(denovo_file)) {
  cat("  HOMER de novo results already exist, skipping...\n")
  denovo_motifs <- read_homer(denovo_file)
  cat(sprintf("  De novo motifs: %d\n", length(denovo_motifs)))
} else if (homer_available) {
  dir.create(homer_denovo_dir, showWarnings = FALSE, recursive = TRUE)
  n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "4"))
  homer_cmd <- sprintf("findMotifs.pl %s fasta %s -len 8,10,12 -p %d -S 25 2>&1",
                        prom_fasta, homer_denovo_dir, n_cores)
  cat(sprintf("  %d promoters, %d cores\n", length(prom_seqs), n_cores))
  system(homer_cmd, intern = TRUE)

  if (file.exists(denovo_file)) {
    denovo_motifs <- read_homer(denovo_file)
    cat(sprintf("  De novo motifs: %d\n", length(denovo_motifs)))
  }
  known_file <- file.path(homer_denovo_dir, "knownResults.txt")
  if (file.exists(known_file)) {
    top_known <- fread(known_file, sep = "\t")[1:min(50, .N)]
    fwrite(top_known, file.path(data_dir, "homer_top_known_enriched.tsv"), sep = "\t")
  }
} else {
  cat("  HOMER not available, skipping de novo.\n")
}
rm(prom_seqs); gc(verbose = FALSE)

# --- D. JASPAR motif scanning (per region type, batched) ---
# Skip if output already exists (scanning takes hours; de novo figs are fast)
existing_hits <- file.path(data_dir, "motif_hits_extended.tsv.gz")
if (file.exists(existing_hits)) {
  cat("\nJASPAR scan output exists, skipping to de novo characterization...\n")
  cat(sprintf("  Existing: %s (%.1f MB)\n", existing_hits, file.size(existing_hits) / 1e6))

  # Still need motif_meta for downstream — already loaded above
  # Load summary data for the final print
  hits_dt <- fread(existing_hits)
  motif_summary <- hits_dt[, .(n_hits = .N, n_genes = uniqueN(gene_id)),
    by = .(motif_id, motif_name, tf_class, region_type)]
  region_overview <- hits_dt[, .(total_hits = .N, unique_motifs = uniqueN(motif_id),
    unique_genes = uniqueN(gene_id)), by = .(region_type)]
  cat("\nRegion overview:\n"); print(region_overview)
  rm(hits_dt, motif_summary, region_overview); gc(verbose = FALSE)

} else {

cat("\nJASPAR motif scanning on extended regions...\n")
cat(sprintf("  %d motifs x %d regions\n", length(jaspar_pfms), length(regions)))

# -- D0. Quick diagnostic: verify matchMotifs works --
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
  cat("  [DIAGNOSTIC] WARNING: 0 test hits! Trying with p.cutoff=1e-4...\n")
  test_result2 <- matchMotifs(jaspar_pfms[1:5], test_seqs, out = "positions", p.cutoff = 1e-4)
  test_nhits2 <- sum(sapply(seq_along(test_result2), function(i) {
    x <- test_result2[[i]]
    if (is(x, "IRangesList")) sum(lengths(x)) else length(unlist(x))
  }))
  cat(sprintf("  [DIAGNOSTIC] p=1e-4 hits: %d\n", test_nhits2))
}
rm(test_reg, test_seqs, test_result); gc(verbose = FALSE)

# -- D1. Detect result structure and build extractor --
# matchMotifs with DNAStringSet returns list of IRangesList or GRangesList depending on version
extract_hits <- function(pos_result, motif_idx, reg_lookup, rtype, motif_meta) {
  x <- pos_result[[motif_idx]]

  # Handle both IRangesList and GRanges/GRangesList
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
  } else {
    return(NULL)
  }

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
    score = NA_real_, stringsAsFactors = FALSE
  )
}

# -- D2. Scan per region type (all chromosomes batched) --
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
    reg_strand = as.character(strand(type_regions)),
    stringsAsFactors = FALSE)

  # Strategy: chunk SEQUENCES (not motifs) and use out="matches" (boolean, fast).
  # matchMotifs with out="positions" hangs even with 50 motifs on 25k seqs.
  # out="matches" returns a sparse boolean matrix (motif present yes/no per gene).
  # For batch09 we need gene_id x motif_id presence — positions not required.
  SEQ_CHUNK <- 2000L
  n_motifs <- length(jaspar_pfms)
  n_seq_chunks <- ceiling(n_reg / SEQ_CHUNK)
  cat(sprintf("  [%s] Scanning %d motifs x %d seqs in %d seq-chunks of %d (out='matches')...\n",
              rtype, n_motifs, n_reg, n_seq_chunks, SEQ_CHUNK))
  flush.console()
  t0 <- proc.time()[3]

  all_chunk_hits <- list()
  for (si in seq(1, n_reg, by = SEQ_CHUNK)) {
    sj <- min(si + SEQ_CHUNK - 1L, n_reg)
    cat(sprintf("    [%s] seqs %d-%d / %d ...", rtype, si, sj, n_reg))
    flush.console()
    ct0 <- proc.time()[3]

    chunk_seqs <- type_seqs[si:sj]
    chunk_lookup <- reg_lookup[si:sj, ]

    match_result <- tryCatch(
      matchMotifs(jaspar_pfms, chunk_seqs, out = "matches", p.cutoff = 5e-5),
      error = function(e) {
        cat(sprintf(" ERROR: %s\n", e$message))
        NULL
      }
    )

    if (!is.null(match_result)) {
      # motifMatches() returns a logical matrix: seqs x motifs
      mat <- motifMatches(match_result)
      if (any(mat)) {
        hits_idx <- which(mat, arr.ind = TRUE)  # row=seq, col=motif
        hit_df <- data.frame(
          chr = chunk_lookup$reg_chr[hits_idx[, 1]],
          start = chunk_lookup$reg_start[hits_idx[, 1]],
          end = chunk_lookup$reg_end[hits_idx[, 1]],
          strand = chunk_lookup$reg_strand[hits_idx[, 1]],
          motif_id = motif_meta$motif_id[hits_idx[, 2]],
          motif_name = motif_meta$motif_name[hits_idx[, 2]],
          tf_class = motif_meta$tf_class[hits_idx[, 2]],
          gene_id = chunk_lookup$gene_id[hits_idx[, 1]],
          region_type = rtype,
          score = NA_real_,
          stringsAsFactors = FALSE
        )
        all_chunk_hits[[length(all_chunk_hits) + 1L]] <- hit_df
        rm(hit_df, hits_idx)
      }
      rm(mat, match_result)
    }
    rm(chunk_seqs, chunk_lookup)

    ct_elapsed <- proc.time()[3] - ct0
    nhits_so_far <- sum(sapply(all_chunk_hits, nrow))
    cat(sprintf(" %.0fs, %s hits total\n", ct_elapsed, format(nhits_so_far, big.mark = ",")))
    flush.console()
    gc(verbose = FALSE)
  }

  scan_time <- proc.time()[3] - t0
  cat(sprintf("  [%s] All seq chunks done (%.0fs = %.1f min)\n", rtype, scan_time, scan_time/60))
  flush.console()

  dt <- rbindlist(all_chunk_hits, fill = TRUE)
  rtype_hits <- if (is.null(dt) || nrow(dt) == 0) 0L else nrow(dt)
  rm(all_chunk_hits, type_seqs); gc(verbose = FALSE)

  if (rtype_hits > 0) {
    chunk_file <- file.path(data_dir, sprintf("_chunk_%s.tsv", rtype))
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

cat(sprintf("\nAll scanning done: %s total hits in %.1f min\n",
            format(total_hits, big.mark = ","), (proc.time()[3] - scan_t0)/60))

# --- E. Combine and annotate ---
cat("\nCombining all hits...\n")
hits_list <- lapply(all_chunk_files, fread, sep = "\t")
hits_dt <- rbindlist(hits_list, fill = TRUE)
rm(hits_list); gc(verbose = FALSE)
cat(sprintf("Total JASPAR hits: %s\n", format(nrow(hits_dt), big.mark = ",")))

# Distance to TSS
tss_idx <- match(hits_dt$gene_id, tss_lookup$gene_id)
hit_tss <- tss_lookup$tss[tss_idx]
hit_str <- tss_lookup$strand[tss_idx]
hit_mid <- (hits_dt$start + hits_dt$end) %/% 2L
hits_dt[, dist_to_tss := ifelse(hit_str == "+", hit_mid - hit_tss, hit_tss - hit_mid)]
rm(hit_tss, hit_str, hit_mid, tss_idx); gc(verbose = FALSE)

file.remove(all_chunk_files)

# --- F. Save output ---
out_hits <- file.path(data_dir, "motif_hits_extended.tsv.gz")
fwrite(hits_dt, out_hits, sep = "\t", compress = "gzip")
cat(sprintf("Saved: %s (%s rows)\n", out_hits, format(nrow(hits_dt), big.mark = ",")))

motif_summary <- hits_dt[, .(n_hits = .N, n_genes = uniqueN(gene_id),
  mean_score = mean(score, na.rm = TRUE), mean_dist = mean(dist_to_tss, na.rm = TRUE)),
  by = .(motif_id, motif_name, tf_class, region_type)]
save_data(motif_summary, BATCH_DIR, "motif_summary_by_region")

region_overview <- hits_dt[, .(total_hits = .N, unique_motifs = uniqueN(motif_id),
  unique_genes = uniqueN(gene_id), mean_hits_per_gene = .N / uniqueN(gene_id)),
  by = .(region_type)]
save_data(region_overview, BATCH_DIR, "motif_region_overview")
save_data(motif_meta, BATCH_DIR, "motif_metadata")

cat("\nRegion overview:\n"); print(region_overview)
cat(sprintf("\nUnique motifs: %d / %d\n", uniqueN(hits_dt$motif_id), nrow(motif_meta)))
cat(sprintf("Unique genes:  %d / %d\n", uniqueN(hits_dt$gene_id), length(genes)))

rm(hits_dt, motif_summary, region_overview); gc(verbose = FALSE)

} # end else (JASPAR scan)

# =============================================================================
# G. DE NOVO MOTIF CHARACTERIZATION FIGURES
# =============================================================================
cat("\n=== De novo motif characterization ===\n")

library(ggplot2)
library(ggseqlogo)

denovo_file <- file.path(BATCH_DIR, "data/homer_denovo/homerMotifs.all.motifs")
nonred_file <- file.path(BATCH_DIR, "data/homer_denovo/nonRedundant.motifs")

# Use non-redundant set for primary analysis
motif_file <- if (file.exists(nonred_file)) nonred_file else denovo_file

if (file.exists(motif_file)) {
  denovo <- read_homer(motif_file)
  cat(sprintf("De novo motifs loaded: %d\n", length(denovo)))

  # Extract consensus + enrichment stats from HOMER header
  denovo_info <- data.table(
    idx = seq_along(denovo),
    consensus = vapply(denovo, function(m) m@consensus, character(1)),
    name = vapply(denovo, function(m) m@name, character(1)),
    width = vapply(denovo, function(m) ncol(m@motif), integer(1))
  )

  # Parse enrichment from name field (HOMER encodes stats in name)
  # Best guess match is in altname for nonredundant
  denovo_info[, best_guess := vapply(denovo, function(m) {
    nm <- m@altname
    if (is.null(nm) || length(nm) == 0 || is.na(nm)) return(NA_character_)
    # Extract BestGuess: field
    bg <- regmatches(nm, regexpr("BestGuess:[^,)]+", nm))
    if (length(bg) > 0) sub("BestGuess:", "", bg) else NA_character_
  }, character(1))]

  # Parse target/background percentages from extrainfo
  denovo_info[, target_pct := vapply(denovo, function(m) {
    ei <- m@extrainfo
    if (length(ei) == 0) return(NA_real_)
    tp <- regmatches(paste(ei, collapse = " "), regexpr("T:[0-9.]+\\([0-9.]+%\\)", paste(ei, collapse = " ")))
    if (length(tp) > 0) {
      as.numeric(sub(".*\\(([0-9.]+)%\\).*", "\\1", tp))
    } else NA_real_
  }, numeric(1))]

  denovo_info[, bg_pct := vapply(denovo, function(m) {
    ei <- m@extrainfo
    if (length(ei) == 0) return(NA_real_)
    bp <- regmatches(paste(ei, collapse = " "), regexpr("B:[0-9.]+\\([0-9.]+%\\)", paste(ei, collapse = " ")))
    if (length(bp) > 0) {
      as.numeric(sub(".*\\(([0-9.]+)%\\).*", "\\1", bp))
    } else NA_real_
  }, numeric(1))]

  denovo_info[, enrichment := target_pct / pmax(bg_pct, 0.01)]

  save_data(denovo_info, BATCH_DIR, "denovo_motif_summary")

  # --- Fig 1.5A: Sequence logos of top de novo motifs ---
  cat("Fig 1.5A: De novo motif sequence logos...\n")

  # Get PPMs for ggseqlogo
  n_show <- min(12, length(denovo))
  ppm_list <- lapply(seq_len(n_show), function(i) {
    m <- denovo[[i]]@motif
    # Ensure it's a proper PPM (rows = ACGT, cols = positions)
    if (nrow(m) == 4) t(m) else m
  })
  names(ppm_list) <- denovo_info$consensus[1:n_show]

  # Create multi-panel logo plot
  p_logos <- ggseqlogo(ppm_list, method = "bits", ncol = 3) +
    theme(strip.text = element_text(size = 8, face = "bold")) +
    ggtitle(sprintf("Top %d de novo motifs discovered in D. laeve promoters", n_show))
  save_fig(p_logos, BATCH_DIR, "fig15a_denovo_sequence_logos", w = 14, h = ceiling(n_show / 3) * 3)

  # --- Fig 1.5B: Enrichment bar plot ---
  cat("Fig 1.5B: Motif enrichment barplot...\n")

  plot_b <- denovo_info[!is.na(target_pct)][order(-enrichment)][1:min(20, .N)]
  plot_b[, label := ifelse(!is.na(best_guess) & best_guess != "",
                           paste0(consensus, "\n(", best_guess, ")"),
                           consensus)]
  plot_b[, label := factor(label, levels = rev(label))]

  p15b <- ggplot(plot_b, aes(x = label, y = enrichment)) +
    geom_col(aes(fill = enrichment), width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%→%.1f%%", bg_pct, target_pct)),
              hjust = -0.1, size = 3) +
    scale_fill_viridis_c(option = "magma", guide = "none") +
    coord_flip() +
    labs(x = NULL, y = "Enrichment fold (target / background)",
         title = "De novo motif enrichment in D. laeve promoters",
         subtitle = "HOMER fasta mode | Background = shuffled promoter sequences") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))
  save_fig(p15b, BATCH_DIR, "fig15b_denovo_enrichment_bar", w = 11, h = 8)

  # --- Fig 1.5C: Motif similarity dendrogram (de novo vs JASPAR best matches) ---
  cat("Fig 1.5C: Motif similarity dendrogram...\n")

  # Compare de novo motifs to each other using universalmotif
  if (length(denovo) >= 3) {
    # Convert to universalmotif list
    denovo_um <- denovo[1:min(20, length(denovo))]
    for (i in seq_along(denovo_um)) {
      denovo_um[[i]]@name <- denovo_info$consensus[i]
    }

    # Pairwise comparison
    comp <- compare_motifs(denovo_um, method = "PCC", min.overlap = 4)

    # Convert similarity to distance
    dist_mat <- as.dist(1 - comp)
    hc <- hclust(dist_mat, method = "average")

    # Plot dendrogram
    png(file.path(BATCH_DIR, "figures/fig15c_denovo_similarity_dendrogram.png"),
        width = 10, height = 6, units = "in", res = 300)
    par(mar = c(8, 4, 3, 1))
    plot(hc, hang = -1, main = "De novo motif similarity (1 - PCC)",
         xlab = "", sub = "", cex = 0.8)
    dev.off()
    cairo_pdf(file.path(BATCH_DIR, "figures/fig15c_denovo_similarity_dendrogram.pdf"),
              width = 10, height = 6)
    par(mar = c(8, 4, 3, 1))
    plot(hc, hang = -1, main = "De novo motif similarity (1 - PCC)",
         xlab = "", sub = "", cex = 0.8)
    dev.off()
    cat("  Saved fig15c\n")
  }

  # --- Fig 1.5D: Positional metaplot (motif position relative to TSS) ---
  cat("Fig 1.5D: Positional metaplot (de novo motifs relative to TSS)...\n")

  # Scan promoter sequences for top de novo motifs to get positions
  prom_regions <- regions[regions$region_type == "promoter"]
  prom_seqs <- getSeq(genome, prom_regions)
  names(prom_seqs) <- prom_regions$gene_id

  # Convert top de novo motifs to PWMs for scanning
  n_scan <- min(6, length(denovo))
  denovo_pwm <- lapply(seq_len(n_scan), function(i) {
    m <- denovo[[i]]
    # Convert to TFBSTools PWMatrix for motifmatchr
    pfm_mat <- m@motif
    if (nrow(pfm_mat) != 4) pfm_mat <- t(pfm_mat)
    rownames(pfm_mat) <- c("A", "C", "G", "T")
    # Scale to counts (motifmatchr needs PFMatrix)
    count_mat <- round(pfm_mat * 100)
    PFMatrix(ID = paste0("denovo_", i), name = denovo_info$consensus[i],
             profileMatrix = count_mat)
  })
  denovo_pfm_list <- do.call(PFMatrixList, denovo_pwm)

  # Scan with positions
  pos_result <- tryCatch(
    matchMotifs(denovo_pfm_list, prom_seqs, out = "positions", p.cutoff = 5e-5),
    error = function(e) { cat("  Position scan error:", e$message, "\n"); NULL }
  )

  if (!is.null(pos_result)) {
    # Build position table relative to TSS
    pos_list <- list()
    prom_strand <- as.character(strand(prom_regions))

    for (mi in seq_len(n_scan)) {
      x <- pos_result[[mi]]
      if (is(x, "IRangesList")) {
        ir <- unlist(x)
        if (length(ir) == 0) next
        hit_genes <- names(ir)
        rel_pos <- (start(ir) + end(ir)) / 2
      } else if (is(x, "GRanges") || is(x, "GRangesList")) {
        x <- unlist(x)
        if (length(x) == 0) next
        hit_genes <- if (!is.null(names(x))) names(x) else as.character(seqnames(x))
        rel_pos <- (start(x) + end(x)) / 2
      } else next

      # Convert relative position in 2kb promoter to distance from TSS
      # Promoter regions are 2kb upstream: position 1 = -2000, position 2000 = -1
      gene_idx <- match(hit_genes, prom_regions$gene_id)
      valid <- !is.na(gene_idx)
      if (!any(valid)) next

      gs <- prom_strand[gene_idx[valid]]
      dist_tss <- ifelse(gs == "+", rel_pos[valid] - 2000, -(rel_pos[valid] - 1))

      pos_list[[mi]] <- data.table(
        motif = denovo_info$consensus[mi],
        dist_to_tss = dist_tss
      )
    }

    if (length(pos_list) > 0) {
      pos_dt <- rbindlist(pos_list)

      p15d <- ggplot(pos_dt, aes(x = dist_to_tss, color = motif)) +
        geom_density(linewidth = 0.8, alpha = 0.7) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
        annotate("text", x = 0, y = Inf, vjust = 2, label = "TSS", size = 3) +
        scale_x_continuous(limits = c(-2000, 0)) +
        labs(x = "Distance to TSS (bp)", y = "Density",
             title = "De novo motif positions relative to TSS",
             subtitle = "D. laeve promoters (2kb upstream)",
             color = "Motif") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "right")
      save_fig(p15d, BATCH_DIR, "fig15d_denovo_tss_metaplot", w = 10, h = 6)
    }
  }
  rm(prom_seqs); gc(verbose = FALSE)

  # --- Fig 1.5E: Co-occurrence heatmap (which de novo motifs appear together) ---
  cat("Fig 1.5E: Co-occurrence heatmap...\n")

  # Re-scan promoters for presence/absence (use matches mode)
  prom_regions2 <- regions[regions$region_type == "promoter"]
  prom_seqs2 <- getSeq(genome, prom_regions2)
  names(prom_seqs2) <- prom_regions2$gene_id

  n_cooc <- min(15, length(denovo))
  denovo_cooc_pfms <- lapply(seq_len(n_cooc), function(i) {
    m <- denovo[[i]]
    pfm_mat <- m@motif
    if (nrow(pfm_mat) != 4) pfm_mat <- t(pfm_mat)
    rownames(pfm_mat) <- c("A", "C", "G", "T")
    count_mat <- round(pfm_mat * 100)
    PFMatrix(ID = paste0("dn_", i), name = denovo_info$consensus[i],
             profileMatrix = count_mat)
  })
  denovo_cooc_list <- do.call(PFMatrixList, denovo_cooc_pfms)

  # Scan in chunks (same pattern as JASPAR scan)
  SEQ_CHUNK2 <- 5000L
  n_prom <- length(prom_seqs2)
  presence_mat <- matrix(FALSE, nrow = n_prom, ncol = n_cooc)

  for (si in seq(1, n_prom, by = SEQ_CHUNK2)) {
    sj <- min(si + SEQ_CHUNK2 - 1L, n_prom)
    chunk_seqs <- prom_seqs2[si:sj]
    match_res <- tryCatch(
      matchMotifs(denovo_cooc_list, chunk_seqs, out = "matches", p.cutoff = 5e-5),
      error = function(e) NULL
    )
    if (!is.null(match_res)) {
      presence_mat[si:sj, ] <- as.matrix(motifMatches(match_res))
    }
  }
  colnames(presence_mat) <- denovo_info$consensus[1:n_cooc]

  # Compute co-occurrence (Jaccard index)
  cooc_mat <- matrix(0, n_cooc, n_cooc)
  for (i in seq_len(n_cooc)) {
    for (j in seq_len(n_cooc)) {
      a <- presence_mat[, i]
      b <- presence_mat[, j]
      intersection <- sum(a & b)
      union_ab <- sum(a | b)
      cooc_mat[i, j] <- if (union_ab > 0) intersection / union_ab else 0
    }
  }
  rownames(cooc_mat) <- colnames(cooc_mat) <- denovo_info$consensus[1:n_cooc]

  # Plot with ComplexHeatmap
  library(ComplexHeatmap)
  library(circlize)

  ht <- Heatmap(cooc_mat,
    name = "Jaccard",
    col = colorRamp2(c(0, 0.15, 0.5), c("white", "#FDE725", "#440154")),
    cluster_rows = TRUE, cluster_columns = TRUE,
    row_names_gp = gpar(fontsize = 8),
    column_names_gp = gpar(fontsize = 8),
    column_names_rot = 45,
    column_title = "De novo motif co-occurrence in D. laeve promoters (Jaccard index)",
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (i != j && cooc_mat[i, j] > 0.1) {
        grid.text(sprintf("%.2f", cooc_mat[i, j]), x, y, gp = gpar(fontsize = 6))
      }
    })

  png(file.path(BATCH_DIR, "figures/fig15e_denovo_cooccurrence_heatmap.png"),
      width = 9, height = 8, units = "in", res = 300)
  draw(ht)
  dev.off()
  cairo_pdf(file.path(BATCH_DIR, "figures/fig15e_denovo_cooccurrence_heatmap.pdf"),
            width = 9, height = 8)
  draw(ht)
  dev.off()
  cat("  Saved fig15e\n")

  save_data(as.data.table(cooc_mat, keep.rownames = "motif"), BATCH_DIR, "denovo_cooccurrence_jaccard")
  rm(prom_seqs2, presence_mat, cooc_mat); gc(verbose = FALSE)

  # --- Fig 1.5F: JASPAR similarity of de novo motifs ---
  cat("Fig 1.5F: De novo vs JASPAR similarity...\n")

  # Compare each de novo motif to all JASPAR motifs
  n_compare <- min(12, length(denovo))
  jaspar_um <- convert_motifs(jaspar_pfms, class = "universalmotif-universalmotif")

  best_matches <- data.table(
    denovo_motif = character(0),
    jaspar_match = character(0),
    jaspar_id = character(0),
    tf_class = character(0),
    similarity = numeric(0)
  )

  for (i in seq_len(n_compare)) {
    dm <- denovo[[i]]
    dm@name <- denovo_info$consensus[i]

    # Compare to all JASPAR
    sims <- tryCatch({
      compare_motifs(c(list(dm), jaspar_um), method = "PCC", min.overlap = 4)
    }, error = function(e) NULL)

    if (!is.null(sims) && nrow(sims) > 1) {
      # First row/col is the de novo motif
      scores <- sims[1, -1]  # all JASPAR similarities to this de novo
      top5_idx <- order(scores, decreasing = TRUE)[1:min(5, length(scores))]

      for (j in top5_idx) {
        best_matches <- rbindlist(list(best_matches, data.table(
          denovo_motif = denovo_info$consensus[i],
          jaspar_match = motif_meta$motif_name[j],
          jaspar_id = motif_meta$motif_id[j],
          tf_class = motif_meta$tf_class[j],
          similarity = scores[j]
        )))
      }
    }
    if (i %% 3 == 0) cat(sprintf("  %d/%d compared\n", i, n_compare))
  }

  if (nrow(best_matches) > 0) {
    save_data(best_matches, BATCH_DIR, "denovo_jaspar_similarity")

    # Plot top match per de novo motif
    top_per_denovo <- best_matches[, .SD[which.max(similarity)], by = denovo_motif]
    top_per_denovo[, label := paste0(jaspar_match, " (", tf_class, ")")]
    top_per_denovo <- top_per_denovo[order(-similarity)]

    p15f <- ggplot(top_per_denovo, aes(x = reorder(denovo_motif, similarity),
                                        y = similarity, fill = similarity)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = label), hjust = -0.05, size = 3) +
      scale_fill_viridis_c(option = "plasma", guide = "none") +
      coord_flip() +
      scale_y_continuous(limits = c(0, 1.3)) +
      labs(x = "De novo motif", y = "PCC similarity to best JASPAR match",
           title = "De novo motifs: best JASPAR 2024 matches",
           subtitle = "Pearson correlation of position weight matrices") +
      theme_minimal(base_size = 11)
    save_fig(p15f, BATCH_DIR, "fig15f_denovo_jaspar_similarity", w = 12, h = 7)
  }

  # --- Fig 1.5G: Known motif enrichment from HOMER (top 30) ---
  cat("Fig 1.5G: Known motif enrichment from HOMER...\n")
  known_file <- file.path(BATCH_DIR, "data/homer_denovo/knownResults.txt")
  if (file.exists(known_file)) {
    known <- fread(known_file, sep = "\t")
    # Clean column names
    setnames(known, 1, "motif_name")
    setnames(known, 2, "consensus")

    # Parse target and background percentages
    if (ncol(known) >= 8) {
      setnames(known, 7, "target_pct_str")
      setnames(known, 9, "bg_pct_str")
      known[, target_pct := as.numeric(gsub("%", "", target_pct_str))]
      known[, bg_pct := as.numeric(gsub("%", "", bg_pct_str))]
      known[, fold := target_pct / pmax(bg_pct, 0.01)]
      known[, log_p := as.numeric(known[[4]])]

      top_known <- known[order(log_p)][1:min(30, .N)]
      # Clean motif names (remove database info)
      top_known[, short_name := sub("/.*", "", motif_name)]
      top_known[, short_name := factor(short_name, levels = rev(short_name))]

      p15g <- ggplot(top_known, aes(x = short_name, y = fold, fill = -log_p)) +
        geom_col(width = 0.7) +
        scale_fill_viridis_c(name = "-log(p)", option = "inferno") +
        coord_flip() +
        labs(x = NULL, y = "Enrichment fold (target / background)",
             title = "Known TF motif enrichment in D. laeve promoters",
             subtitle = "HOMER findMotifs.pl | Top 30 by p-value") +
        theme_minimal(base_size = 10) +
        theme(axis.text.y = element_text(size = 7))
      save_fig(p15g, BATCH_DIR, "fig15g_known_motif_enrichment", w = 11, h = 9)
    }
  }

  cat("De novo motif characterization complete.\n")
} else {
  cat("No HOMER de novo motif file found — skipping characterization.\n")
}

cat("\n=== Batch 1.5 complete ===\n")
