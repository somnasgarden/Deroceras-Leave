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

options(stringsAsFactors = FALSE, scipen = 999)

# -- Cache dir for JASPAR/BiocFileCache --
cache_home <- Sys.getenv("XDG_CACHE_HOME", "~/.cache")
Sys.setenv(XDG_CACHE_HOME = cache_home)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(rtracklayer)
  library(Biostrings)
  library(GenomeInfoDb)
  library(TFBSTools)
  library(JASPAR2024)
  library(motifmatchr)
  library(universalmotif)
})

# -- Paths --
cluster_root <- Sys.getenv("CLUSTER_ROOT", getwd())
cache_dir    <- file.path(cluster_root, "genome/cache")
out_dir      <- file.path(cluster_root, "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

keep_chr <- paste0("chr", 1:31)

# Cluster data paths (fallbacks if RDS caches don't exist)
CLUSTER_GFF <- "/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff"
GENOME_FASTA <- Sys.getenv("GENOME_FASTA", "")

cat("=== Expanded Motif Annotation (Cluster) ===\n")
cat("Cluster root:", cluster_root, "\n")
cat("Cache dir:", cache_dir, "\n")
cat("Output dir:", out_dir, "\n")
cat("Available RAM:", system("free -h | grep Mem | awk '{print $2}'", intern = TRUE), "\n\n")


# =============================================================================
# A. LOAD GENOME + GFF, DEFINE EXTENDED REGIONS
# =============================================================================
cat("=== A. Loading data and defining regions ===\n")

# -- A1. Load genome --
# Primary: BSgenome package installed on cluster R environment
# Fallback: RDS cache or FASTA file
genome_loaded <- FALSE

# Try BSgenome package first (this is how the cluster R env has the genome)
tryCatch({
  suppressPackageStartupMessages(library(BSgenome.Dlaeve.NCBI.dlgm))
  bsg <- BSgenome.Dlaeve.NCBI.dlgm
  cat("Genome loaded from BSgenome.Dlaeve.NCBI.dlgm\n")

  # Extract chr1-31 as DNAStringSet for motifmatchr compatibility
  genome <- getSeq(bsg, keep_chr)
  names(genome) <- keep_chr
  chr_lengths <- setNames(seqlengths(bsg)[keep_chr], keep_chr)

  # Cache as RDS for future runs (faster reload)
  genome_rds <- file.path(cache_dir, "genome_chr1_31.rds")
  if (!file.exists(genome_rds)) {
    cat("Caching genome as RDS for faster future loads...\n")
    saveRDS(genome, genome_rds)
  }
  genome_loaded <- TRUE
}, error = function(e) {
  cat("BSgenome package not available:", e$message, "\n")
})

# Fallback: RDS cache
if (!genome_loaded) {
  genome_rds <- file.path(cache_dir, "genome_chr1_31.rds")
  if (file.exists(genome_rds)) {
    cat("Loading genome from RDS cache...\n")
    genome <- readRDS(genome_rds)
    chr_lengths <- setNames(width(genome), names(genome))
    genome_loaded <- TRUE
  }
}

# Fallback: FASTA file
if (!genome_loaded) {
  fasta_candidates <- c(
    GENOME_FASTA,
    file.path(cluster_root, "genome/derLaeGenome_chr1_31.fasta"),
    "/mnt/data/alfredvar/30-Genoma/derLaeGenome_chr1_31.fasta",
    "/mnt/data/alfredvar/30-Genoma/derLaeGenome_namesDlasi_v2.fasta"
  )
  fasta_path <- ""
  for (f in fasta_candidates) {
    if (f != "" && file.exists(f)) { fasta_path <- f; break }
  }
  if (fasta_path == "") stop("No genome found. Install BSgenome.Dlaeve.NCBI.dlgm or set GENOME_FASTA.")

  cat("Parsing genome FASTA:", fasta_path, "\n")
  genome <- readDNAStringSet(fasta_path)
  genome <- genome[names(genome) %in% keep_chr]
  chr_lengths <- setNames(width(genome), names(genome))

  genome_rds <- file.path(cache_dir, "genome_chr1_31.rds")
  cat("Caching genome RDS for future runs...\n")
  saveRDS(genome, genome_rds)
}

cat(sprintf("Genome: %d chromosomes, %s bp total\n",
            length(genome), format(sum(chr_lengths), big.mark = ",")))

# -- A2. Load GFF --
gff_rds <- file.path(cache_dir, "gff_chr1_31.rds")
if (file.exists(gff_rds)) {
  cat("Loading GFF from RDS cache...\n")
  gff <- readRDS(gff_rds)
} else if (file.exists(CLUSTER_GFF)) {
  cat("Parsing GFF from cluster path...\n")
  gff <- import(CLUSTER_GFF)
  gff <- gff[seqnames(gff) %in% keep_chr]
  seqlevels(gff) <- keep_chr
  saveRDS(gff, gff_rds)
} else {
  stop("No GFF found. Place gff_chr1_31.rds in cache or check CLUSTER_GFF path.")
}

genes <- gff[gff$type == "gene"]
cat("Genes:", length(genes), "\n")

# -- A3. Define extended regions --
regions_rds <- file.path(cache_dir, "extended_regions_10kb.rds")
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
homer_genome_fa <- file.path(cluster_root, "genome", "derLaeGenome_chr1_31.fasta")
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
cat("Scanning per chromosome × per region type to manage memory.\n")
cat(sprintf("  %d motifs × %d regions across %d chromosomes\n",
            length(jaspar_pfms), length(regions), length(keep_chr)))

# Strategy: for each region type, for each chromosome:
#   1. Extract sequences for that chr+region
#   2. motifmatchR(out="positions")
#   3. Convert to genomic coordinates
#   4. Save chunk
#   5. Free memory

region_types <- c("promoter", "upstream_distal", "gene_body", "downstream")
all_chunk_files <- c()
total_hits <- 0L
scan_start_global <- proc.time()

for (rtype in region_types) {
  cat(sprintf("\n--- Scanning region: %s ---\n", rtype))
  type_regions <- regions[regions$region_type == rtype]

  for (chr_name in keep_chr) {
    chr_reg <- type_regions[seqnames(type_regions) == chr_name]
    if (length(chr_reg) == 0) next

    # Extract sequences
    chr_seqs <- getSeq(genome, chr_reg)
    names(chr_seqs) <- chr_reg$gene_id

    if (length(chr_seqs) == 0) next

    # Scan with motifmatchR
    tryCatch({
      pos_result <- matchMotifs(
        jaspar_pfms, chr_seqs,
        out = "positions", p.cutoff = 5e-5
      )
      pos_list <- motifPositions(pos_result)
      rm(pos_result); gc(verbose = FALSE)

      # Convert to genomic coordinates
      reg_lookup <- data.frame(
        gene_id  = chr_reg$gene_id,
        reg_chr  = as.character(seqnames(chr_reg)),
        reg_start = start(chr_reg),
        reg_end   = end(chr_reg),
        reg_strand = as.character(strand(chr_reg)),
        stringsAsFactors = FALSE
      )

      chunk_hits <- vector("list", length(pos_list))
      for (i in seq_along(pos_list)) {
        gr <- pos_list[[i]]
        if (length(gr) == 0) next

        hit_genes <- as.character(seqnames(gr))
        ridx      <- match(hit_genes, reg_lookup$gene_id)
        if (all(is.na(ridx))) next

        rel_s <- start(gr)
        rel_e <- end(gr)
        h_str <- as.character(strand(gr))

        rs <- reg_lookup$reg_strand[ridx]
        ps <- reg_lookup$reg_start[ridx]
        pe <- reg_lookup$reg_end[ridx]

        gen_s <- ifelse(rs == "+", ps + rel_s - 1L, pe - rel_e + 1L)
        gen_e <- ifelse(rs == "+", ps + rel_e - 1L, pe - rel_s + 1L)
        gen_str <- ifelse(h_str == "*", "*",
                          ifelse(rs == "+", h_str,
                                 ifelse(h_str == "+", "-", "+")))

        score_vals <- if ("score" %in% colnames(mcols(gr))) {
          mcols(gr)$score
        } else rep(NA_real_, length(gr))

        chunk_hits[[i]] <- data.frame(
          chr        = reg_lookup$reg_chr[ridx],
          start      = gen_s,
          end        = gen_e,
          strand     = gen_str,
          motif_id   = motif_meta$motif_id[i],
          motif_name = motif_meta$motif_name[i],
          tf_class   = motif_meta$tf_class[i],
          gene_id    = hit_genes,
          region_type = rtype,
          score      = score_vals,
          source     = "JASPAR",
          stringsAsFactors = FALSE
        )
      }

      dt <- rbindlist(chunk_hits, fill = TRUE)
      rm(chunk_hits, pos_list); gc(verbose = FALSE)

      if (nrow(dt) > 0) {
        chunk_file <- file.path(out_dir,
          sprintf("_chunk_%s_%s.tsv", rtype, chr_name))
        fwrite(dt, chunk_file, sep = "\t")
        all_chunk_files <- c(all_chunk_files, chunk_file)
        total_hits <- total_hits + nrow(dt)
      }
      rm(dt)

    }, error = function(e) {
      cat(sprintf("  WARNING: Error on %s/%s: %s\n", rtype, chr_name, e$message))
    })

    rm(chr_seqs); gc(verbose = FALSE)
  }

  cat(sprintf("  %s complete — running total: %s hits\n",
              rtype, format(total_hits, big.mark = ",")))
}

scan_time_global <- (proc.time() - scan_start_global)[3]
cat(sprintf("\nTotal JASPAR scanning time: %.1f minutes\n", scan_time_global / 60))
cat(sprintf("Total JASPAR hits: %s\n", format(total_hits, big.mark = ",")))


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
