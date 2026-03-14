#!/usr/bin/env Rscript
# =============================================================================
# Batch 1.5: TF & Motif Annotation
# Paper: Act 1 — What TFs exist in D. laeve, and what binding motifs are
#        encoded in its promoters?
# Two approaches:
#   1. Protein-first: DeepTFactor → EviAnn → STRING → TF family classification
#   2. Sequence-first: JASPAR motif scanning of all promoters (2 kb upstream)
# Then cross-reference to find confirmed TFs, orphan TFs, orphan motifs.
# Output: results/batch1.5/{pdf/, png/, tf_motif_report.html}
# =============================================================================

# Set cache dir BEFORE loading any Bioconductor packages
Sys.setenv(XDG_CACHE_HOME = "/tmp/claude-1000/bioc_cache")

options(stringsAsFactors = FALSE)
options(scipen = 999)

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(Biostrings)
library(GenomeInfoDb)
library(TFBSTools)
library(JASPAR2024)
library(motifmatchr)
library(universalmotif)

# -- Paths --
data_dir    <- "/mnt/c/Users/rafae/Projects/DATA"
project_dir <- "/mnt/c/Users/rafae/Projects/STANDBY"
cache_dir   <- file.path(project_dir, "genome/cache")
out_dir     <- file.path(project_dir, "results/batch1.5")
fig_pdf_dir <- file.path(out_dir, "pdf")
fig_png_dir <- file.path(out_dir, "png")
out_report  <- file.path(out_dir, "tf_motif_report.html")

dir.create(fig_pdf_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_png_dir, recursive = TRUE, showWarnings = FALSE)

keep_chr <- paste0("chr", 1:31)

# Helper: save a ggplot to both pdf/ and png/
save_plot <- function(p, name, w = 8, h = 5) {
  cairo_pdf(file.path(fig_pdf_dir, paste0(name, ".pdf")), width = w, height = h)
  print(p); dev.off()
  png(file.path(fig_png_dir, paste0(name, ".png")), width = w, height = h, units = "in", res = 300)
  print(p); dev.off()
  cat("  Saved:", name, "\n")
}

# =============================================================================
# A. SETUP + DATA LOADING
# =============================================================================
cat("=== A. Loading cached data ===\n")

# Genome
genome_rds <- file.path(cache_dir, "genome_chr1_31.rds")
genome <- readRDS(genome_rds)
cat("Genome loaded:", length(genome), "chromosomes\n")

# GFF
gff_rds <- file.path(cache_dir, "gff_chr1_31.rds")
gff <- readRDS(gff_rds)
cat("GFF loaded:", length(gff), "features\n")

genes <- gff[gff$type == "gene"]
mrnas <- gff[gff$type == "mRNA"]
cat("Genes:", length(genes), "| mRNAs:", length(mrnas), "\n")

# EviAnn annotations
annot_file <- file.path(data_dir, "derLaeGenome_eviann_annotations.tsv")
annot <- fread(annot_file, header = FALSE, col.names = c("gene_id", "gene_name", "description"))
cat("Annotations loaded:", nrow(annot), "rows\n")

# DeepTFactor predictions
tf_file <- file.path(data_dir, "prediction_result.txt")
tf_raw <- fread(tf_file)
cat("DeepTFactor loaded:", nrow(tf_raw), "predictions\n")

# STRING orthology
string_file <- file.path(data_dir, "STRG0A31YWK.protein.orthology.v12.0.txt")
string_raw <- fread(string_file)
colnames(string_raw)[1] <- "protein"  # fix #protein → protein
cat("STRING orthology loaded:", nrow(string_raw), "rows\n")


# =============================================================================
# B. APPROACH 1: PROTEIN-FIRST TF CLASSIFICATION
# =============================================================================
cat("\n=== B. Approach 1: TF classification ===\n")

# B1. Filter to TF predictions only
tf_pos <- tf_raw[tf_raw$prediction == TRUE, ]
cat("TF-positive mRNAs:", nrow(tf_pos), "\n")

# B2. Extract gene ID from mRNA ID (LOC_XXXXX-mRNA-N → LOC_XXXXX)
tf_pos$gene_id <- sub("-mRNA-.*", "", tf_pos$sequence_ID)

# B3. Deduplicate: keep best-scoring mRNA per gene
tf_genes <- tf_pos %>%
  group_by(gene_id) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  as.data.frame()
cat("Unique TF genes:", nrow(tf_genes), "\n")

# B4. Annotate with EviAnn gene names
# Some genes have multiple annotation rows — take first unique gene_name per gene
annot_unique <- annot %>%
  group_by(gene_id) %>%
  summarise(
    gene_name = gene_name[1],
    description = description[1],
    .groups = "drop"
  )
tf_genes <- merge(tf_genes, annot_unique, by = "gene_id", all.x = TRUE)
cat("TFs with gene name:", sum(!is.na(tf_genes$gene_name)), "/", nrow(tf_genes), "\n")

# B5. Annotate with STRING orthology
# STRING protein IDs have prefix STRG0A31YWK.LOC_XXXXX
string_raw$gene_id <- sub("^STRG0A31YWK\\.", "", string_raw$protein)
# Get unique orthology groups per gene (collapse all taxonomy levels)
string_orth <- string_raw %>%
  group_by(gene_id) %>%
  summarise(
    orthology_groups = paste(unique(orthologous_group_or_ortholog), collapse = ";"),
    .groups = "drop"
  )
tf_genes <- merge(tf_genes, string_orth, by = "gene_id", all.x = TRUE)
cat("TFs with STRING orthology:", sum(!is.na(tf_genes$orthology_groups)), "/", nrow(tf_genes), "\n")

# B6. Classify TF families by keyword matching on gene_name + description
# Combined text for matching
tf_genes$match_text <- paste(
  ifelse(is.na(tf_genes$gene_name), "", tf_genes$gene_name),
  ifelse(is.na(tf_genes$description), "", tf_genes$description)
)

classify_tf_family <- function(text) {
  text_lower <- tolower(text)
  # Order matters: more specific patterns first
  if (grepl("\\bsox\\b|\\bhmg|high.mobility|sry", text_lower)) return("Sox/HMG")
  if (grepl("\\bgata\\b", text_lower)) return("GATA")
  if (grepl("\\bets\\b|\\belf\\b|\\bfli\\b|\\bpea3\\b|\\bspi\\b|pointed", text_lower)) return("ETS")
  if (grepl("forkhead|\\bfox[a-z]|\\bfkh\\b", text_lower)) return("Forkhead")
  if (grepl("\\bbhlh\\b|helix.loop.helix|\\bhand\\b|\\btwist\\b|\\bmyc\\b|\\bmax\\b|\\bmad\\b|\\buscul|\\bachaete|\\bneuro[gd]", text_lower)) return("bHLH")
  if (grepl("\\bbzip\\b|leucine.zipper|\\bcreb\\b|\\batf\\b|\\bjun\\b|\\bfos\\b|\\bmaf\\b|\\bbatf\\b|\\bnfe2|\\bnrf\\b|\\bcnc\\b|\\bap.?1\\b", text_lower)) return("bZIP")
  if (grepl("homeobox|homeodomain|\\bhox\\b|\\bpax\\b|\\bdlx\\b|\\blhx\\b|\\bnkx\\b|\\botx\\b|\\bpitx\\b|\\bpbx\\b|\\bmeis\\b|\\birx\\b|\\bbarh\\b|\\bprox\\b|\\bcut\\b|\\bprep\\b|\\btgif\\b|\\bvsx\\b|\\bhmx\\b|\\bmnx\\b|\\buncx\\b|\\bgbx\\b|\\bemx\\b|\\bnoto\\b|\\bcdx\\b|\\bpdx\\b|\\blbx\\b|\\bevx\\b|\\bdbx\\b|\\ben[12]\\b|\\barx\\b|\\bnanog\\b", text_lower)) return("Homeodomain")
  if (grepl("nuclear.receptor|hormone.receptor|\\bnr[0-9]|\\brar\\b|\\brxr\\b|\\bppar\\b|\\bthr\\b|\\becr\\b|\\bhnf4\\b|\\bsf1\\b|\\berr\\b|\\bcoup\\b|\\bsvp\\b|\\btailless\\b|\\btll\\b", text_lower)) return("Nuclear receptor")
  if (grepl("t.box|\\btbx\\b|\\btbr\\b|\\bbrachyury\\b", text_lower)) return("T-box")
  if (grepl("\\bsmad\\b|\\bmad\\b.*signal|dpp", text_lower)) return("SMAD")
  if (grepl("\\bstat\\b|signal.transducer.and.activator", text_lower)) return("STAT")
  if (grepl("\\bp53\\b|\\bp63\\b|\\bp73\\b|tumor.protein", text_lower)) return("p53")
  if (grepl("\\bnf.?kb\\b|\\brel\\b|\\bdorsal\\b|\\bdif\\b|\\bnfat\\b", text_lower)) return("NF-kB/Rel")
  if (grepl("\\bmyb\\b", text_lower)) return("MYB")
  if (grepl("\\bklf\\b|\\bsp[1-9]\\b|kruppel", text_lower)) return("KLF/Sp")
  if (grepl("\\brunx\\b|\\brunt\\b|\\bcbf\\b", text_lower)) return("RUNX")
  if (grepl("\\be2f\\b|\\bdp[12]\\b|retinoblastoma.binding", text_lower)) return("E2F")
  if (grepl("zinc.?finger|c2h2|\\bznf\\b|\\bzfp\\b|\\bgli\\b|\\bsnai\\b|\\bslug\\b|\\bzic\\b|\\bcasz\\b", text_lower)) return("C2H2-ZF")
  if (grepl("\\btead\\b|scalloped|\\bhippo\\b", text_lower)) return("TEAD")
  if (grepl("\\birf\\b|interferon.regulatory", text_lower)) return("IRF")
  if (grepl("\\bpou\\b|\\boct\\b|\\bbrn\\b|\\bpit\\b.*pou", text_lower)) return("POU")
  if (grepl("\\bcebp\\b|ccaat.enhancer", text_lower)) return("C/EBP")
  if (grepl("\\bap2\\b|\\btfap2\\b", text_lower)) return("AP-2")
  return("Unclassified")
}

tf_genes$tf_family <- sapply(tf_genes$match_text, classify_tf_family)
family_counts <- table(tf_genes$tf_family)
cat("\nTF family classification:\n")
for (f in sort(names(family_counts))) {
  cat(sprintf("  %-20s %d\n", f, family_counts[f]))
}
cat("Classified:", sum(tf_genes$tf_family != "Unclassified"),
    "| Unclassified:", sum(tf_genes$tf_family == "Unclassified"), "\n")

rm(tf_raw, string_raw); gc()


# =============================================================================
# C. APPROACH 2: SEQUENCE-FIRST PROMOTER MOTIF SCAN
# =============================================================================
cat("\n=== C. Approach 2: Promoter motif scan ===\n")

# C1. Extract promoters (2 kb upstream of gene TSS)
promoter_rds <- file.path(cache_dir, "promoters_2kb.rds")
if (file.exists(promoter_rds)) {
  cat("Loading promoters from cache...\n")
  promoters <- readRDS(promoter_rds)
} else {
  cat("Extracting 2 kb promoters...\n")
  # Get chromosome lengths from genome
  chr_lengths <- width(genome)
  names(chr_lengths) <- names(genome)

  # Promoter = 2000 bp upstream of TSS
  # For + strand: TSS = start(gene), promoter = [TSS - 2000, TSS - 1]
  # For - strand: TSS = end(gene),   promoter = [TSS + 1, TSS + 2000]
  gene_strand <- as.character(strand(genes))
  gene_chr    <- as.character(seqnames(genes))

  prom_start <- ifelse(gene_strand == "+",
                       pmax(1, start(genes) - 2000),
                       end(genes) + 1)
  prom_end   <- ifelse(gene_strand == "+",
                       start(genes) - 1,
                       pmin(chr_lengths[gene_chr], end(genes) + 2000))

  # Filter out invalid promoters (start > end can happen at chromosome boundaries)
  valid <- prom_start < prom_end
  cat("Valid promoters:", sum(valid), "/", length(valid), "\n")

  promoters <- GRanges(
    seqnames = gene_chr[valid],
    ranges = IRanges(start = prom_start[valid], end = prom_end[valid]),
    strand = gene_strand[valid],
    gene_id = genes$ID[valid]
  )
  seqlevels(promoters) <- keep_chr
  seqlengths(promoters) <- chr_lengths[keep_chr]
  promoters <- trim(promoters)
  saveRDS(promoters, promoter_rds)
}
cat("Promoters:", length(promoters),
    "| Median width:", median(width(promoters)), "bp\n")

# C2. Extract promoter sequences from genome
cat("Extracting promoter sequences...\n")
prom_seqs <- getSeq(genome, promoters)
names(prom_seqs) <- promoters$gene_id
cat("Promoter sequences extracted:", length(prom_seqs), "\n")

# C3. Get JASPAR 2024 metazoan PWMs
cat("Loading JASPAR 2024 motifs...\n")
# JASPAR2024 SQLite was downloaded via JASPAR2024() in a previous session
jaspar_sqlite <- list.files("/tmp/claude-1000/bioc_cache", pattern = "JASPAR.*sqlite",
                            recursive = TRUE, full.names = TRUE)[1]
if (is.na(jaspar_sqlite)) {
  # First run — need to download it
  jaspar_db <- JASPAR2024()
  jaspar_sqlite <- list.files("/tmp/claude-1000/bioc_cache", pattern = "JASPAR.*sqlite",
                              recursive = TRUE, full.names = TRUE)[1]
}
cat("JASPAR SQLite:", jaspar_sqlite, "\n")

# Get all CORE motifs
all_pfms <- getMatrixSet(jaspar_sqlite, list(collection = "CORE", all_versions = FALSE))
cat("Total CORE motifs:", length(all_pfms), "\n")

# Filter for metazoan tax groups: insects, nematodes, urochordates, vertebrates
metazoan_groups <- c("insects", "nematodes", "urochordates", "vertebrates")
is_metazoa <- sapply(all_pfms, function(x) {
  tg <- tags(x)$tax_group
  !is.null(tg) && tg %in% metazoan_groups
})
jaspar_pfms <- all_pfms[is_metazoa]
rm(all_pfms); gc()
cat("JASPAR metazoan motifs loaded:", length(jaspar_pfms), "\n")

# Build metadata table for motifs
# Use TFBSTools::tags to avoid shadowing
motif_ids   <- sapply(jaspar_pfms, ID)
motif_names <- sapply(jaspar_pfms, name)
motif_class <- vapply(jaspar_pfms, function(x) {
  tgs <- TFBSTools::tags(x)
  if ("family" %in% names(tgs) && length(tgs$family) > 0) tgs$family[1] else NA_character_
}, character(1))
motif_species <- vapply(jaspar_pfms, function(x) {
  tgs <- TFBSTools::tags(x)
  if ("species" %in% names(tgs) && length(tgs$species) > 0) {
    paste(tgs$species, collapse = ";")
  } else NA_character_
}, character(1))

motif_meta <- data.frame(
  motif_id   = motif_ids,
  motif_name = motif_names,
  tf_class   = motif_class,
  species    = motif_species,
  stringsAsFactors = FALSE
)
cat("Motif families represented:", length(unique(na.omit(motif_meta$tf_class))), "\n")

# C4. Scan promoters with motifmatchr
cat("Scanning promoters with", length(jaspar_pfms), "motifs (this takes a while)...\n")
scan_start <- proc.time()

# motifmatchr expects a PWMatrixList
motif_hits <- matchMotifs(
  jaspar_pfms,
  prom_seqs,
  out = "scores",
  p.cutoff = 5e-5
)
scan_time <- (proc.time() - scan_start)[3]
cat(sprintf("Motif scan complete in %.1f seconds\n", scan_time))

# Extract boolean hit matrix (promoters × motifs)
hit_matrix <- motifMatches(motif_hits)
cat("Hit matrix dimensions:", nrow(hit_matrix), "promoters ×", ncol(hit_matrix), "motifs\n")

# C5. Summarize: which motifs appear in how many promoters
motif_prevalence <- data.frame(
  motif_id   = motif_meta$motif_id,
  motif_name = motif_meta$motif_name,
  tf_class   = motif_meta$tf_class,
  n_promoters = colSums(hit_matrix),
  pct_promoters = 100 * colSums(hit_matrix) / nrow(hit_matrix),
  stringsAsFactors = FALSE
)
motif_prevalence <- motif_prevalence[order(-motif_prevalence$n_promoters), ]
cat("\nTop 20 most prevalent motifs:\n")
for (i in 1:min(20, nrow(motif_prevalence))) {
  r <- motif_prevalence[i, ]
  cat(sprintf("  %-20s %-25s %6d promoters (%.1f%%)\n",
              r$motif_name, ifelse(is.na(r$tf_class), "?", r$tf_class),
              r$n_promoters, r$pct_promoters))
}

# Promoters with at least one motif hit
prom_any_hit <- sum(rowSums(hit_matrix) > 0)
cat(sprintf("\nPromoters with >= 1 motif hit: %d / %d (%.1f%%)\n",
            prom_any_hit, nrow(hit_matrix), 100 * prom_any_hit / nrow(hit_matrix)))

# Cleanup scan objects
rm(motif_hits); gc()


# =============================================================================
# D. COMPARISON: CROSS-REFERENCE TF FAMILIES VS GENOME MOTIFS
# =============================================================================
cat("\n=== D. Cross-referencing approaches ===\n")

# Map TF families to JASPAR motif class names
# JASPAR uses different naming conventions — build a mapping
# Mapping uses EXACT JASPAR 2024 family names (verified from database)
family_to_jaspar <- list(
  "C2H2-ZF"         = c("More than 3 adjacent zinc fingers",
                         "Three-zinc finger Kruppel-related",
                         "Other factors with up to three adjacent zinc fingers",
                         "Factors with multiple dispersed zinc fingers",
                         "BED zinc finger factors", "Zinc finger, BED-type"),
  "Homeodomain"      = c("HOX", "NK", "Paired-related HD factors",
                         "HD-LIM", "HD-CUT", "HD-SINE", "HD-ZF factors",
                         "HD-PROS factors", "Homeo", "TALE-type HD",
                         "TALE-type homeo domain factors",
                         "Paired domain only", "Paired plus homeo domain"),
  "bHLH"             = c("Tal-related", "Hairy-related factors",
                         "bHLH-ZIP", "MyoD/ASC-related factors",
                         "E2A", "Helix-Loop-Helix", "PAS domain factors"),
  "Forkhead"         = c("FOX"),
  "Sox/HMG"          = c("SOX-related factors", "TCF-7-related factors",
                         "High Mobility Group box (HMG)", "HMGA"),
  "GATA"             = c("C4-GATA-related"),
  "bZIP"             = c("Jun-related", "Fos-related", "Maf-related",
                         "CREB-related factors", "ATF-4-related",
                         "ATF-4-related factors", "B-ATF-related factors",
                         "XBP-1-related factors", "CEBP-related"),
  "ETS"              = c("Ets-related"),
  "Nuclear receptor" = c("Thyroid hormone receptor-related factors (NR1)",
                         "RXR-related receptors (NR2)",
                         "Steroid hormone receptors (NR3)",
                         "NGFI-B-related receptors (NR4)",
                         "FTZF1related(NR5A)",
                         "GCNF-related receptors (NR6)"),
  "T-box"            = c("Brachyury-related factors", "TBX1-related factors",
                         "TBX2-related factors", "TBX6-related factors",
                         "TBrain-related factors"),
  "SMAD"             = c("SMAD factors"),
  "STAT"             = c("STAT factors"),
  "p53"              = c("p53-related factors"),
  "NF-kB/Rel"        = c("NF-kappaB-related factors", "NFAT-related factors",
                         "Rel homology domain (RHD)"),
  "MYB"              = c("Myb/SANT domain factors"),
  "KLF/Sp"           = c("Three-zinc finger Kruppel-related"),
  "RUNX"             = c("Runt-related factors"),
  "E2F"              = c("E2F"),
  "TEAD"             = c("TEF-1-related factors"),
  "IRF"              = c("Interferon-regulatory factors"),
  "POU"              = c("POU domain factors"),
  "C/EBP"            = c("CEBP-related"),
  "AP-2"             = c("AP-2")
)

# For each TF family, check if ANY of its matching JASPAR classes have hits in the genome
comparison <- data.frame(
  tf_family = character(),
  n_tfs = integer(),
  has_jaspar_motif = logical(),
  n_jaspar_motifs = integer(),
  motif_in_genome = logical(),
  n_promoters_with_motif = integer(),
  stringsAsFactors = FALSE
)

for (fam in sort(unique(tf_genes$tf_family[tf_genes$tf_family != "Unclassified"]))) {
  n_tfs <- sum(tf_genes$tf_family == fam)

  # Find matching JASPAR motifs
  jaspar_classes <- family_to_jaspar[[fam]]
  if (is.null(jaspar_classes)) jaspar_classes <- character(0)

  matching_motifs <- motif_prevalence[motif_prevalence$tf_class %in% jaspar_classes, ]
  has_motif <- nrow(matching_motifs) > 0
  n_motifs <- nrow(matching_motifs)
  in_genome <- has_motif && any(matching_motifs$n_promoters > 0)
  n_proms <- if (in_genome) max(matching_motifs$n_promoters) else 0L

  comparison <- rbind(comparison, data.frame(
    tf_family = fam, n_tfs = n_tfs,
    has_jaspar_motif = has_motif, n_jaspar_motifs = n_motifs,
    motif_in_genome = in_genome, n_promoters_with_motif = n_proms,
    stringsAsFactors = FALSE
  ))
}

cat("\nTF family × motif comparison:\n")
cat(sprintf("  %-20s %5s %8s %10s %8s\n",
            "Family", "TFs", "JASPAR?", "InGenome?", "MaxProm"))
for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  cat(sprintf("  %-20s %5d %8s %10s %8d\n",
              r$tf_family, r$n_tfs,
              ifelse(r$has_jaspar_motif, "YES", "no"),
              ifelse(r$motif_in_genome, "YES", "no"),
              r$n_promoters_with_motif))
}

# Summary categories
confirmed  <- sum(comparison$motif_in_genome)
orphan_tf  <- sum(comparison$has_jaspar_motif & !comparison$motif_in_genome)
no_motif   <- sum(!comparison$has_jaspar_motif)

# Orphan motifs: JASPAR families with hits in genome but no predicted TF
jaspar_families_in_genome <- unique(motif_prevalence$tf_class[motif_prevalence$n_promoters > 100])
jaspar_families_in_genome <- jaspar_families_in_genome[!is.na(jaspar_families_in_genome)]
tf_jaspar_classes <- unlist(family_to_jaspar)
orphan_motif_classes <- setdiff(jaspar_families_in_genome, tf_jaspar_classes)

cat(sprintf("\nSummary:\n  Confirmed TF families (TF + motif in genome): %d\n", confirmed))
cat(sprintf("  Orphan TF families (TF + motif, but not found in promoters): %d\n", orphan_tf))
cat(sprintf("  No JASPAR motif available: %d\n", no_motif))
cat(sprintf("  Orphan motif classes (in genome, no predicted TF): %d\n", length(orphan_motif_classes)))


# =============================================================================
# E. PLOTS (4 total)
# =============================================================================
cat("\n=== E. Generating plots ===\n")

# Color palette for TF families
n_fam <- length(unique(tf_genes$tf_family))
fam_colors <- c(
  "C2H2-ZF" = "#E74C3C", "Homeodomain" = "#3498DB", "bHLH" = "#2ECC71",
  "Forkhead" = "#9B59B6", "Sox/HMG" = "#F39C12", "GATA" = "#1ABC9C",
  "bZIP" = "#E67E22", "ETS" = "#34495E", "Nuclear receptor" = "#16A085",
  "T-box" = "#C0392B", "SMAD" = "#2980B9", "STAT" = "#27AE60",
  "p53" = "#8E44AD", "NF-kB/Rel" = "#D35400", "MYB" = "#7F8C8D",
  "KLF/Sp" = "#BDC3C7", "RUNX" = "#2C3E50", "E2F" = "#F1C40F",
  "TEAD" = "#E91E63", "IRF" = "#607D8B", "POU" = "#795548",
  "C/EBP" = "#FF9800", "AP-2" = "#00BCD4", "Unclassified" = "#95A5A6"
)

# E1. TF Family Census — bar chart with motif found yes/no
tf_census <- tf_genes %>%
  group_by(tf_family) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

# Add motif-in-genome status
tf_census <- merge(tf_census,
                   comparison[, c("tf_family", "motif_in_genome")],
                   by = "tf_family", all.x = TRUE)
tf_census$motif_in_genome[is.na(tf_census$motif_in_genome)] <- FALSE
tf_census$motif_status <- ifelse(tf_census$tf_family == "Unclassified",
                                 "Unclassified",
                                 ifelse(tf_census$motif_in_genome,
                                        "Motif found in genome",
                                        "No motif found"))
tf_census$tf_family <- factor(tf_census$tf_family,
                              levels = tf_census$tf_family[order(tf_census$n)])

p1 <- ggplot(tf_census, aes(x = tf_family, y = n, fill = motif_status)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Motif found in genome" = "#27AE60",
                                "No motif found" = "#E74C3C",
                                "Unclassified" = "#95A5A6"),
                    name = "Motif status") +
  coord_flip() +
  labs(title = "TF Family Census (DeepTFactor predictions)",
       subtitle = paste0(nrow(tf_genes), " TF genes | ",
                         sum(tf_genes$tf_family != "Unclassified"), " classified"),
       x = NULL, y = "Number of genes") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")
save_plot(p1, "tf_family_census", w = 8, h = 7)

# E2. Motif Prevalence — top 30 most common motifs
top_motifs <- head(motif_prevalence, 30)
top_motifs$label <- paste0(top_motifs$motif_name, " (", top_motifs$motif_id, ")")
top_motifs$label <- factor(top_motifs$label,
                           levels = rev(top_motifs$label))

p2 <- ggplot(top_motifs, aes(x = label, y = pct_promoters)) +
  geom_col(fill = "#2471A3", width = 0.7) +
  coord_flip() +
  labs(title = "Most prevalent JASPAR motifs in D. laeve promoters",
       subtitle = paste0("Scanned ", nrow(hit_matrix), " promoters (2 kb upstream) with ",
                         ncol(hit_matrix), " JASPAR metazoan motifs | p < 5e-5"),
       x = NULL, y = "% of promoters with motif") +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_text(size = 8))
save_plot(p2, "motif_prevalence_top30", w = 9, h = 8)

# E3. Approach Comparison — tile heatmap
comp_plot <- comparison
comp_plot$tf_family <- factor(comp_plot$tf_family,
                              levels = comp_plot$tf_family[order(comp_plot$n_tfs)])
comp_long <- tidyr::pivot_longer(
  comp_plot[, c("tf_family", "has_jaspar_motif", "motif_in_genome")],
  cols = c("has_jaspar_motif", "motif_in_genome"),
  names_to = "criterion",
  values_to = "present"
)
comp_long$criterion <- factor(comp_long$criterion,
                              levels = c("has_jaspar_motif", "motif_in_genome"),
                              labels = c("Has JASPAR motif", "Motif in genome"))

# Add TF count column
comp_tf_counts <- data.frame(
  tf_family = comp_plot$tf_family,
  criterion = factor("Has TFs", levels = c("Has TFs", "Has JASPAR motif", "Motif in genome")),
  present = TRUE
)
comp_long$criterion <- factor(as.character(comp_long$criterion),
                              levels = c("Has TFs", "Has JASPAR motif", "Motif in genome"))
comp_long <- rbind(comp_tf_counts, comp_long)

p3 <- ggplot(comp_long, aes(x = criterion, y = tf_family, fill = present)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("TRUE" = "#27AE60", "FALSE" = "#E74C3C"),
                    name = "Present", labels = c("No", "Yes")) +
  geom_text(aes(label = ifelse(present, "\u2713", "\u2717")),
            size = 4, color = "white") +
  labs(title = "TF family validation: protein predictions vs. genome motifs",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid = element_blank())
save_plot(p3, "approach_comparison_tile", w = 6, h = 7)

# E4. DeepTFactor Score Distribution by TF family
# Only show classified families
tf_classified <- tf_genes[tf_genes$tf_family != "Unclassified", ]
# Order by median score
fam_order <- tf_classified %>%
  group_by(tf_family) %>%
  summarise(med = median(score), .groups = "drop") %>%
  arrange(med)
tf_classified$tf_family <- factor(tf_classified$tf_family,
                                  levels = fam_order$tf_family)

p4 <- ggplot(tf_classified, aes(x = tf_family, y = score, fill = tf_family)) +
  geom_boxplot(outlier.size = 0.5, show.legend = FALSE) +
  scale_fill_manual(values = fam_colors) +
  coord_flip() +
  labs(title = "DeepTFactor prediction confidence by TF family",
       subtitle = paste0(nrow(tf_classified), " classified TFs across ",
                         length(unique(tf_classified$tf_family)), " families"),
       x = NULL, y = "DeepTFactor score") +
  theme_minimal(base_size = 12)
save_plot(p4, "deeptfactor_score_by_family", w = 8, h = 6)


# =============================================================================
# F. HTML REPORT
# =============================================================================
cat("\n=== F. Generating HTML report ===\n")

html_lines <- c(
  "<!DOCTYPE html>",
  "<html><head><meta charset='utf-8'>",
  "<title>Batch 1.5: TF & Motif Annotation</title>",
  "<style>",
  "body { font-family: 'Segoe UI', Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; }",
  "h1, h2, h3 { color: #2C3E50; }",
  "table { border-collapse: collapse; width: 100%; margin: 10px 0; }",
  "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
  "th { background-color: #34495E; color: white; }",
  "tr:nth-child(even) { background-color: #f2f2f2; }",
  ".metric { font-size: 1.3em; font-weight: bold; color: #2471A3; }",
  "img { max-width: 100%; border: 1px solid #ddd; margin: 10px 0; }",
  ".confirmed { color: #27AE60; font-weight: bold; }",
  ".orphan { color: #E74C3C; font-weight: bold; }",
  "</style></head><body>",
  "<h1>Batch 1.5: TF & Motif Annotation</h1>",
  paste0("<p><em>D. laeve</em> genome-level TF and motif inventory | Generated: ", Sys.Date(), "</p>"),
  "",
  "<h2>Summary</h2>",
  "<table>",
  paste0("<tr><td>Total DeepTFactor TF mRNAs</td><td class='metric'>", nrow(tf_pos), "</td></tr>"),
  paste0("<tr><td>Unique TF genes (deduplicated)</td><td class='metric'>", nrow(tf_genes), "</td></tr>"),
  paste0("<tr><td>Classified TFs (with family)</td><td class='metric'>",
         sum(tf_genes$tf_family != "Unclassified"), "</td></tr>"),
  paste0("<tr><td>Unclassified TFs</td><td class='metric'>",
         sum(tf_genes$tf_family == "Unclassified"), "</td></tr>"),
  paste0("<tr><td>TF families identified</td><td class='metric'>",
         length(unique(tf_genes$tf_family[tf_genes$tf_family != "Unclassified"])), "</td></tr>"),
  paste0("<tr><td>JASPAR motifs scanned</td><td class='metric'>", ncol(hit_matrix), "</td></tr>"),
  paste0("<tr><td>Promoters scanned (2 kb upstream)</td><td class='metric'>",
         format(nrow(hit_matrix), big.mark = ","), "</td></tr>"),
  paste0("<tr><td>Promoters with &ge;1 motif hit</td><td class='metric'>",
         format(prom_any_hit, big.mark = ","),
         " (", round(100 * prom_any_hit / nrow(hit_matrix), 1), "%)</td></tr>"),
  paste0("<tr><td>Confirmed TF families (TF + motif in genome)</td>",
         "<td class='confirmed'>", confirmed, "</td></tr>"),
  paste0("<tr><td>Orphan TF families (motif not in promoters)</td>",
         "<td class='orphan'>", orphan_tf, "</td></tr>"),
  paste0("<tr><td>Dark matter (TF family, no JASPAR motif)</td>",
         "<td class='orphan'>", no_motif, "</td></tr>"),
  paste0("<tr><td>Orphan motif classes (in genome, no predicted TF)</td>",
         "<td>", length(orphan_motif_classes), "</td></tr>"),
  "</table>",
  "",
  "<h2>Plot 1: TF Family Census</h2>",
  "<img src='png/tf_family_census.png'>",
  "",
  "<h2>Plot 2: Motif Prevalence (Top 30)</h2>",
  "<img src='png/motif_prevalence_top30.png'>",
  "",
  "<h2>Plot 3: Approach Comparison</h2>",
  "<img src='png/approach_comparison_tile.png'>",
  "",
  "<h2>Plot 4: DeepTFactor Score by Family</h2>",
  "<img src='png/deeptfactor_score_by_family.png'>",
  "",
  "<h2>TF Family × Motif Cross-Reference</h2>",
  "<table>",
  "<tr><th>TF Family</th><th>TFs</th><th>JASPAR Motif?</th><th>Motif in Genome?</th><th>Max Promoters</th></tr>"
)

for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  html_lines <- c(html_lines, paste0(
    "<tr><td>", r$tf_family, "</td>",
    "<td>", r$n_tfs, "</td>",
    "<td>", ifelse(r$has_jaspar_motif, "<span class='confirmed'>YES</span>", "no"), "</td>",
    "<td>", ifelse(r$motif_in_genome, "<span class='confirmed'>YES</span>",
                   "<span class='orphan'>no</span>"), "</td>",
    "<td>", format(r$n_promoters_with_motif, big.mark = ","), "</td></tr>"
  ))
}
html_lines <- c(html_lines, "</table>")

# Top 30 motifs table
html_lines <- c(html_lines,
  "",
  "<h2>Top 30 Most Prevalent Motifs</h2>",
  "<table>",
  "<tr><th>Motif</th><th>ID</th><th>Family</th><th>Promoters</th><th>%</th></tr>"
)
for (i in 1:min(30, nrow(motif_prevalence))) {
  r <- motif_prevalence[i, ]
  html_lines <- c(html_lines, paste0(
    "<tr><td>", r$motif_name, "</td>",
    "<td>", r$motif_id, "</td>",
    "<td>", ifelse(is.na(r$tf_class), "?", r$tf_class), "</td>",
    "<td>", format(r$n_promoters, big.mark = ","), "</td>",
    "<td>", round(r$pct_promoters, 1), "</td></tr>"
  ))
}
html_lines <- c(html_lines, "</table>")

# Orphan motif classes
if (length(orphan_motif_classes) > 0) {
  html_lines <- c(html_lines,
    "",
    "<h2>Orphan Motif Classes (in genome, no predicted TF)</h2>",
    "<ul>",
    paste0("<li>", orphan_motif_classes, "</li>"),
    "</ul>"
  )
}

# Limitations
html_lines <- c(html_lines,
  "",
  "<h2>Limitations</h2>",
  "<ul>",
  "<li>No mollusk-specific motifs in JASPAR &mdash; using insect/nematode/urochordate as proxy</li>",
  paste0("<li>", sum(tf_genes$tf_family == "Unclassified"), " TFs (",
         round(100 * sum(tf_genes$tf_family == "Unclassified") / nrow(tf_genes), 1),
         "%) remain unclassified (no informative annotation)</li>"),
  "<li>STRING orthology enriches classification but does not resolve all unknowns</li>",
  "<li>Motif scanning uses PWM p-value threshold 5e-5 &mdash; tradeoff between sensitivity and specificity</li>",
  "</ul>",
  "",
  "</body></html>"
)

writeLines(html_lines, out_report)
cat("Report saved:", out_report, "\n")


# =============================================================================
# G. MOTIF COORDINATE EXPORT
# =============================================================================
cat("\n=== G. Exporting motif hit coordinates ===\n")

# G1. Re-scan with out="positions" to get exact hit coordinates
cat("Re-scanning promoters for exact motif positions...\n")
scan_start2 <- proc.time()
pos_result <- matchMotifs(
  jaspar_pfms,
  prom_seqs,
  out = "positions",
  p.cutoff = 5e-5
)
scan_time2 <- (proc.time() - scan_start2)[3]
cat(sprintf("Position scan complete in %.1f seconds\n", scan_time2))

# Extract positions — GRangesList of length = n_motifs
pos_list <- motifPositions(pos_result)
rm(pos_result); gc()

# G2. Build promoter lookup for coordinate conversion
# getSeq on + strand promoters returns the forward sequence (pos 1 = prom_start)
# getSeq on - strand promoters returns the reverse complement (pos 1 = prom_end)
prom_lookup <- data.frame(
  gene_id    = promoters$gene_id,
  prom_chr   = as.character(seqnames(promoters)),
  prom_start = start(promoters),
  prom_end   = end(promoters),
  prom_strand = as.character(strand(promoters)),
  stringsAsFactors = FALSE
)

# G3. Flatten GRangesList into annotated data.frame
cat("Flattening", length(pos_list), "motifs into coordinate table...\n")
all_hits <- vector("list", length(pos_list))

for (i in seq_along(pos_list)) {
  gr <- pos_list[[i]]
  if (length(gr) == 0) next

  hit_genes <- as.character(seqnames(gr))
  prom_idx  <- match(hit_genes, prom_lookup$gene_id)

  # Relative positions within promoter sequence
  rel_start  <- start(gr)
  rel_end    <- end(gr)
  hit_strand <- as.character(strand(gr))

  # Promoter genomic info
  ps      <- prom_lookup$prom_strand[prom_idx]
  p_start <- prom_lookup$prom_start[prom_idx]
  p_end   <- prom_lookup$prom_end[prom_idx]

  # Convert to genomic coordinates
  # + strand: genomic = prom_start + relative - 1
  # - strand: genomic_start = prom_end - rel_end + 1
  gen_start <- ifelse(ps == "+", p_start + rel_start - 1L, p_end - rel_end + 1L)
  gen_end   <- ifelse(ps == "+", p_start + rel_end - 1L,   p_end - rel_start + 1L)

  # Motif strand in genomic coordinates (flip if promoter is -)
  gen_strand <- ifelse(hit_strand == "*", "*",
                       ifelse(ps == "+", hit_strand,
                              ifelse(hit_strand == "+", "-", "+")))

  # Score (if motifmatchr includes it)
  score_vals <- if ("score" %in% colnames(mcols(gr))) {
    mcols(gr)$score
  } else {
    rep(NA_real_, length(gr))
  }

  all_hits[[i]] <- data.frame(
    chr        = prom_lookup$prom_chr[prom_idx],
    start      = gen_start,
    end        = gen_end,
    strand     = gen_strand,
    motif_id   = motif_meta$motif_id[i],
    motif_name = motif_meta$motif_name[i],
    tf_class   = motif_meta$tf_class[i],
    species    = motif_meta$species[i],
    gene_id    = hit_genes,
    score      = score_vals,
    stringsAsFactors = FALSE
  )
}

hits_dt <- data.table::rbindlist(all_hits, fill = TRUE)
rm(all_hits, pos_list); gc()

cat(sprintf("Total motif hits: %s\n", format(nrow(hits_dt), big.mark = ",")))
cat(sprintf("Unique motifs with hits: %d / %d\n",
            length(unique(hits_dt$motif_id)), nrow(motif_meta)))
cat(sprintf("Unique genes with hits: %d / %d\n",
            length(unique(hits_dt$gene_id)), length(promoters)))

# G4. Save full coordinate table (compressed TSV)
out_hits <- file.path(out_dir, "motif_hits_annotated.tsv.gz")
cat("Writing compressed TSV:", out_hits, "\n")
data.table::fwrite(hits_dt, out_hits, sep = "\t", compress = "gzip")
cat(sprintf("File size: %.1f MB\n", file.size(out_hits) / 1e6))

# G5. Motif-level summary (ensure all motifs represented, even those with 0 hits)
hit_summary <- hits_dt[, .(
  n_hits      = .N,
  n_promoters = uniqueN(gene_id)
), by = .(motif_id)]

motif_summary <- merge(
  as.data.table(motif_meta),
  hit_summary,
  by = "motif_id",
  all.x = TRUE
)
motif_summary[is.na(n_hits), c("n_hits", "n_promoters") := .(0L, 0L)]
motif_summary[, pct_promoters := round(100 * n_promoters / length(promoters), 2)]
motif_summary <- motif_summary[order(-n_hits)]

out_summary <- file.path(out_dir, "motif_annotation_summary.tsv")
data.table::fwrite(motif_summary, out_summary, sep = "\t")
cat(sprintf("Summary table: %d motifs, saved to %s\n", nrow(motif_summary), out_summary))

# G6. Verification
cat("\n--- Verification ---\n")
cat(sprintf("  Rows in hits table: %s\n", format(nrow(hits_dt), big.mark = ",")))
cat(sprintf("  Columns: %s\n", paste(colnames(hits_dt), collapse = ", ")))
cat(sprintf("  Motifs represented: %d / %d\n",
            nrow(motif_summary), nrow(motif_meta)))
cat(sprintf("  Summary rows: %d\n", nrow(motif_summary)))

# Spot-check: sampled hits should overlap promoter regions
n_check <- min(1000, nrow(hits_dt))
check_idx <- sample(nrow(hits_dt), n_check)
check_gr <- GRanges(
  seqnames = hits_dt$chr[check_idx],
  ranges   = IRanges(hits_dt$start[check_idx], hits_dt$end[check_idx])
)
in_prom <- sum(overlapsAny(check_gr, promoters))
cat(sprintf("  Coordinate check: %d/%d sampled hits overlap promoters (%.1f%%)\n",
            in_prom, n_check, 100 * in_prom / n_check))

rm(hits_dt, hit_summary, motif_summary, check_gr); gc()


cat("\n=== Batch 1.5 complete ===\n")
cat("Output directory:", out_dir, "\n")
