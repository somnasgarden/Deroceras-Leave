#!/usr/bin/env Rscript
# =============================================================================
# Prepare extended gene regions (10kb upstream + gene body + 10kb downstream)
# =============================================================================
# Lightweight — only needs GFF cache, no genome loading.
# Creates a GRanges object with 4 region types per gene:
#   promoter        : TSS - 2kb   to TSS - 1      (strand-aware)
#   upstream_distal : TSS - 10kb  to TSS - 2kb - 1 (strand-aware)
#   gene_body       : gene start  to gene end
#   downstream      : TES + 1     to TES + 10kb    (strand-aware)
#
# Output: genome/cache/extended_regions_10kb.rds
#         cluster/genome/cache/extended_regions_10kb.rds (if cluster/ exists)
# =============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

library(GenomicRanges)
library(Biostrings)

keep_chr <- paste0("chr", 1:31)

# -- Paths --
project_dir  <- "/mnt/c/Users/rafae/Projects/STANDBY"
cache_dir    <- file.path(project_dir, "genome/cache")
genome_fasta <- file.path(project_dir, "genome/derLaeGenome_chr1_31.fasta")

# -- Load GFF --
gff <- readRDS(file.path(cache_dir, "gff_chr1_31.rds"))
genes <- gff[gff$type == "gene"]
cat("Genes:", length(genes), "\n")

# -- Chromosome lengths --
# Try FASTA index, then Biostrings index, then GFF-derived
if (file.exists(paste0(genome_fasta, ".fai"))) {
  fai <- read.delim(paste0(genome_fasta, ".fai"), header = FALSE,
                    col.names = c("chr", "length", "offset", "linebases", "linewidth"))
  chr_lengths <- setNames(as.integer(fai$length), fai$chr)
  chr_lengths <- chr_lengths[keep_chr]
  cat("Chromosome lengths from FASTA index\n")
} else if (file.exists(genome_fasta)) {
  idx <- fasta.index(genome_fasta)
  chr_lengths <- setNames(as.integer(idx$seqlength), idx$desc)
  chr_lengths <- chr_lengths[names(chr_lengths) %in% keep_chr]
  cat("Chromosome lengths from Biostrings index\n")
} else {
  # Fallback: max coordinate per chromosome from GFF (slightly imprecise)
  cat("WARNING: No genome FASTA found — using GFF-derived chromosome lengths\n")
  chr_lengths <- tapply(end(gff), as.character(seqnames(gff)), max)
  chr_lengths <- as.integer(chr_lengths[keep_chr])
  names(chr_lengths) <- keep_chr
}
cat("Chromosomes with lengths:", length(chr_lengths), "\n")

# -- Gene info --
gene_strand <- as.character(strand(genes))
gene_chr    <- as.character(seqnames(genes))
gene_ids    <- genes$ID

# TSS and TES (strand-aware)
# + strand: TSS = start, TES = end
# - strand: TSS = end,   TES = start
tss <- ifelse(gene_strand == "+", start(genes), end(genes))
tes <- ifelse(gene_strand == "+", end(genes), start(genes))

chrlen <- as.integer(chr_lengths[gene_chr])

# -- Region 1: Promoter (2kb upstream of TSS) --
prom_s <- ifelse(gene_strand == "+",
                 pmax(1L, tss - 2000L),
                 tss + 1L)
prom_e <- ifelse(gene_strand == "+",
                 tss - 1L,
                 pmin(chrlen, tss + 2000L))

# -- Region 2: Upstream distal (2-10kb upstream, beyond promoter) --
up_s <- ifelse(gene_strand == "+",
               pmax(1L, tss - 10000L),
               tss + 2001L)
up_e <- ifelse(gene_strand == "+",
               pmax(1L, tss - 2001L),
               pmin(chrlen, tss + 10000L))

# -- Region 3: Gene body --
gb_s <- start(genes)
gb_e <- end(genes)

# -- Region 4: Downstream (10kb past TES) --
dn_s <- ifelse(gene_strand == "+",
               tes + 1L,
               pmax(1L, tes - 10000L))
dn_e <- ifelse(gene_strand == "+",
               pmin(chrlen, tes + 10000L),
               tes - 1L)

# -- Build GRanges (filter invalid: start > end at chr boundaries) --
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

regions <- c(
  build_gr(prom_s, prom_e, "promoter"),
  build_gr(up_s, up_e, "upstream_distal"),
  build_gr(gb_s, gb_e, "gene_body"),
  build_gr(dn_s, dn_e, "downstream")
)

seqlevels(regions) <- keep_chr
seqlengths(regions) <- chr_lengths[keep_chr]
regions <- trim(regions)

# -- Summary --
cat("\nExtended regions (10kb flanks):\n")
for (rt in c("promoter", "upstream_distal", "gene_body", "downstream")) {
  n <- sum(regions$region_type == rt)
  med_w <- median(width(regions[regions$region_type == rt]))
  cat(sprintf("  %-18s %6d regions  median width %6d bp\n", rt, n, med_w))
}
cat(sprintf("  %-18s %6d\n", "TOTAL", length(regions)))

# -- Save --
out_file <- file.path(cache_dir, "extended_regions_10kb.rds")
saveRDS(regions, out_file)
cat("\nSaved:", out_file, "\n")

# Copy to cluster cache if available
cluster_cache <- file.path(project_dir, "cluster/genome/cache")
if (dir.exists(cluster_cache)) {
  file.copy(out_file, file.path(cluster_cache, "extended_regions_10kb.rds"),
            overwrite = TRUE)
  cat("Copied to:", cluster_cache, "\n")
}

cat("Done!\n")
