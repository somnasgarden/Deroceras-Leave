#!/usr/bin/env Rscript
# =============================================================================
# Annotate existing whole-genome motif scan chunk files
# Runs steps 4-7 from whole_genome_motif_scan.R on pre-scanned data
# =============================================================================

source("methylation_pipeline/_config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

out_dir <- file.path(PROJECT_DIR, "cluster/results/whole_genome")
cat("=== Annotating whole-genome motif hits ===\n")
cat("Chunk dir:", out_dir, "\n\n")
t0 <- proc.time()[3]

# =============================================================================
# 1. Load chunk files
# =============================================================================
chunk_files <- list.files(out_dir, pattern = "^_chunk_chr.*\\.tsv$", full.names = TRUE)
cat(sprintf("Found %d chunk files\n", length(chunk_files)))
stopifnot(length(chunk_files) > 0)

all_hits <- rbindlist(lapply(chunk_files, fread, sep = "\t"))
cat(sprintf("Combined: %s hits\n", format(nrow(all_hits), big.mark = ",")))

# =============================================================================
# 2. Load annotation data
# =============================================================================
cat("\nLoading GFF + TEs...\n")
gff <- load_gff()
genes <- gff[gff$type == "gene"]
genes <- genes[as.character(seqnames(genes)) %in% keep_chr]
exons <- gff[gff$type == "exon"]
exons <- exons[as.character(seqnames(exons)) %in% keep_chr]

genome <- load_genome()
chr_lengths <- setNames(width(genome), names(genome))
rm(genome); gc(verbose = FALSE)

# Promoters
gene_strand <- as.character(strand(genes))
gene_tss <- ifelse(gene_strand == "+", start(genes), end(genes))
prom_start <- ifelse(gene_strand == "+", pmax(1L, gene_tss - 2000L), gene_tss + 1L)
prom_end   <- ifelse(gene_strand == "+", gene_tss - 1L,
                     pmin(as.integer(chr_lengths[as.character(seqnames(genes))]), gene_tss + 2000L))
promoters <- GRanges(seqnames = seqnames(genes),
                     ranges = IRanges(start = prom_start, end = prom_end),
                     gene_id = genes$ID)

# TSS lookup
tss_lookup <- data.table(gene_id = genes$ID, chr = as.character(seqnames(genes)),
                         tss = gene_tss, strand = gene_strand)

# TEs
te <- load_te()
te_gr <- te$te_gr
if (is.null(te_gr) || length(te_gr) == 0) {
  cat("  WARNING: TE data empty — skipping TE annotation\n")
  te_gr <- GRanges()
}
cat(sprintf("  Genes: %d | Exons: %d | TEs: %d\n", length(genes), length(exons), length(te_gr)))

# =============================================================================
# 3. Annotate
# =============================================================================
cat("\nAnnotating...\n")

hit_gr <- GRanges(seqnames = all_hits$chr,
                  ranges = IRanges(start = all_hits$start, end = all_hits$end))

cat("  Overlapping with promoters...\n")
ov_prom <- overlapsAny(hit_gr, promoters)
cat("  Overlapping with exons...\n")
ov_exon <- overlapsAny(hit_gr, exons)
cat("  Overlapping with genes...\n")
ov_gene <- overlapsAny(hit_gr, genes)
cat("  Overlapping with TEs...\n")
ov_te   <- overlapsAny(hit_gr, te_gr)

all_hits[, region := fifelse(ov_prom, "Promoter",
                    fifelse(ov_exon, "Exon",
                    fifelse(ov_gene, "Intron", "Intergenic")))]
all_hits[, in_te := ov_te]

# Nearest gene + distance to TSS
cat("  Computing distance to nearest gene TSS...\n")
hit_mid <- (all_hits$start + all_hits$end) %/% 2L
hit_gr_mid <- GRanges(seqnames = all_hits$chr,
                      ranges = IRanges(start = hit_mid, width = 1))
gene_gr <- GRanges(seqnames = tss_lookup$chr,
                   ranges = IRanges(start = tss_lookup$tss, width = 1),
                   gene_id = tss_lookup$gene_id,
                   gene_strand = tss_lookup$strand)
nearest_idx <- nearest(hit_gr_mid, gene_gr)

all_hits[, nearest_gene := gene_gr$gene_id[nearest_idx]]
all_hits[, nearest_gene_strand := gene_gr$gene_strand[nearest_idx]]
all_hits[, nearest_tss := start(gene_gr)[nearest_idx]]
all_hits[, dist_to_nearest_tss := fifelse(nearest_gene_strand == "+",
                                          hit_mid - nearest_tss,
                                          nearest_tss - hit_mid)]
all_hits[, c("nearest_gene_strand", "nearest_tss") := NULL]

rm(hit_gr, hit_gr_mid, gene_gr, ov_prom, ov_exon, ov_gene, ov_te, nearest_idx, hit_mid)
gc(verbose = FALSE)

# =============================================================================
# 4. Save output
# =============================================================================
cat("\nSaving...\n")
out_file <- file.path(out_dir, "whole_genome_motif_hits.tsv.gz")
fwrite(all_hits, out_file, sep = "\t", compress = "gzip")
cat(sprintf("Saved: %s (%s hits, %.1f MB)\n",
            out_file, format(nrow(all_hits), big.mark = ","),
            file.size(out_file) / 1e6))

# Summaries
region_summary <- all_hits[, .(total_hits = .N, unique_motifs = uniqueN(motif_id),
                                in_te_pct = 100 * sum(in_te) / .N), by = region]
print(region_summary)
fwrite(region_summary, file.path(out_dir, "region_summary.tsv"), sep = "\t")

motif_summary <- all_hits[, .(total_hits = .N, n_promoter = sum(region == "Promoter"),
  n_exon = sum(region == "Exon"), n_intron = sum(region == "Intron"),
  n_intergenic = sum(region == "Intergenic"), n_in_te = sum(in_te)), by = .(motif_id, motif_name, tf_class)]
fwrite(motif_summary, file.path(out_dir, "motif_summary_whole_genome.tsv"), sep = "\t")

elapsed <- (proc.time()[3] - t0) / 60
cat(sprintf("\n=== Annotation complete (%.1f min) ===\n", elapsed))
