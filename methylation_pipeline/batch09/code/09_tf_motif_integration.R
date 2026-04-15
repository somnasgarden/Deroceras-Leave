#!/usr/bin/env Rscript
# =============================================================================
# Batch 09: TF + GENIE3 + Motif Integration
# Question: Do TF binding sites near DMPs explain regulatory connections?
#           Do methylation-concordant-with-expression genes share TF logic?
# Output: data/ (TSV) + figures/ (PNG+PDF)
# Requires: Batch 06 DMPs/DMRs, Batch 07 concordant genes, Batch 1.5 motif hits,
#           DeepTFactor, GENIE3, BSseq cache
# =============================================================================

source("methylation_pipeline/_config.R")
t0 <- proc.time()

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(GenomicRanges)
  library(IRanges)
})

BATCH_DIR <- file.path(PIPE_DIR, "batch09")
cat("=== Batch 09: TF + Motif + GENIE3 Integration ===\n\n")

unlink(file.path(BATCH_DIR, "figures"), recursive = TRUE)
unlink(file.path(BATCH_DIR, "data"), recursive = TRUE)
dir.create(file.path(BATCH_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BATCH_DIR, "data"), showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# Methylation-sensitivity TF list (Yin et al. 2017, Science)
# Curated by JASPAR tf_class + symbol overrides.
# Class-level priors reflect Yin's dominant SELEX patterns:
#   Homeodomain               -> MethylPlus  (often prefer methylated CpG)
#   bHLH (E-box CACGTG)       -> MethylMinus (CpG inside core; methylation disrupts)
#   bZIP (CRE TGACGTCA)       -> MethylMinus
#   ETS (GGAA, often CCGGAA)  -> MethylMinus
#   C2H2 ZF                   -> Mixed/Unknown (KLF/SP1 = MethylMinus; KRAB-ZFP = MethylPlus)
#   Nuclear receptor          -> LittleEffect
# Specific symbol overrides take precedence over class.
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

# Symbol-level overrides from Yin 2017 supplementary tables (high-confidence calls).
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
  yin_call = c("MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus",
               "MethylMinus","MethylMinus","MethylPlus","MethylPlus",
               "MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus",
               "MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus",
               "MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus",
               "MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus",
               "MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus",
               "MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus","MethylMinus",
               "MethylMinus","MethylMinus","MethylMinus",
               "MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus",
               "MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus",
               "MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus",
               "MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus","MethylPlus",
               "LittleEffect","LittleEffect","LittleEffect")
)

assign_methyl_sens <- function(motif_name, tf_class) {
  base <- methyl_sens_class(tf_class)
  sym  <- toupper(sub("\\(.*$", "", as.character(motif_name)))
  sym  <- sub("::.*$", "", sym)  # heterodimers like JUN::FOS -> JUN
  ov   <- yin_overrides[match(sym, yin_overrides$symbol), yin_call]
  ifelse(!is.na(ov), ov, base)
}

# =============================================================================
# Step 1: Load all inputs
# =============================================================================
cat("[1/8] Loading inputs...\n")

dmp_raw <- fread(file.path(PIPE_DIR, "batch06/data/dmps_annotated.tsv"))
dmp <- unique(dmp_raw, by = c("chr", "pos"))
cat(sprintf("  DMPs: %s unique positions\n", format(nrow(dmp), big.mark = ",")))

dmr <- fread(file.path(PIPE_DIR, "batch06/data/dmrs_annotated.tsv"))
dmr <- unique(dmr, by = c("chr", "start", "end"))
cat(sprintf("  DMRs: %s unique regions\n", format(nrow(dmr), big.mark = ",")))

tf_pred <- fread(OG$deeptf)
tf_pred[, gene_id := sub("-mRNA-.*$", "", sequence_ID)]
tf_genes <- unique(tf_pred[prediction == TRUE, gene_id])
cat(sprintf("  Predicted TFs: %d\n", length(tf_genes)))

# Annotation (headerless: gene_id, gene_name, description)
annot_dt <- if (file.exists(OG$annot)) {
  a <- fread(OG$annot, header = FALSE, col.names = c("gene_id", "gene_name", "description"))
  a[gene_name == "", gene_name := gene_id]
  a[, .SD[1], by = gene_id][, .(gene_id, gene_name)]
} else NULL
cat(sprintf("  Annotations: %s rows\n", if (is.null(annot_dt)) "0" else format(nrow(annot_dt), big.mark = ",")))

# Motif hits (promoter + upstream_distal)
motif_file <- file.path(PIPE_DIR, "batch1.5/data/motif_hits_extended.tsv.gz")
if (!file.exists(motif_file)) stop("Motif hits not found: ", motif_file)
cat("  Loading motif hits (promoter + upstream_distal)...\n")
motif <- fread(cmd = paste0("zcat ", shQuote(motif_file),
                            " | awk -F'\\t' 'NR==1 || $9==\"promoter\" || $9==\"upstream_distal\"'"))
cat(sprintf("  Motif hits: %s rows\n", format(nrow(motif), big.mark = ",")))

# Tag every motif with methylation sensitivity
motif[, methyl_sens := assign_methyl_sens(motif_name, tf_class)]

gff <- load_gff()
genes <- gff[gff$type == "gene"]
all_gene_ids <- genes$ID

# Concordant genes from batch07
concord_dmp <- fread(file.path(PIPE_DIR, "batch07/data/concordant_dmp_genes.tsv"))
concord_dmr <- fread(file.path(PIPE_DIR, "batch07/data/concordant_dmr_genes.tsv"))
concord_dmp_set <- unique(concord_dmp[concordant == TRUE, gene_id])
concord_dmr_set <- unique(concord_dmr[concordant == TRUE, gene_id])
cat(sprintf("  Concordant DMP genes: %d  DMR genes: %d\n",
            length(concord_dmp_set), length(concord_dmr_set)))

# Helper: attach gene_name
attach_name <- function(dt, key = "nearest_gene", out = "nearest_gene_name") {
  if (is.null(annot_dt)) { dt[[out]] <- dt[[key]]; return(dt) }
  m <- annot_dt[match(dt[[key]], annot_dt$gene_id), gene_name]
  m[is.na(m)] <- dt[[key]][is.na(m)]
  dt[[out]] <- m
  dt
}

# =============================================================================
# Step 2: Background universe = CpGs in motif-scannable regions
# Fix for Fisher artifact: motifs only exist in promoter+upstream_distal,
# so universe must be restricted to CpGs in those regions, not genome-wide.
# =============================================================================
cat("\n[2/8] Building region-restricted background universe...\n")

if (file.exists(CACHE$bsseq)) {
  bs <- readRDS(CACHE$bsseq)
  bs_gr <- granges(bs)
  cat(sprintf("  BSseq CpGs: %s\n", format(length(bs_gr), big.mark = ",")))
} else {
  stop("CACHE$bsseq not found; required for background fix.")
}

# Build motif-scannable region GRanges = promoter (2kb upstream of TSS) + 5kb upstream distal
gene_gr <- genes
prom_gr <- promoters(gene_gr, upstream = 2000, downstream = 0)
upstream_distal <- flank(gene_gr, width = 5000, start = TRUE)
upstream_distal <- shift(upstream_distal, ifelse(strand(upstream_distal) == "-", 2000, -2000))
scan_gr <- reduce(c(prom_gr, upstream_distal), ignore.strand = TRUE)
scan_gr <- scan_gr[seqnames(scan_gr) %in% keep_chr]

n_region_cpgs <- countOverlaps(scan_gr, bs_gr) |> sum()
cat(sprintf("  CpGs in motif-scannable regions: %s\n", format(n_region_cpgs, big.mark = ",")))

# Restrict DMPs to scannable regions for fair comparison
dmp_gr_full <- GRanges(dmp$chr, IRanges(dmp$pos, dmp$pos))
in_scan <- overlapsAny(dmp_gr_full, scan_gr)
dmp_scan <- dmp[in_scan]
cat(sprintf("  DMPs in scannable regions: %s / %s (%.1f%%)\n",
            format(nrow(dmp_scan), big.mark = ","),
            format(nrow(dmp), big.mark = ","),
            100 * nrow(dmp_scan) / nrow(dmp)))

# =============================================================================
# Step 3: DMP-motif overlap (with cpg_in_motif strict flag + gene names)
# =============================================================================
cat("\n[3/8] DMP-motif overlap (DMP +/-100bp; strict cpg_in_motif flag)...\n")

dmp_gr <- GRanges(dmp_scan$chr,
                  IRanges(pmax(1L, dmp_scan$pos - 100L), dmp_scan$pos + 100L))
mcols(dmp_gr)$dmp_idx <- seq_len(nrow(dmp_scan))
motif_gr <- GRanges(motif$chr, IRanges(motif$start, motif$end))
mcols(motif_gr)$motif_idx <- seq_len(nrow(motif))

ov <- findOverlaps(dmp_gr, motif_gr)
qh <- queryHits(ov); sh <- subjectHits(ov)
cat(sprintf("  Overlaps (within 100bp): %s\n", format(length(ov), big.mark = ",")))

dmp_motif <- data.table(
  chr           = dmp_scan$chr[qh],
  pos           = dmp_scan$pos[qh],
  diff          = dmp_scan$diff[qh],
  direction     = dmp_scan$direction[qh],
  annotation    = dmp_scan$annotation[qh],
  nearest_gene  = dmp_scan$nearest_gene[qh],
  motif_start   = motif$start[sh],
  motif_end     = motif$end[sh],
  motif_id      = motif$motif_id[sh],
  motif_name    = motif$motif_name[sh],
  tf_class      = motif$tf_class[sh],
  methyl_sens   = motif$methyl_sens[sh],
  motif_gene    = motif$gene_id[sh],
  region_type   = motif$region_type[sh],
  score         = motif$score[sh]
)
dmp_motif[, cpg_in_motif := pos >= motif_start & pos <= motif_end]
dmp_motif <- attach_name(dmp_motif, "nearest_gene", "nearest_gene_name")
save_data(dmp_motif, BATCH_DIR, "dmp_motif_overlaps")

n_dmp_with_motif <- uniqueN(dmp_motif[, paste(chr, pos)])
n_dmp_in_motif   <- uniqueN(dmp_motif[cpg_in_motif == TRUE, paste(chr, pos)])
cat(sprintf("  DMPs with any motif overlap (+/-100bp): %s\n", format(n_dmp_with_motif, big.mark = ",")))
cat(sprintf("  DMPs strictly inside a motif:           %s\n", format(n_dmp_in_motif, big.mark = ",")))

# =============================================================================
# Step 4: Motif enrichment at DMPs (Fisher's exact, region-restricted universe)
# =============================================================================
cat("\n[4/8] Motif enrichment at DMPs (region-restricted background)...\n")

dmp_per_motif <- dmp_motif[, .(n_dmp_hits = uniqueN(paste(chr, pos))),
                           by = .(motif_id, motif_name, tf_class, methyl_sens)]
motif_summary <- motif[, .(n_total_hits = .N),
                       by = .(motif_id, motif_name, tf_class, methyl_sens)]
enrich <- merge(motif_summary, dmp_per_motif,
                by = c("motif_id","motif_name","tf_class","methyl_sens"), all.x = TRUE)
enrich[is.na(n_dmp_hits), n_dmp_hits := 0L]

n_dmp_total <- nrow(dmp_scan)  # restricted to scannable regions

fisher_results <- enrich[, {
  a <- n_dmp_hits
  b <- n_dmp_total - a
  c <- max(0L, n_total_hits - a)
  d <- max(0L, n_region_cpgs - n_dmp_total - c)
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
  list(odds_ratio = as.numeric(ft$estimate), pvalue = ft$p.value,
       n_dmp_hits = a, n_total_hits = n_total_hits)
}, by = .(motif_id, motif_name, tf_class, methyl_sens)]
fisher_results[, fdr := p.adjust(pvalue, method = "BH")]
fisher_results[, log2OR := log2(pmax(odds_ratio, 1e-10))]
setorder(fisher_results, fdr, -odds_ratio)
save_data(fisher_results, BATCH_DIR, "motif_enrichment_at_dmps")
cat(sprintf("  Motifs tested: %d  significant (FDR<0.05): %d\n",
            nrow(fisher_results), sum(fisher_results$fdr < 0.05)))

# Top gene name per motif (by |meth diff|) for figure labels
top_gene_per_motif <- dmp_motif[order(-abs(diff)),
                                .(top_gene = first(nearest_gene_name)),
                                by = motif_id]
fisher_results <- merge(fisher_results, top_gene_per_motif, by = "motif_id", all.x = TRUE)

# =============================================================================
# Step 5: DMR-motif overlap + direction split
# =============================================================================
cat("\n[5/8] DMR-motif overlap (split by direction)...\n")

dmr_gr <- GRanges(dmr$chr, IRanges(dmr$start, dmr$end))
ov_dmr <- findOverlaps(dmr_gr, motif_gr)
qhd <- queryHits(ov_dmr); shd <- subjectHits(ov_dmr)

dmr_motif <- data.table(
  chr          = dmr$chr[qhd],
  start        = dmr$start[qhd],
  end          = dmr$end[qhd],
  direction    = dmr$direction[qhd],
  nearest_gene = dmr$nearest_gene[qhd],
  motif_id     = motif$motif_id[shd],
  motif_name   = motif$motif_name[shd],
  tf_class     = motif$tf_class[shd],
  methyl_sens  = motif$methyl_sens[shd],
  motif_gene   = motif$gene_id[shd],
  region_type  = motif$region_type[shd]
)
dmr_motif <- attach_name(dmr_motif, "nearest_gene", "nearest_gene_name")
save_data(dmr_motif, BATCH_DIR, "dmr_motif_overlaps")

dmr_fisher_split <- function(direction_label) {
  d <- dmr_motif[direction == direction_label]
  if (nrow(d) == 0) return(data.table())
  per <- d[, .(n_dmr_hits = uniqueN(paste(chr, start, end))),
           by = .(motif_id, motif_name, tf_class, methyl_sens)]
  m <- merge(motif_summary, per,
             by = c("motif_id","motif_name","tf_class","methyl_sens"), all.x = TRUE)
  m[is.na(n_dmr_hits), n_dmr_hits := 0L]
  n_dmr_total <- uniqueN(d[, paste(chr, start, end)])
  res <- m[, {
    a <- n_dmr_hits
    b <- n_dmr_total - a
    c <- max(0L, n_total_hits - a)
    d_ <- max(0L, n_region_cpgs - n_dmr_total - c)
    ft <- fisher.test(matrix(c(a, b, c, d_), nrow = 2), alternative = "greater")
    list(odds_ratio = as.numeric(ft$estimate), pvalue = ft$p.value,
         n_dmr_hits = a, n_total_hits = n_total_hits)
  }, by = .(motif_id, motif_name, tf_class, methyl_sens)]
  res[, fdr := p.adjust(pvalue, method = "BH")]
  res[, direction := direction_label]
  setorder(res, fdr, -odds_ratio)
  res
}
dmr_hyper <- dmr_fisher_split("Hyper")
dmr_hypo  <- dmr_fisher_split("Hypo")
dmr_fisher_all <- rbind(dmr_hyper, dmr_hypo)
save_data(dmr_fisher_all, BATCH_DIR, "motif_enrichment_at_dmrs_split")
cat(sprintf("  Hyper-DMR sig motifs: %d  |  Hypo-DMR sig motifs: %d\n",
            sum(dmr_hyper$fdr < 0.05), sum(dmr_hypo$fdr < 0.05)))

# =============================================================================
# Step 6: Concordant gene motif scan (your asks 2 + 3)
# =============================================================================
cat("\n[6/8] Concordant-gene motif enrichment + per-locus shortlist...\n")

# 6a. Motif enrichment restricted to concordant DMP genes
concord_motif <- motif[gene_id %in% concord_dmp_set]
cat(sprintf("  Motif hits in concordant DMP genes: %s\n",
            format(nrow(concord_motif), big.mark = ",")))

bg_motif <- motif  # universe = all motif-region motifs
concord_per <- concord_motif[, .(n_concord = .N),
                             by = .(motif_id, motif_name, tf_class, methyl_sens)]
bg_per <- bg_motif[, .(n_bg = .N),
                   by = .(motif_id, motif_name, tf_class, methyl_sens)]
cm <- merge(bg_per, concord_per,
            by = c("motif_id","motif_name","tf_class","methyl_sens"), all.x = TRUE)
cm[is.na(n_concord), n_concord := 0L]

n_concord_genes <- length(concord_dmp_set)
n_bg_genes <- uniqueN(motif$gene_id)

concord_enrich <- cm[, {
  a <- n_concord
  b <- n_bg - a
  c <- n_concord_genes
  d <- n_bg_genes - n_concord_genes
  # Hits per gene basis
  ft <- fisher.test(matrix(c(a, max(0L, b), c, max(0L, d)), nrow = 2),
                    alternative = "greater")
  list(odds_ratio = as.numeric(ft$estimate), pvalue = ft$p.value,
       n_concord = a, n_bg = n_bg)
}, by = .(motif_id, motif_name, tf_class, methyl_sens)]
concord_enrich[, fdr := p.adjust(pvalue, method = "BH")]
setorder(concord_enrich, fdr, -odds_ratio)
save_data(concord_enrich, BATCH_DIR, "concordant_motif_enrichment")
cat(sprintf("  Concordant-gene motifs sig (FDR<0.05): %d\n",
            sum(concord_enrich$fdr < 0.05)))

# 6b. Per-locus DMP-motif shortlist for concordant genes (the cloning table)
concord_loci <- dmp_motif[nearest_gene %in% concord_dmp_set]
# Merge in quadrant + LFC from batch07
concord_meta <- concord_dmp[concordant == TRUE,
                            .(gene_id, log2FoldChange, padj, quadrant, n_dmp)]
concord_loci <- merge(concord_loci, concord_meta,
                      by.x = "nearest_gene", by.y = "gene_id", all.x = TRUE)
setcolorder(concord_loci, c("nearest_gene","nearest_gene_name","quadrant",
                            "chr","pos","diff","direction","log2FoldChange","padj",
                            "motif_id","motif_name","tf_class","methyl_sens",
                            "cpg_in_motif","score","region_type"))
concord_loci[, abs_diff := abs(diff)]
setorder(concord_loci, quadrant, -abs_diff)
concord_loci[, abs_diff := NULL]
save_data(concord_loci, BATCH_DIR, "concordant_dmp_motif_loci")
cat(sprintf("  Concordant DMP-motif loci: %s rows (%d genes)\n",
            format(nrow(concord_loci), big.mark = ","),
            uniqueN(concord_loci$nearest_gene)))

# =============================================================================
# Step 7: TF enrichment + GENIE3 integration
# =============================================================================
cat("\n[7/8] TF enrichment + GENIE3 integration...\n")

dmp_genes <- unique(dmp$nearest_gene)
tf_with_dmp    <- sum(tf_genes %in% dmp_genes)
tf_without_dmp <- length(tf_genes) - tf_with_dmp
nontf_with_dmp <- length(setdiff(dmp_genes, tf_genes))
nontf_without_dmp <- length(all_gene_ids) - length(tf_genes) - nontf_with_dmp
ft_tf <- fisher.test(matrix(c(tf_with_dmp, tf_without_dmp,
                              nontf_with_dmp, nontf_without_dmp), nrow = 2))
cat(sprintf("  TF DMP enrichment: OR=%.2f, p=%s\n",
            ft_tf$estimate, format(ft_tf$p.value, digits = 3)))
save_data(data.table(
  test = "TF_enrichment_at_DMPs",
  OR = as.numeric(ft_tf$estimate), pvalue = ft_tf$p.value,
  tf_with_dmp = tf_with_dmp, tf_total = length(tf_genes)
), BATCH_DIR, "tf_dmp_enrichment")

expr <- NULL
if (file.exists(CACHE$transcriptome)) {
  txome <- readRDS(CACHE$transcriptome)
  expr <- as.data.table(txome$res_tail)
  if (!"gene_id" %in% names(expr)) expr[, gene_id := rownames(txome$res_tail)]
}

genie3_tf <- NULL
if (file.exists(OG$genie3)) {
  genie3 <- fread(OG$genie3)
  if (ncol(genie3) == 3) setnames(genie3, c("regulatoryGene","targetGene","weight"))
  if (nrow(genie3) > 500000) { setorder(genie3, -weight); genie3 <- genie3[1:500000] }
  genie3_tf <- genie3[regulatoryGene %in% tf_genes]
  target_with_motif <- unique(motif$gene_id)
  dmp_gene_set <- unique(dmp$nearest_gene)
  dmr_gene_set <- unique(dmr$nearest_gene)
  de_gene_set <- if (!is.null(expr)) unique(expr[!is.na(padj) & padj < 0.05, gene_id]) else character(0)
  genie3_tf[, has_motif := targetGene %in% target_with_motif]
  genie3_tf[, has_dmp   := targetGene %in% dmp_gene_set]
  genie3_tf[, has_dmr   := targetGene %in% dmr_gene_set]
  genie3_tf[, is_de     := targetGene %in% de_gene_set]
  genie3_tf[, full_evidence     := has_motif & has_dmp & is_de]
  genie3_tf[, full_evidence_dmr := has_motif & has_dmr & is_de]
  genie3_tf[, full_evidence_any := has_motif & (has_dmp | has_dmr) & is_de]
  if (!is.null(annot_dt)) {
    annot_unique <- unique(annot_dt, by = "gene_id")
    genie3_tf <- merge(genie3_tf, annot_unique, by.x = "regulatoryGene", by.y = "gene_id", all.x = TRUE)
    setnames(genie3_tf, "gene_name", "tf_name")
    genie3_tf <- merge(genie3_tf, annot_unique, by.x = "targetGene", by.y = "gene_id", all.x = TRUE)
    setnames(genie3_tf, "gene_name", "target_name")
  }
  if (!is.null(expr)) {
    expr_unique <- unique(expr[, .(gene_id, log2FoldChange, padj)], by = "gene_id")
    genie3_tf <- merge(genie3_tf, expr_unique,
                       by.x = "targetGene", by.y = "gene_id", all.x = TRUE)
  }
  save_data(genie3_tf, BATCH_DIR, "tf_target_evidence")
  cat(sprintf("  GENIE3 TF edges: %s  full_evidence: %d\n",
              format(nrow(genie3_tf), big.mark = ","), sum(genie3_tf$full_evidence)))
}

# =============================================================================
# Step 8: GO enrichment on concordant + motif-hit shortlist
# =============================================================================
cat("\n[8/8] GO enrichment on concordant DMP-motif shortlist...\n")
go_res_all <- NULL
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    file.exists(OG$string_enrich)) {
  string_enr <- fread(OG$string_enrich)
  # Use STRING TERM2GENE for clusterProfiler::enricher (org-agnostic)
  if ("#string_protein_id" %in% names(string_enr)) setnames(string_enr, "#string_protein_id", "string_id")
  # Map STRING IDs to gene_id via OG$string if available
  if (file.exists(OG$string)) {
    str_map <- fread(OG$string)
    if ("#string_protein_id" %in% names(str_map)) setnames(str_map, "#string_protein_id", "string_id")
    nm <- names(str_map)
    gid_col <- nm[grep("gene", nm, ignore.case = TRUE)][1]
    if (!is.na(gid_col)) {
      str_map <- str_map[, .(string_id, gene_id = get(gid_col))]
      string_enr <- merge(string_enr, str_map, by = "string_id", allow.cartesian = TRUE)
      term2gene <- unique(string_enr[, .(term, gene_id)])
      shortlist_genes <- unique(concord_loci$nearest_gene)
      bg <- unique(motif$gene_id)
      go_res_all <- tryCatch({
        clusterProfiler::enricher(gene = shortlist_genes, universe = bg,
                                  TERM2GENE = term2gene, pvalueCutoff = 0.1,
                                  qvalueCutoff = 0.2, minGSSize = 5)
      }, error = function(e) { cat("  enricher failed:", conditionMessage(e), "\n"); NULL })
      if (!is.null(go_res_all) && nrow(as.data.frame(go_res_all)) > 0) {
        save_data(as.data.table(go_res_all), BATCH_DIR, "concordant_shortlist_GO")
        cat(sprintf("  GO terms enriched: %d\n", nrow(as.data.frame(go_res_all))))
      } else cat("  No GO terms enriched.\n")
    }
  }
} else cat("  clusterProfiler or STRING enrichment not available — skipped.\n")

# =============================================================================
# Figures
# =============================================================================
cat("\n--- Figures ---\n")

ms_colors <- c(MethylPlus = "#C0392B", MethylMinus = "#2471A3",
               LittleEffect = "grey60", Mixed = "#7D6608", Unknown = "grey80")

# 9a: Top enriched motifs at DMPs (with top target gene name)
top30 <- fisher_results[fdr < 0.05][order(fdr)][1:min(30, sum(fisher_results$fdr < 0.05))]
if (nrow(top30) > 0) {
  top30[, label := paste0(motif_name, " - ", ifelse(is.na(top_gene), "", top_gene))]
  top30[, neg_log10_fdr := -log10(pmax(fdr, 1e-300))]
  top30[, label := make.unique(label)]   # avoid duplicate factor level error
  top30[, label := factor(label, levels = rev(label))]
  fig9a <- ggplot(top30, aes(x = neg_log10_fdr, y = label, fill = methyl_sens)) +
    geom_col() +
    scale_fill_manual(values = ms_colors, name = "Yin 2017\nclass") +
    labs(x = expression(-log[10](FDR)), y = NULL,
         title = "Top enriched TF motifs at DMPs (with top target gene)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_text(size = 8))
  save_fig(fig9a, BATCH_DIR, "fig9a_top_enriched_motifs_dmps", w = 11, h = 9)
}

# 9b: Volcano (gene labels for top hits)
volc <- copy(fisher_results)
volc[, neg_log10_p := pmin(-log10(pmax(pvalue, 1e-300)), 50)]
volc[, sig := fdr < 0.05]
top_label <- volc[sig == TRUE][order(-abs(log2OR))][1:min(15, sum(volc$sig))]
fig9b <- ggplot(volc, aes(x = log2OR, y = neg_log10_p)) +
  geom_point(aes(color = methyl_sens), alpha = 0.6, size = 1.6) +
  scale_color_manual(values = ms_colors, name = "Methyl\nsensitivity") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  labs(x = expression(log[2](Odds~Ratio)),
       y = expression(-log[10](p)~(capped~at~50)),
       title = "Motif enrichment at DMPs (region-restricted background)") +
  theme_minimal(base_size = 11)
if (nrow(top_label) > 0) {
  top_label[, lab := paste0(motif_name, "\n", ifelse(is.na(top_gene), "", top_gene))]
  fig9b <- fig9b + ggrepel::geom_text_repel(data = top_label, aes(label = lab),
                                            size = 2.4, max.overlaps = 25, seed = 42,
                                            min.segment.length = 0, box.padding = 0.5)
}
save_fig(fig9b, BATCH_DIR, "fig9b_motif_enrichment_volcano", w = 12, h = 8)

# 9c: TF class distribution (faceted by methyl_sens)
sig_motifs <- fisher_results[fdr < 0.05]
if (nrow(sig_motifs) > 0) {
  cc <- sig_motifs[, .N, by = .(tf_class, methyl_sens)][order(-N)]
  fig9c <- ggplot(cc, aes(x = reorder(tf_class, N), y = N, fill = methyl_sens)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = ms_colors) +
    labs(x = NULL, y = "Enriched motifs (FDR<0.05)",
         title = "TF classes enriched at DMPs, by methylation sensitivity") +
    theme_minimal(base_size = 10)
  save_fig(fig9c, BATCH_DIR, "fig9c_tf_class_methyl_sens", w = 11, h = 7)
}

# 9d: TF-target evidence (DMP vs DMR side-by-side)
if (!is.null(genie3_tf) && nrow(genie3_tf) > 0) {
  total <- nrow(genie3_tf)
  evidence_dt <- data.table(
    category = rep(c("Motif in promoter", "Meth change at target",
                     "Target is DE", "Full evidence"), each = 2),
    track    = rep(c("DMP", "DMR"), times = 4),
    count    = c(sum(genie3_tf$has_motif),       sum(genie3_tf$has_motif),
                 sum(genie3_tf$has_dmp),         sum(genie3_tf$has_dmr),
                 sum(genie3_tf$is_de),           sum(genie3_tf$is_de),
                 sum(genie3_tf$full_evidence),   sum(genie3_tf$full_evidence_dmr)),
    total = total)
  evidence_dt[, pct := 100 * count / total]
  evidence_dt[, category := factor(category, levels = rev(unique(category)))]
  fig9d <- ggplot(evidence_dt, aes(x = category, y = count, fill = track)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_text(aes(label = sprintf("%s (%.1f%%)", format(count, big.mark = ","), pct)),
              position = position_dodge(width = 0.75),
              hjust = -0.08, size = 3, color = "grey20") +
    scale_fill_manual(values = c(DMP = "#2471A3", DMR = "#C0392B"), name = NULL) +
    coord_flip() + scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
    labs(x = NULL, y = "TF-target edges",
         title = sprintf("Evidence for %s GENIE3 TF-target edges (DMP vs DMR)",
                         format(total, big.mark = ","))) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")
  save_fig(fig9d, BATCH_DIR, "fig9d_tf_target_evidence_summary", w = 10, h = 5)
  cat(sprintf("  full_evidence DMP: %d | DMR: %d | any: %d\n",
              sum(genie3_tf$full_evidence),
              sum(genie3_tf$full_evidence_dmr),
              sum(genie3_tf$full_evidence_any)))
}

# 9e: Top TF->target pairs with full evidence (gene names already used)
if (!is.null(genie3_tf) && sum(genie3_tf$full_evidence, na.rm = TRUE) > 0) {
  full_ev <- genie3_tf[full_evidence == TRUE]
  setorder(full_ev, -weight)
  top20_ev <- full_ev[1:min(20, nrow(full_ev))]
  top20_ev[, tf_label  := ifelse(is.na(tf_name) | tf_name == "", regulatoryGene, tf_name)]
  top20_ev[, tgt_label := ifelse(is.na(target_name) | target_name == "", targetGene, target_name)]
  top20_ev[, pair := factor(paste0(tf_label, " -> ", tgt_label),
                            levels = rev(paste0(tf_label, " -> ", tgt_label)))]
  fig9e <- ggplot(top20_ev, aes(x = weight, y = pair)) +
    geom_point(aes(color = log2FoldChange, size = abs(log2FoldChange))) +
    scale_color_gradient2(low = COLORS$direction["Hypo"], mid = "grey80",
                          high = COLORS$direction["Hyper"], midpoint = 0,
                          name = expression(log[2]~FC)) +
    scale_size_continuous(range = c(2, 6), name = "|log2FC|") +
    labs(x = "GENIE3 weight", y = NULL,
         title = "Top TF -> target pairs with full evidence") +
    theme_minimal(base_size = 11)
  save_fig(fig9e, BATCH_DIR, "fig9e_top_tf_target_pairs", w = 10, h = 7)
}

# 9f: DMR motif enrichment paired (Hyper vs Hypo) — with top DMR gene
top_dmr_split <- dmr_fisher_all[fdr < 0.05][order(fdr)]
if (nrow(top_dmr_split) > 0) {
  top_gene_per_dmr_motif <- dmr_motif[, .(top_gene = first(nearest_gene_name)),
                                      by = .(motif_id, direction)]
  top_dmr_split <- merge(top_dmr_split, top_gene_per_dmr_motif,
                         by = c("motif_id", "direction"), all.x = TRUE)
  top_motifs_each <- top_dmr_split[, head(.SD, 15), by = direction]
  top_motifs_each[, label := paste0(motif_name, " - ", ifelse(is.na(top_gene), "", top_gene))]
  top_motifs_each[, neg_log10_fdr := -log10(pmax(fdr, 1e-300))]
  fig9f <- ggplot(top_motifs_each, aes(x = neg_log10_fdr, y = reorder(label, neg_log10_fdr),
                                       fill = methyl_sens)) +
    geom_col() + facet_wrap(~ direction, scales = "free_y") +
    scale_fill_manual(values = ms_colors, name = "Yin 2017") +
    labs(x = expression(-log[10](FDR)), y = NULL,
         title = "DMR motif enrichment, split by direction") +
    theme_minimal(base_size = 10)
  save_fig(fig9f, BATCH_DIR, "fig9f_dmr_motif_enrichment_split", w = 12, h = 7)
}

# 9g: Concordant-gene motif enrichment with gene names
top_concord <- concord_enrich[fdr < 0.05][order(fdr)][1:min(30, sum(concord_enrich$fdr < 0.05))]
if (nrow(top_concord) > 0) {
  # Top gene per motif among concordant set
  top_gene_concord <- concord_loci[order(-abs(diff)),
                                   .(top_gene = first(nearest_gene_name)),
                                   by = motif_id]
  top_concord <- merge(top_concord, top_gene_concord, by = "motif_id", all.x = TRUE)
  top_concord[, label := paste0(motif_name, " - ", ifelse(is.na(top_gene), "", top_gene))]
  top_concord[, neg_log10_fdr := -log10(pmax(fdr, 1e-300))]
  top_concord[, label := make.unique(label)]   # avoid duplicate factor level error
  top_concord[, label := factor(label, levels = rev(label))]
  fig9g <- ggplot(top_concord, aes(x = neg_log10_fdr, y = label, fill = methyl_sens)) +
    geom_col() +
    scale_fill_manual(values = ms_colors) +
    labs(x = expression(-log[10](FDR)), y = NULL,
         title = "TF motifs enriched in methylation-concordant-with-expression genes") +
    theme_minimal(base_size = 10)
  save_fig(fig9g, BATCH_DIR, "fig9g_concordant_motif_enrichment", w = 11, h = 9)
}

# 9h: Per-locus shortlist heatmap (gene x TF, fill = meth_diff direction)
if (nrow(concord_loci) > 0) {
  hm <- concord_loci[cpg_in_motif == TRUE]
  if (nrow(hm) > 0) {
    top_genes <- hm[, .(score = sum(abs(diff))), by = nearest_gene_name][
      order(-score)][1:min(40, .N), nearest_gene_name]
    top_tfs <- hm[, .N, by = motif_name][order(-N)][1:min(25, .N), motif_name]
    hm_sub <- hm[nearest_gene_name %in% top_genes & motif_name %in% top_tfs]
    if (nrow(hm_sub) > 0) {
      hm_agg <- hm_sub[, .(meth_diff = mean(diff), quadrant = first(quadrant)),
                       by = .(nearest_gene_name, motif_name, methyl_sens)]
      fig9h <- ggplot(hm_agg, aes(x = motif_name, y = nearest_gene_name, fill = meth_diff)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = COLORS$direction["Hypo"], mid = "grey95",
                             high = COLORS$direction["Hyper"], midpoint = 0,
                             name = "Δβ") +
        labs(x = NULL, y = NULL,
             title = "Concordant DMP loci strictly inside TF motifs (top genes x top TFs)") +
        theme_minimal(base_size = 9) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      save_fig(fig9h, BATCH_DIR, "fig9h_concordant_locus_heatmap", w = 12, h = 10)
    }
  }
}

# 9i: Methylation sensitivity composition of significantly enriched motifs
if (sum(fisher_results$fdr < 0.05) > 0) {
  ms_dt <- fisher_results[fdr < 0.05, .N, by = methyl_sens]
  ms_dt[, pct := 100 * N / sum(N)]
  fig9i <- ggplot(ms_dt, aes(x = "", y = N, fill = methyl_sens)) +
    geom_col(width = 1, color = "white") + coord_polar(theta = "y") +
    scale_fill_manual(values = ms_colors) +
    labs(title = "Methylation sensitivity of TFs enriched at DMPs",
         x = NULL, y = NULL, fill = "Yin 2017") +
    theme_void(base_size = 11) +
    geom_text(aes(label = sprintf("%s\n%.0f%%", methyl_sens, pct)),
              position = position_stack(vjust = 0.5), size = 3)
  save_fig(fig9i, BATCH_DIR, "fig9i_methyl_sensitivity_pie", w = 7, h = 7)
}

# 9j: GO enrichment of concordant shortlist
if (!is.null(go_res_all) && nrow(as.data.frame(go_res_all)) > 0) {
  go_dt <- as.data.table(go_res_all)
  go_top <- go_dt[order(p.adjust)][1:min(20, .N)]
  go_top[, neg_log10_q := -log10(pmax(p.adjust, 1e-300))]
  # Map geneID strings to top 3 gene names per term for the label
  if (!is.null(annot_dt) && "geneID" %in% names(go_top)) {
    name_lookup <- setNames(annot_dt$gene_name, annot_dt$gene_id)
    go_top[, gene_label := sapply(strsplit(geneID, "/"), function(g) {
      nm <- name_lookup[g]; nm[is.na(nm)] <- g[is.na(nm)]
      paste(head(nm, 3), collapse = ", ")
    })]
    go_top[, label := paste0(Description, " (", gene_label, ")")]
  } else go_top[, label := Description]
  fig9j <- ggplot(go_top, aes(x = neg_log10_q, y = reorder(label, neg_log10_q))) +
    geom_col(fill = "#2471A3") +
    labs(x = expression(-log[10](q)), y = NULL,
         title = "GO/STRING enrichment: concordant DMP-motif shortlist (top genes shown)") +
    theme_minimal(base_size = 9)
  save_fig(fig9j, BATCH_DIR, "fig9j_concordant_shortlist_GO", w = 11, h = 7)
}

# =============================================================================
# Summary
# =============================================================================
cat("\n--- Batch 09 Summary ---\n")
cat(sprintf("DMPs in scannable regions: %s / %s\n",
            format(nrow(dmp_scan), big.mark = ","), format(nrow(dmp), big.mark = ",")))
cat(sprintf("Region-restricted background CpGs: %s\n", format(n_region_cpgs, big.mark = ",")))
cat(sprintf("Motifs sig at DMPs (FDR<0.05): %d\n", sum(fisher_results$fdr < 0.05)))
cat(sprintf("Concordant-gene motifs sig: %d\n", sum(concord_enrich$fdr < 0.05)))
cat(sprintf("Concordant DMP-motif loci (cpg_in_motif): %d\n",
            sum(concord_loci$cpg_in_motif, na.rm = TRUE)))
cat(sprintf("Hyper-DMR sig motifs: %d  Hypo-DMR sig motifs: %d\n",
            sum(dmr_hyper$fdr < 0.05), sum(dmr_hypo$fdr < 0.05)))
elapsed <- (proc.time() - t0)[["elapsed"]]
cat(sprintf("\n=== Batch 09 complete (%.1f s) ===\n", elapsed))
