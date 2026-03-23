#!/usr/bin/env Rscript
# =============================================================================
# Batch 06: Differential Methylation (DMP/DMR)
# Question: What changes during regeneration?
# Output: data/ (DMPs, DMRs, enrichment) + figures/ (8 plots)
# Requires: BSseq + DMLtest cache
# =============================================================================

source("methylation_pipeline/_config.R")

library(bsseq)
library(DSS)
library(GenomicRanges)
library(rtracklayer)
library(data.table)
library(ggplot2)
library(dplyr)
library(scales)
library(pheatmap)

BATCH_DIR <- file.path(PIPE_DIR, "batch06")
cat("=== Batch 06: Differential Methylation ===\n\n")

# --- Load cached DMLtest or run it ---
if (file.exists(CACHE$dmltest)) {
  cat("Loading cached DMLtest...\n")
  dml_test <- readRDS(CACHE$dmltest)
} else {
  cat("No DMLtest cache. Loading BSseq and running DMLtest (this takes hours)...\n")
  bs_obj <- readRDS(CACHE$bsseq)
  dml_test <- DMLtest(bs_obj, group1 = c("C1","C2"), group2 = c("A1","A2"), smoothing = TRUE)
  saveRDS(dml_test, CACHE$dmltest)
  rm(bs_obj); gc(verbose = FALSE)
}
cat(sprintf("DMLtest: %s sites\n", format(nrow(dml_test), big.mark = ",")))

# --- Call DMPs + DMRs ---
dmp <- callDML(dml_test, p.threshold = 0.05, delta = 0.1)
dmr <- callDMR(dml_test, p.threshold = 0.05, delta = 0.1, minlen = 50, minCG = 3)
if (is.null(dmr)) dmr <- data.frame(chr=character(), start=integer(), end=integer(),
                                      length=integer(), nCG=integer(), areaStat=numeric())

n_dmp <- nrow(dmp); n_dmr <- nrow(dmr)
n_hyper <- sum(dmp$diff > 0); n_hypo <- sum(dmp$diff < 0)
cat(sprintf("DMPs: %s (hyper: %s, hypo: %s)\n", format(n_dmp, big.mark=","),
            format(n_hyper, big.mark=","), format(n_hypo, big.mark=",")))
cat(sprintf("DMRs: %s\n", format(n_dmr, big.mark=",")))

# --- Annotate ---
gff <- load_gff(); genes <- gff[gff$type=="gene"]; exons <- gff[gff$type=="exon"]
promoters <- if(file.exists(CACHE$promoters)) readRDS(CACHE$promoters) else trim(promoters(genes, upstream=2000, downstream=0))

dmp$annotation <- annotate_regions(dmp$chr, dmp$pos, promoters, exons, genes)
dmp$direction <- ifelse(dmp$diff > 0, "Hyper", "Hypo")

# TE overlap (separate column)
te_data <- load_te()
dmp_gr <- GRanges(seqnames=dmp$chr, ranges=IRanges(start=dmp$pos, width=1))
dmp$in_te <- overlapsAny(dmp_gr, te_data$te_gr)

# Nearest gene
nearest_idx <- nearest(dmp_gr, genes)
dmp$nearest_gene <- ifelse(!is.na(nearest_idx), genes$ID[nearest_idx], NA)
dmp$dist_to_gene <- mcols(distanceToNearest(dmp_gr, genes))$distance

save_data(dmp, BATCH_DIR, "dmps_annotated")
save_data(as.data.frame(dmr), BATCH_DIR, "dmrs")

# --- Enrichment vs background ---
bs_obj <- readRDS(CACHE$bsseq)
all_anno <- annotate_regions(as.character(seqnames(bs_obj)), start(bs_obj), promoters, exons, genes)
bg_counts <- as.data.frame(table(all_anno)); colnames(bg_counts) <- c("Region","BG_Count")
dmp_counts <- as.data.frame(table(dmp$annotation)); colnames(dmp_counts) <- c("Region","DMP_Count")
enrichment <- merge(dmp_counts, bg_counts, by="Region")
enrichment$DMP_Pct <- round(100*enrichment$DMP_Count/sum(enrichment$DMP_Count),1)
enrichment$BG_Pct <- round(100*enrichment$BG_Count/sum(enrichment$BG_Count),1)
enrichment$Fold <- round(enrichment$DMP_Pct/enrichment$BG_Pct, 2)
enrichment$pvalue <- sapply(1:nrow(enrichment), function(i) {
  fisher.test(matrix(c(enrichment$DMP_Count[i], n_dmp-enrichment$DMP_Count[i],
                        enrichment$BG_Count[i], sum(bg_counts$BG_Count)-enrichment$BG_Count[i]), nrow=2))$p.value
})
enrichment$sig <- ifelse(enrichment$pvalue<0.001,"***",ifelse(enrichment$pvalue<0.01,"**",ifelse(enrichment$pvalue<0.05,"*","ns")))
save_data(enrichment, BATCH_DIR, "dmp_region_enrichment")
rm(bs_obj); gc(verbose=FALSE)

# --- Figures ---
# (A) PCA - load BSseq for this
bs_obj <- readRDS(CACHE$bsseq)
set.seed(42); pca_idx <- sort(sample(nrow(bs_obj), min(500000, nrow(bs_obj))))
beta_pca <- getMeth(bs_obj[pca_idx,], type="raw")
beta_pca <- beta_pca[complete.cases(beta_pca),]
pca <- prcomp(t(beta_pca), center=TRUE, scale.=FALSE)
ve <- summary(pca)$importance[2,1:2]*100
pca_df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Sample=rownames(pca$x),
                     Condition=c("Control","Control","Amputated","Amputated"))
p6a <- ggplot(pca_df, aes(x=PC1,y=PC2,color=Condition,label=Sample)) +
  geom_point(size=5) + geom_text(vjust=-1.2,size=4,show.legend=FALSE) +
  scale_color_manual(values=COLORS$condition) +
  labs(title="PCA of methylation profiles", x=sprintf("PC1 (%.1f%%)",ve[1]), y=sprintf("PC2 (%.1f%%)",ve[2])) +
  theme_minimal(base_size=13)
save_fig(p6a, BATCH_DIR, "fig6a_pca", w=8, h=6)
rm(beta_pca, pca, bs_obj); gc(verbose=FALSE)

# (D) Volcano
p6d <- ggplot(data.frame(diff=dmp$diff, pval=dmp$pval, dir=dmp$direction),
              aes(x=diff, y=-log10(pval), color=dir)) +
  geom_point(alpha=0.4, size=0.5) +
  scale_color_manual(values=COLORS$direction) +
  geom_vline(xintercept=c(-0.1,0.1), linetype="dashed", color="gray40") +
  labs(title="DMP Volcano Plot", x="Methylation Difference", y=expression(-log[10](p))) +
  theme_minimal(base_size=12)
save_fig(p6d, BATCH_DIR, "fig6d_dmp_volcano", w=10, h=7)

# (F) Annotation bar
p6f <- ggplot(as.data.frame(table(dmp$annotation)), aes(x=Var1, y=Freq, fill=Var1)) +
  geom_col() + scale_fill_manual(values=COLORS$region, guide="none") +
  labs(title="DMP genomic distribution", x=NULL, y="Count") +
  theme_minimal(base_size=12) + theme(axis.text.x=element_text(angle=45,hjust=1))
save_fig(p6f, BATCH_DIR, "fig6f_dmp_annotation", w=7, h=5)

# (G) Enrichment
enrich_long <- rbind(
  data.frame(Region=enrichment$Region, Pct=enrichment$DMP_Pct, Source="DMPs"),
  data.frame(Region=enrichment$Region, Pct=enrichment$BG_Pct, Source="All CpGs"))
p6g <- ggplot(enrich_long, aes(x=Region, y=Pct, fill=Source)) +
  geom_col(position="dodge", width=0.7) +
  scale_fill_manual(values=c("DMPs"="#E74C3C","All CpGs"="#95A5A6")) +
  geom_text(data=enrichment, aes(x=Region, y=pmax(DMP_Pct,BG_Pct)+1.5,
            label=paste0(Fold,"x ",sig)), inherit.aes=FALSE, size=3.5) +
  labs(title="DMP enrichment by region", x=NULL, y="% of positions") +
  theme_minimal(base_size=12) + theme(axis.text.x=element_text(angle=45,hjust=1))
save_fig(p6g, BATCH_DIR, "fig6g_dmp_enrichment", w=8, h=6)

rm(dml_test); gc(verbose=FALSE)
cat("\n=== Batch 06 complete ===\n")
