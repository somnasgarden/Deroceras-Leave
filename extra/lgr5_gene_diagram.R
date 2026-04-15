#!/usr/bin/env Rscript
# =============================================================================
# LGR5 Gene Diagram — Differential methylation landscape
#
# NOT part of the pipeline — standalone exploratory figure.
# Shows gene structure, methylation levels, DMP positions, and TFBS motifs
# with DNA-binding domain (DBD) classification and cross-phyla conservation.
#
# Output: extra/lgr5_gene_diagram.png + .pdf
# =============================================================================

source("methylation_pipeline/_config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(bsseq)
  library(ggplot2)
  library(patchwork)
})

out_dir <- file.path(PROJECT_DIR, "extra")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== LGR5 Gene Diagram ===\n")

# =============================================================================
# 1. WHICH TF DBD FAMILIES HAVE CROSS-PHYLA CONSERVATION?
# =============================================================================
# The JASPAR motifs are derived from vertebrate/insect/nematode SELEX/ChIP.
# For a mollusk, only DBD families with CONSERVED DNA-recognition specificity
# are trustworthy. The key principle: the DBD fold determines which bases are
# contacted, and ancient DBD folds contact the same bases across metazoa.
#
# HIGH CONFIDENCE (DBD conserved across Bilateria, >500 My):
#   - Homeodomain (HOX, POU, PAX, NKX, LHX, etc.) — helix-turn-helix, recognition
#     helix contacts major groove at same positions. TFBS = same in fly, worm, human.
#   - bHLH (MYC, TWIST, HAND, ATOH) — E-box (CANNTG) universal across metazoa.
#   - bZIP (CREB, JUN, FOS, ATF, MAF) — leucine zipper, AP-1 site conserved.
#   - Forkhead/FOX — winged-helix domain, same recognition across phyla.
#   - ETS (ELK, FLI, SPI1) — winged-helix-turn-helix, core GGA(A/T) conserved.
#   - Nuclear receptors (RXRA, RARA, Nr2e3) — zinc finger + P-box, conserved.
#   - Runt/RUNX — immunoglobulin-fold, TGTGGT conserved.
#   - p53 — immunoglobulin-fold, conserved in bilaterians.
#   - Sox/HMG — HMG box, minor groove binder, AACAAT conserved.
#   - MADS/MEF2 — conserved in metazoa (but border case).
#   - T-box (TBX, Brachyury) — conserved DNA-binding.
#   - TCF/LEF (TCF7L2) — HMG box (Wnt pathway), conserved.
#
# LOW CONFIDENCE (lineage-specific or highly variable DBD):
#   - C2H2 zinc fingers (>3 adjacent) — KRAB/ZNF expansion is vertebrate-specific,
#     individual finger combinations are not conserved. Cannot trust specific motifs.
#   - KRAB-ZNF — entirely vertebrate-specific.
#   - Zinc finger "other" — too diverse.
#
# MODERATE (conserved fold but variable specificity):
#   - C2H2 zinc fingers (≤3 fingers, e.g., SP1/KLF, EGR) — Sp/KLF family is
#     ancient (GC-box), but many C2H2 are not.
#   - GATA (2-finger) — conserved in metazoa.

# DBD confidence lookup (maps JASPAR tf_class to conservation level)
dbd_confidence <- data.table(
  tf_class_pattern = c(
    "HOX|Homeo|homeo|POU|PAX|NK|LIM|TALE|Paired|CUT",
    "Helix-Loop-Helix|bHLH|Tal-related",
    "Leucine zipper|bZIP|Jun-related|Maf-related",
    "Forkhead|FOX|Fox",
    "Tryptophan|ETS|Ets",
    "receptor|NR1|NR2|NR3|NR4|Steroid",
    "Runt",
    "p53",
    "HMG|Sox|TCF-7|LEF",
    "T-box",
    "MADS",
    "GATA",
    "Rel|NF-kB",
    "More than 3 adjacent zinc|multiple dispersed zinc|KRAB",
    "C2H2|Zinc finger"
  ),
  dbd_family = c(
    "Homeodomain", "bHLH", "bZIP", "Forkhead", "ETS",
    "Nuclear receptor", "Runt", "p53", "HMG/Sox/TCF", "T-box",
    "MADS", "GATA", "Rel/NF-kB",
    "C2H2-ZF (multi)", "C2H2-ZF (other)"
  ),
  conservation = c(
    "High", "High", "High", "High", "High",
    "High", "High", "High", "High", "High",
    "Moderate", "Moderate", "Moderate",
    "Low", "Low"
  ),
  yin_class = c(
    "MethylPlus", "MethylMinus", "MethylMinus", "MethylMinus", "MethylMinus",
    "LittleEffect", "LittleEffect", "MethylMinus", "LittleEffect", "Mixed",
    "Unknown", "Mixed", "LittleEffect",
    "Mixed", "Mixed"
  )
)

assign_dbd <- function(tf_class_str) {
  for (i in seq_len(nrow(dbd_confidence))) {
    if (grepl(dbd_confidence$tf_class_pattern[i], tf_class_str, ignore.case = TRUE)) {
      return(dbd_confidence[i])
    }
  }
  data.table(tf_class_pattern = NA, dbd_family = "Unknown", conservation = "Unknown", yin_class = "Unknown")
}

# =============================================================================
# 2. LOAD DATA
# =============================================================================
cat("Loading data...\n")

# BSseq
bs <- readRDS(CACHE$bsseq)

# GFF
gff <- load_gff()
genes <- gff[gff$type == "gene"]
exons <- gff[gff$type == "exon"]

# DMPs
dmps <- fread(file.path(PIPE_DIR, "batch06/data/dmps_annotated.tsv"))
dmps <- unique(dmps, by = c("chr", "pos"))  # deduplicate

# Whole-genome motif hits (exact positions)
wg_chunk <- file.path(PROJECT_DIR, "cluster/results/whole_genome/_chunk_chr14.tsv")
has_wg <- file.exists(wg_chunk)
if (has_wg) {
  wg_motifs <- fread(wg_chunk)
  cat(sprintf("  Whole-genome motif hits on chr14: %s\n", format(nrow(wg_motifs), big.mark = ",")))
}

# =============================================================================
# 3. DEFINE THE TWO LGR5 LOCI
# =============================================================================

lgr5_loci <- list(
  list(
    gene_id = "LOC_00009623",
    gene_name = "LGR5",
    chr = "chr14",
    gene_start = 43985725, gene_end = 43991921, strand = "+",
    exons = data.table(start = c(43985725, 43987962, 43988940, 43990375),
                       end   = c(43986681, 43988074, 43989044, 43991921)),
    # Show wider region: gene + 3kb flanks
    view_start = 43982000, view_end = 43995000,
    dmp_positions = c(43986630),
    title = "LGR5 (LOC_00009623) — exonic DMP"
  ),
  list(
    gene_id = "LOC_00009602",
    gene_name = "LGR5",
    chr = "chr14",
    gene_start = 43335997, gene_end = 43341949, strand = "-",
    exons = data.table(start = c(43335997, 43340254, 43341846),
                       end   = c(43339083, 43340358, 43341949)),
    # Intergenic DMPs are 45kb upstream — show both gene and DMP region
    view_start = 43288000, view_end = 43345000,
    dmp_positions = c(43290398, 43290584),
    title = "LGR5 (LOC_00009602) — intergenic DMPs (45kb upstream)"
  )
)

# =============================================================================
# 4. BUILD FIGURE FOR EACH LOCUS
# =============================================================================

for (locus in lgr5_loci) {
  cat(sprintf("\n--- %s ---\n", locus$title))

  # --- Methylation data ---
  idx <- which(as.character(seqnames(bs)) == locus$chr &
               start(bs) >= locus$view_start &
               start(bs) <= locus$view_end)
  beta <- getMeth(bs[idx,], type = "raw")
  meth_dt <- data.table(
    pos = start(bs)[idx],
    C1 = beta[,1], C2 = beta[,2], A1 = beta[,3], A2 = beta[,4]
  )
  meth_dt[, ctrl := rowMeans(cbind(C1, C2), na.rm = TRUE)]
  meth_dt[, ampu := rowMeans(cbind(A1, A2), na.rm = TRUE)]
  meth_dt[, diff := ampu - ctrl]
  meth_dt[, is_dmp := pos %in% locus$dmp_positions]

  # --- DMP details ---
  dmp_dt <- dmps[chr == locus$chr & pos %in% locus$dmp_positions]
  cat(sprintf("  DMPs: %d\n", nrow(dmp_dt)))

  # --- Motifs at DMP positions ---
  motif_at_dmps <- NULL
  if (has_wg) {
    motif_list <- list()
    for (dmp_pos in locus$dmp_positions) {
      hits <- wg_motifs[start <= dmp_pos + 10 & end >= dmp_pos - 10]
      if (nrow(hits) > 0) {
        hits[, dmp_pos := dmp_pos]
        # Assign DBD family and conservation
        hits[, c("dbd_family", "conservation", "yin_class") := {
          res <- rbindlist(lapply(tf_class, function(tc) assign_dbd(tc)))
          list(res$dbd_family, res$conservation, res$yin_class)
        }]
        motif_list[[length(motif_list) + 1]] <- hits
      }
    }
    if (length(motif_list) > 0) {
      motif_at_dmps <- rbindlist(motif_list)
      # Keep only HIGH confidence DBD families
      motif_high <- motif_at_dmps[conservation == "High"]
      cat(sprintf("  Motifs at DMPs (total): %d | High-confidence DBD: %d\n",
                  nrow(motif_at_dmps), nrow(motif_high)))
      if (nrow(motif_high) > 0) {
        cat("  High-confidence motifs:\n")
        print(motif_high[, .(motif_name, tf_class, dbd_family, yin_class)])
      }
    }
  }

  # =================================================================
  # PANEL A: Gene structure
  # =================================================================
  gene_df <- data.frame(
    xmin = locus$gene_start, xmax = locus$gene_end,
    ymin = -0.15, ymax = 0.15
  )
  exon_df <- data.frame(
    xmin = locus$exons$start, xmax = locus$exons$end,
    ymin = -0.3, ymax = 0.3
  )

  # DMP markers
  dmp_marks <- data.frame(x = locus$dmp_positions, y = 0.5)

  # Motif annotations at DMPs (high confidence only)
  motif_labels <- NULL
  if (!is.null(motif_at_dmps)) {
    motif_high <- motif_at_dmps[conservation == "High"]
    if (nrow(motif_high) > 0) {
      # Deduplicate by DBD family per DMP
      motif_summary <- motif_high[, .(
        motifs = paste(unique(motif_name), collapse = ", "),
        n = .N
      ), by = .(dmp_pos, dbd_family, yin_class)]
      motif_labels <- motif_summary
    }
  }

  pa <- ggplot() +
    # Intron line
    geom_segment(aes(x = locus$gene_start, xend = locus$gene_end, y = 0, yend = 0),
                 linewidth = 1, color = "gray40") +
    # Gene body
    geom_rect(data = gene_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "gray80", color = "gray40") +
    # Exons
    geom_rect(data = exon_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#2C3E50", color = "#2C3E50") +
    # DMP positions
    geom_segment(data = dmp_marks, aes(x = x, xend = x, y = -0.6, yend = 0.5),
                 color = "#E74C3C", linewidth = 0.8, linetype = "dashed") +
    geom_point(data = dmp_marks, aes(x = x, y = 0.5),
               color = "#E74C3C", size = 3, shape = 25, fill = "#E74C3C") +
    # Strand arrow
    annotate("text", x = locus$gene_start - (locus$view_end - locus$view_start) * 0.02,
             y = 0, label = if (locus$strand == "+") "\u2192" else "\u2190",
             size = 6, color = "gray40") +
    # Gene name
    annotate("text", x = (locus$gene_start + locus$gene_end) / 2, y = -0.5,
             label = paste0(locus$gene_name, " (", locus$gene_id, ")"),
             size = 4, fontface = "italic") +
    scale_x_continuous(limits = c(locus$view_start, locus$view_end),
                       labels = function(x) paste0(round(x/1e6, 3), " Mb")) +
    ylim(-0.8, 1.5) +
    labs(x = NULL, y = NULL, title = locus$title) +
    theme_minimal(base_size = 11) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
          panel.grid.minor.x = element_blank())

  # Add motif labels at DMPs
  if (!is.null(motif_labels) && nrow(motif_labels) > 0) {
    # Stagger labels vertically
    motif_labels[, y_pos := 0.7 + (seq_len(.N) - 1) * 0.15, by = dmp_pos]
    # Color by Yin class
    yin_colors <- c("MethylPlus" = "#E74C3C", "MethylMinus" = "#3498DB",
                    "LittleEffect" = "#95A5A6", "Mixed" = "#F39C12", "Unknown" = "#BDC3C7")
    pa <- pa +
      geom_label(data = motif_labels,
                 aes(x = dmp_pos, y = y_pos,
                     label = paste0(dbd_family, ": ", motifs),
                     fill = yin_class),
                 size = 2.5, hjust = 0, nudge_x = (locus$view_end - locus$view_start) * 0.01,
                 label.size = 0.3, alpha = 0.85, show.legend = TRUE) +
      scale_fill_manual(values = yin_colors, name = "Methylation\nsensitivity")
  }

  # =================================================================
  # PANEL B: Methylation level (control vs amputated)
  # =================================================================
  meth_long <- rbindlist(list(
    meth_dt[, .(pos, beta = ctrl, condition = "Control")],
    meth_dt[, .(pos, beta = ampu, condition = "Amputated")]
  ))

  pb <- ggplot(meth_long, aes(x = pos, y = beta * 100, color = condition)) +
    geom_point(size = 0.5, alpha = 0.4) +
    geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 1) +
    # Mark DMPs
    geom_vline(xintercept = locus$dmp_positions, color = "#E74C3C",
               linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(values = c("Control" = "#2471A3", "Amputated" = "#C0392B")) +
    scale_x_continuous(limits = c(locus$view_start, locus$view_end),
                       labels = function(x) paste0(round(x/1e6, 3), " Mb")) +
    labs(x = NULL, y = "Methylation (%)", color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "top")

  # =================================================================
  # PANEL C: Differential methylation (amputated - control)
  # =================================================================
  pc <- ggplot(meth_dt, aes(x = pos, y = diff * 100)) +
    geom_point(aes(color = is_dmp), size = ifelse(meth_dt$is_dmp, 3, 0.5),
               alpha = ifelse(meth_dt$is_dmp, 1, 0.3)) +
    geom_smooth(method = "loess", span = 0.15, se = TRUE, color = "gray30", linewidth = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_vline(xintercept = locus$dmp_positions, color = "#E74C3C",
               linetype = "dashed", linewidth = 0.5) +
    scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "#E74C3C"),
                       labels = c("CpG", "DMP"), name = NULL) +
    scale_x_continuous(limits = c(locus$view_start, locus$view_end),
                       labels = function(x) paste0(round(x/1e6, 3), " Mb")) +
    labs(x = paste0(locus$chr, " position"), y = "Delta methylation (%)\n(Amputated - Control)") +
    theme_minimal(base_size = 11)

  # =================================================================
  # COMBINE
  # =================================================================
  combined <- pa / pb / pc + plot_layout(heights = c(2, 1.5, 1.5))

  # Adjust height based on motif labels
  fig_h <- if (!is.null(motif_labels) && nrow(motif_labels) > 3) 14 else 12

  ggsave(file.path(out_dir, sprintf("lgr5_%s_diagram.png", locus$gene_id)),
         combined, width = 14, height = fig_h, dpi = 300)
  ggsave(file.path(out_dir, sprintf("lgr5_%s_diagram.pdf", locus$gene_id)),
         combined, width = 14, height = fig_h, device = cairo_pdf)
  cat(sprintf("  Saved: lgr5_%s_diagram.png/pdf\n", locus$gene_id))
}

# =============================================================================
# 5. MOTIF DBD CONSERVATION TABLE (for paper methods)
# =============================================================================
cat("\n=== TF DBD Conservation Table ===\n")
cat("Families with HIGH cross-phyla conservation (safe for mollusk):\n\n")
cat(sprintf("%-20s %-15s %s\n", "DBD Family", "Yin 2017", "Conservation basis"))
cat(sprintf("%-20s %-15s %s\n", "---", "---", "---"))
high_conf <- dbd_confidence[conservation == "High"]
notes <- c(
  "Homeodomain"      = "HTH fold, recognition helix contacts identical bases in all bilaterians",
  "bHLH"             = "E-box (CANNTG) universal; TWIST, HAND, MYC conserved >600 My",
  "bZIP"             = "Leucine zipper, AP-1 (TGA[CG]TCA) identical fly→human",
  "Forkhead"         = "Winged-helix, core TRTTKRY conserved across metazoa",
  "ETS"              = "Winged-HTH, GGA(A/T) core conserved; present in all metazoa",
  "Nuclear receptor" = "Zn finger DBD + P-box, hormone response elements conserved",
  "Runt"             = "Ig-fold DBD, TGTGGT conserved; Drosophila Lozenge = human RUNX",
  "p53"              = "Ig-fold, conserved in bilaterians (fly p53 binds same motif)",
  "HMG/Sox/TCF"      = "HMG box, minor groove binder, AACAAT; TCF/LEF = Wnt pathway",
  "T-box"            = "T-domain, half-palindrome conserved across deuterostomes"
)
for (i in seq_len(nrow(high_conf))) {
  fam <- high_conf$dbd_family[i]
  note <- if (fam %in% names(notes)) notes[fam] else ""
  cat(sprintf("%-20s %-15s %s\n", fam, high_conf$yin_class[i], note))
}

cat("\n\nFamilies with LOW cross-phyla conservation (DO NOT trust for mollusk):\n")
cat("  - C2H2 zinc fingers (>3 adjacent): KRAB expansion is vertebrate-specific\n")
cat("  - Zinc fingers 'multiple dispersed': finger combinations are lineage-specific\n")
cat("  - Any TF class = 'Unknown' or empty in JASPAR\n")

# =============================================================================
# 6. SUMMARY TABLE: LGR5 DMPs with motif annotations
# =============================================================================
cat("\n=== LGR5 DMP Motif Summary ===\n")

if (has_wg) {
  all_lgr5_motifs <- list()
  for (dmp_pos in c(43290398, 43290584, 43986630)) {
    hits <- wg_motifs[start <= dmp_pos + 10 & end >= dmp_pos - 10]
    if (nrow(hits) > 0) {
      hits[, dmp_pos := dmp_pos]
      hits[, c("dbd_family", "conservation", "yin_class") := {
        res <- rbindlist(lapply(tf_class, function(tc) assign_dbd(tc)))
        list(res$dbd_family, res$conservation, res$yin_class)
      }]
      all_lgr5_motifs[[length(all_lgr5_motifs) + 1]] <- hits
    }
  }
  all_motifs <- rbindlist(all_lgr5_motifs)

  # Print sorted by conservation
  setorder(all_motifs, dmp_pos, -conservation, dbd_family)
  summary_table <- all_motifs[, .(dmp_pos, motif_name, tf_class, dbd_family, conservation, yin_class)]
  cat("\nAll motifs overlapping LGR5 DMPs (HIGH confidence marked):\n\n")
  for (dp in unique(summary_table$dmp_pos)) {
    cat(sprintf("--- DMP at chr14:%d ---\n", dp))
    sub <- summary_table[dmp_pos == dp]
    for (j in seq_len(nrow(sub))) {
      marker <- if (sub$conservation[j] == "High") " ***" else ""
      cat(sprintf("  %-12s %-40s %-20s %-15s%s\n",
                  sub$motif_name[j], sub$tf_class[j],
                  sub$dbd_family[j], sub$yin_class[j], marker))
    }
  }

  fwrite(summary_table, file.path(out_dir, "lgr5_dmp_motif_table.tsv"), sep = "\t")
  cat(sprintf("\nSaved: %s\n", file.path(out_dir, "lgr5_dmp_motif_table.tsv")))
}

cat("\n=== Done ===\n")
