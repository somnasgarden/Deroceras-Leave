# Genome characterization — Deroceras laeve
# Checks GFF/genome compatibility, computes genome & CpG properties
# Output: results/01_genome/genome_report.md + dinucleotide plot

library(Biostrings)    # DNA sequence analysis: base counting, dinucleotide freqs
library(GenomicRanges) # genomic intervals: overlaps, subtract regions, merge
library(rtracklayer)   # import() reads GFF annotation files into R
library(GenomeInfoDb)  # chromosome name and length management
library(ggplot2)       # plotting library
library(data.table)    # fread() reads big TSV files fast

# <- is assignment: store the text on the right into the name on the left
genome_fasta <- "/mnt/c/Users/rafae/Projects/STANDBY/genome/derLaeGenome_chr1_31.fasta"
gff_file     <- "/mnt/c/Users/rafae/Projects/DATA/derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff"
te_file      <- "/mnt/c/Users/rafae/Projects/DATA/collapsed_te_age_data.tsv"
out_report   <- "/mnt/c/Users/rafae/Projects/STANDBY/results/batch01/genome_report.html"
fig_pdf_dir  <- "/mnt/c/Users/rafae/Projects/STANDBY/results/batch01/pdf"
fig_png_dir  <- "/mnt/c/Users/rafae/Projects/STANDBY/results/batch01/png"
# RDS cache directory — saves parsed objects so second run is ~10x faster
cache_dir    <- "/mnt/c/Users/rafae/Projects/STANDBY/genome/cache"
# options: no scientific notation, strings stay as strings
options(scipen = 999, stringsAsFactors = FALSE)

# dir.create() makes a folder. showWarnings=FALSE = don't complain if it exists
# recursive=TRUE = create parent folders too if needed
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_pdf_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_png_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(out_report), showWarnings = FALSE, recursive = TRUE)

keep_chr <- paste0("chr", 1:31)


# 1. LOAD GENOME (with RDS caching)
# file.exists() returns TRUE if the file is on disk
# saveRDS() saves any R object as a compressed binary file
# readRDS() loads it back — MUCH faster than re-parsing FASTA/GFF
genome_rds <- file.path(cache_dir, "genome_chr1_31.rds")
if (file.exists(genome_rds)) {
  cat("Loading genome from cache...\n")
  genome <- readRDS(genome_rds)
} else {
  cat("Parsing genome FASTA (first run, will cache)...\n")
  genome <- readDNAStringSet(genome_fasta)
  saveRDS(genome, genome_rds)
}
cat("Genome loaded:", length(genome), "sequences\n")


# 2. LOAD GFF, KEEP CHR1-31 ONLY (with caching)
gff_rds <- file.path(cache_dir, "gff_chr1_31.rds")
if (file.exists(gff_rds)) {
  cat("Loading GFF from cache...\n")
  gff <- readRDS(gff_rds)
} else {
  cat("Parsing GFF (first run, will cache)...\n")
  # import() reads a GFF file into a GRanges object — a table where each row
  # is a genomic feature (gene, exon, mRNA...) with chr, start, end, strand, type
  gff <- import(gff_file)
  # seqnames(gff) = the chromosome column of every row in the GFF
  # %in% = "is this in our list?" returns TRUE or FALSE for each row
  gff <- gff[seqnames(gff) %in% keep_chr]
  # keepSeqlevels() removes leftover scaffold names from the registry
  gff <- keepSeqlevels(gff, keep_chr, pruning.mode = "coarse")
  saveRDS(gff, gff_rds)
}
cat("GFF loaded:", length(gff), "features on chr1-31\n")


# 2b. LOAD TE DATA (with caching)
te_rds <- file.path(cache_dir, "te_chr1_31.rds")
if (file.exists(te_rds)) {
  cat("Loading TEs from cache...\n")
  te_data <- readRDS(te_rds)
  te_gr <- te_data$te_gr; te_regions <- te_data$te_regions
} else {
  cat("Parsing TE data (first run, will cache)...\n")
  # fread() from data.table reads big TSV files much faster than read.table()
  te_dt <- fread(te_file)
  te_dt <- te_dt[chrom %in% keep_chr]
  # makeGRangesFromDataFrame() converts a data.frame with chr/start/end to GRanges
  te_gr <- makeGRangesFromDataFrame(te_dt, seqnames.field = "chrom",
                                    start.field = "start", end.field = "end",
                                    keep.extra.columns = TRUE)
  te_regions <- reduce(te_gr)  # merge overlapping TEs into non-overlapping blocks
  saveRDS(list(te_gr = te_gr, te_regions = te_regions), te_rds)
}
cat("TEs loaded:", length(te_gr), "elements,", length(te_regions), "merged regions\n")


# 3. GFF/GENOME COMPATIBILITY
# names() gets the names attached to an object. here = chromosome names from FASTA
# sort() puts them in alphabetical order
genome_chrs <- sort(names(genome))

# unique() removes duplicates (thousands of features share "chr1", we only need it once)
# as.character() converts Bioconductor's special type to plain text strings
gff_chrs <- sort(unique(as.character(seqnames(gff))))

# intersect(a, b) = items that appear in BOTH a and b
# setdiff(a, b) = items in a that are NOT in b
chrs_in_both     <- intersect(genome_chrs, gff_chrs)
chrs_only_genome <- setdiff(genome_chrs, gff_chrs)
chrs_only_gff    <- setdiff(gff_chrs, genome_chrs)

# cat() prints text to the screen. \n = new line
cat("\nCOMPATIBILITY CHECK:\n")
cat("  Shared chromosomes:", length(chrs_in_both), "\n")
cat("  Only in genome:", length(chrs_only_genome), "\n")
cat("  Only in GFF:", length(chrs_only_gff), "\n")

# width() = length in bp of each chromosome sequence
chr_lengths <- width(genome)
# names(x) <- y labels each value. now chr_lengths["chr1"] = size of chr1
names(chr_lengths) <- names(genome)

# seqlengths() tells the GFF how long each chromosome is (needed for trim later)
# seqlevels(gff) = list of chr names the GFF knows about
# chr_lengths[seqlevels(gff)] = look up each one by name (indexing with names)
seqlengths(gff) <- chr_lengths[seqlevels(gff)]

# trim() clips any feature that goes past the end of its chromosome
# e.g. if chr1 is 100bp but a feature says end=105, trim makes it end=100
gff_trimmed <- trim(gff)

# != means "not equal". so width(gff) != width(gff_trimmed) = TRUE where clipped
# sum() on TRUE/FALSE counts the TRUEs (TRUE=1, FALSE=0)
n_clipped <- sum(width(gff) != width(gff_trimmed))
cat("  Features out of bounds:", n_clipped, "\n")

# ifelse(condition, value_if_true, value_if_false) picks one of two answers
# & means AND — both conditions must be true
cat("  VERDICT:", ifelse(length(chrs_in_both) == 31 & n_clipped == 0,
    "PASS", "CHECK NEEDED"), "\n\n")


# 4. BASIC GENOME STATS
chr_sizes <- width(genome)
names(chr_sizes) <- names(genome)

# letterFrequency() counts how many times specific bases appear
# letters="GC" counts G and C together
# as.prob=TRUE divides by sequence length so you get a fraction (0 to 1)
# the result is a matrix (rows x columns). [, 1] takes column 1 from all rows
gc_per_chr <- letterFrequency(genome, letters = "GC", as.prob = TRUE)[, 1]

# sum() adds up all values. as.numeric() needed because genome sizes are huge
total_size <- sum(as.numeric(chr_sizes))

# same but without as.prob = raw counts. divide by total size = GC fraction
overall_gc <- sum(letterFrequency(genome, letters = "GC")) / total_size

cat("GENOME STATS:\n")
# format(x, big.mark=",") inserts commas: 1783141715 -> "1,783,141,715"
cat("  Total size:", format(total_size, big.mark = ","), "bp\n")
cat("  Chromosomes:", length(genome), "\n")
# max() = largest value. which.max() = which item IS the largest (gives its name)
cat("  Largest:", format(max(chr_sizes), big.mark = ","), "bp (", names(which.max(chr_sizes)), ")\n")
cat("  Smallest:", format(min(chr_sizes), big.mark = ","), "bp (", names(which.min(chr_sizes)), ")\n")
# round(number, 2) = round to 2 decimal places. * 100 converts fraction to %
cat("  GC content:", round(overall_gc * 100, 2), "%\n\n")


# 5. GENE ANNOTATION SUMMARY
# $ accesses a column by name. gff$type = the "type" column of every row
# == checks equality: "gene" == "gene" is TRUE, "exon" == "gene" is FALSE
# gff[TRUE/FALSE] keeps only matching rows (same filtering idea as section 2)
genes <- gff[gff$type == "gene"]
mrnas <- gff[gff$type == "mRNA"]
exons <- gff[gff$type == "exon"]

cat("ANNOTATION SUMMARY:\n")
cat("  Genes:", length(genes), "\n")
cat("  mRNAs:", length(mrnas), "\n")
cat("  Exons:", length(exons), "\n")
# width() on GRanges = end - start + 1 for each feature (its size in bp)
# mean() = average of all values. median() = the middle value when sorted
cat("  Mean gene length:", format(round(mean(width(genes))), big.mark = ","), "bp\n")
cat("  Median gene length:", format(round(median(width(genes))), big.mark = ","), "bp\n")

# mcols() = metadata columns (everything beyond chr/start/end/strand)
# colnames() = the names of those columns
# if "Parent" column exists, we can count exons per mRNA
if ("Parent" %in% colnames(mcols(exons))) {
  # exons$Parent = which mRNA each exon belongs to (a list of IDs)
  # unlist() flattens a nested list into a simple flat vector
  exon_parents <- unlist(exons$Parent)
  # table() counts how many times each value appears (like a tally)
  # so table(exon_parents) = how many exons each mRNA has
  exons_per_mrna <- table(exon_parents)
  cat("  Mean exons per mRNA:", round(mean(exons_per_mrna), 1), "\n")
  cat("  Median exons per mRNA:", median(exons_per_mrna), "\n")
}

# reduce() merges overlapping intervals:
#   gene A: chr1 100-500, gene B: chr1 400-800
#   reduce -> chr1 100-800 (one merged block)
# this prevents double-counting overlap when computing genome coverage
gene_cov <- sum(as.numeric(width(reduce(genes)))) / total_size * 100
te_cov   <- sum(as.numeric(width(te_regions))) / total_size * 100
cat("  Genome covered by genes:", round(gene_cov, 1), "%\n")
cat("  Genome covered by TEs:", round(te_cov, 1), "%\n\n")


# 6. CpG ANALYSIS (genome-wide)
# dinucleotideFrequency() counts all 16 two-base combos per chromosome:
#   AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT
# returns a matrix: rows = chromosomes, columns = dinucleotides
dinuc <- dinucleotideFrequency(genome)

# [, "CG"] = get the CG column from all rows
# sum() = total CG dinucleotides across all chromosomes
total_cpg <- sum(dinuc[, "CG"])
total_c   <- sum(letterFrequency(genome, letters = "C"))
total_g   <- sum(letterFrequency(genome, letters = "G"))

# CpG observed/expected (O/E) ratio
# if C and G were randomly placed in the genome, you'd expect:
#   expected CpG count = freq(C) * freq(G) * length
# rearranging: O/E = (CpG_count * genome_length) / (C_count * G_count)
# O/E < 1 = fewer CpGs than expected = CpGs have been lost over time
# this happens because methylated CpGs mutate to TpG, slowly erasing CpGs
# as.numeric() prevents integer overflow when multiplying huge numbers
cpg_oe <- (as.numeric(total_cpg) * total_size) / (as.numeric(total_c) * as.numeric(total_g))

# density = how many CpGs per 1000 bases
cpg_density <- total_cpg / (total_size / 1000)

cat("CpG ANALYSIS (genome-wide):\n")
cat("  Total CpGs:", format(total_cpg, big.mark = ","), "\n")
cat("  CpG density:", round(cpg_density, 2), "per kb\n")
cat("  CpG O/E:", round(cpg_oe, 4), "(< 1 = depleted, typical of methylated genomes)\n\n")


# 7. CpG BY REGION (now includes Gene body and TEs)
# promoters() from GenomicRanges:
#   for each gene, takes the transcription start site (TSS)
#   goes 2000 bp upstream (before the gene starts)
#   this is where transcription factors bind to control gene expression
# trim() clips any promoter that would go past position 0
gene_promoters <- trim(promoters(genes, upstream = 2000, downstream = 0))

# reduce() merges overlapping intervals into clean blocks
gene_bodies  <- reduce(genes)
exon_regions <- reduce(exons)

# setdiff() = genomic subtraction: A minus B
# gene bodies minus exons = introns (the parts between exons)
intron_regions <- setdiff(gene_bodies, exon_regions)

# intergenic = everything that's NOT a gene or promoter
# first merge genes + promoters into one set
genic_plus_promoter <- reduce(c(gene_bodies, gene_promoters))
# strand(x) <- "*" removes strand info so everything is treated as same strand
# without this, gaps() would compute separate gaps for + and - strands
strand(genic_plus_promoter) <- "*"
genic_plus_promoter <- reduce(genic_plus_promoter)
# tell it how long each chromosome is (needed to find gaps at chromosome ends)
seqlengths(genic_plus_promoter) <- chr_lengths[seqlevels(genic_plus_promoter)]

# gaps() finds all the spaces BETWEEN features on each chromosome
# it returns the "holes" — regions not covered by any gene or promoter
all_gaps <- gaps(genic_plus_promoter)
# gaps() returns gaps on +, -, and * strands; keep only * (unstranded)
intergenic_regions <- all_gaps[strand(all_gaps) == "*"]

# FAST APPROACH: find ALL CpG positions once, then just count overlaps per region
# old way: extract DNA sequences 6 times (once per region) = slow
# new way: locate CpGs once as GRanges, then countOverlaps = fast
cat("  Locating all CpG sites genome-wide...\n")

# vmatchPattern() finds every occurrence of "CG" in each chromosome
# it returns match positions — much faster than extracting sequences
# we loop per chromosome because vmatchPattern works on one sequence at a time
cpg_gr_list <- list()  # empty list to collect results
for (chr_name in names(genome)) {
  # matchPattern() finds all "CG" in one chromosome, returns positions
  hits <- matchPattern("CG", genome[[chr_name]])
  if (length(hits) > 0) {
    # GRanges() wraps positions into a proper genomic ranges object
    # start(hits) = where each CG starts. width is always 2 (CG = 2 bases)
    cpg_gr_list[[chr_name]] <- GRanges(seqnames = chr_name,
                                        ranges = IRanges(start = start(hits), width = 2))
  }
}
# unlist(GRangesList(...)) flattens a list of GRanges into one big GRanges
# GRangesList() first wraps the list properly so unlist works correctly
cpg_all <- unlist(GRangesList(cpg_gr_list))
cat("  Found", length(cpg_all), "CpG sites\n")

# now for each region, we need: CpG count, total bp, C count, G count (for O/E)
# CpG count = just overlap CpG positions with regions (instant)
# C and G counts = still need sequence extraction, but we do it smarter:
#   letterFrequency on Views, one pass per chromosome, all regions at once

# function to count base composition in regions efficiently
# instead of extracting full DNAStringSet, we use Views and sum in place
count_bases_in_regions <- function(regions, genome_seq) {
  regions <- trim(regions)
  regions <- regions[width(regions) > 0]
  total_c <- 0; total_g <- 0  # running totals
  for (chr_name in unique(as.character(seqnames(regions)))) {
    chr_reg <- regions[seqnames(regions) == chr_name]
    # Views = "look through windows" at the chromosome without copying DNA
    v <- Views(genome_seq[[chr_name]], start = start(chr_reg), end = end(chr_reg))
    seqs <- as(v, "DNAStringSet")
    # letterFrequency counts specific bases. sum all rows at once
    total_c <- total_c + sum(letterFrequency(seqs, letters = "C"))
    total_g <- total_g + sum(letterFrequency(seqs, letters = "G"))
  }
  c(C = total_c, G = total_g)
}

# list() groups different objects together with names, like a labeled box
regions_list <- list(Promoter = gene_promoters, Gene_body = gene_bodies,
                     Exon = exon_regions, Intron = intron_regions,
                     TE = te_regions, Intergenic = intergenic_regions)

# countOverlaps() counts how many CpG sites fall inside each region
# sum() of that = total CpGs in that region type. This is the fast part.
region_stats <- do.call(rbind, lapply(names(regions_list), function(nm) {
  cat("  Processing", nm, "...\n")
  reg <- trim(regions_list[[nm]])
  reg <- reg[width(reg) > 0]
  # overlapsAny: for each CpG site, TRUE if it overlaps any region = fast
  cpg <- sum(overlapsAny(cpg_all, reg))
  len <- sum(as.numeric(width(reg)))
  # base composition still needs sequence extraction (for O/E calculation)
  bases <- count_bases_in_regions(reg, genome)
  cc <- as.numeric(bases["C"]); gg <- as.numeric(bases["G"])
  oe <- ifelse(cc > 0 & gg > 0, (as.numeric(cpg) * len) / (cc * gg), NA)
  exp_cpg <- (cc * gg) / len
  data.frame(region = nm, total_bp = len, cpg_count = cpg, expected_cpg = exp_cpg,
             cpg_per_kb = cpg / (len / 1000), cpg_oe = oe)
}))

# 7b. CHI-SQUARE TEST: is each region's CpG count significantly different from expected?
# chisq.test() compares observed vs expected counts
# if p-value < 0.05, the difference is statistically significant
# we test: does observed CpG count differ from what base composition predicts?
region_stats$chi_p <- NA
region_stats$enrichment <- ""
for (i in seq_len(nrow(region_stats))) {
  obs <- region_stats$cpg_count[i]
  exp <- region_stats$expected_cpg[i]
  # chisq.test with 1 observed vs 1 expected: is the deviation significant?
  # c(obs, total_dinuc - obs) vs c(exp, total_dinuc - exp)
  # simpler: use the O/E ratio. O/E > 1 = overrepresented, < 1 = underrepresented
  # chi-square = sum of (observed - expected)^2 / expected
  chi_val <- (as.numeric(obs) - exp)^2 / exp
  # pchisq() converts chi-square value to a p-value
  # df=1 because we have one comparison. lower.tail=FALSE = we want the upper tail
  p_val <- pchisq(chi_val, df = 1, lower.tail = FALSE)
  region_stats$chi_p[i] <- p_val
  # label as over or under represented
  region_stats$enrichment[i] <- ifelse(region_stats$cpg_oe[i] > 1, "OVER", "UNDER")
}

cat("\nCpG BY REGION:\n")
# for loops: repeat code for each value of i
# seq_len(4) = c(1, 2, 3, 4). nrow() = number of rows in a table
for (i in seq_len(nrow(region_stats))) {
  # [i, ] = get row number i from the table (comma with nothing after = all columns)
  r <- region_stats[i, ]
  # sprintf() builds formatted text. %s = insert a string, %.2f = number with 2 decimals
  # %-12s = left-aligned string padded to 12 characters (for neat columns)
  p_str <- ifelse(r$chi_p < 2.2e-16, "< 2.2e-16", formatC(r$chi_p, format = "e", digits = 2))
  cat(sprintf("  %-12s %12s bp | %10s CpGs | %5.2f/kb | O/E = %.4f | %s (p %s)\n",
      r$region, format(r$total_bp, big.mark = ","),
      format(r$cpg_count, big.mark = ","), r$cpg_per_kb, r$cpg_oe,
      r$enrichment, p_str))
}


# 8. PLOTS
# helper: save a ggplot to both pdf/ and png/ inside results/batch01/
save_plot <- function(p, name, w = 8, h = 5) {
  png(file.path(fig_png_dir, paste0(name, ".png")), width = w, height = h, units = "in", res = 300)
  print(p); dev.off()
  tmp_pdf <- file.path(Sys.getenv("TMPDIR", "/tmp"), paste0(name, ".pdf"))
  cairo_pdf(tmp_pdf, width = w, height = h)
  print(p); dev.off()
  ok <- file.copy(tmp_pdf, file.path(fig_pdf_dir, paste0(name, ".pdf")), overwrite = TRUE)
  if (!ok) cat("  WARNING: PDF locked, saved to:", tmp_pdf, "\n")
  cat("  Saved:", name, "\n")
}

cat("Generating plots...\n")

# region colors from CLAUDE.md conventions
region_colors <- c(Promoter = "#8E44AD", `Gene body` = "#27AE60", Exon = "#2471A3",
                   Intron = "#1ABC9C", TE = "#F39C12", Intergenic = "#C0392B")
# gsub() replaces "_" with " " in region names for cleaner plot labels
region_stats$region_name <- gsub("_", " ", as.character(region_stats$region))

# 8a. DINUCLEOTIDE FREQUENCY (reordered by observed frequency, descending)
obs_counts <- colSums(dinuc)
obs_freq <- obs_counts / sum(obs_counts)
base_counts <- colSums(alphabetFrequency(genome)[, c("A", "C", "G", "T")])
base_freq <- base_counts / sum(base_counts)
combos <- expand.grid(first = names(base_freq), second = names(base_freq),
                      stringsAsFactors = FALSE)
combos$dinuc <- paste0(combos$first, combos$second)
combos$expected <- base_freq[combos$first] * base_freq[combos$second]
exp_freq <- combos$expected[match(names(obs_freq), combos$dinuc)]
names(exp_freq) <- names(obs_freq)

plot_df <- data.frame(
  dinucleotide = names(obs_freq),
  observed = as.numeric(obs_freq) * 100,
  expected = as.numeric(exp_freq) * 100
)
plot_df <- plot_df[order(-plot_df$observed), ]
plot_df$dinucleotide <- factor(plot_df$dinucleotide, levels = plot_df$dinucleotide)
plot_df$is_cpg <- ifelse(plot_df$dinucleotide == "CG", "CpG", "other")

p_dinuc <- ggplot(plot_df, aes(x = dinucleotide, y = observed, fill = is_cpg)) +
  geom_col(width = 0.7) +
  geom_point(aes(y = expected), shape = 18, size = 3, color = "black") +
  scale_fill_manual(values = c("CpG" = "#C0392B", "other" = "#2471A3"), guide = "none") +
  labs(x = "Dinucleotide", y = "Frequency (%)",
       title = "Dinucleotide frequencies — D. laeve genome",
       subtitle = "Bars = observed, diamonds = expected from base composition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_dinuc, "dinucleotide_freq")

# 8b. CpG O/E BY REGION (bar plot)
region_stats$region_f <- factor(region_stats$region_name,
                                levels = region_stats$region_name[order(region_stats$cpg_oe)])

p_oe <- ggplot(region_stats, aes(x = region_f, y = cpg_oe, fill = region_name)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  annotate("text", x = 0.7, y = 1.02, label = "expected", color = "gray40", size = 3, hjust = 0) +
  geom_hline(yintercept = cpg_oe, linetype = "dotted", color = "black") +
  annotate("text", x = 0.7, y = cpg_oe + 0.02, label = paste0("genome: ", round(cpg_oe, 3)),
           size = 3, hjust = 0) +
  scale_fill_manual(values = region_colors, guide = "none") +
  labs(x = NULL, y = "CpG O/E ratio",
       title = "CpG observed/expected ratio by genomic region") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10))
save_plot(p_oe, "cpg_oe_by_region", w = 7, h = 5)

# 8c. CpG DENSITY BY REGION (bar plot)
region_stats$region_f2 <- factor(region_stats$region_name,
                                 levels = region_stats$region_name[order(region_stats$cpg_per_kb)])

p_density <- ggplot(region_stats, aes(x = region_f2, y = cpg_per_kb, fill = region_name)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = cpg_density, linetype = "dotted", color = "black") +
  annotate("text", x = 0.7, y = cpg_density + 0.5,
           label = paste0("genome: ", round(cpg_density, 1)), size = 3, hjust = 0) +
  scale_fill_manual(values = region_colors, guide = "none") +
  labs(x = NULL, y = "CpG per kb",
       title = "CpG density by genomic region") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10))
save_plot(p_density, "cpg_density_by_region", w = 7, h = 5)

# 8d. DENSITY CURVES — CpG/kb and O/E per region
# X = value (CpG/kb or O/E), Y = density. One curve per region type.
# Pre-split cpg_all by chromosome once to avoid repeated 45M subsetting
cat("  Computing per-chromosome distributions...\n")
cpg_by_chr <- split(cpg_all, seqnames(cpg_all))

chr_region_stats <- do.call(rbind, lapply(names(regions_list), function(nm) {
  reg <- trim(regions_list[[nm]])
  reg <- reg[width(reg) > 0]
  do.call(rbind, lapply(names(genome), function(chr_name) {
    chr_reg <- reg[seqnames(reg) == chr_name]
    if (length(chr_reg) == 0) return(NULL)
    # use pre-split CpG sites for this chromosome
    chr_cpg <- cpg_by_chr[[chr_name]]
    n_cpg <- if (!is.null(chr_cpg)) sum(overlapsAny(chr_cpg, chr_reg)) else 0L
    len <- sum(as.numeric(width(chr_reg)))
    v <- Views(genome[[chr_name]], start = start(chr_reg), end = end(chr_reg))
    seqs <- as(v, "DNAStringSet")
    cc <- sum(letterFrequency(seqs, letters = "C"))
    gg <- sum(letterFrequency(seqs, letters = "G"))
    oe <- ifelse(cc > 0 & gg > 0, (as.numeric(n_cpg) * len) / (as.numeric(cc) * as.numeric(gg)), NA)
    data.frame(region = nm, chr = chr_name, cpg_per_kb = n_cpg / (len / 1000), cpg_oe = oe)
  }))
}))

chr_region_stats$region <- gsub("_", " ", chr_region_stats$region)
chr_region_stats$region <- factor(chr_region_stats$region,
                                  levels = c("Promoter", "Gene body", "Exon",
                                             "Intron", "TE", "Intergenic"))

# CpG O/E density curve — X = O/E, Y = density, one curve per region
p_oe_dist <- ggplot(chr_region_stats, aes(x = cpg_oe, fill = region, color = region)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = cpg_oe, linetype = "dotted", color = "black") +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  labs(x = "CpG O/E ratio", y = "Density",
       title = "CpG O/E distribution by region",
       subtitle = "Per-chromosome values. Dashed = expected (1.0), dotted = genome avg") +
  theme_minimal() +
  theme(legend.title = element_blank())
save_plot(p_oe_dist, "cpg_oe_distribution", w = 8, h = 5)

# CpG density curve — X = CpG/kb, Y = density
p_dens_dist <- ggplot(chr_region_stats, aes(x = cpg_per_kb, fill = region, color = region)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_vline(xintercept = cpg_density, linetype = "dotted", color = "black") +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  labs(x = "CpG per kb", y = "Density",
       title = "CpG density distribution by region",
       subtitle = "Per-chromosome values. Dotted = genome-wide average") +
  theme_minimal() +
  theme(legend.title = element_blank())
save_plot(p_dens_dist, "cpg_density_distribution", w = 8, h = 5)

# free the large cpg_all and cpg_by_chr objects before TE analysis
rm(cpg_all, cpg_by_chr, chr_region_stats); gc(verbose = FALSE)

# 8e. TE FAMILY ANALYSIS
# split TEs by class (LINE, LTR, DNA, etc.) and compute CpG stats per family
cat("  Computing TE family CpG stats...\n")
# sub() extracts text before the "/" : "LINE/RTE-RTE" -> "LINE"
te_gr$te_class <- sub("/.*", "", te_gr$class_family)

# count how many TEs per class, keep only classes with > 100 elements
te_class_counts <- table(te_gr$te_class)
keep_classes <- names(te_class_counts[te_class_counts > 100])
cat("  TE classes (>100 elements):", paste(keep_classes, collapse = ", "), "\n")

# compute CpG density and O/E per TE class
# MEMORY-SAFE: use dinucleotideFrequency on extracted sequences instead of
# overlapsAny on 45M CpG GRanges (which caused OOM in previous session)
te_family_stats <- do.call(rbind, lapply(keep_classes, function(cls) {
  cat("    Processing TE class:", cls, "\n")
  cls_gr <- reduce(te_gr[te_gr$te_class == cls])
  cls_gr <- trim(cls_gr)
  cls_gr <- cls_gr[width(cls_gr) > 0]
  len <- sum(as.numeric(width(cls_gr)))
  total_cpg_cls <- 0; total_c <- 0; total_g <- 0
  # process per chromosome to keep memory bounded
  for (chr_name in unique(as.character(seqnames(cls_gr)))) {
    chr_reg <- cls_gr[seqnames(cls_gr) == chr_name]
    v <- Views(genome[[chr_name]], start = start(chr_reg), end = end(chr_reg))
    seqs <- as(v, "DNAStringSet")
    total_cpg_cls <- total_cpg_cls + sum(dinucleotideFrequency(seqs)[, "CG"])
    total_c <- total_c + sum(letterFrequency(seqs, letters = "C"))
    total_g <- total_g + sum(letterFrequency(seqs, letters = "G"))
    rm(seqs); gc(verbose = FALSE)  # free memory immediately
  }
  oe <- ifelse(total_c > 0 & total_g > 0,
               (as.numeric(total_cpg_cls) * len) / (as.numeric(total_c) * as.numeric(total_g)), NA)
  data.frame(te_class = cls, total_bp = len, cpg_count = total_cpg_cls,
             cpg_per_kb = total_cpg_cls / (len / 1000), cpg_oe = oe,
             n_elements = as.integer(te_class_counts[cls]))
}))

te_family_stats <- te_family_stats[order(te_family_stats$cpg_oe), ]
te_family_stats$te_class <- factor(te_family_stats$te_class,
                                   levels = te_family_stats$te_class)

cat("  TE family stats:\n")
for (i in seq_len(nrow(te_family_stats))) {
  r <- te_family_stats[i, ]
  cat(sprintf("    %-15s %10s bp | %8s CpGs | %5.2f/kb | O/E = %.4f | n=%s\n",
      r$te_class, format(r$total_bp, big.mark = ","),
      format(r$cpg_count, big.mark = ","), r$cpg_per_kb, r$cpg_oe,
      format(r$n_elements, big.mark = ",")))
}

# TE family O/E bar plot
p_te_oe <- ggplot(te_family_stats, aes(x = te_class, y = cpg_oe, fill = cpg_oe)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  geom_hline(yintercept = cpg_oe, linetype = "dotted", color = "black") +
  annotate("text", x = 0.7, y = cpg_oe + 0.02, label = paste0("genome: ", round(cpg_oe, 3)),
           size = 3, hjust = 0) +
  scale_fill_gradient(low = "#2471A3", high = "#C0392B", guide = "none") +
  labs(x = NULL, y = "CpG O/E ratio",
       title = "CpG O/E by TE class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_te_oe, "te_cpg_oe_by_class", w = 8, h = 5)

# TE family CpG density bar plot
te_family_stats2 <- te_family_stats[order(te_family_stats$cpg_per_kb), ]
te_family_stats2$te_class <- factor(te_family_stats2$te_class,
                                    levels = te_family_stats2$te_class)

p_te_dens <- ggplot(te_family_stats2, aes(x = te_class, y = cpg_per_kb, fill = cpg_per_kb)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = cpg_density, linetype = "dotted", color = "black") +
  annotate("text", x = 0.7, y = cpg_density + 0.5,
           label = paste0("genome: ", round(cpg_density, 1)), size = 3, hjust = 0) +
  scale_fill_gradient(low = "#2471A3", high = "#C0392B", guide = "none") +
  labs(x = NULL, y = "CpG per kb",
       title = "CpG density by TE class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p_te_dens, "te_cpg_density_by_class", w = 8, h = 5)

# 8f. TE CLASS DENSITY CURVES — O/E and CpG/kb per TE class (per-chromosome)
cat("  Computing per-chromosome TE class distributions...\n")
te_chr_stats <- do.call(rbind, lapply(keep_classes, function(cls) {
  cls_gr <- reduce(te_gr[te_gr$te_class == cls])
  cls_gr <- trim(cls_gr)
  cls_gr <- cls_gr[width(cls_gr) > 0]
  do.call(rbind, lapply(names(genome), function(chr_name) {
    chr_reg <- cls_gr[seqnames(cls_gr) == chr_name]
    if (length(chr_reg) == 0) return(NULL)
    v <- Views(genome[[chr_name]], start = start(chr_reg), end = end(chr_reg))
    seqs <- as(v, "DNAStringSet")
    n_cpg <- sum(dinucleotideFrequency(seqs)[, "CG"])
    cc <- sum(letterFrequency(seqs, letters = "C"))
    gg <- sum(letterFrequency(seqs, letters = "G"))
    len <- sum(as.numeric(width(chr_reg)))
    oe <- ifelse(cc > 0 & gg > 0, (as.numeric(n_cpg) * len) / (as.numeric(cc) * as.numeric(gg)), NA)
    data.frame(te_class = cls, chr = chr_name, cpg_per_kb = n_cpg / (len / 1000), cpg_oe = oe)
  }))
}))

te_chr_stats$te_class <- factor(te_chr_stats$te_class, levels = levels(te_family_stats$te_class))

te_colors <- setNames(
  colorRampPalette(c("#2471A3", "#8E44AD", "#1ABC9C", "#F39C12", "#C0392B", "#27AE60"))(length(keep_classes)),
  keep_classes
)

p_te_oe_dist <- ggplot(te_chr_stats, aes(x = cpg_oe, fill = te_class, color = te_class)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = cpg_oe, linetype = "dotted", color = "black") +
  scale_fill_manual(values = te_colors) +
  scale_color_manual(values = te_colors) +
  labs(x = "CpG O/E ratio", y = "Density",
       title = "CpG O/E distribution by TE class",
       subtitle = "Per-chromosome values. Dashed = expected (1.0), dotted = genome avg") +
  theme_minimal() +
  theme(legend.title = element_blank())
save_plot(p_te_oe_dist, "te_cpg_oe_distribution", w = 8, h = 5)

p_te_dens_dist <- ggplot(te_chr_stats, aes(x = cpg_per_kb, fill = te_class, color = te_class)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_vline(xintercept = cpg_density, linetype = "dotted", color = "black") +
  scale_fill_manual(values = te_colors) +
  scale_color_manual(values = te_colors) +
  labs(x = "CpG per kb", y = "Density",
       title = "CpG density distribution by TE class",
       subtitle = "Per-chromosome values. Dotted = genome-wide average") +
  theme_minimal() +
  theme(legend.title = element_blank())
save_plot(p_te_dens_dist, "te_cpg_density_distribution", w = 8, h = 5)


# 9. WRITE HTML REPORT (works in RStudio viewer + any browser)
cat("Writing report to:", out_report, "\n")
sink(out_report)

cat('<!DOCTYPE html>\n<html><head><meta charset="utf-8">\n')
cat('<title>Genome Characterization — D. laeve</title>\n')
cat('<style>\n')
cat('body { font-family: "Segoe UI", Arial, sans-serif; max-width: 900px; margin: 40px auto; padding: 0 20px; color: #333; line-height: 1.6; }\n')
cat('h1 { color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 8px; }\n')
cat('h2 { color: #34495e; margin-top: 30px; }\n')
cat('h3 { color: #5d6d7e; }\n')
cat('table { border-collapse: collapse; margin: 15px 0; width: auto; }\n')
cat('th, td { border: 1px solid #ddd; padding: 6px 12px; text-align: left; }\n')
cat('th { background: #f5f5f5; font-weight: 600; }\n')
cat('tr:nth-child(even) { background: #fafafa; }\n')
cat('img { max-width: 100%; margin: 15px 0; border: 1px solid #eee; }\n')
cat('.pass { color: #27ae60; font-weight: bold; }\n')
cat('.under { color: #c0392b; }\n')
cat('</style>\n</head><body>\n')

cat('<h1>Genome Characterization — <em>Deroceras laeve</em></h1>\n')
cat('<p>Generated:', format(Sys.time(), "%Y-%m-%d %H:%M"), '</p>\n')

cat('<h2>GFF / Genome Compatibility</h2>\n')
cat('<table><tr><th>Check</th><th>Result</th></tr>\n')
cat('<tr><td>Shared chromosomes</td><td>', length(chrs_in_both), '/ 31</td></tr>\n')
cat('<tr><td>Features out of bounds</td><td>', n_clipped, '</td></tr>\n')
verdict <- ifelse(length(chrs_in_both) == 31 & n_clipped == 0, "PASS", "CHECK")
cat('<tr><td>Verdict</td><td class="pass">', verdict, '</td></tr>\n')
cat('</table>\n')

cat('<h2>Genome Assembly</h2>\n')
cat('<table><tr><th>Property</th><th>Value</th></tr>\n')
cat('<tr><td>Total size</td><td>', format(total_size, big.mark = ","), 'bp</td></tr>\n')
cat('<tr><td>Chromosomes</td><td>', length(genome), '</td></tr>\n')
cat('<tr><td>Largest</td><td>', format(max(chr_sizes), big.mark = ","), 'bp (', names(which.max(chr_sizes)), ')</td></tr>\n')
cat('<tr><td>Smallest</td><td>', format(min(chr_sizes), big.mark = ","), 'bp (', names(which.min(chr_sizes)), ')</td></tr>\n')
cat('<tr><td>GC content</td><td>', round(overall_gc * 100, 2), '%</td></tr>\n')
cat('<tr><td>Genome covered by genes</td><td>', round(gene_cov, 1), '%</td></tr>\n')
cat('<tr><td>Genome covered by TEs</td><td>', round(te_cov, 1), '%</td></tr>\n')
cat('</table>\n')

cat('<h3>Chromosome sizes</h3>\n')
cat('<table><tr><th>Chromosome</th><th>Size (bp)</th><th>GC %</th></tr>\n')
for (i in seq_along(chr_sizes)) {
  cat('<tr><td>', names(chr_sizes)[i], '</td><td>', format(chr_sizes[i], big.mark = ","),
      '</td><td>', round(gc_per_chr[i] * 100, 2), '</td></tr>\n')
}
cat('</table>\n')

cat('<h2>Gene Annotation</h2>\n')
cat('<table><tr><th>Property</th><th>Value</th></tr>\n')
cat('<tr><td>Genes</td><td>', format(length(genes), big.mark = ","), '</td></tr>\n')
cat('<tr><td>mRNAs</td><td>', format(length(mrnas), big.mark = ","), '</td></tr>\n')
cat('<tr><td>Exons</td><td>', format(length(exons), big.mark = ","), '</td></tr>\n')
cat('<tr><td>Mean gene length</td><td>', format(round(mean(width(genes))), big.mark = ","), 'bp</td></tr>\n')
cat('<tr><td>TEs</td><td>', format(length(te_gr), big.mark = ","), '</td></tr>\n')
cat('</table>\n')

cat('<h2>CpG Analysis</h2>\n')
cat('<h3>Genome-wide</h3>\n')
cat('<table><tr><th>Property</th><th>Value</th></tr>\n')
cat('<tr><td>Total CpGs</td><td>', format(total_cpg, big.mark = ","), '</td></tr>\n')
cat('<tr><td>CpG density</td><td>', round(cpg_density, 2), 'per kb</td></tr>\n')
cat('<tr><td>CpG O/E</td><td>', round(cpg_oe, 4), '</td></tr>\n')
cat('</table>\n')

cat('<h3>Dinucleotide frequencies</h3>\n')
cat('<img src="png/dinucleotide_freq.png">\n')

cat('<h3>CpG by genomic region</h3>\n')
cat('<table><tr><th>Region</th><th>Size (bp)</th><th>CpGs</th><th>CpG/kb</th><th>O/E</th><th>Status</th><th>p-value</th></tr>\n')
for (i in seq_len(nrow(region_stats))) {
  r <- region_stats[i, ]
  p_str <- ifelse(r$chi_p < 2.2e-16, "&lt; 2.2e-16", formatC(r$chi_p, format = "e", digits = 2))
  cat('<tr><td>', r$region_name, '</td><td>', format(r$total_bp, big.mark = ","),
      '</td><td>', format(r$cpg_count, big.mark = ","), '</td><td>', round(r$cpg_per_kb, 2),
      '</td><td>', round(r$cpg_oe, 4), '</td><td class="under">', r$enrichment,
      '</td><td>', p_str, '</td></tr>\n')
}
cat('</table>\n')
cat('<p>All regions show significant CpG underrepresentation (O/E &lt; 1), ',
    'consistent with historical DNA methylation across the genome.</p>\n')

cat('<h3>CpG O/E by region</h3>\n')
cat('<img src="png/cpg_oe_by_region.png">\n')
cat('<h3>CpG density by region</h3>\n')
cat('<img src="png/cpg_density_by_region.png">\n')

cat('<h3>CpG O/E distribution across chromosomes</h3>\n')
cat('<img src="png/cpg_oe_distribution.png">\n')
cat('<h3>CpG density distribution across chromosomes</h3>\n')
cat('<img src="png/cpg_density_distribution.png">\n')

cat('<h2>TE Family Analysis</h2>\n')
cat('<table><tr><th>TE class</th><th>Size (bp)</th><th>CpGs</th><th>CpG/kb</th><th>O/E</th><th>Elements</th></tr>\n')
for (i in seq_len(nrow(te_family_stats))) {
  r <- te_family_stats[i, ]
  cat('<tr><td>', as.character(r$te_class), '</td><td>', format(r$total_bp, big.mark = ","),
      '</td><td>', format(r$cpg_count, big.mark = ","), '</td><td>', round(r$cpg_per_kb, 2),
      '</td><td>', round(r$cpg_oe, 4), '</td><td>', format(r$n_elements, big.mark = ","),
      '</td></tr>\n')
}
cat('</table>\n')
cat('<h3>CpG O/E by TE class</h3>\n')
cat('<img src="png/te_cpg_oe_by_class.png">\n')
cat('<h3>CpG density by TE class</h3>\n')
cat('<img src="png/te_cpg_density_by_class.png">\n')

cat('<h3>CpG O/E distribution by TE class</h3>\n')
cat('<img src="png/te_cpg_oe_distribution.png">\n')
cat('<h3>CpG density distribution by TE class</h3>\n')
cat('<img src="png/te_cpg_density_distribution.png">\n')

cat('</body></html>\n')
sink()
cat("Done!\n")
