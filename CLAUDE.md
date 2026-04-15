# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

Epigenetic regulation of tail regeneration in *Deroceras laeve* (gastropod mollusk). WGBS + RNA-seq → Science Advances paper: **does DNA methylation regulate gene expression during invertebrate regeneration?**

All analysis code in R. The `methylation_pipeline/` folder is supplementary code for the paper — it must produce all data and figures used in the manuscript.

## Organism

- **Species**: *Deroceras laeve* (gray field slug, gastropod mollusk)
- **Non-model invertebrate** — no UCSC, no Ensembl, no CpG islands, no mammalian assumptions
- **Genome**: 1.78 Gb, 31 chromosomes (chr1–chr31), 42.99% GC. ALWAYS filter HiC_scaffold_*
- **Annotation**: EviAnn GFF
- **Key**: Has DNMT1 (maintenance) but NO DNMT3 (de novo). Mosaic gene body methylation.

## Experimental design

- **WGBS**: 2 control (C1, C2) vs 2 amputated (A1, A2) tail. Bismark CpG_report.txt. DSS for DM.
- **RNA-seq**: Multi-tissue (tail, eye, bodywall, head, juvenile, ovotestis), multi-condition (control, amputated, irradiated, fungicide). DESeq2 + WGCNA.
- **Strand collapse required**: CpG reports have both strands. Collapse before any analysis.

## Key numbers (strand-collapsed, verified)

| Metric | Value |
|--------|-------|
| Genome CpGs (chr1-31) | 45,535,149 |
| CpG O/E | 0.553 |
| CpGs tested (cov>=5, all 4 samples) | 25,439,068 |
| DMPs (FDR<0.05, \|diff\|>10%) | 17,729 (8,285 hyper + 9,444 hypo) |
| DMRs | 1,462 (690 hyper + 772 hypo) |
| Promoter classes | 85.1% ICP, 0.1% HCP, 14.9% LCP |
| Predicted TFs | 2,217 genes (DeepTFactor) |
| JASPAR TFBS motifs | 1,362 metazoan PWMs |
| Baseline gene body meth vs expression | rho = 0.407 (p < 10^-80) |
| Differential meth vs expression MI | 0.014 bits (p = 0.496, n.s.) |

## Repository structure

```
methylation_pipeline/       ← MAIN: supplementary code for paper (reviewers see this)
  _config.R                 — shared config: paths, colors, helpers (see below)
  run_pipeline.slurm        — SLURM runner (4 CPUs, 64 GB): 0.5→1.5→02→03→...→10
  generate_report.R         — builds pipeline_report.html from all batch figures (run from repo root)

  batch0.5/ "What is the clean, filtered expression dataset?"
            Transcriptome Processing — all 44 samples, 6 tissues, 4 conditions.
            Within-group outlier detection (centroid distance + PCA). DESeq2 per-tissue
            (not one big model). Gene filter: >= 10 counts in >= min_group_size.
            LFC shrinkage: apeglm. Caches CACHE$transcriptome for all downstream batches.
            → 5+ TSV, 4 figures (dendrogram, PCA all, PCA tail, correlation heatmap)

  batch01/  "What does the D. laeve CpG landscape look like?"
            Genomic CpG Annotation — dinucleotide freq, CpG density + O/E + GC%
            distributions by region & TE class
            → 3 TSV, 7 figures (a-g)

  batch1.5/ "Where are transcription factor binding sites genome-wide?"
            TFBS Motif Annotation — JASPAR 2024 (1362 PWMs) + HOMER de novo, genome-wide scan.
            matchMotifs chunked (50 motifs/batch) to avoid MOODS hang on large inputs.
            → motif hit table for batch09. Requires 64 GB RAM.

  batch02/  "Does D. laeve have a functional methylation toolkit?"
            Methylation Machinery — 26 genes (18 present + 8 absent), expression across 6 tissues
            → 1 TSV, 3 figures (a-c)

  batch03/  "What does the promoter CpG landscape look like?"
            Promoter CpG Classification — Weber HCP/ICP/LCP via O/E + GC% thresholds,
            HCP gene list with function annotations
            → 5 TSV, 5 figures (a-e)

  batch04/  "Where is methylation? Does gene body meth correlate with expression?"
            Baseline Methylation Landscape — strand-collapse BSseq, metagene profile,
            gene body meth vs expression (Spearman by region + decile),
            genome-wide 1Mb windows, single-exon metagene, expanded region correlations
            → 10 TSV, 10 figures (a-j,l) + 27 metagene panels (3 bin sizes × 9 categories)

  batch05/  "Is TE methylation related to evolutionary age?"
            TE Methylation + Age — ridgeplots by class, intronic vs intergenic, Kimura age,
            DMP/DMR overlap, Cohen's d effect sizes
            → 7 TSV, 6 figures (a,b,e,f,g,j)

  batch06/  "What methylation changes during regeneration?"
            Differential Methylation — DMP/DMR calling (DSS), PCA, Manhattan, volcano
            (labeled + unlabeled), annotation, enrichment (DMP + DMR Fisher's),
            heatmaps with gene names, top genes barplots (DMP + DMR), GO/KEGG/Reactome
            → 12+ TSV, 15+ figures

  batch07/  "Does differential methylation predict expression changes?"
            Differential Methylation vs Expression — DMP and DMR treated SEPARATELY.
            Per-region correlation (Spearman, basic + expanded), DE rate, dose-response.
            No MI (that's batch 10). No module analysis (that's batch 08).
            → 12 TSV, 11 figures (a-k)

  batch08/  "Are methylation changes concentrated in specific expression modules?"
            Methylation of WGCNA Co-expression Modules — multi-tissue (6 tissues),
            outlier detection, module-trait correlations, DMP burden + Fisher, hub genes, GO
            → 10 TSV, 10 figures (a-j)

  batch09/  "Do TF binding sites near DMPs explain regulatory connections?"
            TF + GENIE3 + Motif Integration — DeepTFactor enrichment, GENIE3 network,
            motif overlap with DMPs
            → 8 TSV, 9 figures (a-i)

  batch10/  "Does regeneration increase or decrease methylation entropy?"
            Entropy Analysis — Level 1 binary Shannon + Wilcoxon, per-region delta,
            Fang framework: sliding-window entropy regions, hidden layer,
            per-gene delta entropy, CpG density stratification, direction-specific.
            10b: Level 2 per-read NME from BAMs (Xie 4-CpG windows), cluster only
            → 14+ TSV, 14 figures (a-n)

  batch11/  "Does TF binding enforce ordered methylation at DE genes?"
            Motif x Per-read Entropy x Expression — closes the loop
            TF motif → low NME (ordered) → DMP → DE gene → phenotype.
            Loads batch10b NME windows + batch1.5 genome-wide motif hits +
            CACHE$transcriptome. Min reads/sample = 8. Tests:
            (A) beta-matched in-motif vs out-of-motif baseline NME (Wilcoxon)
            (B) delta NME at TFBS stratified by Yin 2017 methyl-sensitivity
            (C) directed-mechanism shortlist: in_motif & ctrl_nme<0.3 & has_dmp,
                joined to DESeq2 res_tail → cloning candidates
            Read coverage diagnostics per TF class (sanity).
            → 5 TSV, 6 figures (a-f) [NOT YET RUN — blocked on batch10b]

cluster/scripts/            — R + SLURM for HPC jobs
genome/cache/               — RDS caches (genome, GFF, BSseq, DMLtest, TE, promoters, extended)
papers/                     — Reference PDFs + paper drafts
PROTOCOL.md                 — Internal protocol (not for reviewers)
```

- **No HTML in batch scripts.** Batches produce data (TSV) + figures (PNG+PDF) only.
- Each batch sources `methylation_pipeline/_config.R` for shared config.
- Each batch cleans old output (`unlink()`) before regenerating.

## Architecture: `_config.R`

Every batch script starts with `source("methylation_pipeline/_config.R")`. It provides:

- **Environment auto-detection**: `ON_CLUSTER <- dir.exists("/mnt/data/alfredvar")` — switches all paths automatically.
- **`OG` list**: Named list of all input data file paths (CpG reports, GFF, TE age, counts, DeepTFactor, STRING, GENIE3). Access as `OG$gff`, `OG$cpg_C1`, etc.
- **`CACHE` list**: Named list of cached RDS file paths (`CACHE$genome`, `CACHE$bsseq`, `CACHE$dmltest`, `CACHE$extended`, etc.).
- **Loaders**: `load_genome()`, `load_gff()`, `load_te()` — load from cache if available, otherwise build from source and cache.
- **`annotate_regions(chr, pos, promoters, exons, genes)`**: Assigns mutually exclusive labels (Promoter > Exon > Intron > Intergenic). TE is NOT included here.
- **`save_fig(p, batch_dir, name, w, h)`**: Saves both PNG (300 dpi) and PDF (cairo_pdf). Always use this.
- **`save_data(df, batch_dir, name)`**: Writes TSV via `data.table::fwrite()`.
- **`COLORS` list**: Standardized palettes for direction, region, TE class, condition.
- **`keep_chr`**: Always `paste0("chr", 1:31)`.

## Batch dependencies

`run_pipeline.slurm` runs: **0.5 → 1.5 → 02 → 03 → 04 → 05 → 06 → 07 → 08 → 09 → 10 → 11** (batch01 pre-run separately or cached; batch10b BAM-based NME submitted separately). Key data flow:
- **batch0.5** → `CACHE$transcriptome` (filtered expression, apeglm LFC) → used by batch04, batch07, batch08, batch10
- **batch01** → genome CpG stats, cached genome → used by batch03, batch1.5, batch10
- **batch1.5** → genome-wide TFBS motif hits (`CACHE$extended`) → used by batch09
- **batch03** → promoter classifications (HCP/ICP/LCP) → used by batch04, batch06
- **batch04** → BSseq object (`CACHE$bsseq`), baseline methylation → used by batch05, batch06, batch07, batch10
- **batch05** → TE methylation (reads batch06 DMPs/DMRs conditionally if available)
- **batch06** → DMP/DMR lists (`dmps_annotated.tsv`, `dmrs_annotated.tsv`) → used by batch05, batch07, batch08, batch09, batch10
- **batch07** → decoupling stats (reads batch08 WGCNA modules if available for module-level correlation)
- **batch08** → WGCNA module assignments + DMP burden → used by batch07, batch09
- **batch09** → TF-DMP enrichment + GENIE3 network (reads batch1.5 motif hits)
- **batch10b** → `batch10/data/perread_nme_windows.tsv` (per-read 4-CpG NME from BAMs) → used by batch11
- **batch11** → reads batch10b NME windows + batch1.5 motif hits + batch06 DMPs + `CACHE$transcriptome`; outputs directed-mechanism shortlist

## R package dependencies

Core: `data.table`, `GenomicRanges`, `IRanges`, `rtracklayer`, `BSgenome`, `DSS`, `ggplot2`, `patchwork`
Batch-specific: `DESeq2` (batch02, batch04, batch07, batch08), `WGCNA` (batch08), `TFBSTools`/`JASPAR2024`/`motifmatchr`/`universalmotif` (batch1.5), `clusterProfiler` (batch06, batch08), `Rsamtools` (batch10b), `ggridges` (batch05), `pheatmap` (batch02, batch06, batch07, batch08)

## Data paths

### Local (Windows)

`C:/Users/rafae/Projects/DATA/` — raw data (READ-ONLY)

| File | What |
|------|------|
| `{C1,C2,A1,A2}.CpG_report.txt` | Bismark CpG reports (~3 GB each) |
| `derLaeGenome_namesDlasi_v2.fasta.functional_note.pseudo_label.gff` | EviAnn GFF (145 MB) |
| `derLaeGenome_eviann_annotations.tsv` | Gene names/descriptions |
| `collapsed_te_age_data.tsv` | TE age (Kimura divergence, 271 MB) |
| `counts_HTseq_EviAnn/` | HTSeq count matrices |
| `prediction_result.txt` | DeepTFactor TF predictions (`prediction` column is logical TRUE/FALSE) |
| `genie3_top500k.tsv` | GENIE3 regulatory network (top 500K edges) |
| `STRG0A31YWK.protein.orthology.v12.0.txt` | STRING orthology (`#string_protein_id` has leading `#`) |
| `protein.enrichment.terms.v12.0.txt` | STRING enrichment terms (GO/KEGG/Reactome) |

### Cluster (HPC)

`/mnt/data/alfredvar/` — shared data root

| File | Cluster path |
|------|-------------|
| CpG reports (.gz) | `/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/` |
| GFF | `/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/` |
| TE age | `/mnt/data/alfredvar/30-Genoma/32-Repeats/age_of_transposons/collapsed_te_age_data.tsv` |
| HTSeq counts | `/mnt/data/alfredvar/jmiranda/20-Transcriptomic_Bulk/25-metaAnalysisTranscriptome/counts_HTseq_EviAnn/` |
| GENIE3 (full) | `/mnt/data/alfredvar/wgutierrez/genie3_2/genie3_all_links.tsv` |
| STRING enrichment | `/mnt/data/alfredvar/rlopezt/Metilacion/protein.enrichment.terms.v12.0.txt` |
| Bismark BAMs | `/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/08_deduplication/` |

**Repo on HPC**: `/mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave/`
**Partition**: `defq` | **RAM**: 1 TB | **HOMER**: installed at `cluster/homer/`

### Cached objects (genome/cache/)

| File | What | Size |
|------|------|------|
| `genome_chr1_31.rds` | Genome DNAStringSet (from BSgenome.Dlaeve.NCBI.dlgm) | 365 MB |
| `gff_chr1_31.rds` | GFF GRanges (chr1-31) | 12 MB |
| `te_chr1_31.rds` | TE GRanges + reduced regions | 59 MB |
| `promoters_2kb.rds` | 2kb upstream promoter GRanges | 178 KB |
| `bsseq_tutorial.rds` | BSseq (25.4M CpGs × 4 samples, strand-collapsed, cov>=5) | 170 MB |
| `dmltest_full.rds` | DMLtest results (smoothing=TRUE, all 25.4M sites) | 1.4 GB |
| `extended_regions_10kb.rds` | Per-gene extended regions (promoter, upstream_distal, gene_body, downstream) for motif scanning | ~5 MB |

## Batch-specific implementation notes

### Batch 1.5: Chunked motif scanning
`matchMotifs()` from `motifmatchr` hangs when given all 1362 JASPAR motifs at once — the MOODS engine builds a combined automaton that exceeds memory. **Fixed by chunking: 50 motifs per batch**, with progress logging. Each chunk scans, extracts hits, frees memory. Gene body regions (~500 MB sequence) are the slowest.

### Batch 01: GC% added to region and TE class stats
Region stats and TE class stats now include `gc_pct` column alongside CpG O/E and density.

### Batch 02: Absent genes shown
Presence/absence heatmap now includes 8 **absent** genes (DNMT3A, DNMT3B, DNMT3L, DNMT2, TET1, MBD1, MBD3, MECP2) alongside the 18 present genes. Categories include "Writer (de novo)" to highlight the missing DNMT3 family.

### Batch 04: BSseq construction + expanded analysis
If `CACHE$bsseq` doesn't exist, batch04 builds it: reads all 4 CpG reports, strand-collapses (minus-strand pos -1, sum counts), filters cov≥5, intersects common sites across all samples. Metagene profile uses 20 bins + 5kb flanks. Expression from `CACHE$transcriptome` (fallback: DESeq2 tail controls). Correlation by region (promoter, exon, intron, downstream) and by expression decile. Also: genome-wide 1 Mb windows (fig4i), CpG region pie (fig4j), single-exon metagene (fig4k), first/last exon + first intron expanded correlation (fig4l).

### Batch 05: DMP/DMR overlap conditional
Loads batch06 DMPs/DMRs if available for TE overlap analysis. Cohen's d effect sizes per TE class. Age conversion: 1% Kimura ≈ 2.27 My. Uses GAM instead of loess for age plots (all data, no subsampling).

### Batch 06: Gene names + DMR parity
EviAnn annotation (`OG$annot`) merged into DMPs and DMRs — `gene_name` column flows to all downstream. DMP heatmap shows gene name row labels. DMP volcano has labeled version (top 30, `ggrepel`). DMR volcano uses `-log10(p)` from areaStat: `Z = areaStat/sqrt(nCG)` → `p = 2*pnorm(-|Z|)`. DMR enrichment by region (Fisher's exact, like DMP version). DMR top genes barplot. GO/KEGG/Reactome via `clusterProfiler` + STRING enrichment terms. Background annotation chunked in 5M blocks. Hyper vs Hypo GO BP run separately.

### Batch 07: Window-based gene assignment + concordant gene lists
**Gene assignment via 5kb windows**: DMPs/DMRs assigned to ALL genes whose extended window they fall within (promoter 2kb + distal upstream 5kb + gene body + downstream 5kb). A DMP between two genes can be assigned to BOTH (e.g., KANK1/sox19a locus). Expression from `CACHE$transcriptome` (apeglm shrinkage, T1S5 removed). Expanded annotation for first/last exon, intron-level correlation. Dose-response: DMP count bins vs DE rate. **Exports concordant gene lists**: genes with both significant methylation (>=3 DMPs or DMR) AND significant expression change (padj<0.05), classified into quadrants (Hyper_Up, Hyper_Down, Hypo_Up, Hypo_Down). These are cloning targets for batch09 motif analysis.

### Batch 08: WGCNA outlier detection
Uses **within-tissue centroid distance**: for each tissue group (≥3 samples), computes Euclidean distance to tissue centroid. Samples with distance > 2× median within-tissue distance are flagged as outliers and removed before module detection. Power = 12. Multi-tissue: tail (control/amputated), eye (control/amputated), bodywall (control/irradiated/fungicide), head, juvenile, ovotestis (control). Fisher's exact test for DMP enrichment per module. Hub genes identified. GO BP enrichment per module.

### Batch 10 + 10b: Entropy analysis
- **Level 1** (batch 10): Per-site binary Shannon entropy from BSseq: `H = -(beta * log2(beta) + (1-beta) * log2(1-beta))`. Per-region delta, Wilcoxon paired test, mutual information (meth diff vs LFC).
- **Level 2** (batch 10b, cluster only): Per-read NME from Bismark BAMs. Xie et al. 2011 method: 4-CpG sliding windows, epiallele frequencies, NME = H/4. Requires 64+ GB RAM, separate submission. BAMs at `/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/08_deduplication/`. **Note:** C1 BAM is unsorted — needs `samtools sort` before indexing.
- **Critical control**: Compare DMPs to **beta-matched non-DMPs**, not all CpGs (sites near 0/1 have low entropy by definition). Also stratify by CpG density.
- **Drift vs reprogramming test**: DMPs from low-entropy baseline → directed reprogramming (rejecting null for DNMT3-absent organism). High entropy → stochastic drift.
- Full scientific rationale: see PROTOCOL.md.

## Cross-species validation

See PROTOCOL.md for the full 25-point table mapping comparisons to responsible batches. Key species: Nautilus pompilius (Wu 2025), annelids (Guynes 2024), Weber/Schubeler promoter architecture.

## Critical rules

1. **Strand collapse.** CpG_report.txt has both strands. Minus-strand pos shifted by -1, counts summed. One measurement per CpG dinucleotide.
2. **NOT a mammal.** No CpG islands, no promoter silencing, no DNMT3.
3. **Filter scaffolds.** `keep_chr <- paste0("chr", 1:31)` always.
4. **TE is NOT a functional region.** Report as separate boolean column, not in pie charts.
5. **GenomicRanges for annotation.** NOT ChIPseeker. Priority: Promoter > Exon > Intron > Intergenic.
6. **TFBS motifs only.** JASPAR 2024 CORE metazoan. All are transcription factor binding sites.
7. **Let data speak.** Code asks questions, doesn't confirm assumptions.
8. **No HTML in batch scripts.** Data + figures only.
9. **PDF + PNG for every figure.** `cairo_pdf` + `png(res=300)`.
10. **All numbers from running code.** Never paste from comments.
11. **Vectorize everything.** Never use for-loops on large data. Use `findOverlaps()`, `data.table` joins/grouped ops, and bulk GRanges construction instead of row-by-row iteration.
12. **`dos2unix` before cluster runs.** Scripts edited on Windows have `\r\n` line endings that break `#!/usr/bin/env Rscript` and SLURM. Always run `dos2unix` on `.R`, `.slurm`, and `_config.R` before `sbatch`.
13. **Never subsample data.** Use ALL data points for every analysis and plot. If a method is too slow (loess O(n^2)), use a faster method (GAM, binned means) — never drop data to make a method work.
14. **HTML report must be self-contained.** `generate_report.R` embeds all figures as base64 data URIs inside the HTML. No external file path references — the HTML must work when moved or shared standalone.

## Running locally (Windows)

```bash
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" methylation_pipeline/batchNN/code/NN_script.R
```

BSgenome package: `BSgenome.Dlaeve.NCBI.dlgm` (installed on Windows R)

## Running on cluster (HPC)

### Full pipeline (batch 1.5 through 10)
```bash
cd Deroceras-Leave
git pull
dos2unix methylation_pipeline/*.slurm methylation_pipeline/_config.R
sbatch methylation_pipeline/run_pipeline.slurm   # runs 0.5→1.5→02→03→...→10 sequentially
```
**Note:** Pipeline requests 64 GB (batch 1.5 motif scanning needs it). Batch 01 is NOT in the runner (genome cache must pre-exist). Batch 10b (per-read NME) requires separate submission.

### HTML figure report
```bash
Rscript generate_report.R   # produces pipeline_report.html at repo root
```
Scans all `methylation_pipeline/batch*/figures/*.png` and builds a styled HTML with analysis notes per figure. Re-run after pipeline completes.

`_config.R` auto-detects cluster vs local — all paths in `OG` list switch automatically.

**NOTE:** CpG reports are `.txt.gz` on cluster (compressed). `data.table::fread()` reads gzipped files natively.

### Individual batches
```bash
cd Deroceras-Leave
Rscript methylation_pipeline/batch03/code/03_promoter_cpg_classification.R
```

### Prerequisite: genome cache
If `genome/cache/genome_chr1_31.rds` doesn't exist, batch01 will generate it from BSgenome on first run. Or run separately:
```bash
sbatch cluster/scripts/generate_genome_cache.slurm
```

### Other cluster jobs
```bash
cd Deroceras-Leave/cluster
sbatch scripts/dmltest_full.slurm            # DMLtest with smoothing (32 GB)
sbatch scripts/batch1.5_expanded_motif.slurm # TFBS genome-wide positions (64 GB)
```

## TODO

### Folder renaming (deferred)
- [ ] Rename batch folders: `batch01` → `Batch 01 — Genomic CpG Landscape`, etc. Requires updating ~36 path references across all scripts + run_pipeline.slurm + CLAUDE.md.

### Artifact notes
- Old Metilacion/ DMPs (18,755) used **non-strand-collapsed data** — pipeline uses correct strand-collapsed (17,729)
- Old MXT analysis correlation sign flip vs pipeline is noise around zero (both ~0), not a real discrepancy
- ChIPseeker annotation differs from GenomicRanges — we use GenomicRanges (standard, no external DB needed)

## Git

- Remote: `https://github.com/somnasgarden/Deroceras-Leave.git`
- Author: `GIT_AUTHOR_NAME="rlopezt" GIT_AUTHOR_EMAIL="rafae@users.noreply.github.com"`
- HPC clone: `/mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave/`
