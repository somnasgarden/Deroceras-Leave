# Methylation Pipeline Protocol

**Internal document — not for reviewers.** This guides development of the `methylation_pipeline/` code.

The `methylation_pipeline/` folder IS the supplementary code for the paper. Reviewers will read it. It must:
- Produce ALL data tables and figures used in the paper
- Be self-contained and reproducible (source `_config.R`, use cached data)
- Generate NO HTML — only TSV data + PNG/PDF figures
- Use strand-collapsed CpG data throughout
- Use GenomicRanges for annotation (NOT ChIPseeker)
- Let the data answer questions — don't hardcode expected results

## Overview

This pipeline analyzes whole-genome bisulfite sequencing (WGBS) data from *Deroceras laeve* tail regeneration. Each batch is a self-contained analysis that produces data (TSV) and publication-quality figures (PNG + PDF).

## Execution order

Batches must run in order. Each depends on outputs from previous batches.

```
Batch 01 (genome CpGs) → Batch 02 (machinery) → Batch 03 (Weber promoters)
                                                          ↓
Batch 1.5 (TFBS motifs) ──────────────→ Batch 04 (baseline meth+expr) → Batch 05 (TE+age)
                                                          ↓
                                          Batch 06 (DMP/DMR) → Batch 07 (decoupling)
                                                          ↓
                                          Batch 08 (WGCNA) → Batch 09 (TF+GENIE3) → Batch 10 (entropy)

Data flow:
  OG Data (CpG_report.txt, GFF, HTSeq counts) → genome/cache/*.rds
                                                      ↓
                                              _config.R (shared paths)
                                                      ↓
                                              batchNN/code/*.R
                                                      ↓
                                        batchNN/data/*.tsv  +  batchNN/figures/*.{png,pdf}
```

## Running a batch

```bash
cd STANDBY
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" methylation_pipeline/batchNN/code/NN_script_name.R
```

All scripts source `methylation_pipeline/_config.R` for shared paths, colors, and helper functions.

## Key procedures

### Strand collapse (CRITICAL)
CpG_report.txt files report both strands separately. We collapse them:
- Minus-strand CpG at position X is the same dinucleotide as plus-strand at X-1
- Shift minus-strand positions by -1, then sum methylated/unmethylated counts
- This produces **one measurement per CpG dinucleotide**
- Result: 25,439,068 strand-collapsed CpGs (from 45,535,149 genome CpGs)

### Functional region annotation (NOT ChIPseeker)
Mutually exclusive categories using GenomicRanges:
- Priority: Promoter > Exon > Intron > Intergenic
- TE overlap reported as separate boolean column (TEs overlap with all regions)
- Enrichment tested with Fisher's exact test against all tested CpGs as background

### Weber promoter classification
- Window: TSS +/- 1 kb (Weber et al., 2007)
- HCP: CpG O/E >= 0.75 AND GC% >= 55%
- LCP: CpG O/E < 0.48
- ICP: everything else
- D. laeve result: 91.6% ICP (no CpG islands)

### DSS differential methylation
- DMLtest with smoothing=TRUE (Bayesian shrinkage)
- DMPs: FDR < 0.05, |methylation difference| > 10%
- DMRs: >= 3 CpGs, >= 50 bp
- Strand-collapsed CpG data input

### TFBS annotation
- Database: JASPAR 2024 CORE, metazoan motifs (1,362 PWMs)
- All motifs are transcription factor binding sites
- Scan scope: promoters (2kb upstream) locally; full genome on cluster
- Tool: motifmatchr with p-cutoff = 5e-5
- Output must include: motif name, consensus sequence, chr, start, end, nearest gene, distance to TSS

### Mutual information
- 2x2 contingency table: methylation direction x expression direction
- MI computed from joint vs marginal probabilities
- Significance: 10,000 permutations shuffling expression labels
- Reports MI in bits with permutation p-value

## Batch descriptions

### Batch 01: Genomic CpG Annotation
- **Input**: Genome FASTA (from BSgenome), GFF, TE annotations
- **Output**: Region CpG stats, TE class stats, genome summary
- **Figures**: Dinucleotide freq, CpG O/E by region, CpG density by region, per-chr distributions, TE class O/E

### Batch 02: Methylation Machinery
- **Input**: GFF, HTSeq count matrices
- **Output**: Methylation toolkit expression table
- **Figures**: Presence/absence grid, expression boxplot, z-scored heatmap
- **Key finding**: DNMT1 present, no DNMT3

### Batch 03: Promoter CpG Classification
- **Input**: Genome, GFF (gene TSS positions)
- **Output**: Per-promoter O/E + GC% + Weber class
- **Figures**: O/E histogram, O/E vs GC% scatter, class proportions, faceted histograms
- **Key finding**: 91.6% ICP, no CpG islands

### Batch 04: Baseline Methylation Landscape
- **Input**: BSseq cache, DESeq2 counts
- **Output**: Per-sample stats, region methylation, gene body meth vs expression
- **Figures**: Metagene profile, meth vs expression scatter, expression by meth quartile, density per sample
- **Key finding**: Gene body meth positively correlated with expression (rho = 0.307)

### Batch 05: TE Methylation + Age
- **Input**: BSseq cache, TE age data (Kimura divergence)
- **Output**: Per-TE methylation by age, intronic vs intergenic comparison
- **Figures**: Age vs methylation by TE class, intronic vs intergenic ridge plots
- **Key finding**: Younger TEs more methylated, intronic > intergenic, no change upon amputation

### Batch 06: Differential Methylation
- **Input**: BSseq cache (or DMLtest cache)
- **Output**: Annotated DMPs, DMRs, enrichment table
- **Figures**: PCA, correlation heatmap, chromosome methylation, DMP volcano, DMR volcano, annotation bar, enrichment bar
- **Key finding**: 17,729 DMPs, 1,462 DMRs, exon enriched 1.95x

### Batch 07: Methylation-Expression Decoupling
- **Input**: Batch 06 DMPs, DESeq2 DE results
- **Output**: Gene-level meth vs expression, MI statistics, dose-response table
- **Figures**: Delta-meth vs delta-expression scatter, dose-response bar, contingency heatmap
- **Key finding**: MI = 0.000087 bits (p = 0.55), rho = -0.021 (p = 0.24)

### Batch 08: WGCNA + Methylation
- **Input**: Multi-tissue HTSeq counts, Batch 06 DMPs
- **Output**: Module assignments, module DMP burden
- **Figures**: Module DMP burden bar, within-module correlations, module enrichment
- **Key finding**: Simpson's paradox — module-level effects cancel globally

### Batch 09: TF + Motif + GENIE3 Integration
- **Input**: DeepTFactor predictions, JASPAR motif hits, GENIE3 network, Batch 06 DMPs, Batch 07 meth-expression table
- **Output**: TF-DMP enrichment, DMP-TFBS overlap, TF network nodes, motif-near-DMP-near-DE-gene list
- **Figures**: TF census, motif prevalence, cross-reference, Sox19a locus
- **Key finding**: TFs enriched for DMPs (OR = 1.30) but none DE
- **Critical check**: For genes where a TFBS has a DMP nearby AND the gene is DE, cross-reference with the Batch 07 methylation-expression table. This is the subset where methylation could actually be functional — if a TFBS is differentially methylated AND the target gene changes expression, that's evidence of regulatory methylation. Compare this subset to the global MI = 0 finding.

### Batch 10: Entropy Analysis
- **Input**: BSseq cache
- **Output**: Per-CpG entropy, entropy by region
- **Figures**: Entropy distribution (control vs amputated), per-region entropy, entropy vs expression
- **Key question**: Does regeneration increase or decrease methylation entropy?

## Data flow

```
OG Data (CpG reports, GFF, counts) → genome/cache/ (RDS objects)
                                          ↓
                                    Batch scripts
                                          ↓
                               batchNN/data/ (TSV)
                               batchNN/figures/ (PNG + PDF)
```

## Key numbers (strand-collapsed, verified)

| Metric | Value |
|--------|-------|
| Genome size | 1,783,141,715 bp |
| Chromosomes | 31 (chr1-31, scaffolds excluded) |
| GC content | 42.99% |
| Total CpGs | 45,535,149 |
| CpG O/E | 0.553 |
| CpGs tested (cov >= 5, all samples) | 25,439,068 |
| DMPs (FDR < 0.05, |diff| > 10%) | 17,729 |
| DMRs | 1,462 |
| Predicted TFs | 2,217 genes |
| JASPAR motifs scanned | 1,362 |
| Promoter classes | 91.6% ICP, 6.5% HCP, 1.9% LCP |

## Expected biological findings (to be verified by code)

These are expected from the dlaeve project, but each batch must independently confirm or reject them:

1. **MI = 0** between methylation and expression direction. Methylation does not predict expression.
2. **Intergenic/intronic DMPs dominate** (~65% of all DMPs). Not promoter-driven.
3. **TF methylation paradox**: TFs enriched for DMPs but none are DE.
4. **Module-specific structure**: Global correlation ~0 but within-module correlations can be strong (Simpson's paradox).
5. **91.6% ICP promoters**: CpG island model is inapplicable.
6. **Binary switch**: No dose-response between DMP count and DE probability.
7. **TE methylation unchanged**: 53% of DMPs overlap TEs but TE methylation itself doesn't change.
8. **Gene body meth ~ expression**: Baseline positive correlation (rho ~0.3).

## Comparison with dlaeve-regeneration-epigenomics project

A previous analysis (`C:/Users/rafae/Projects/dlaeve-regeneration-epigenomics/`) ran the same data but with two methodological issues:
1. **No strand collapse** — CpG_report.txt strands kept separate, inflating CpG count (24.8M vs our 25.4M collapsed)
2. **ChIPseeker for annotation** — wrong tool for WGBS (designed for ChIP-seq peaks)

Their numbers (18,754 DMPs, 1,424 DMRs) are close but include strand duplicates. Our strand-collapsed numbers (17,729 DMPs, 1,462 DMRs) are cleaner. The biological conclusions are the same — the bugs don't change the story.

Their CLAUDE.md has useful elements: dependency graph, 21-script table, key findings list. But their pipeline is harder to reproduce (WSL-only, hardcoded cluster paths, HTML mixed with analysis).

## Reference papers (in papers/ folder)

- **Weber et al. 2007 (ng1990)**: Promoter HCP/ICP/LCP classification. TSS +/- 1kb window.
- **Schubeler 2015 (nature14192)**: Review. HCPs unmethylated by default. Methylation at HCPs is consequence not cause.
- **Guynes et al. 2024 (Genome Biology)**: Annelid methylomes. GbM ancestral. Methylation erodes with age. 17-21% CpG meth.
- **Wu et al. 2025 (BMC Genomics)**: Nautilus pompilius. Fellow mollusk. CpG meth 21.2%. Has DNMT1+DNMT3. GbM rho=0.32.
- **Results_code_guide.pdf**: Master plan defining what each batch must produce.

## Critical rules

1. **ALL numbers from running code.** Never paste from comments.
2. **NOT a mammal.** No CpG islands, no promoter silencing, no DNMT3.
3. **ALL code in R.** Except HOMER (Perl) for de novo motifs.
4. **Filter scaffolds.** `keep_chr <- paste0("chr", 1:31)` always.
5. **Strand collapse before analysis.** Two strand positions = one CpG.
6. **TE is NOT a region.** It's a separate annotation layer.
7. **Let data speak.** Code asks questions, doesn't confirm assumptions.
8. **No HTML in batch scripts.** Batches produce data + figures only.
9. **TFBS motifs only.** All JASPAR motifs are TF binding sites. Filter anything else.
10. **Motif output must include:** motif name, consensus sequence, chr, start, end, nearest gene, distance to gene.
