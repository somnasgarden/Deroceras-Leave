# Methylation Pipeline Protocol

**Internal document for the research team and any Claude agent working on this project.** Not for reviewers. The `methylation_pipeline/` folder is what reviewers see — it must be clean, reproducible code that produces all paper figures and data.

---

## Project goal

Test whether DNA methylation regulates gene expression during tail regeneration in *Deroceras laeve*. The answer from the dlaeve project (which ran the same data with a strand-collapse bug) is: **no, at least not globally.** MI = 0, rho = 0. But module-specific effects exist (Simpson's paradox). This pipeline must independently verify or refute that finding using strand-collapsed data.

---

## The organism

*Deroceras laeve* (gray field slug, Stylommatophora, Gastropoda, Mollusca). Pulmonate gastropod capable of tail regeneration with blastema formation by 2 dpa.

**Genome paper:** Miranda-Rodriguez et al. 2025, G3 (jkaf164). Assembly: Verkko (PacBio HiFi + ONT + Hi-C). 2.05 Gb total, **1.78 Gb on 31 chromosomes**. Unplaced scaffolds (~200 Mb) are repeats + bacterial contamination — ALWAYS EXCLUDED (`keep_chr <- paste0("chr", 1:31)`). BUSCO 93.5% (metazoa). 24,253 protein-coding genes (EviAnn). Repeats: ~40% interspersed (LINEs dominant, RTE family), 72% total with tandem repeats. 178 miRNA loci, 170 piRNA clusters.

**Key biological features:**
- Has DNMT1 (maintenance methyltransferase) but **NO DNMT3** (no de novo methyltransferase)
- Gene body methylation pattern (exons > introns > promoters > intergenic)
- 91.6% ICP promoters — NO CpG islands
- Mosaic methylome typical of invertebrates (not global like mammals)

---

## Experimental design

### WGBS (methylation)
- 2 control (C1, C2) vs 2 amputated (A1, A2) tail samples, 2 dpa
- Bismark alignment → CpG_report.txt (7 columns: chr, pos, strand, count_meth, count_unmeth, context, trinuc)
- **STRAND COLLAPSE MANDATORY:** Both strands report the same CpG dinucleotide. Minus-strand position shifted by -1, counts summed. Result: 25,439,068 strand-collapsed CpGs with cov >= 5 in all 4 samples.
- DSS for differential methylation (NOT MethylKit). DMLtest with smoothing=TRUE.

### RNA-seq (expression)
- Multi-tissue: tail, eye, bodywall, head, juvenile, ovotestis
- Multi-condition: control, amputated, irradiated, fungicide
- Tail contrast for DE: C1S1-C4S4 (control) vs T2S6-T4S8 (amputated). T1S5 removed as WGCNA outlier.
- DESeq2 for DE. WGCNA for co-expression modules.

---

## Pipeline structure

```
methylation_pipeline/
  _config.R               — shared paths, colors, helpers (source this in every batch)
  batch01/ {code/,data/,figures/}  — Genomic CpG annotation
  batch02/ {code/,data/,figures/}  — Methylation machinery (18 genes, DNMT1 only)
  batch03/ {code/,data/,figures/}  — Promoter Weber classification (HCP/ICP/LCP)
  batch04/ {code/,data/,figures/}  — Baseline methylation + gene body meth vs expression
  batch05/ {code/,data/,figures/}  — TE methylation + evolutionary age (Kimura)
  batch06/ {code/,data/,figures/}  — Differential methylation (DMP/DMR)
  batch07/ {code/,data/,figures/}  — Methylation-expression decoupling (MI, correlation)
  batch08/ {code/,data/,figures/}  — WGCNA modules + methylation enrichment
  batch09/ {code/,data/,figures/}  — TF + GENIE3 + motif integration
  batch10/ {code/,data/,figures/}  — Entropy analysis
```

### Dependency graph

```
Batch 01 (genome CpGs) → Batch 02 (machinery) → Batch 03 (Weber promoters)
                                                          ↓
Batch 1.5 (TFBS motifs) ──────────────→ Batch 04 (baseline meth+expr) → Batch 05 (TE+age)
                                                          ↓
                                          Batch 06 (DMP/DMR) → Batch 07 (decoupling)
                                                          ↓
                                          Batch 08 (WGCNA) → Batch 09 (TF+GENIE3) → Batch 10 (entropy)
```

### Rules for batch scripts
1. **No HTML.** Batches produce TSV data + PNG/PDF figures only. Reviewers read the code.
2. **Source `_config.R`.** Every batch starts with `source("methylation_pipeline/_config.R")`.
3. **Save with helpers.** Use `save_fig(p, BATCH_DIR, "name")` and `save_data(df, BATCH_DIR, "name")`.
4. **Let data speak.** Code asks questions and reports answers. Don't hardcode expected results.
5. **Strand-collapsed data throughout.** The BSseq cache (`bsseq_tutorial.rds`) is already strand-collapsed.

---

## Batch details

### Batch 01: Genomic CpG Annotation
- **Question:** What does the D. laeve CpG landscape look like?
- **Input:** Genome (BSgenome.Dlaeve.NCBI.dlgm or genome_chr1_31.rds), GFF, TE data
- **Output:** Region CpG stats, TE class stats, genome summary
- **Figures:** (A) Dinucleotide freq, (B) CpG O/E by region, (C) CpG density by region, (D) Per-chr distributions, (E) TE O/E by class, (F) TE per-chr distributions
- **Method:** CpG positions found with `matchPattern("CG", genome)`. O/E = (n_CpG × len) / (n_C × n_G). Per-chr distributions filtered to >= 10 kb per region (removes outliers from small intergenic or SINE regions on some chromosomes).

### Batch 02: Methylation Machinery
- **Question:** Does D. laeve have a functional methylation toolkit?
- **Input:** GFF, HTSeq count matrices, EviAnn annotations
- **Output:** 18-gene expression table
- **Figures:** (A) Presence/absence, (B) Expression boxplot, (C) Z-scored heatmap
- **Key finding expected:** DNMT1 present, NO DNMT3. All 18 genes expressed. DNMT1 highest in ovotestis.

### Batch 03: Promoter CpG Classification (Weber)
- **Question:** Does D. laeve have CpG islands?
- **Input:** Genome, GFF (gene TSS)
- **Output:** Per-promoter O/E, GC%, Weber class
- **Figures:** (A) O/E histogram with thresholds, (B) O/E vs GC% scatter, (C) Class proportions, (D) Faceted histograms by class
- **Method:** Weber et al. 2007 criteria on TSS +/- 1 kb window:
  - HCP: CpG O/E >= 0.75 AND GC% >= 55%
  - LCP: CpG O/E < 0.48
  - ICP: everything else
- **Key finding expected:** 91.6% ICP, no CpG islands. Totally unlike mammals (72% HCP).

### Batch 04: Baseline Methylation Landscape
- **Question:** Where is methylation? Does gene body meth correlate with expression?
- **Input:** BSseq cache, DESeq2 counts (tail controls)
- **Output:** Per-sample stats, region methylation, gene body meth vs expression
- **Figures:** (A) Metagene profile (upstream → gene body → downstream), (B) Meth vs expression scatter, (C) Expression by meth quartile, (D) Methylation density per sample
- **Method:** Mean beta per gene body (control samples) vs VST expression. Spearman correlation.
- **Key finding expected:** rho ~ 0.3, exon-peaked methylation, bimodal beta distribution.
- **Reference plots:** Nautilus pompilius Fig 2B (Wu et al. 2025), Annelid metagene Fig 3 (Guynes et al. 2024).

### Batch 05: TE Methylation + Age
- **Question:** Is TE methylation related to evolutionary age? Intronic vs intergenic?
- **Input:** BSseq cache, TE age data (Kimura divergence)
- **Output:** Per-TE methylation, intronic vs intergenic comparison
- **Figures:** (A) Age vs methylation by TE class, (B) Intronic vs intergenic ridge plots
- **Method:** Kimura divergence → millions of years (1% ≈ 2.27 My). Per-TE mean beta (control). Cohen's d for change upon amputation.
- **Key finding expected:** Younger TEs more methylated (passive gene body meth, not active silencing). TE meth unchanged by amputation (Cohen's d < 0.01). Without DNMT3, no de novo TE silencing.

### Batch 06: Differential Methylation
- **Question:** What changes during regeneration?
- **Input:** BSseq cache, DMLtest cache (or run DMLtest)
- **Output:** Annotated DMPs, DMRs, enrichment table
- **Figures:** (A) PCA, (B) Correlation heatmap, (C) Chromosome methylation, (D) DMP volcano, (E) DMR volcano, (F) Annotation bar, (G) Enrichment vs background, (H) Ontology enrichment
- **Method:** DSS DMLtest with smoothing=TRUE. DMPs: FDR < 0.05, |diff| > 10%. DMRs: >= 3 CpGs, >= 50 bp. Annotation: Promoter > Exon > Intron > Intergenic (GenomicRanges, NOT ChIPseeker). TE as separate boolean. Enrichment: Fisher's exact vs all 25.4M CpGs.
- **Key finding expected:** 17,729 DMPs, 1,462 DMRs. Exon 1.95x enriched, promoter 1.44x. TE not enriched.

### Batch 07: Methylation-Expression Decoupling
- **Question:** Do methylation changes predict expression changes?
- **Input:** Batch 06 DMPs, DESeq2 DE results
- **Output:** Gene-level table, MI statistics, dose-response table
- **Figures:** (A) Delta-meth vs delta-expression scatter, (B) Dose-response bar, (C) Contingency heatmap
- **Method:** 2x2 table (hyper/hypo × up/down). MI from joint vs marginal probabilities. 10,000 permutations. Spearman correlation. Cochran-Armitage trend test for dose-response.
- **Key finding expected:** MI ≈ 0, rho ≈ 0, no dose-response. Binary switch, not rheostat.

### Batch 08: WGCNA + Methylation
- **Question:** Is the global null masking module-specific effects?
- **Input:** Multi-tissue HTSeq counts, Batch 06 DMPs
- **Output:** Module assignments, module DMP burden, within-module correlations
- **Figures:** (A) Module DMP burden, (B) Within-module correlations, (C) Module enrichment
- **Method:** WGCNA signed network, soft-thresholding, dynamic tree cut. Per-module Fisher's test for DMP enrichment. Per-module Spearman meth vs expression.
- **Key finding expected:** 12 modules. Global rho ≈ 0 but within-module |r| up to 0.35-0.81. Simpson's paradox.

### Batch 09: TF + Motif + GENIE3 Integration
- **Question:** Do TFBS near DMPs explain regulatory connections?
- **Input:** DeepTFactor, JASPAR motif hits, GENIE3, Batch 06 DMPs, **Batch 07 meth-expression table**
- **Output:** TF-DMP enrichment, DMP-TFBS overlap, TF network nodes, motif-near-DMP-near-DE-gene list
- **Figures:** (A) TF census, (B) Motif prevalence, (C) Cross-reference, (D) Sox19a locus
- **Method:** Fisher's test for TF DMP enrichment. Overlap DMPs with TFBS coordinates. Cross-reference: for genes where a TFBS has a DMP nearby AND the gene is DE, compare with Batch 07 methylation-expression table. This subset is where methylation could be functional.
- **Critical check:** GENIE3 regulatory edges where the TF has a DMP → does the target change expression? This is the mechanistic test.
- **Key finding expected:** TFs enriched for DMPs (OR ~ 1.3) but none DE. Sox19a: 13 DMPs 4-5 kb upstream in intergenic region.

### Batch 10: Methylation Entropy
- **Question:** Does regeneration increase or decrease methylation entropy?
- **Input:** BSseq cache (for per-site entropy). BAM files needed for per-read entropy (cluster).
- **Output:** Per-CpG entropy, entropy by region, entropy at DMPs vs non-DMPs
- **Figures:** (A) Entropy distribution control vs amputated, (B) Per-region entropy, (C) Entropy at DMPs vs background

#### Two types of methylation entropy

**Per-site (binary) entropy** — what we can compute from CpG_report.txt:
```
H = -(beta * log2(beta) + (1-beta) * log2(1-beta))
```
Maximum at beta = 0.5, zero at 0 or 1. Measures uncertainty about each CpG's state.

**Per-read (epiallele) entropy** — needs BAM files:
```
H = -sum(p_i * log2(p_i)) for all 2^k epiallele patterns in a 4-CpG window
```
Captures co-methylation structure across neighboring CpGs on individual DNA molecules. More informative but requires aligned reads, not just counts.

#### Key question
**Do DMPs come from high-entropy (already disordered) or low-entropy (stably committed) sites?**
- High-entropy DMPs → stochastic drift (noise, passive demethylation during replication)
- Low-entropy DMPs → directed reprogramming (active, targeted methylation changes)
This distinguishes "the methylome is falling apart" from "the methylome is being specifically reprogrammed."

#### Literature
- **No published study has computed methylation entropy in regeneration or invertebrates.** Both would be novel.
- All published entropy papers used BAM files for per-read entropy. Per-site entropy is a valid approximation.
- Check if BAMs exist on cluster at `/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/`

---

## Key numbers (strand-collapsed, verified)

| Metric | Value |
|--------|-------|
| Genome size (chr1-31) | 1,783,141,715 bp |
| GC content | 42.99% |
| Chromosomes | 31 |
| Total CpGs | 45,535,149 |
| CpG O/E | 0.553 |
| CpGs tested (cov>=5, all 4 samples) | 25,439,068 |
| DMPs (FDR<0.05, \|diff\|>10%) | 17,729 |
| Hypermethylated / Hypomethylated | 8,285 / 9,444 |
| DMRs | 1,462 |
| DMR mean length / CpGs | 264 bp / 7.5 |
| DMP enrichment: Exon | 1.95x (p < 0.001) |
| DMP enrichment: Promoter | 1.44x (p < 0.001) |
| DMP enrichment: Intergenic | 0.81x (p < 0.001, depleted) |
| TE overlap of DMPs | 53.4% (no enrichment, fold=1.01) |
| Predicted TFs (DeepTFactor) | 2,217 genes (4,824 isoforms) |
| JASPAR TFBS motifs | 1,362 metazoan PWMs |
| Promoter classes | 91.6% ICP, 6.5% HCP, 1.9% LCP |
| Baseline gene body meth vs expression | rho = 0.307 (p < 10^-80) |
| Differential meth vs expression MI | 0.000087 bits (p = 0.55) |
| Differential meth vs expression rho | -0.021 (p = 0.24) |

---

## Comparison with dlaeve-regeneration-epigenomics project

A previous analysis (`C:/Users/rafae/Projects/dlaeve-regeneration-epigenomics/`) ran the same raw data. Key differences:

| Issue | dlaeve (other) | STANDBY (ours) |
|-------|---------------|----------------|
| Strand collapse | **NO** (bug at line 150 of 01_quality_control.R) | **YES** |
| CpGs | 24,771,877 (inflated by dual strands) | 25,439,068 (true collapsed) |
| DMPs | 18,754 (may include strand duplicates) | 17,729 (clean) |
| DMRs | 1,424 | 1,462 |
| Annotation tool | ChIPseeker (wrong for WGBS) | GenomicRanges (correct) |
| Region categories | Promoter / Gene Body / Intergenic | Promoter / Exon / Intron / Intergenic |
| TE in annotation | Not reported separately | Separate boolean column |
| HTML in code | Mixed with analysis | NOT in batch scripts |

**Biological conclusions are the same** — the strand-collapse bug doesn't change the story. But our numbers are cleaner and more defensible for the paper.

**Useful from dlaeve project:** Pipeline dependency graph, 21 local analysis scripts, WGCNA results, paper drafts (v3, v4), interactive HTML report. These can be referenced but should not be copied without fixing the strand issue.

---

## Bibliography context (key papers, short summaries)

### Genome & organism
- **Miranda-Rodriguez et al. 2025 (G3, jkaf164):** D. laeve genome assembly. Verkko + Hi-C. 1.78 Gb, 31 chr. 24,253 genes. BUSCO 93.5%. ~40% interspersed repeats (LINEs dominant). Scaffolds = contamination/repeats, excluded. This is OUR genome paper.

### Promoter classification
- **Weber et al. 2007 (Nature Genetics, ng1990):** Defined HCP/ICP/LCP classification. TSS +/- 1 kb. HCP: O/E >= 0.75 + GC >= 55%. LCP: O/E < 0.48. ICP: rest. In humans: 72% HCP, 16% ICP, 12% LCP. HCPs mostly unmethylated regardless of expression.
- **Schubeler 2015 (Nature, nature14192):** Review. HCPs protected from methylation by CFP1/KDM2B (CxxC domain proteins recruit H3K4 methyltransferases). Methylation at HCPs is CONSEQUENCE of silencing, not cause. LCP methylation more functionally relevant. Invertebrates show gene body methylation correlated with constitutive expression — opposite of mammalian promoter silencing.

### Comparative invertebrate methylation
- **Guynes et al. 2024 (Genome Biology, s13059-024-03346-z):** Annelid methylomes (3 species). GbM ancestral to bilaterians. Methylation erodes during development + aging (global 5mC decreases from embryo to adult). 17-21% CpG methylation. Gene body methylation correlates with expression + transcriptional stability. Strongest in early development, decays with age. D. gyrociliatus lost methylation entirely despite having DNMT1.
- **Wu et al. 2025 (BMC Genomics, s12864-025-11907-0):** Nautilus pompilius methylome. Fellow mollusk (cephalopod). CpG meth 21.2%. Has BOTH DNMT1 and DNMT3 (unlike D. laeve which lost DNMT3). Gene body meth rho = 0.32 with expression. Tissue-specific genes are hypomethylated. Promoter meth weakly negatively correlated with expression (rho = -0.05).
- **Haidar et al. 2024:** Review of mollusk methylation. DNMT1 conservation, MBD2/3 readers, mosaic gene body pattern. Many gastropods lost DNMT3.

### Methylation entropy
- **Xie et al. 2011 (NAR):** Defined per-read methylation entropy. 4-CpG windows, 2^4 epiallele patterns, Shannon entropy. **Needs BAMs.** The foundational method.
- **Landan et al. 2012 (Nature Genetics):** Epigenetic polymorphism increases in culture/aging. Coined "epipolymorphism" (Simpson's index on epiallele patterns). Stochastic DNMT1 errors accumulate. **Used BAMs.**
- **Landau et al. 2014 (Cancer Cell):** CLL has dramatically increased methylation entropy. PDR (Proportion of Discordant Reads) + Shannon entropy. High entropy = transcriptional variability. **Used BAMs.**
- **Jenkinson et al. 2017 (Nature Genetics):** Ising model (informME). Both stem cells AND differentiated tissues have LOW entropy. Entropy increases with aging + cancer, NOT with differentiation per se. **Used BAMs.**
- **Shao et al. 2014 (BMC Genomics):** Entropy DECREASES at promoters during adipose stem cell → adipocyte differentiation. Entropy INCREASES during reprogramming to iPSC. **Used BAMs.**
- **Vaidya et al. 2023 (Genome Biology):** "Aging but NOT differentiation leads to increasing entropy." Very little entropy change between stem cells and differentiated progeny at same age. **Used BAMs (RRBS).**
- **NO PUBLISHED STUDY on methylation entropy in regeneration or invertebrates.** Both would be novel contributions.

### TFBS motif annotation
- **Hashim et al. 2019 (AJMB-11-130):** Review of 119 motif discovery algorithms. JASPAR highlighted as primary public TFBS database for metazoans. PWM scanning + de novo (MEME/HOMER) are the two complementary approaches for non-model organisms.
- **Benchmarking TFBS paper (papers/batch1.5):** PWM scanning (JASPAR + motifmatchr) competitive with deep learning. For non-model organisms without training data, PWM from cross-species databases is the recommended fallback.

### Methylation-expression relationship
- **Stefansson et al. 2024:** Nanopore in humans showed much of meth-expression correlation driven by sequence variants, not methylation per se. Challenges causality even in mammals.
- **Sarda et al. 2012:** Invertebrate gene body methylation marks housekeeping genes. Bimodal: highly methylated vs unmethylated. Methylated genes = broadly expressed.
- **Zemach et al. 2010:** Foundational paper on gene body methylation across eukaryotes.

---

## Critical rules

1. **Strand collapse before any analysis.** Minus-strand pos - 1, sum counts. One CpG per dinucleotide.
2. **NOT a mammal.** No CpG islands, no promoter silencing, no DNMT3.
3. **Filter scaffolds.** `keep_chr <- paste0("chr", 1:31)` — scaffolds are contamination/repeats.
4. **TE is NOT a functional region.** Separate boolean column, not in pie charts or annotation categories.
5. **GenomicRanges for annotation.** NOT ChIPseeker. Priority: Promoter > Exon > Intron > Intergenic.
6. **TFBS motifs only.** JASPAR 2024 CORE metazoan (1,362 PWMs). All are TF binding sites.
7. **Motif output format:** motif name, consensus sequence, chr, start, end, nearest gene, distance to gene.
8. **Let data speak.** Code asks questions, doesn't confirm assumptions.
9. **No HTML in batch scripts.** Data + figures only. Code is supplementary for reviewers.
10. **PDF + PNG for every figure.** `cairo_pdf` + `png(res=300)`.
11. **All numbers from running code.** Never paste from comments or markdown.
12. **Per-read entropy needs BAMs.** Per-site entropy (from beta values) is the fallback. State which you used.

---

## Cluster information

- **Repo:** `/mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave/`
- **Partition:** `defq` | **RAM:** 1 TB | **HOMER:** installed at `cluster/homer/`
- **BSgenome:** `BSgenome.Dlaeve.NCBI.dlgm` (NOT `ASM5140357v1`)
- **SLURM uses `SLURM_SUBMIT_DIR`** (not `dirname $0` which resolves to /var/spool/)
- **CpG reports on cluster are .txt.gz** (compressed) vs .txt locally
- **genome_chr1_31.rds needs to be generated on cluster** (run `generate_genome_cache.slurm` first)

### Cluster jobs available
| Job | Script | RAM | Status |
|-----|--------|-----|--------|
| Genome cache | `generate_genome_cache.R` | 8 GB | Ready to submit |
| DMLtest (full) | `dmltest_full.R` | 32 GB | Ready to submit |
| Expanded motif | `batch1.5_expanded_motif.R` | 64 GB | Failed (needs genome cache first) |

### Run order on cluster
```bash
cd Deroceras-Leave && git pull
cd cluster && dos2unix scripts/*.slurm scripts/*.sh
sbatch scripts/generate_genome_cache.slurm   # FIRST
# wait, then:
sbatch scripts/dmltest_full.slurm
sbatch scripts/batch1.5_expanded_motif.slurm
```
