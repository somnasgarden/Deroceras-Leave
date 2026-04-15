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
  batch01/  {code/,data/,figures/}  — Genomic CpG annotation
  batch1.5/ {code/,data/,figures/}  — TFBS motif annotation (JASPAR + HOMER, genome-wide)
  batch02/  {code/,data/,figures/}  — Methylation machinery (18 genes, DNMT1 only)
  batch03/  {code/,data/,figures/}  — Promoter Weber classification (HCP/ICP/LCP)
  batch04/ {code/,data/,figures/}  — Baseline methylation + gene body meth vs expression
  batch05/ {code/,data/,figures/}  — TE methylation + evolutionary age (Kimura)
  batch06/ {code/,data/,figures/}  — Differential methylation (DMP/DMR)
  batch07/ {code/,data/,figures/}  — Methylation-expression decoupling (MI, correlation)
  batch08/ {code/,data/,figures/}  — WGCNA modules + methylation enrichment
  batch09/ {code/,data/,figures/}  — TF + GENIE3 + motif integration
  batch10/ {code/,data/,figures/}  — Entropy analysis (Level 1 marginal + Level 2 per-read NME from BAMs as 10b)
  batch11/ {code/,data/,figures/}  — Motif × per-read NME × expression (closes the loop)
```

### Dependency graph

```
Batch 01 (genome CpGs) → Batch 02 (machinery) → Batch 03 (Weber promoters)
                                                          ↓
Batch 1.5 (TFBS motifs) ──────────────→ Batch 04 (baseline meth+expr) → Batch 05 (TE+age)
                                                          ↓
                                          Batch 06 (DMP/DMR) → Batch 07 (decoupling)
                                                          ↓
                                          Batch 08 (WGCNA) → Batch 09 (TF+GENIE3) → Batch 10 (entropy L1)
                                                                                          ↓
                                                                              Batch 10b (per-read NME, BAMs)
                                                                                          ↓
                                                          Batch 11 (motif × NME × expression — directed-mechanism shortlist)
```

### Rules for batch scripts
1. **No HTML.** Batches produce TSV data + PNG/PDF figures only. Reviewers read the code.
2. **Source `_config.R`.** Every batch starts with `source("methylation_pipeline/_config.R")`.
3. **Save with helpers.** Use `save_fig(p, BATCH_DIR, "name")` and `save_data(df, BATCH_DIR, "name")`.
4. **Let data speak.** Code asks questions and reports answers. Don't hardcode expected results or gene names.
5. **Strand-collapsed data throughout.** The BSseq cache (`bsseq_tutorial.rds`) is already strand-collapsed.
6. **Self-contained.** All objects built from OG raw data by the pipeline. No external pre-computed objects.
7. **Batch 04 builds BSseq** from CpG reports if cache doesn't exist. Batch 06 runs DMLtest if cache missing.
8. **WGCNA power = 12.** Filter outlier samples via hierarchical clustering before module detection.
9. **Motif annotation = genome-wide.** Scan all chr1-31, annotate each hit with nearest gene + region type.
10. **Batch 09 integration chain:** GENIE3 (TF→target) + motif (TFBS near gene) + DMP (methylation change at motif) + DE (expression change at target). All 4 must align for a regulatory link.

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
- **Question:** Does regeneration increase or decrease methylation entropy? Are DMPs directed or stochastic?
- **Input:** BSseq cache (Level 1). BAM files for Level 2 at `/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/08_deduplication/` (confirmed: C1=25GB, C2=20GB, A1=26GB, A2=22GB).
- **Output:** Per-CpG entropy, entropy by region, entropy at DMPs vs matched non-DMPs, entropy-expression variability correlation
- **Figures:** (A) Entropy distribution control vs amputated, (B) Per-region entropy, (C) DMP vs matched non-DMP entropy, (D) Baseline entropy vs beta for DMPs and non-DMPs

#### Two levels of entropy analysis (Fang et al. 2023, Feinberg lab)

Fang et al. showed that **entropy carries comparable or more developmental information than mean methylation**. 22-43% of developmental changes are entropy-dominant and **invisible to DMP/DMR analysis**. The field has been systematically missing half the picture by only looking at mean shifts.

**Level 1: Per-site (binary) entropy** — from CpG_report.txt (available now):
```
H = -(beta * log2(beta) + (1-beta) * log2(1-beta))
```
A proxy for NME. Cannot distinguish two loci both at beta = 0.5 where one has a clean 50/50 split (low true entropy) and the other has random patterns (high true entropy).

**Level 2: Per-read NME (normalized methylation entropy)** — from BAM files (available later):
```
NME = -(1/N) * sum(p_i * log2(p_i)) for all 2^N epiallele patterns in a window
```
The real methylation entropy. Captures co-methylation structure. Tools: DMEAS, informME (CPEL pipeline), or custom from Bismark BAMs. **BAMs exist on cluster — to be added.**

**Level 2 unlocks:**
- Entropy-only changes (mean unchanged, entropy changed) — the "hidden layer"
- These are enriched at enhancers and pioneer TF binding sites (KLF4, SOX, GATA in Fang et al.)
- Could find regeneration-associated changes invisible to our 17,729 DMPs

#### Key question: drift vs reprogramming
**Null hypothesis:** Without DNMT3 (de novo methyltransferase), methylation changes should be stochastic. DMPs should come from already-disordered (high-entropy) baseline sites.

- **DMPs at LOW baseline entropy → directed reprogramming** (reject null). Striking for DNMT3-absent organism.
- **DMPs at HIGH baseline entropy → stochastic drift** (fail to reject).

From Shao et al. 2014: lineage commitment sites show **intermediate methylation + low entropy** = coordinated directed change. If D. laeve DMPs show this same signature, that's evidence for conserved programmed reprogramming.

#### CRITICAL CONTROL: matched comparison
Sites near beta = 0 or 1 have low entropy by definition. Comparing DMP entropy to all CpGs just recapitulates the beta distribution. **Must compare DMPs to matched non-DMPs at similar beta values.** For each DMP, sample non-DMPs with similar baseline beta (e.g., within 0.05) and compare entropy distributions. Also stratify by local CpG density (Fang showed entropy inversely related to CpG density, 2.1-fold enrichment, p < 2.2e-16).

#### Entropy and expression variability
Mean methylation correlates with mean expression. Entropy correlates with **expression variability** (cell-to-cell noise). High NME near TSS = high expression variability across cells. This is a fundamentally different axis. Fang et al.: entropy-associated TF motifs include KLF4, SOX, GATA (pluripotency/reprogramming factors). Mean-associated motifs include NFI family (lineage commitment).

#### Literature
- **No published study has computed methylation entropy in regeneration or invertebrates.** Both would be novel.
- All published entropy papers used BAM files for per-read entropy. Per-site entropy is a valid approximation.
- Check if BAMs exist on cluster at `/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/`

### Batch 11: Motif × Per-read NME × Expression
- **Question:** Do TF binding sites enforce ordered methylation patterns at differentially expressed genes? (Closes the directed-mechanism loop.)
- **Input:** `batch10/data/perread_nme_windows.tsv` (4-CpG NME windows from batch10b) + `batch1.5/data/motif_hits_extended.tsv.gz` (75M genome-wide JASPAR hits, all 4 region classes) + `batch06/data/dmps_annotated.tsv` + `CACHE$transcriptome` (DESeq2 control vs amputated tail).
- **Coverage filter:** Minimum 8 reads/sample/window across all 4 samples.
- **Output:** `window_motif_assignments.tsv`, `nme_in_vs_out_motif_test.tsv`, `delta_nme_by_methyl_sens.tsv`, `read_coverage_by_motif_class.tsv`, `directed_mechanism_loci.tsv`.
- **Figures:** (a) β-matched in-motif vs out-of-motif baseline NME, (b) Δ NME at TFBS by Yin 2017 methyl-sensitivity, (c) directed-loci volcano (log2FC vs −ctrl_nme, colored by Δ NME), (d) reads per window in vs out of motif, (e) reads per TF class top 20, (f) NME vs read depth sanity.

#### The logic chain
A 4-CpG window with low NME (≈0) means every read at that locus carries the same methylation pattern. That uniformity must be **enforced** by something — most plausibly a sequence-specific factor. High NME (≈1) means reads disagree → no enforcement → drift-compatible.

If 4-CpG windows that overlap a JASPAR TFBS have **systematically lower NME than β-matched non-motif windows**, that is evidence the bound TF is the enforcement mechanism (directed). Then if those low-NME-in-motif windows belong to genes that are **differentially expressed during regeneration**, the loop closes:

**TF motif → ordered methylation pattern (low NME) → DMP → DE gene → regeneration phenotype**

Converts the small global Δ NME (~0.002) reported by batch10 into a per-locus story: the change is concentrated at TFBS of methylation-sensitive TFs in DE genes, with N specific cloning candidates.

#### Three tests
- **Test A — ordered enforcement.** β-matched (±0.05) Wilcoxon paired test on baseline NME, in-motif vs out-of-motif. β-matching is mandatory — without it the test rediscovers that promoter CpGs have different entropy than gene-body CpGs.
- **Test B — ordering during regeneration.** For in-motif windows, Δ NME stratified by Yin 2017 methyl-sensitivity class, one-sample Wilcoxon vs 0. Expectation: negative.
- **Test C — directed shortlist.** `in_motif & ctrl_nme < 0.30 & has_dmp`, join `CACHE$transcriptome$res_tail`, rank by `|Δ NME| × −log10(padj)`. Mechanistic priors for cloning targets.

#### Implementation notes
- Yin 2017 sensitivity is **inlined from batch09** (~100 hand-curated TF symbols + class-prior fallback). Batch11 does not depend on batch09 having run. The override list is partial — Yin assayed ~542 TFs. Upgrade target: HOCOMOCO v12 annotation (~1300 TFs, includes methyl-aware PWM variants for ~80) or Yin Table S2 directly.
- 4-CpG windows are 50–500 bp; JASPAR PWMs are 6–20 bp. Windows are intentional **supersets** of the motif core — a TF footprint protects ~30–60 bp of flanking CpGs, where the entropy signal lives.
- Whole gene body, not just promoters: input motif file covers 9M promoter + 20M upstream_distal + 24M gene_body + 22M downstream hits.
- Read-coverage diagnostics (figs 11d/e/f) are mandatory: NME is unstable below ~8 reads/window, reviewers will ask whether low NME at TFBS is a coverage artifact.

#### Cluster execution
Batch10b must finish first (`perread_nme_windows.tsv` exists). Then `sbatch cluster/scripts/batch11_motif_entropy.slurm` (4 CPU / 48 GB / 4h on defq). Pure data joins, no BAM parsing. Wired into `run_pipeline.slurm` as the final batch.

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
- **Fang et al. 2023 (NAR, Feinberg lab):** THE key entropy paper. Entropy carries comparable or MORE developmental information than mean methylation. 22-43% of developmental changes are entropy-dominant (invisible to DMPs). Entropy inversely related to CpG density. Entropy-associated TF motifs: KLF4, SOX, GATA (pluripotency factors). Mean-associated: NFI family (lineage commitment). Entropy predicts expression variability, not expression level. Conserved mouse-human. **Used BAMs + informME/CPEL pipeline.** Framework: MML (mean) + NME (entropy) + UC (uncertainty coefficient) per 250bp window.
- **NO PUBLISHED STUDY on methylation entropy in regeneration or invertebrates.** Both would be novel contributions.

### TFBS motif annotation
- **Hashim et al. 2019 (AJMB-11-130):** Review of 119 motif discovery algorithms. JASPAR highlighted as primary public TFBS database for metazoans. PWM scanning + de novo (MEME/HOMER) are the two complementary approaches for non-model organisms.
- **Benchmarking TFBS paper (papers/batch1.5):** PWM scanning (JASPAR + motifmatchr) competitive with deep learning. For non-model organisms without training data, PWM from cross-species databases is the recommended fallback.

### Methylation-expression relationship
- **Stefansson et al. 2024:** Nanopore in humans showed much of meth-expression correlation driven by sequence variants, not methylation per se. Challenges causality even in mammals.
- **Abbott et al. 2024 (Evolutionary Applications, e13662):** Coral gene body methylation changes do NOT correlate with expression on physiological timescales. Meta-analysis of 8 invertebrate taxa confirms. Supports MI=0 finding as expected, not surprising.
- **Sarda et al. 2012:** Invertebrate gene body methylation marks housekeeping genes. Bimodal: highly methylated vs unmethylated. Methylated genes = broadly expressed.
- **Zemach et al. 2010:** Foundational paper on gene body methylation across eukaryotes.

### DNMT3 loss — mechanistic parallels
- **Langer et al. 2025 (J Exp Zool B):** Colorado potato beetle (*L. decemlineata*) also lost DNMT3 but retains DNMT1. H3K36me3 mirrors CpG methylation patterns, suggesting histone marks guide methylation maintenance without de novo enzyme. Dynamic methylation changes between embryo and adult. Closest mechanistic parallel to D. laeve.
- **Haggerty et al. 2021:** DNMT1 has de novo activity at retrotransposons in mouse ESCs, dependent on UHRF1. D. laeve has UHRF1 — could explain how young TEs get methylated.

### Regeneration epigenetics
- **Axolotl epigenetic clocks 2024 (bioRxiv):** Methylomes stable across lifespan but show "structure-specific rejuvenation events upon regeneration." If D. laeve shows entropy decrease, it parallels axolotl across phyla.
- **Planarian regeneration:** No detectable 5mC in Schmidtea mediterranea. Makes D. laeve one of very few invertebrate regeneration systems where WGBS methylation analysis is meaningful.

### Methylation entropy and aging
- **Wang et al. 2024 (Aging Cell):** Long-lived humans suppress methylation entropy increase. Lower entropy = transcriptional noise suppression. Framework for regeneration-as-rejuvenation.
- **Daphnia magna 2025 (Epigenetics & Chromatin):** No age-associated methylation changes. Counterpoint: not all invertebrates show aging-entropy pattern.

---

## Cross-species validation points

Each batch must explicitly test these against our data. These are not assumptions — they are testable predictions from the literature. Report whether *D. laeve* confirms or contradicts each one.

### From Nautilus pompilius (Wu et al. 2025) — same phylum, closest comparison

| # | Point | Batch | How to test |
|---|-------|-------|-------------|
| 1 | Metagene profile: low at promoter → peak in gene body → drop at terminator | 04 | Compute mean beta in windows around TSS (upstream → exons → introns → downstream) |
| 2 | Moderately expressed genes have HIGHEST gene body methylation (not the most expressed) | 04 | Bin genes by expression decile, plot mean gene body beta per decile |
| 3 | Unexpressed genes still have moderate gene body methylation | 04 | Check beta for genes with FPKM/VST ~ 0 |
| 4 | Promoter methylation weakly negatively correlated with expression (Nautilus rho = -0.05) | 04 | Compute Spearman(promoter_beta, expression) |
| 5 | First exon methylation negatively correlated with expression | 04 | Compute Spearman(first_exon_beta, expression) |
| 6 | Tissue-specific genes (Tau > 0.8) are hypomethylated vs housekeeping | 04 | Compute Tau index from multi-tissue RNA-seq, compare gene body beta |
| 7 | Nautilus has DNMT1+DNMT3, D. laeve has DNMT1 only. Any pattern difference? | 02,04 | Qualitative comparison of metagene profiles and correlations |

### From Annelids (Guynes et al. 2024) — spiralian sister group

| # | Point | Batch | How to test |
|---|-------|-------|-------------|
| 8 | Exon vs intron methylation: which is higher? (varies by species) | 04 | Compare mean beta exon vs intron |
| 9 | Methylation erodes from embryo → adult. Is regenerating tissue different from control? | 06 | Compare global mean beta: control vs amputated (even small difference matters) |
| 10 | GbM-expression correlation strongest in early development. Where does our adult rho fit? | 04 | Report our rho = 0.307 in context of their developmental series |
| 11 | Promoter always hypomethylated regardless of stage | 04 | Confirm promoter mean beta << gene body in both conditions |
| 12 | D. gyrociliatus lost methylation despite DNMT1. D. laeve kept it. What's different? | 02 | Compare DNMT toolkit, note presence of UHRF1, TET, etc. |

### From Weber et al. 2007 / Schubeler 2015 — promoter architecture

| # | Point | Batch | How to test |
|---|-------|-------|-------------|
| 13 | HCP promoters unmethylated regardless of expression | 03,04 | Mean beta of our 6.5% HCP promoters — are they unmethylated? |
| 14 | ICP promoters are targets for de novo methylation in mammals | 03,04 | Mean beta of our 91.6% ICP promoters — methylated or not? |
| 15 | LCP promoters methylated but methylation doesn't silence them | 03,04 | Expression of genes with LCP promoters vs their methylation |
| 16 | H3K4me protects CpG islands from methylation. Without islands, what protects D. laeve promoters? | 03 | Describe the structural absence; note no chromatin data (ATAC-seq needed) |
| 17 | Weak CpG islands are predisposed to de novo methylation (Weber) — do ICP promoters gain methylation in amputation? | 03,06 | Compare promoter methylation between control and amputated for each Weber class |
| 18 | Gene body methylation suppresses intragenic spurious transcription (Schubeler) | 04 | Can we detect cryptic TSS in unmethylated gene bodies? (needs RNA-seq, may not be testable) |
| 19 | Methylation at CG-poor regulatory regions occurs DOWNSTREAM of TF binding (Schubeler 2015) | 09 | If TFBS overlaps DMP, is the TF expressed? If TF binds first then methylation follows, we expect TF expression to precede methylation change |
| 20 | DNMT1 alone can perform de novo methylation at retrotransposons (Haggerty 2021) | 05 | Are young TEs methylated despite no DNMT3? If yes, DNMT1 may be doing de novo |

### From entropy literature — novel territory

| # | Point | Batch | How to test |
|---|-------|-------|-------------|
| 21 | DMPs from high-entropy (disordered) or low-entropy (committed) baseline sites? | 10 | Compare baseline per-site entropy at DMP positions vs non-DMP positions |
| 22 | Global methylation entropy change upon amputation | 10 | Paired Wilcoxon test: control entropy vs amputated entropy (per-site) |
| 23 | Entropy by region (promoter vs gene body vs intergenic) | 10 | Mean entropy per region, both conditions |
| 24 | Intergenic DMP clusters (Sox19a) — high or low entropy at baseline? | 10 | Extract Sox19a upstream region, compute entropy |
| 25 | Aging increases entropy (Jenkinson 2017). Regeneration = rejuvenation? Does entropy decrease? | 10 | If amputated entropy < control → dedifferentiation looks like rejuvenation |

### Summary: what each batch must cross-validate

| Batch | Points to test |
|-------|---------------|
| 02 | 7, 12 |
| 03 | 13, 14, 15, 16, 17 |
| 04 | 1, 2, 3, 4, 5, 6, 8, 10, 11, 18 |
| 05 | 20 |
| 06 | 9, 17 |
| 09 | 19 |
| 10 | 21, 22, 23, 24, 25 |

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

### BAM files (for Level 2 per-read entropy)
BAMs confirmed at `/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/08_deduplication/`:
- C1_paired_bismark_bt2_pe.deduplicated.bam (25 GB)
- C2_paired_bismark_bt2_pe.deduplicated.bam (20 GB)
- A1_paired_bismark_bt2_pe.deduplicated.bam (26 GB)
- A2_paired_bismark_bt2_pe.deduplicated.bam (22 GB)

### Running the pipeline
All code is self-contained in `methylation_pipeline/`. From repo root:
```bash
sbatch methylation_pipeline/run_pipeline.slurm   # runs batch 01 → 1.5 → 02-10
```
The pipeline auto-generates all caches (genome, GFF, TE, BSseq, DMLtest) from raw data on first run. Subsequent runs use caches. Requires 64GB RAM, 8 CPUs.

### Cluster scripts (legacy, now integrated into pipeline)
`cluster/scripts/` contains standalone versions that also source `_config.R`. These were used for initial cache generation but the pipeline is now self-contained.
