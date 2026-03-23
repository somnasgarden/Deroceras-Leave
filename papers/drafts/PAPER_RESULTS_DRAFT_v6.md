# Results

Each section corresponds to a batch in `methylation_pipeline/`. All numbers come from strand-collapsed CpG data (25,439,068 sites, coverage >= 5 in all 4 samples). Figures reference the outputs of each batch's code.

---

## 1. Genomic CpG annotation (Batch 1)

The *D. laeve* genome (1.78 Gb, 31 chromosomes, 42.99% GC) contains 45,535,149 CpG dinucleotides (25.5 per kb). The genome-wide CpG observed/expected ratio is 0.553 (computed from actual base composition: A = 28.5%, C = 21.5%, G = 21.5%, T = 28.5%), indicating moderate CpG depletion consistent with historical methylation-driven deamination. This depletion is substantially less severe than in humans (O/E = 0.25) but comparable to other invertebrate mollusks and annelids (Guynes et al., 2024).

CpG dinucleotides are the most underrepresented of all 16 dinucleotides relative to base composition expectation (**Fig. 1A**). CpG O/E varies by genomic region: transposable elements show the highest ratio (0.66), followed by promoters (0.64), intergenic regions (0.56), introns (0.53), and exons (0.50). Among TE classes, SINEs have the highest O/E (0.76) and LTRs the lowest (0.59), suggesting differential historical methylation pressure across repeat families (**Fig. 1B-D**).

**Figures:**
- (A) Dinucleotide frequency: observed vs expected from base composition. CpG highlighted.
- (B) CpG O/E by genomic region (bar chart).
- (C) CpG O/E by TE class (bar chart).
- (D) Per-chromosome CpG O/E density distributions by region and TE class.

---

## 2. Methylation machinery (Batch 2)

We identified 18 methylation-related genes in the *D. laeve* genome: a single DNMT1 homologue (maintenance methyltransferase), TET2 and TET3 (active demethylation), MBD2 and MBD4 (methyl-CpG readers), two UHRF1 copies (DNMT1 recruiters), base excision repair components (APEX1, TDG), GADD45 family members, chromatin remodelers (HELLS), and HDAC2. No DNMT3 homologue was identified, indicating that *D. laeve* lacks canonical *de novo* methylation capacity. This contrasts with the cephalopod *Nautilus pompilius*, which retains both DNMT1 and DNMT3 (Wu et al., 2025).

All 18 genes are expressed in control tail tissue. DNMT1 shows highest expression in ovotestis, consistent with a role in maintaining methylation during germ cell divisions. Expression across six tissues (tail, eye, bodywall, head, juvenile, ovotestis) confirms a functional methylation toolkit (**Fig. 2A-C**).

**Figures:**
- (A) Presence/absence of methylation toolkit genes (genome annotation).
- (B) Expression boxplot by tissue (VST-normalized, controls only).
- (C) Z-scored expression heatmap across tissues.

---

## 3. Promoter CpG architecture (Batch 3)

We classified 25,413 protein-coding gene promoters (TSS +/- 1 kb) using the Weber criteria (Weber et al., 2007). Unlike the bimodal distribution characteristic of mammalian genomes, *D. laeve* promoters display a unimodal, right-skewed CpG O/E distribution. The vast majority of promoters (91.6%) fall into the intermediate-CpG (ICP) class, with only 6.5% high-CpG (HCP) and 1.9% low-CpG (LCP). In mammals, proportions are approximately 72% HCP, 16% ICP, and 12% LCP.

The dominance of ICP promoters means: (1) canonical CpG islands are absent from the *D. laeve* genome; (2) the CpG island model of promoter regulation -- where methylation of HCP promoters stably silences genes -- does not apply; and (3) CpG island "shores" (Irizarry et al., 2009) cannot exist without islands (**Fig. 3A-D**).

**Figures:**
- (A) CpG O/E histogram of all 25,413 promoters (with Weber thresholds marked).
- (B) CpG O/E vs GC% scatter plot (Weber classification overlay).
- (C) Promoter class proportions (pie or bar: HCP / ICP / LCP).
- (D) Faceted CpG O/E histograms by class.

---

## 4. Baseline methylation landscape and gene body methylation vs expression (Batch 4)

WGBS of control and amputated tail tissue yielded methylation data for 25,439,068 strand-collapsed CpG sites (coverage >= 5x in all four samples: C1, C2, A1, A2). The global methylation distribution is bimodal: ~80% of CpGs are unmethylated (beta < 0.2), ~16% are highly methylated (beta > 0.8), and ~4% are intermediate. Mean genome-wide methylation is 18%.

Methylation varies by genomic region in a pattern consistent with other invertebrates: exons have the highest mean methylation (beta = 0.437), followed by introns (0.229), promoters (0.177), and intergenic regions (0.024). This exon-peaked gradient is distinct from the mammalian pattern and consistent with annelid (Guynes et al., 2024) and cephalopod (Wu et al., 2025) methylomes.

Baseline gene body methylation is positively correlated with gene expression (Spearman rho = 0.307, p < 10^-80), confirming the canonical invertebrate pattern. The relationship holds for exons, introns, and gene body as a whole, with exonic methylation showing the strongest association (**Fig. 4A-D**).

**Figures:**
- (A) Metagene methylation profile (upstream -> promoter -> exons -> introns -> downstream), comparing control vs amputated.
- (B) Gene body methylation vs expression scatter (Spearman correlation).
- (C) Expression by methylation quartile (boxplots), separately for promoter, gene body, downstream.
- (D) Global methylation distribution per sample (density plots).

---

## 5. Transposable element methylation and evolutionary age (Batch 5)

Transposable elements constitute 41.2% of the *D. laeve* genome (3,127,881 elements in 1,401,889 merged regions). TE methylation shows an inverse relationship with evolutionary age (Kimura 2-parameter divergence): younger TEs tend to be more methylated. Among TE classes, SINEs show the highest mean methylation and DNA transposons the lowest.

TEs within genes (intronic) are more methylated than intergenic TEs, consistent with passive acquisition of gene body methylation rather than active TE-targeted silencing. This is expected given the absence of DNMT3 in *D. laeve*.

TE methylation does not change between control and amputated tissue (all Cohen's d < 0.01), despite 53.4% of DMPs overlapping TE sequences. TEs appear to serve as genomic substrate for regulatory elements rather than being direct targets of methylation-based silencing (**Fig. 5A-B**).

**Figures:**
- (A) Kimura divergence vs methylation by TE class (scatter/line plot, age in millions of years).
- (B) Intronic vs intergenic TE methylation comparison (ridge plots showing bimodal distribution).

---

## 6. Differential methylation upon tail amputation (Batch 6)

Using DSS with Bayesian smoothing on strand-collapsed data, we identified 17,729 differentially methylated positions (DMPs; FDR < 0.05, |methylation difference| > 10%) and 1,462 differentially methylated regions (DMRs; >= 3 CpGs, >= 50 bp) between control and amputated tail tissue.

DMPs are slightly hypo-dominant: 9,444 (53.3%) hypomethylated vs 8,285 (46.7%) hypermethylated. DMRs: 772 (52.8%) hypo vs 690 (47.2%) hyper. Mean DMR length: 264 bp, mean CpGs per DMR: 7.5. Sample correlations exceed 0.985 across all four samples.

DMPs are enriched in exons (1.95x vs genome background, Fisher's p < 0.001) and promoters (1.44x, p < 0.001), depleted in intergenic regions (0.81x, p < 0.001), and at baseline in introns (1.00x, ns). TE overlap shows no enrichment (53.4% vs 53.0% background, p = 0.38) (**Fig. 6A-H**).

**Figures:**
- (A) PCA of methylation profiles (4 samples).
- (B) Spearman correlation heatmap.
- (C) Per-chromosome mean methylation (control vs amputated).
- (D) DMP volcano plot.
- (E) DMR volcano plot.
- (F) DMP/DMR genomic region annotation (bar chart, functional regions only).
- (G) DMP enrichment vs background (fold change + Fisher's test).
- (H) GO/KEGG enrichment of genes near DMPs/DMRs (if available; note: non-model organism, limited annotation).

---

## 7. Differential methylation is decoupled from differential expression (Batch 7)

To test whether methylation changes predict expression changes, we computed mutual information (MI) between methylation direction (hyper/hypo) and expression direction (up/down) for all genes harboring at least one DMP.

The 2x2 contingency table:

|                | Down-regulated | Up-regulated |
|----------------|---------------|-------------|
| Hypermethylated | 752            | 581          |
| Hypomethylated  | 984            | 795          |

MI = 0.000087 bits (maximum = 1 bit). Permutation test (10,000 iterations): p = 0.55. The global Spearman correlation between methylation change and expression change is rho = -0.021 (p = 0.24, n = 3,112 genes), indistinguishable from zero.

Methylation operates as a binary switch: DE rate does not increase with DMP count (Cochran-Armitage z = -1.21, p = 0.23). The 60 largest DMRs (>= 20 CpGs) overlap zero differentially expressed genes (**Fig. 7A-C**).

**Figures:**
- (A) Delta-methylation vs delta-expression scatter (3,112 genes, with rho and MI annotated).
- (B) DMR-focused: top DMRs mapped to nearest genes with expression change.
- (C) Heatmap of top DMP-bearing genes with expression data.

---

## 8. Co-expression modules reveal module-specific methylation (Batch 8)

WGCNA on multi-tissue RNA-seq identified 12 co-expression modules (62 to 3,546 genes). The proportion of genes with DMPs varies from 10.0% to 22.6% across modules. All modules are hypo-dominant.

Within-module methylation-expression correlations range from rho = -0.346 (magenta) to rho = +0.064 (tan). The global null correlation (rho = -0.021) arises not because methylation has no effect, but because positive and negative module-level effects partially cancel -- a form of Simpson's paradox.

Specific modules are enriched for methylation changes beyond what is expected by chance. These modules are associated with distinct biological functions (to be determined by enrichment analysis) (**Fig. 8A-C**).

**Figures:**
- (A) Module DMP burden (% genes with DMPs per module).
- (B) Within-module methylation-expression correlation (forest plot or bar chart).
- (C) Module enrichment analysis (GO/functional terms for enriched modules).

---

## 9. Transcription factor binding sites and regulatory network integration (Batch 9)

We identified 2,217 transcription factor genes using DeepTFactor (from 73,425 proteins). TFs are enriched for DMPs (Fisher's OR = 1.30, p = 0.0002). In silico TFBS annotation using 1,362 JASPAR 2024 metazoan motifs identified binding sites across promoters, with 19 confirmed TF families (TF protein + binding motif both present in genome) and 19 orphan motif classes (motif present but no predicted TF).

Integration with GENIE3 regulatory network predictions tests whether methylation changes at TF loci or their binding sites propagate through regulatory networks. The Sox19a locus illustrates the intergenic regulatory mechanism: 13 hypermethylated DMPs cluster 4-5 kb upstream of the TSS in an intergenic region, sparing the gene body. Whether these DMP clusters coincide with predicted TFBS will determine if they mark functional regulatory elements (**Fig. 9**).

**Figures:**
- (A) TF family census (DeepTFactor + classification).
- (B) TFBS motif prevalence (top 30 motifs across promoters).
- (C) TF-motif cross-reference (confirmed vs orphan).
- (D) Sox19a locus: DMP positions + TFBS predictions + gene structure.

---

## 10. Methylation entropy during regeneration (Batch 10)

To test whether regeneration induces a global shift in methylation organization, we compute per-CpG Shannon entropy across samples. The question: does amputation drive the methylome toward higher entropy (disorder, dedifferentiation) or lower entropy (directed reprogramming)?

[Results pending -- entropy analysis to be completed.]

**Figures:**
- (A) Global entropy distribution: control vs amputated.
- (B) Per-region entropy changes.
- (C) Entropy vs expression relationship.
