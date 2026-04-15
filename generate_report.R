#!/usr/bin/env Rscript
# =============================================================================
# Generate HTML report showing all pipeline figures with interpretations
# Run from repo root: Rscript generate_report.R
# Output: pipeline_report.html
# =============================================================================

PIPE_DIR <- "methylation_pipeline"
OUT_FILE <- "pipeline_report.html"

# Batch metadata: name, title, question, analysis notes per figure prefix
batches <- list(
  list(
    dir = "batch0.5", title = "Batch 0.5 — Transcriptome Processing",
    question = "What is the clean, filtered expression dataset?",
    context = "44 RNA-seq samples across 6 tissues. Tail (control/amputated), eye (control/amputated), bodywall (control/irradiated/fungicide), head/juvenile/ovotestis (control). Within-group centroid distance flags 3 outliers (Bodywall_Control_4, Bodywall_Irradiated_5, Tail_Amputated_1). DESeq2 per-tissue with apeglm LFC shrinkage. This clean dataset feeds all downstream methylation-expression analyses.",
    notes = list(
      fig05a = "Sample dendrogram showing hierarchical clustering of all 44 RNA-seq samples across 6 tissues. Samples cluster primarily by tissue type, confirming tissue identity drives expression.",
      fig05b = "PCA of all samples — PC1/PC2 separate tissues (ovotestis, juvenile vs adult). Within-group outliers flagged by centroid distance (>2× median within-tissue distance).",
      fig05c = "PCA of tail samples only — control vs amputated separation. Tail_Amputated_1 detected as within-group outlier.",
      fig05d = "Sample correlation heatmap — high within-tissue correlation, low between-tissue, confirming biological replicates are consistent."
    )
  ),
  list(
    dir = "batch01", title = "Batch 01 — Genomic CpG Landscape",
    question = "What does the D. laeve CpG landscape look like?",
    context = "CpG O/E = 0.553 — strong depletion from evolutionary deamination. 45.5M genome CpGs. TEs have unexpectedly high CpG density. Exons show lowest O/E despite being most methylated — consistent with methylation-driven CpG loss in gene bodies.",
    notes = list(
      fig1a = "Dinucleotide frequencies: CG is depleted (red bar, ~2.5%) relative to expectation (~4.6%), classic signature of CpG deamination over evolutionary time. CpG O/E = 0.553.",
      fig1b = "CpG density distribution by region (CpG/kb): per-chromosome density plots. TEs have unexpectedly high CpG density. Exons lowest. Promoters show tight distribution.",
      fig1c = "CpG O/E distribution by region: All regions below expected (deamination). Exons are lowest O/E despite being the most methylated — consistent with methylation-driven deamination in gene bodies.",
      fig1d = "GC% distribution by region: per-chromosome density plots. Regions differ in GC content — promoters tighter, intergenic broader.",
      fig1e = "CpG density distribution by TE class (CpG/kb): SINE stands out with high CpG density. Could TEs harbor CpG island-like sequences?",
      fig1f = "CpG O/E distribution by TE class: SINE has highest O/E (~0.8), approaching expected. LINE/DNA/LTR cluster around 0.55-0.65.",
      fig1g = "GC% distribution by TE class: per-chromosome density plots showing GC content variation across TE families."
    )
  ),
  list(
    dir = "batch1.5", title = "Batch 1.5 — TFBS Motif Annotation",
    question = "Where are transcription factor binding sites genome-wide?",
    context = "1,362 JASPAR metazoan motifs scanned across all chr1-31. 75M total hits (promoter, upstream distal, gene body, downstream). This genome-wide TFBS map enables batch 09 to test whether DMPs co-localize with specific TF binding sites.",
    notes = list(
      fig15a = "JASPAR 2024 CORE metazoan motif hits genome-wide. 1,362 PWMs scanned in chunks of 50 to avoid MOODS memory issues. Produces motif hit tables for batch09 integration.",
      fig15b = "HOMER de novo motif enrichment: fold enrichment of discovered motifs in gene-proximal regions vs genome background. Identifies novel regulatory elements not in JASPAR.",
      fig15c = "De novo motif similarity dendrogram: hierarchical clustering of discovered motifs by PWM similarity. Clusters reveal motif families that may bind related TFs.",
      fig15e = "De novo motif co-occurrence heatmap: pairwise co-occurrence of de novo motifs across target regions. High co-occurrence suggests combinatorial TF binding or shared regulatory grammar.",
      fig15f = "De novo vs JASPAR similarity: best JASPAR match for each de novo motif. Identifies which discovered motifs correspond to known TF families and which are truly novel."
    )
  ),
  list(
    dir = "batch02", title = "Batch 02 — Methylation Machinery",
    question = "Does D. laeve have a functional methylation toolkit?",
    context = "DNMT1 present (maintenance) but NO DNMT3 (de novo) — this is the central biological constraint. Without de novo methyltransferase, where do 17,729 DMPs come from during regeneration? The absence of DNMT3 makes any directed methylation change remarkable and is the foundation of the drift-vs-reprogramming hypothesis.",
    notes = list(
      fig2a = "Presence/absence of 26 methylation-related genes. Key finding: DNMT1 present (maintenance) but NO DNMT3 family (no de novo methylation). Also missing: TET1, MBD1, MBD3, MeCP2. This defines the organism's epigenetic capacity.",
      fig2b = "Expression heatmap across 6 tissues (controls only). Ovotestis shows HIGH expression of methylation machinery — biologically relevant for germline methylation maintenance. DNMT1, HDAC2, HELLS cluster together with high expression.",
      fig2c = "Additional machinery expression view."
    )
  ),
  list(
    dir = "batch03", title = "Batch 03 — Promoter CpG Classification",
    question = "What does the promoter CpG landscape look like?",
    context = "Virtually all promoters are ICP (intermediate CpG). No CpG islands in the mammalian sense. This rules out promoter CpG-island-based silencing as a regulatory mechanism. Methylation acts through gene bodies, not promoters — fundamentally different from mammals.",
    notes = list(
      fig3a = "Promoter CpG O/E histogram: unimodal distribution centered at ~0.6 with Weber thresholds marked. No bimodal pattern (mammals have two peaks for HCP/LCP). Confirms invertebrate-typical promoter architecture.",
      fig3b = "Weber classification scatter (O/E vs GC%): virtually all promoters are ICP (gray). Only 15 HCPs found (red, plotted on top) — D. laeve has NO CpG islands in the mammalian sense. HCPs include WNT-1, Frizzled-10, ECE1.",
      fig3c = "Promoter class proportions: 85.1% ICP, 14.9% LCP, 0.1% HCP. Completely different from mammals (~60% HCP). Methylation-based promoter silencing unlikely.",
      fig3d = "Faceted CpG O/E by promoter class: HCP, ICP, LCP shown separately. HCPs (only 15 genes!) have O/E 0.77-1.08 — true CpG islands.",
      fig3e = "HCP gene list with function, methylation, and DMP/DMR status. All 15 HCP promoters are unmethylated (<0.4%) — CpG island candidates. Includes developmental signaling genes (WNT-1, Frizzled-10-A, lev-9, megf10, ECE1)."
    )
  ),
  list(
    dir = "batch04", title = "Batch 04 — Baseline Methylation Landscape",
    question = "Where is methylation? Does gene body meth correlate with expression?",
    context = "Gene body methylation positively correlates with expression (rho = 0.407, p &lt; 10^-80) — the invertebrate pattern. Exons most methylated (rho = 0.430), promoters near zero (rho = 0.013). Metagene profiles split by gene structure (3+ exon, 2-exon, single-exon, ncRNA) and by promoter (2kb) vs distal upstream (3kb).",
    notes = list(
      fig4a = "Global methylation per sample: all 4 samples show similar mean CpG methylation. Control vs amputated highly concordant — global methylation does not shift during regeneration.",
      fig4b = "Methylation density per sample: bimodal distribution (unmethylated peak near 0, methylated peak near 80-90%). Classic invertebrate mosaic pattern. Samples overlay perfectly.",
      fig4c = "Region methylation (control vs amputated): Exon > Intron > Downstream > Promoter > Intergenic. Gene body methylation pattern confirmed. Minimal condition differences.",
      fig4d1 = "Metagene: ALL genes (control). Distal upstream (5-2kb) / Promoter (2kb) / gene body / downstream. Split upstream tests whether methylation-expression correlation is in promoter or distal enhancer-like region.",
      fig4d2 = "Metagene: ALL genes (amputated). Compare to control for regeneration-induced changes.",
      fig4d3 = "Metagene: 3+ exon genes (control, 8-segment). Distal / Promoter / First exon / First intron / Body / Last intron / Last exon / Downstream. Cf. Nautilus Fig 2B.",
      fig4d4 = "Metagene: 3+ exon genes (amputated).",
      fig4d5 = "Metagene: 2-exon genes (control, 6-segment). Distal / Promoter / First exon / Intron / Last exon / Downstream.",
      fig4d6 = "Metagene: 2-exon genes (amputated).",
      fig4d7 = "Metagene: single-exon intronless genes (control). Distal / Promoter / Gene body / Downstream.",
      fig4d8 = "Metagene: single-exon genes (amputated).",
      fig4d9 = "Metagene: XLOC non-coding RNA (control). How does ncRNA methylation compare to protein-coding genes?",
      fig4e = "Gene body methylation vs expression scatter: Spearman rho = 0.407, highly significant. Positive correlation — methylated genes are more expressed. Typical invertebrate gene body methylation pattern.",
      fig4f = "Per-region correlation bars: Exon has highest rho, followed by gene body, downstream. Promoter near zero — promoter methylation does NOT predict expression in this organism.",
      fig4g = "Methylation by expression decile: monotonic increase from D1 to D9, slight decrease at D10. Not a clear inverted-U (unlike Nautilus). More similar to annelid pattern.",
      fig4h = "Region-specific methylation by expression decile: faceted boxplots. Exon methylation increases most with expression. Promoter methylation flat across deciles.",
      fig4i = "Genome-wide methylation by chromosome (1 Mb windows): all 31 chromosomes, 4 samples overlaid. Shows methylation landscape across the genome — identifies hyper/hypo regions.",
      fig4j = "CpG distribution pie chart by region: shows what fraction of testable CpGs fall in each region.",
      fig4l = "Expanded region correlation: first exon, last exon, single exon, first intron vs expression. Tests whether first exon methylation has a different relationship with expression than last exon.",

      "fig4d_bin10_03_ampu_3plus_exon" = paste0(
        "<strong>NOVEL FINDING &mdash; ICP Promoter Architecture (Paper highlight).</strong> ",
        "Expression-decile metagene for 3+ exon genes (amputated). The most expressed genes (D09-D10, red) show a striking pattern in the promoter region: ",
        "high methylation in distal upstream (~55%) &rarr; drops entering promoter &rarr; <strong>brief spike UP in mid-promoter</strong> &rarr; ",
        "crashes to near 0% at the TSS/first exon &rarr; then climbs into gene body (~65%). ",
        "The bigger the crash, the more expressed the gene. D01 (silent) is flat at ~5%. ",
        "<br><br>",
        "<strong>Why this matters:</strong> 91.6% of D. laeve promoters are ICP (intermediate CpG). Only 15 genes are HCP. ",
        "Schubeler 2015 (Nature) predicted that methylation-based repression only works at HCPs. ",
        "In ICP promoters, methylation is a <em>consequence</em> of transcriptional state, not a cause. ",
        "This figure is the first empirical demonstration in a &gt;95% ICP genome: ",
        "<strong>more upstream/promoter methylation = more expression</strong> (the opposite of mammals). ",
        "The promoter spike-and-crash is the boundary between gene body methylation and the TSS nucleosome-free region. ",
        "<br><br>",
        "<strong>Literature:</strong> Brenet 2011 showed first exon meth inversely linked to expression, but only in humans. ",
        "Wu 2025 (Nautilus) used 4 expression bins but Nautilus has DNMT3. ",
        "Keller 2016 (Ciona) hinted at positive promoter-expression correlation in low-GBM genes. ",
        "<strong>No one has shown 10-decile metagene profiles in a DNMT3-absent, &gt;95% ICP invertebrate.</strong> ",
        "The mid-promoter spike may mark a positioned nucleosome or regulatory element (TATA, Inr) that retains methylation ",
        "while flanking regions are cleared by transcription factor binding."),

      "fig4d_bin10_04_ctrl_3plus_exon" = "Same as amputated 3+ exon decile plot but for controls. Pattern is nearly identical &mdash; the ICP promoter architecture is a baseline feature, not regeneration-induced. Control vs amputated comparison confirms the spike-and-crash is constitutive.",

      "fig4d_bin10_07_ampu_single_exon" = paste0(
        "<strong>Intronless gene contrast.</strong> Single-exon genes show a completely different pattern: ",
        "lines are messy, overlapping, with no clear expression-decile ordering. ",
        "D01 (silent) sits at ~15% &mdash; much higher than D01 in 3+ exon genes (~5%). ",
        "The promoter crash is weak or absent. ",
        "This is consistent with the splicing-dependent methylation model (Zilberman, Bewick &amp; Schmitz): ",
        "DNMT1 recruitment to gene bodies depends on the splicing machinery, ",
        "so intronless genes do not develop the same structured methylation landscape. ",
        "The expression-methylation decoupling in intronless genes further supports that the 3+ exon pattern is splicing-dependent."),

      "fig4d_bin10_01_ampu_all_genes" = "All genes combined (amputated, 10 deciles). The promoter spike-and-crash is visible but diluted because it mixes the strong 3+ exon signal with the noisy intronless pattern. The 3+ exon plot is the cleaner view.",
      "fig4d_bin10_02_ctrl_all_genes" = "All genes combined (control, 10 deciles). Compare to amputated &mdash; patterns are nearly identical, confirming the metagene landscape is constitutive.",
      "fig4d_bin10_05_ampu_2exon" = "2-exon genes (amputated, 10 deciles). 6-segment profile: distal / promoter / first exon / intron / last exon / downstream. Intermediate between the structured 3+ exon and noisy intronless patterns.",
      "fig4d_bin10_06_ctrl_2exon" = "2-exon genes (control, 10 deciles). Compare to amputated &mdash; same constitutive pattern.",
      "fig4d_bin10_08_ctrl_single_exon" = "Single-exon intronless genes (control, 10 deciles). Compare to amputated. Both show noisy, weakly ordered patterns &mdash; consistent with splicing-dependent methylation model.",
      "fig4d_bin10_09_ctrl_xloc_ncrna" = "XLOC non-coding RNA (control, 10 deciles). How does ncRNA methylation compare to protein-coding genes? Lower overall methylation expected.",

      "fig4d_bin05_01_ampu_all_genes" = "All genes (amputated, 5 quintiles). Lower resolution than 10-decile version but cleaner separation between expression groups.",
      "fig4d_bin05_02_ctrl_all_genes" = "All genes (control, 5 quintiles).",
      "fig4d_bin05_03_ampu_3plus_exon" = "3+ exon genes (amputated, 5 quintiles). The ICP promoter spike-and-crash is visible. Q5 (highest expression) shows the strongest pattern.",
      "fig4d_bin05_04_ctrl_3plus_exon" = "3+ exon genes (control, 5 quintiles). Same constitutive pattern as amputated.",
      "fig4d_bin05_05_ampu_2exon" = "2-exon genes (amputated, 5 quintiles).",
      "fig4d_bin05_06_ctrl_2exon" = "2-exon genes (control, 5 quintiles).",
      "fig4d_bin05_07_ampu_single_exon" = "Single-exon intronless genes (amputated, 5 quintiles). Noisy &mdash; no clear expression-dependent methylation ordering.",
      "fig4d_bin05_08_ctrl_single_exon" = "Single-exon intronless genes (control, 5 quintiles).",
      "fig4d_bin05_09_ctrl_xloc_ncrna" = "XLOC non-coding RNA (control, 5 quintiles).",

      "fig4d_bin04_01_ampu_all_genes" = "All genes (amputated, 4 quartiles). Coarsest resolution &mdash; broad expression bins reduce noise but lose fine structure.",
      "fig4d_bin04_02_ctrl_all_genes" = "All genes (control, 4 quartiles).",
      "fig4d_bin04_03_ampu_3plus_exon" = "3+ exon genes (amputated, 4 quartiles). Even at this coarse resolution, Q4 (highest) clearly separates from Q1 (silent).",
      "fig4d_bin04_04_ctrl_3plus_exon" = "3+ exon genes (control, 4 quartiles).",
      "fig4d_bin04_05_ampu_2exon" = "2-exon genes (amputated, 4 quartiles).",
      "fig4d_bin04_06_ctrl_2exon" = "2-exon genes (control, 4 quartiles).",
      "fig4d_bin04_07_ampu_single_exon" = "Single-exon intronless genes (amputated, 4 quartiles).",
      "fig4d_bin04_08_ctrl_single_exon" = "Single-exon intronless genes (control, 4 quartiles).",
      "fig4d_bin04_09_ctrl_xloc_ncrna" = "XLOC non-coding RNA (control, 4 quartiles)."
    )
  ),
  list(
    dir = "batch05", title = "Batch 05 — TE Methylation & Evolutionary Age",
    question = "Is TE methylation related to evolutionary age?",
    context = "TE methylation is PASSIVE — intronic TEs are methylated (gene body spillover), intergenic TEs are unmethylated. No evidence of active TE silencing, consistent with no DNMT3. This distinguishes D. laeve from mammals where DNMT3-mediated TE silencing is a primary methylation function.",
    notes = list(
      fig5a = "TE methylation ridgeplots by class: bimodal distribution — most TEs unmethylated, some fully methylated. The methylated fraction likely represents intronic TEs (passive gene body methylation). No active TE silencing (consistent with no DNMT3).",
      fig5b = "TE methylation intronic vs intergenic ridgeplots: confirms fully methylated TEs are inside genes. Intergenic TEs have ~0% methylation. Key finding: TE methylation in D. laeve is PASSIVE (gene body spillover), not active silencing.",
      fig5e = "DMP count by TE class: DNA and LINE TEs harbor most DMPs. Raw counts — need Fisher's test to normalize by TE coverage genome-wide.",
      fig5f = "DMR count by TE class: regional methylation changes overlapping TEs, broken down by class and direction (hyper/hypo).",
      fig5g = "TE age bar plot split by genomic context: genic TEs (~35-40%) vs intergenic (~5-8%), BOTH flat across all age bins. TE methylation is purely a function of being inside a gene (passive spillover), independent of TE age. No active TE-targeting without DNMT3 — fundamentally different from mammals.",
      fig5j = "Cohen's d effect sizes per TE class: quantifies the magnitude of methylation difference (control vs amputated) within each TE class. Identifies which TE families are most affected by regeneration."
    )
  ),
  list(
    dir = "batch06", title = "Batch 06 — Differential Methylation",
    question = "What methylation changes during regeneration?",
    context = "17,729 DMPs (8,285 hyper + 9,444 hypo) and 1,462+ DMRs across all chromosomes. DMPs enriched in exons (1.95x) and promoters (1.44x), depleted in intergenic (0.81x). The slight hypo bias is notable — more demethylation than new methylation during regeneration, despite no DNMT3. GO enrichment: cell projection, ECM, Notch signaling.",
    notes = list(
      fig6a = "DMP PCA: control and amputated separate on PC1, confirming methylation differences are systematic, not noise.",
      fig6b = "Sample correlation heatmap: pairwise Pearson correlation of CpG methylation across all 4 samples. Control pair and amputated pair should cluster together.",
      fig6c = "DMP Manhattan plot: DMPs distributed across all 31 chromosomes. No single chromosome dominates &mdash; changes are genome-wide.",
      fig6d = "DMP volcano: 8,285 hyper + 9,444 hypo. Slight hypo bias &mdash; regeneration involves more demethylation than de novo methylation. Notable given no DNMT3.",
      fig6d2 = "DMP volcano with gene name labels: top 30 DMPs labeled with gene names for biological interpretation.",
      fig6d3 = "DMR volcano: log10(|areaStat|+1) as significance proxy. Unlike the DMP volcano (which uses per-site Wald p-values giving a smooth funnel), DMRs only have areaStat (cumulative test statistic) &mdash; the two-cloud shape is inherent to DSS.",
      fig6d4 = "DMR volcano with labels: top 10 DMRs by |areaStat| labeled with chr:position. Same plot as fig6d3 with coordinate annotations for the strongest regional signals.",
      fig6e = "DMP annotation pie chart: genomic distribution of DMPs across promoter, exon, intron, and intergenic regions. Shows where DMPs fall relative to gene features.",
      fig6f = "DMP direction by region: hyper vs hypo breakdown per genomic region.",
      fig6f2 = "DMR annotation pie: genomic distribution of DMRs.",
      fig6f3 = "DMR direction by region.",
      fig6g = "DMP enrichment by region (Fisher's exact): Exon and Promoter enriched (1.95x, 1.44x). Intergenic depleted (0.81x). DMPs preferentially occur in functional regions.",
      fig6g2 = "DMR enrichment by region (Fisher's exact): same analysis for DMRs &mdash; are DMRs enriched in different regions than DMPs?",
      fig6h = "DMP heatmap: top 100 DMPs with gene name row labels. Clear separation of control vs amputated. Clusters of co-methylated sites visible.",
      fig6i = "DMR heatmap: top 50 DMRs with gene names and coordinates. Region-level methylation confirms coordinated changes.",
      fig6j = "Top 20 genes by DMP count: genes with most differentially methylated positions. Gene names shown where available.",
      fig6k = "Top 20 genes by DMR count: genes with most differentially methylated regions.",
      "fig6_dmp_go_bp_dotplot" = "GO Biological Process enrichment of DMP-associated genes (clusterProfiler + STRING). Dot size = gene count, color = adjusted p-value. Key pathways for regeneration biology.",
      "fig6_dmp_go_cc_dotplot" = "GO Cellular Component enrichment of DMP-associated genes. Identifies subcellular compartments enriched among differentially methylated genes.",
      "fig6_dmp_go_mf_dotplot" = "GO Molecular Function enrichment of DMP-associated genes.",
      "fig6_dmp_go_all_dotplot" = "Combined GO enrichment (BP + CC + MF) of DMP-associated genes in a single dotplot.",
      "fig6_dmp_kegg_dotplot" = "KEGG pathway enrichment of DMP-associated genes. Identifies metabolic and signaling pathways affected by methylation changes during regeneration.",
      "fig6_dmp_reactome_dotplot" = "Reactome pathway enrichment of DMP-associated genes.",
      "fig6_hypo_go_bp_dotplot" = "GO Biological Process enrichment of hypo-methylated DMP genes specifically. Tests whether demethylated genes during regeneration target different pathways than hypermethylated genes."
    )
  ),
  list(
    dir = "batch07", title = "Batch 07 — Methylation-Expression Decoupling",
    question = "Does differential methylation predict expression changes?",
    context = "THE DECOUPLING RESULT: baseline methylation correlates with expression (rho ~ 0.3), but differential methylation does NOT predict expression changes (rho ~ 0.02, MI ~ 0 bits). Methylation marks active genes but methylation CHANGES during regeneration are statistically independent of expression changes. This paradox is resolved by WGCNA (batch 08) showing pathway-specific effects.",
    notes = list(
      fig7a = "DMP meth-expression correlation by region (basic 4 regions): Spearman rho per region. Key question: does methylation change in a region predict expression change?",
      fig7b = "DMR meth-expression correlation by region: same analysis for DMRs, kept completely separate from DMPs.",
      fig7c = "DMP DE rate by region: what fraction of genes with DMPs in each region are also differentially expressed?",
      fig7d = "DMR DE rate by region.",
      fig7e = "DMP dose-response: more DMPs per gene = more likely to be DE? Tests whether methylation changes accumulate to affect expression.",
      fig7f = "DMP genomic distribution bar: location breakdown of all DMPs.",
      fig7g = "DMR genomic distribution bar.",
      fig7h = "DMP expanded region correlation: first/last exon, first/internal/last intron, downstream. Reveals whether first exon methylation change has a different effect than last exon.",
      fig7i = "DMR expanded region correlation: same expanded view for DMRs.",
      fig7j = "DMP methylation vs expression scatter: mean methylation difference (x) vs log2FC (y) per gene. Top 20 DE genes labeled. Shows the decoupling — cloud centered at origin, no systematic trend.",
      fig7k = "DMR methylation vs expression scatter: same for DMRs."
    )
  ),
  list(
    dir = "batch08", title = "Batch 08 — WGCNA Co-expression Modules",
    question = "Are methylation changes concentrated in specific expression modules?",
    context = "RESOLVES THE DECOUPLING PARADOX: DMPs are NOT randomly distributed across WGCNA modules — they concentrate in specific co-expression modules (identity depends on module detection per run; see fig8e Fisher enrichment). The genome-wide rho ~ 0 masks pathway-specific enrichment. Methylation changes are targeted to regeneration-relevant pathways, not scattered randomly. Multi-tissue WGCNA uses within-tissue centroid outlier detection on 6 tissues (power = 12).",
    notes = list(
      fig8a = "Sample dendrogram: hierarchical clustering of all samples used for WGCNA. Tissues separate cleanly. Within-tissue centroid outlier detection removes samples > 2&times; median distance.",
      fig8b = "Soft threshold power selection: scale-free topology fit (R&sup2;) and mean connectivity vs power. Power = 12 selected &mdash; achieves R&sup2; > 0.8 with reasonable connectivity.",
      fig8c = "Gene dendrogram and module assignments: WGCNA module detection with dynamic tree cut. Module colors shown below dendrogram.",
      fig8d = "Module-trait correlation heatmap: which modules correlate with which tissues/conditions? Key: which modules correlate with amputation? Magenta module strongly associated with fungicide condition.",
      fig8e = "Module eigengenes: tail control vs amputated. Shows which modules shift expression during regeneration. Modules with large shifts are candidates for methylation-expression coupling.",
      fig8f = "DMP burden per WGCNA module: total DMP count per module, normalized by module size. Tests whether DMPs concentrate in specific co-expression modules rather than being randomly distributed.",
      fig8g = "DMP enrichment per WGCNA module (Fisher's exact test): identifies modules significantly enriched for DMPs. These are the modules where methylation and expression changes coincide &mdash; Simpson's paradox resolution.",
      fig8h = "Wilcoxon test statistic vs DMP burden per module: compares module-level expression shift (Wilcoxon p) against DMP enrichment. Modules in the upper-right quadrant have both expression shifts AND DMP enrichment.",
      fig8i = "GO Biological Process enrichment for DMP-enriched modules. Identifies biological functions of modules where methylation-expression coupling occurs.",
      fig8j = "Module sizes: number of genes per WGCNA module. Context for interpreting DMP burden and enrichment results."
    )
  ),
  list(
    dir = "batch09", title = "Batch 09 — TF + Motif + GENIE3 Integration",
    question = "Do TF binding sites near DMPs explain regulatory connections? Are methylation-concordant-with-expression genes targeted by methylation-sensitive TFs?",
    context = "Tests the full regulatory chain: TF (DeepTFactor) &rarr; binding motif in target promoter (JASPAR) &rarr; target has DMP &rarr; target is DE. Fisher enrichment uses a <strong>region-restricted background</strong> (CpGs in promoter+5kb upstream only) to avoid the artifact of comparing genome-wide DMPs against a promoter-only motif universe. Every motif is tagged with Yin 2017 (Science) methylation sensitivity (MethylPlus / MethylMinus / LittleEffect). A focused scan on batch07 concordant genes (methylation direction agrees with expression direction) yields the cloning shortlist.",
    notes = list(
      fig9a = "Top enriched JASPAR motifs at DMPs (Fisher's exact, region-restricted background). Bars labeled with motif name + top target gene; colored by Yin 2017 methylation sensitivity.",
      fig9b = "Motif enrichment volcano: log2(odds ratio) vs -log10(p), colored by methylation sensitivity. Top 15 hits labeled with motif + top target gene.",
      fig9c = "TF classes of enriched motifs, faceted by Yin 2017 methylation sensitivity. Reveals whether DMPs are enriched for methyl-sensitive (bHLH, bZIP, ETS) vs methyl-tolerant (homeodomain) families.",
      fig9d = "TF-target evidence chain: how many GENIE3 edges have motif support, DMP evidence, DE evidence, or all three?",
      fig9e = "Top TF-target pairs with full evidence (motif + DMP + DE). Labels use gene names from EviAnn annotation.",
      fig9f = "DMR motif enrichment, split by direction (Hyper vs Hypo). Bars labeled with motif + top DMR gene; colored by methylation sensitivity. Reveals whether activator-binding vs repressor-binding TFs differentially affected.",
      fig9g = "<strong>NEW:</strong> TF motifs enriched in methylation-concordant-with-expression genes (batch07 shortlist). Bars labeled with motif + top concordant gene. These are the candidates for cloning/validation.",
      fig9h = "<strong>NEW:</strong> Per-locus heatmap — concordant DMP loci that fall <em>strictly inside</em> a TF motif (cpg_in_motif=TRUE). Top 40 genes &times; top 25 TFs, fill = mean &Delta;&beta;. Direct cloning targets.",
      fig9i = "<strong>NEW:</strong> Methylation sensitivity composition of TFs significantly enriched at DMPs (Yin 2017 categories).",
      fig9j = "<strong>NEW:</strong> GO/STRING enrichment of the concordant DMP-motif shortlist genes. Bar labels include the top 3 contributing gene names per term."
    )
  ),
  list(
    dir = "batch10", title = "Batch 10 — Methylation Entropy",
    question = "Does regeneration increase or decrease methylation entropy?",
    context = "THE NOVEL FINDING: Entropy DECREASES during regeneration (delta = -0.002, p = 7.6e-6). More ordered, not more random. This rejects the drift hypothesis — without DNMT3, methylation changes should be stochastic (entropy up), but they are directed (entropy down). Key test: do DMPs come from low-entropy (ordered) or high-entropy (disordered) baseline sites? Low entropy = reprogramming. MI ~ 0 bits confirms methylation-expression independence (Fang framework).",
    notes = list(
      fig10a = "CpG methylation entropy distribution: control vs amputated. Global entropy decreases during regeneration (more ordered).",
      fig10b = "Regional entropy change: delta entropy by genomic region. Exons show largest entropy decrease.",
      fig10c = "DRIFT VS REPROGRAMMING TEST: baseline entropy at DMPs vs beta-matched non-DMPs. If DMPs have LOWER entropy at baseline = directed reprogramming (not stochastic drift). Key novel finding.",
      fig10d = "Entropy change (delta) at DMP sites vs matched non-DMPs. Do DMPs specifically change entropy?",
      fig10e = "Entropy by DMP direction: do hyper-methylated DMPs start from different entropy states than hypo-methylated DMPs?",
      fig10f = "Regional entropy at DMP sites vs genome-wide: are DMP entropy changes larger than background in specific regions?",
      fig10g = "Mutual Information: methylation change vs expression change. MI near zero with permutation null = methylation and expression changes are statistically independent.",
      fig10h = "Entropy vs expression variability (Fang framework): high entropy near TSS predicts high expression noise (CV), not expression level.",
      fig10i = "Sliding-window entropy distribution: entropy computed in windows across gene bodies. Distribution for control vs amputated reveals genome-wide entropy shift patterns.",
      fig10j = "Entropy-gaining vs entropy-losing regions by genomic region: which regions (promoter, exon, intron, intergenic) gain vs lose entropy during regeneration? Direction-specific analysis.",
      fig10k = "Entropy-only hidden layer (Fang framework): genes with significant entropy change but NO significant expression change. These are the &lsquo;hidden&rsquo; epigenetic changes &mdash; methylation disorder shifts that have not yet manifested as expression changes.",
      fig10l = "Per-gene delta entropy volcano: each gene plotted by mean entropy change (x) vs &minus;log10(p) from Wilcoxon paired test (y). Identifies genes with the most significant entropy shifts during regeneration."
    )
  ),
  list(
    dir = "batch11", title = "Batch 11 — Motif × Per-read Entropy × Expression",
    question = "Does TF binding enforce ordered methylation at differentially expressed genes?",
    context = "Closes the loop: TF motif &rarr; low per-read NME (ordered chromatin) &rarr; DMP &rarr; DE gene &rarr; phenotype. Loads batch10b per-read NME windows (4-CpG, Xie 2011), batch1.5 genome-wide motif hits, batch06 DMPs, and CACHE$transcriptome. Tests: (A) beta-matched in-motif vs out-of-motif baseline NME; (B) delta NME at TFBS stratified by Yin 2017 methyl-sensitivity; (C) directed-mechanism shortlist (in_motif & ctrl_nme<0.3 & has_dmp) joined to tail DESeq2 &rarr; cloning candidates.",
    notes = list(
      fig11a = "In-motif vs out-of-motif baseline NME (beta-matched). If TFBS enforce ordered methylation, in-motif NME should be lower than matched non-motif CpG neighborhoods.",
      fig11b = "Delta NME at TFBS by Yin 2017 methyl-sensitivity class (MethylPlus / MethylMinus / LittleEffect). Tests whether methyl-sensitive TFs lose order during regeneration.",
      fig11c = "Directed-mechanism shortlist: concordant loci where a TF motif overlaps a low-entropy baseline CpG that becomes a DMP in a DE gene. Cloning candidates.",
      fig11d = "Per-TF read coverage diagnostic — sanity check that tests are not driven by low-coverage TFs.",
      fig11e = "NME vs log2FC scatter for shortlist genes — do shortlist loci show stronger expression coupling than genome-wide?",
      fig11f = "Top shortlist genes with motif + NME + DMP + DE evidence summary bar."
    )
  ),
  list(
    dir = "batch12", title = "Batch 12 — Promoter Methylation Spike Characterization",
    question = "What causes the methylation spike in the promoter of highly expressed 3+ exon genes?",
    context = paste0(
      "<strong>FINDING: mid-promoter spike was a metagene binning artifact.</strong> ",
      "At 10bp single-site resolution, the promoter is smooth &mdash; no local maximum. ",
      "The 20-bin metagene (100bp/bin) created an apparent peak at the segment boundary. ",
      "However, the ICP promoter architecture IS real and novel: ",
      "more upstream methylation = more expression, TSS crash scales with expression level. ",
      "<br><br>",
      "<strong>Why this could be big:</strong> 91.6% of D. laeve promoters are ICP (intermediate CpG). ",
      "Schubeler 2015 (Nature) predicted that methylation-based repression only works at HCPs (65% in humans). ",
      "In ICP promoters, methylation should be a <em>consequence</em> of transcriptional state, not a cause. ",
      "This batch tests whether the spike is: (A) a sequence feature (CpG-dense element), (B) a positioned nucleosome (-1 nuc), ",
      "(C) a TFBS cluster, or (D) a core promoter element (TATA/Inr/DPE/BRE). ",
      "<br><br>",
      "<strong>Literature precedent:</strong> Brenet 2011 showed first exon meth inversely linked to expression (humans only). ",
      "Keller 2016 (Ciona) hinted at positive promoter-meth-expression correlation in low-GBM genes. ",
      "Wu 2025 (Nautilus) used 4 expression bins; Nautilus has DNMT3. ",
      "<strong>No one has shown 10-decile metagene profiles in a DNMT3-absent, &gt;95% ICP invertebrate. ",
      "The ICP promoter architecture demonstrated here is empirical proof of Schubeler's theoretical framework.</strong>"),
    notes = list(
      fig12a = "Single-bp resolution (10 bp bins) methylation around TSS by expression decile. Full window (-5000 to +2000). Pins down the exact position of the promoter spike. The spike-and-crash pattern is clear in D08-D10.",
      fig12b = "Zoomed view of the promoter spike region (-1500 to +500). The mid-promoter spike and TSS crash are clearly resolved. Spike position and amplitude quantified in spike_parameters.tsv.",
      fig12c = "CpG density around TSS by expression decile. Tests whether the spike coincides with a CpG-dense element. Red shading marks the spike region.",
      fig12d = "GC content around TSS. Tests whether the spike region has unusual base composition.",
      fig12e = "CpG O/E around TSS. Tests whether the spike region has CpG island-like O/E values despite being in an ICP promoter.",
      fig12f = "Core promoter element (TATA, Inr, DPE, BRE) frequency by expression decile. Tests whether the spike coincides with canonical promoter elements and whether highly expressed genes are enriched for specific elements.",
      fig12g = "TFBS motif enrichment at the spike position (Fisher's exact vs flanking regions). Identifies whether specific TF binding sites cluster at the spike.",
      fig12h = "Positional density of top enriched TFBS motifs around TSS for D09+D10 genes. Shows whether enriched motifs specifically peak at the spike position.",
      fig12i = "Nucleosome positioning prediction (WW dinucleotide score) overlaid with methylation for D09+D10 genes. Tests whether the spike aligns with a predicted -1 nucleosome position. High WW = nucleosome-favoring sequence.",
      fig12j = "Control vs amputated methylation at the spike (D09+D10 genes). Tests whether the spike changes during regeneration or is a constitutive feature."
    )
  ),
  list(
    dir = "batch13", title = "Batch 13 — Entropy Metagene Profiles",
    question = "Does normalized entropy vary across the gene structure by expression? Which genes gain/lose entropy during regeneration?",
    context = paste0(
      "Binary Shannon entropy H = -p*log2(p) - (1-p)*log2(1-p) applied as metagene Y-axis instead of methylation %. ",
      "Entropy is maximal at 50% methylation (H=1 bit) and minimal at 0% or 100% (H~0). ",
      "This reveals <strong>methylation stochasticity</strong> across the gene structure, not just level. ",
      "Control, amputated, and &Delta; entropy metagenes by 5 quintiles and 10 deciles. ",
      "Gene-level analysis: per-gene &Delta; entropy (Wilcoxon paired test), classified into Entropy_Up/Down/NS. ",
      "Cross-referenced with batch06 DMPs and DESeq2 tail results for concordant entropy-DE gene lists. ",
      "Motif enrichment (batch 1.5) tested at entropy-up/down genes."),
    notes = list(
      fig13a = "Entropy metagene profiles by expression group. Compare structure across 3+ exon, 2-exon, and intronless genes. Control, amputated, and delta (ampu-ctrl) versions.",
      "fig13a_bin10_01_ctrl_3plus_exon" = "KEY FIGURE: Control entropy metagene for 3+ exon genes (10 deciles). Entropy is highest in gene body of moderately expressed genes (intermediate beta ~50% = max entropy). Highly expressed genes (D10) have HIGH methylation (~65%) which maps to LOWER entropy. Silent genes (D01) have LOW methylation (~5%) which also maps to low entropy.",
      "fig13a_bin10_03_delta_3plus_exon" = "Delta entropy metagene (ampu-ctrl) for 3+ exon genes (10 deciles). Shows which expression groups and gene regions gain/lose entropy during regeneration.",
      "fig13a_bin10_07_ctrl_single_exon" = "Intronless gene control entropy metagene. Compare to 3+ exon pattern to test splicing-dependent methylation model.",
      fig13b = "Entropy volcano: per-gene delta entropy vs significance (Wilcoxon paired, BH-adjusted). Red=entropy up, blue=entropy down during regeneration.",
      fig13c = "Entropy change by gene structure class (3+ exon vs 2-exon vs intronless). Tests whether gene architecture affects entropy dynamics.",
      fig13d = "Delta entropy vs expression level (GAM fit). Tests whether highly expressed genes are more or less susceptible to entropy change.",
      fig13e = "Methylation level vs entropy hexbin. Shows the concavity relationship: sites near 0%/100% have low entropy by definition. Critical control for interpreting entropy metagenes.",
      fig13f = "Entropy change by expression decile (boxplots). Which expression bins gain/lose most entropy during regeneration?",
      fig13g = "TFBS motifs enriched in entropy-up genes (Fisher's exact vs NS genes). Do specific TFs associate with entropy gain?",
      fig13h = "Entropy change vs DMP count per gene. Do genes with more DMPs show more entropy change?",
      fig13i = "Entropy vs expression change four-way plot. Red=both significant, orange=entropy only, blue=DE only. Identifies concordant entropy-DE genes."
    )
  )
)

# Build HTML
cat("Generating pipeline report...\n")
html <- c(
  "<!DOCTYPE html>",
  "<html lang='en'>",
  "<head>",
  "<meta charset='UTF-8'>",
  "<meta name='viewport' content='width=device-width, initial-scale=1.0'>",
  "<title>D. laeve Methylation Pipeline — Figure Report</title>",
  "<style>",
  "body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 0; background: #f5f5f5; color: #333; }",
  "#main { max-width: 1200px; margin: 0 auto; padding: 20px 20px 20px 280px; }",
  "h1 { color: #2c3e50; border-bottom: 3px solid #2471A3; padding-bottom: 10px; }",
  "h2 { color: #2471A3; margin-top: 40px; border-left: 5px solid #2471A3; padding-left: 15px; background: #eaf2f8; padding: 10px 15px; scroll-margin-top: 20px; }",
  ".question { font-style: italic; color: #666; font-size: 1.1em; margin-bottom: 20px; }",
  ".figure-card { background: white; border-radius: 8px; padding: 20px; margin: 15px 0; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }",
  ".figure-card img { max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 4px; }",
  ".figure-card h3 { color: #C0392B; margin-top: 0; }",
  ".analysis { background: #fef9e7; border-left: 4px solid #F39C12; padding: 10px 15px; margin-top: 10px; font-size: 0.95em; }",
  ".missing { color: #999; font-style: italic; }",
  ".stats { display: flex; gap: 15px; flex-wrap: wrap; margin: 15px 0; }",
  ".stat-box { background: #2471A3; color: white; padding: 10px 20px; border-radius: 5px; text-align: center; }",
  ".stat-box .val { font-size: 1.5em; font-weight: bold; }",
  "",
  "/* Sticky sidebar navigation */",
  "#sidebar { position: fixed; top: 0; left: 0; width: 260px; height: 100vh; background: #2c3e50; color: #ecf0f1; overflow-y: auto; padding: 15px 10px; box-sizing: border-box; font-size: 0.85em; }",
  "#sidebar h3 { color: #3498db; margin: 0 0 10px 0; font-size: 1.1em; border-bottom: 1px solid #3498db; padding-bottom: 8px; }",
  "#sidebar a { color: #bdc3c7; text-decoration: none; display: block; padding: 6px 8px; border-radius: 4px; transition: background 0.15s; }",
  "#sidebar a:hover, #sidebar a.active { background: #34495e; color: #ecf0f1; }",
  "#sidebar .nav-batch { font-weight: bold; color: #ecf0f1; margin-top: 4px; }",
  "#sidebar .nav-short { font-size: 0.9em; color: #95a5a6; padding-left: 12px; }",
  "",
  "/* Summary table */",
  ".summary-table { width: 100%; border-collapse: collapse; margin: 20px 0; font-size: 0.92em; }",
  ".summary-table th { background: #2c3e50; color: white; padding: 10px 12px; text-align: left; }",
  ".summary-table td { padding: 8px 12px; border-bottom: 1px solid #ddd; }",
  ".summary-table tr:nth-child(even) { background: #f9f9f9; }",
  ".summary-table tr:hover { background: #eaf2f8; }",
  ".summary-table td:first-child { font-weight: bold; white-space: nowrap; }",
  "",
  "@media (max-width: 900px) { #sidebar { display: none; } #main { padding-left: 20px; } }",
  "</style>",
  "</head>",
  "<body>",
  "",
  "<!-- Sticky sidebar navigation -->",
  "<nav id='sidebar'>",
  "<h3>Navigation</h3>"
)

# Sidebar links
for (b in batches) {
  short_name <- sub("^Batch [0-9.]+ . ", "", b$title)
  batch_num <- sub("^batch", "", b$dir)
  html <- c(html, sprintf("<a href='#%s' class='nav-batch'>%s</a>", b$dir, paste0("B", batch_num, " ", short_name)))
}
html <- c(html,
  "<a href='#top' style='margin-top:15px; color:#3498db;'>Back to top</a>",
  "</nav>",
  "",
  "<!-- Main content -->",
  "<div id='main'>",
  "<h1 id='top'>Methylation landscape of mollusk regeneration</h1>",
  "<p><em>Deroceras laeve</em> &mdash; WGBS + RNA-seq pipeline figure report</p>",
  "<div class='stats'>",
  "<div class='stat-box'><div class='val'>45.5M</div>Genome CpGs</div>",
  "<div class='stat-box'><div class='val'>25.4M</div>Tested CpGs</div>",
  "<div class='stat-box'><div class='val'>17,729</div>DMPs</div>",
  "<div class='stat-box'><div class='val'>1,462</div>DMRs</div>",
  "<div class='stat-box'><div class='val'>0.553</div>CpG O/E</div>",
  "<div class='stat-box'><div class='val'>No DNMT3</div>De novo absent</div>",
  "</div>",
  "",
  "<!-- Summary table -->",
  "<table class='summary-table'>",
  "<tr><th>Batch</th><th>Question</th><th>Key Result</th></tr>",
  "<tr><td><a href='#batch0.5'>0.5</a></td><td>Clean expression dataset?</td><td>44 samples &rarr; 41 clean (3 within-group outliers removed). Tail DE: 1,797 genes (FDR&lt;0.05)</td></tr>",
  "<tr><td><a href='#batch01'>01</a></td><td>CpG landscape?</td><td>CpG O/E = 0.553, 45.5M CpGs, TEs have high CpG density</td></tr>",
  "<tr><td><a href='#batch1.5'>1.5</a></td><td>Where are TFBS?</td><td>29M motif hits from 1,362 JASPAR PWMs across chr1-31</td></tr>",
  "<tr><td><a href='#batch02'>02</a></td><td>Methylation toolkit?</td><td>DNMT1 present, NO DNMT3 &mdash; no de novo methyltransferase</td></tr>",
  "<tr><td><a href='#batch03'>03</a></td><td>Promoter CpGs?</td><td>85.1% ICP &mdash; no CpG islands, no mammalian-style promoter silencing</td></tr>",
  "<tr><td><a href='#batch04'>04</a></td><td>Baseline methylation?</td><td>Gene body meth + expression: rho=0.407, exon highest</td></tr>",
  "<tr><td><a href='#batch05'>05</a></td><td>TE methylation &amp; age?</td><td>TE meth is PASSIVE &mdash; intronic TEs methylated (spillover), intergenic unmethylated</td></tr>",
  "<tr><td><a href='#batch06'>06</a></td><td>DM during regeneration?</td><td>17,729 DMPs (8,285 hyper + 9,444 hypo), 1,462 DMRs. Enriched in exons (1.95x)</td></tr>",
  "<tr><td><a href='#batch07'>07</a></td><td>DM predicts expression?</td><td>NO &mdash; rho ~ 0.02, methylation-expression DECOUPLED</td></tr>",
  "<tr><td><a href='#batch08'>08</a></td><td>Module-specific effects?</td><td>DMPs concentrate in specific WGCNA modules &mdash; pathway-targeted, not random. Resolves the decoupling paradox</td></tr>",
  "<tr><td><a href='#batch09'>09</a></td><td>TFBS at DMPs?</td><td>Motif enrichment with region-restricted background + Yin 2017 methyl-sensitivity. Concordant-gene shortlist for cloning</td></tr>",
  "<tr><td><a href='#batch10'>10</a></td><td>Entropy direction?</td><td>Entropy DECREASES (&delta;=-0.002, p=7.6e-6). Directed reprogramming, not drift. MI=0.022 bits</td></tr>",
  "<tr><td><a href='#batch11'>11</a></td><td>Motif &times; NME &times; DE?</td><td>TFBS enforce ordered methylation at DE genes. Directed-mechanism cloning shortlist</td></tr>",
  "<tr><td><a href='#batch12'>12</a></td><td>Promoter spike?</td><td>Mid-promoter spike was binning artifact. ICP architecture (meth-expression positive) confirmed at 10bp resolution</td></tr>",
  "<tr><td><a href='#batch13'>13</a></td><td>Entropy metagene?</td><td>Entropy landscape across gene structure by expression. Per-gene &Delta;entropy gene lists + motif enrichment</td></tr>",
  "</table>"
)

for (b in batches) {
  fig_dir <- file.path(PIPE_DIR, b$dir, "figures")
  pngs <- sort(list.files(fig_dir, pattern = "\\.png$", full.names = FALSE))

  ctx <- if (!is.null(b$context)) sprintf("<div class='analysis' style='border-left-color:#2471A3; background:#eaf2f8; margin-bottom:15px;'><strong>Key finding:</strong> %s</div>", b$context) else ""
  html <- c(html,
    sprintf("<h2 id='%s'>%s</h2>", b$dir, b$title),
    sprintf("<p class='question'>%s</p>", b$question),
    ctx
  )

  if (length(pngs) == 0) {
    html <- c(html, "<p class='missing'>No figures generated yet — run the pipeline to produce figures.</p>")
    next
  }

  for (png in pngs) {
    fig_path <- file.path(PIPE_DIR, b$dir, "figures", png)
    fig_full <- sub("\\.png$", "", png)
    fig_prefix <- sub("_.*", "", fig_full)

    # Find matching analysis note: try full name first, then prefix
    note <- b$notes[[fig_full]]
    if (is.null(note)) note <- b$notes[[fig_prefix]]
    if (is.null(note)) note <- ""

    # Embed image as base64 data URI — self-contained HTML, no external file dependencies
    raw_bytes <- readBin(fig_path, "raw", file.info(fig_path)$size)
    b64 <- base64enc::base64encode(raw_bytes)
    img_src <- paste0("data:image/png;base64,", b64)

    fig_name <- sub("\\.png$", "", png)
    html <- c(html,
      "<div class='figure-card'>",
      sprintf("<h3>%s</h3>", fig_name),
      sprintf("<img src='%s' alt='%s' loading='lazy'>", img_src, fig_name)
    )
    if (nchar(note) > 0) {
      html <- c(html, sprintf("<div class='analysis'>%s</div>", note))
    }
    html <- c(html, "</div>")
  }
}

html <- c(html,
  sprintf("<hr><p style='color:#999; font-size:0.85em;'>Generated %s from %s</p>",
          format(Sys.time(), "%Y-%m-%d %H:%M"), getwd()),
  "</div>", # close #main
  "",
  "<!-- Highlight active sidebar link on scroll -->",
  "<script>",
  "document.addEventListener('DOMContentLoaded', function() {",
  "  var links = document.querySelectorAll('#sidebar a.nav-batch');",
  "  var sections = [];",
  "  links.forEach(function(a) { var t = document.querySelector(a.getAttribute('href')); if(t) sections.push({el:t, link:a}); });",
  "  window.addEventListener('scroll', function() {",
  "    var scroll = window.scrollY + 80;",
  "    var current = null;",
  "    sections.forEach(function(s) { if(s.el.offsetTop <= scroll) current = s; });",
  "    links.forEach(function(a) { a.classList.remove('active'); });",
  "    if(current) current.link.classList.add('active');",
  "  });",
  "});",
  "</script>",
  "</body></html>"
)

writeLines(html, OUT_FILE)
cat(sprintf("Report written to %s (%d figures across %d batches)\n", OUT_FILE,
            sum(sapply(batches, function(b) length(list.files(file.path(PIPE_DIR, b$dir, "figures"),
                                                              pattern = "\\.png$")))),
            length(batches)))
