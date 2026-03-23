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
| Promoter classes | 91.6% ICP, 6.5% HCP, 1.9% LCP |
| Predicted TFs | 2,217 genes (DeepTFactor) |
| JASPAR TFBS motifs | 1,362 metazoan PWMs |
| Baseline gene body meth vs expression | rho = 0.307 (p < 10^-80) |
| Differential meth vs expression MI | 0.000087 bits (p = 0.55) |

## Repository structure

```
STANDBY/
  methylation_pipeline/     ← MAIN: supplementary code for paper (reviewers see this)
    _config.R               — shared paths, colors, helpers
    batch01/ {code/,data/,figures/}  — Genomic CpG annotation
    batch02/ {code/,data/,figures/}  — Methylation machinery
    batch03/ {code/,data/,figures/}  — Promoter Weber classification (HCP/ICP/LCP)
    batch04/ {code/,data/,figures/}  — Baseline methylation + gene body meth vs expression
    batch05/ {code/,data/,figures/}  — TE methylation + evolutionary age
    batch06/ {code/,data/,figures/}  — Differential methylation (DMP/DMR)
    batch07/ {code/,data/,figures/}  — Methylation-expression decoupling (MI, correlation)
    batch08/ {code/,data/,figures/}  — WGCNA modules + methylation enrichment
    batch09/ {code/,data/,figures/}  — TF + GENIE3 + motif integration
    batch10/ {code/,data/,figures/}  — Entropy analysis
  cluster/scripts/          — R + SLURM for HPC jobs
  genome/cache/             — RDS caches (genome, GFF, BSseq, DMLtest, TE, promoters)
  papers/                   — Reference PDFs + paper drafts
  PROTOCOL.md               — Internal protocol (not for reviewers)
```

- **No HTML in batch scripts.** Batches produce data (TSV) + figures (PNG+PDF) only.
- Each batch sources `methylation_pipeline/_config.R` for shared config.

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

### Cluster (HPC)

`/mnt/data/alfredvar/` — shared data root

| File | Cluster path |
|------|-------------|
| CpG reports (.gz) | `/mnt/data/alfredvar/jmiranda/50-Genoma/51-Metilacion/09_methylation_calls/` |
| GFF | `/mnt/data/alfredvar/30-Genoma/31-Alternative_Annotation_EviAnn/` |
| TE age | `/mnt/data/alfredvar/30-Genoma/32-Repeats/age_of_transposons/collapsed_te_age_data.tsv` |
| HTSeq counts | `/mnt/data/alfredvar/jmiranda/20-Transcriptomic_Bulk/25-metaAnalysisTranscriptome/counts_HTseq_EviAnn/` |
| GENIE3 (full) | `/mnt/data/alfredvar/wgutierrez/genie3_2/genie3_all_links.tsv` |

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

## Running locally (Windows)

```bash
"C:\Program Files\R\R-4.5.2\bin\Rscript.exe" methylation_pipeline/batchNN/code/NN_script.R
```

BSgenome package: `BSgenome.Dlaeve.NCBI.dlgm` (installed on Windows R)

## Running on cluster

```bash
cd Deroceras-Leave/cluster
dos2unix scripts/*.slurm scripts/*.sh
sbatch scripts/generate_genome_cache.slurm   # run first if genome_chr1_31.rds missing
sbatch scripts/dmltest_full.slurm            # full DMLtest (32 GB)
sbatch scripts/batch1.5_expanded_motif.slurm # TFBS motif positions (64 GB)
```

SLURM template:
```bash
#!/usr/bin/env bash
#SBATCH --job-name=Rjob
#SBATCH --output=slurm_%j.out
#SBATCH --error=slurm_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
set -euo pipefail
Rscript scripts/my_script.R
```

## Git

- Remote: `https://github.com/somnasgarden/Deroceras-Leave.git`
- Author: `GIT_AUTHOR_NAME="rlopezt" GIT_AUTHOR_EMAIL="rafae@users.noreply.github.com"`
- HPC clone: `/mnt/data/alfredvar/rlopezt/repos/Deroceras-Leave/`
