#!/bin/bash
# Install to /tmp and run the script in a single call
LOG=/mnt/c/Users/rafae/Projects/STANDBY/methylation_tutorials/results/run_log.txt
mkdir -p /mnt/c/Users/rafae/Projects/STANDBY/methylation_tutorials/results
mkdir -p /tmp/claude-1000

echo "=== Starting at $(date) ===" > $LOG

# Install micromamba + R env if needed
if [ ! -f /tmp/claude-1000/mamba/envs/renv/bin/Rscript ]; then
  echo "Installing R environment..." | tee -a $LOG
  mkdir -p /tmp/claude-1000/bin
  if [ ! -f /tmp/claude-1000/bin/micromamba ]; then
    curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C /tmp/claude-1000/ bin/micromamba
  fi
  export MAMBA_ROOT_PREFIX=/tmp/claude-1000/mamba
  /tmp/claude-1000/bin/micromamba create -n renv -y -c conda-forge -c bioconda \
    r-base=4.4 r-ggplot2 r-dplyr r-tidyr r-data.table r-rcolorbrewer r-scales \
    r-ggrepel r-patchwork r-pheatmap r-reshape2 r-stringr r-gridextra r-cairo \
    bioconductor-genomicranges bioconductor-rtracklayer bioconductor-biostrings \
    bioconductor-bsseq bioconductor-dss 2>&1 | tail -3
  echo "Env installed at $(date)" | tee -a $LOG
else
  echo "Env exists" | tee -a $LOG
fi

export PATH="/tmp/claude-1000/mamba/envs/renv/bin:/usr/bin:/bin:$PATH"
export TMPDIR=/tmp/claude-1000
export R_LIBS_SITE="/tmp/claude-1000/mamba/envs/renv/lib/R/library"
export R_LIBS_USER="/tmp/claude-1000/rlibs"

echo "=== Script start at $(date) ===" >> $LOG
Rscript /mnt/c/Users/rafae/Projects/STANDBY/methylation_tutorials/part3_downstream_analysis.R >> $LOG 2>&1
echo "=== Done at $(date), exit=$? ===" | tee -a $LOG
