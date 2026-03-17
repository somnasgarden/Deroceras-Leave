#!/bin/bash
# =============================================================================
# Install HOMER and configure D. laeve custom genome
# =============================================================================
# Idempotent — safe to rerun. Installs to $CLUSTER_ROOT/homer/
#
# Usage (from cluster/ directory):
#   bash scripts/setup_homer.sh
#
# Requires: perl, wget or curl
# =============================================================================

set -euo pipefail

CLUSTER_ROOT="${CLUSTER_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
HOMER_DIR="${CLUSTER_ROOT}/homer"
HOMER_BIN="${HOMER_DIR}/bin"

echo "=== HOMER Setup ==="
echo "Cluster root: ${CLUSTER_ROOT}"
echo "HOMER dir:    ${HOMER_DIR}"

# ---- 1. Install HOMER if not present ----
if [ -x "${HOMER_BIN}/homer" ] || [ -x "${HOMER_BIN}/findMotifs.pl" ]; then
  echo "HOMER already installed at ${HOMER_DIR}"
else
  echo "Installing HOMER..."
  mkdir -p "${HOMER_DIR}"
  cd "${HOMER_DIR}"

  # Download configureHomer.pl
  if command -v wget &>/dev/null; then
    wget -q http://homer.ucsd.edu/homer/configureHomer.pl
  elif command -v curl &>/dev/null; then
    curl -sO http://homer.ucsd.edu/homer/configureHomer.pl
  else
    echo "ERROR: Need wget or curl to download HOMER"
    exit 1
  fi

  # Install HOMER base
  perl configureHomer.pl -install homer

  echo "HOMER installed successfully"
  cd "${CLUSTER_ROOT}"
fi

export PATH="${HOMER_BIN}:${PATH}"
echo "HOMER version: $(perl ${HOMER_BIN}/homer 2>&1 | head -1 || echo 'installed')"

# ---- 2. Configure custom D. laeve genome ----
GENOME_NAME="dlaeve"
HOMER_GENOME_DIR="${HOMER_DIR}/data/genomes/${GENOME_NAME}"

# Find the genome FASTA — the R script exports it from BSgenome.Dlaeve.NCBI.dlgm
# to cluster/genome/derLaeGenome_chr1_31.fasta (chr1-31 only).
# If R hasn't run yet, check common cluster locations.
GENOME_FASTA="${GENOME_FASTA:-}"

if [ -z "${GENOME_FASTA}" ] || [ ! -f "${GENOME_FASTA}" ]; then
  for candidate in \
    "${CLUSTER_ROOT}/genome/derLaeGenome_chr1_31.fasta" \
    "/mnt/data/alfredvar/30-Genoma/derLaeGenome_chr1_31.fasta"; do
    if [ -f "${candidate}" ]; then
      GENOME_FASTA="${candidate}"
      break
    fi
  done
fi

if [ -z "${GENOME_FASTA}" ] || [ ! -f "${GENOME_FASTA}" ]; then
  echo "NOTE: Genome FASTA not found yet."
  echo "  The R script will export it from BSgenome.Dlaeve.NCBI.dlgm on first run."
  echo "  HOMER genome config will be done after the R script creates the FASTA."
  echo "  Skipping HOMER genome setup for now."
  exit 0
fi

echo "Genome FASTA: ${GENOME_FASTA}"

if [ -d "${HOMER_GENOME_DIR}" ] && [ -f "${HOMER_GENOME_DIR}/genome.fa" ]; then
  echo "Custom genome '${GENOME_NAME}' already configured"
else
  echo "Configuring custom genome '${GENOME_NAME}'..."
  mkdir -p "${HOMER_GENOME_DIR}"

  # Filter to chr1-31 only (exclude HiC_scaffold_*)
  echo "Filtering genome to chr1-31..."
  awk 'BEGIN{keep=0} /^>/{keep=0; for(i=1;i<=31;i++){if($0 ~ "^>chr"i"$" || $0 ~ "^>chr"i" "){keep=1;break}}} keep{print}' \
    "${GENOME_FASTA}" > "${HOMER_GENOME_DIR}/genome.fa"

  # If the FASTA is already filtered (no scaffolds), just symlink
  if [ ! -s "${HOMER_GENOME_DIR}/genome.fa" ]; then
    echo "  Filter produced empty file — source may already be filtered. Symlinking."
    rm "${HOMER_GENOME_DIR}/genome.fa"
    ln -sf "$(realpath "${GENOME_FASTA}")" "${HOMER_GENOME_DIR}/genome.fa"
  fi

  # Create chromosome sizes file
  echo "Creating chrom.sizes..."
  if command -v samtools &>/dev/null; then
    samtools faidx "${HOMER_GENOME_DIR}/genome.fa"
    cut -f1,2 "${HOMER_GENOME_DIR}/genome.fa.fai" > "${HOMER_GENOME_DIR}/chrom.sizes"
  else
    # Parse FASTA headers for sizes (slower but no samtools needed)
    awk '/^>/{if(name)print name"\t"len; name=substr($1,2); len=0; next} {len+=length($0)} END{if(name)print name"\t"len}' \
      "${HOMER_GENOME_DIR}/genome.fa" > "${HOMER_GENOME_DIR}/chrom.sizes"
  fi

  # Preparsed sequences for faster HOMER scanning
  echo "Preparsing genome (this takes a few minutes)..."
  mkdir -p "${HOMER_GENOME_DIR}/preparsed"
  perl "${HOMER_BIN}/homerTools" extract "${HOMER_GENOME_DIR}/genome.fa" \
    -o "${HOMER_GENOME_DIR}/preparsed/" 2>/dev/null || true

  echo "Custom genome '${GENOME_NAME}' configured"
fi

echo ""
echo "=== HOMER Setup Complete ==="
echo "Add to PATH: export PATH=${HOMER_BIN}:\$PATH"
echo "Genome name for HOMER commands: ${GENOME_NAME}"
echo "Genome dir: ${HOMER_GENOME_DIR}"
