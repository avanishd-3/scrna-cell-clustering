#!/usr/bin/env bash
# Download the GEO GSE337533 input data used by this pipeline.
#   data/raw/  : 22 gene-expression libraries (barcodes / features / matrix), ~3.2 GB
#   data/vdj/  : 22 BCR filtered_contig_annotations tables (validation only), ~23 MB
# Files already present and passing a gzip integrity check are skipped, so the script can be re-run.
# usage: bash scripts/00_download.sh [data_dir]   (default: data)
set -euo pipefail
cd "$(dirname "$0")/.."

SERIES=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE337nnn/GSE337533/suppl
SAMPLES=https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9857nnn
DATA=${1:-data}
mkdir -p "$DATA/raw" "$DATA/vdj"

fetch() {  # fetch <url> <dest>
  if [ -s "$2" ] && gzip -t "$2" 2>/dev/null; then return 0; fi
  echo "downloading $(basename "$2")"
  curl -sSfL --retry 3 -o "$2.part" "$1" && mv "$2.part" "$2"
  gzip -t "$2"
}

# 1. Gene-expression matrices (series-level supplementary files)
gex=$(curl -sSfL "$SERIES/" | grep -o 'GSE337533_[A-Za-z0-9_]*_\(barcodes\.tsv\|features\.tsv\|matrix\.mtx\)\.gz' | sort -u)
n=0; for f in $gex; do fetch "$SERIES/$f" "$DATA/raw/$f"; n=$((n+1)); done
echo "gene expression: $n files in $DATA/raw (expected 66)"

# 2. BCR contig annotations (sample-level supplementary files)
vdj=$(curl -sSfL "$SERIES/filelist.txt" | awk -F'\t' '$2 ~ /_filtered_contig_annotations\.csv\.gz$/ {print $2}')
n=0; for f in $vdj; do gsm=${f%%_*}; fetch "$SAMPLES/$gsm/suppl/$f" "$DATA/vdj/$f"; n=$((n+1)); done
echo "BCR contigs: $n files in $DATA/vdj (expected 22)"
