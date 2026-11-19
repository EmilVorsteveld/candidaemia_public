#!/usr/bin/env bash
#
# Gene body coverage (5' -> 3') of the Kinnex genome BAMs over housekeeping genes, with RSeQC geneBody_coverage.py.
# Plotted in 01_kinnex_QC.Rmd ("Gene body coverage").
#
# Usage:
#   bash gene_body_coverage.sh            # all BAMs in $BAM_DIR, $JOBS at a time
#   bash gene_body_coverage.sh ALEX1_D0   # one sample
#
# Before the first run, copy the gene file to the server ($BED below). It is the RSeQC hg38 housekeeping set
# (3,802 RefSeq transcripts, BED12, chr-prefixed like the BAMs).
# Output: $OUT/<sample>.geneBodyCoverage.txt (one row, percentiles 1-100). Copy $OUT back to
# data/kinnex/gene_body_coverage/ in the project.
set -euo pipefail

# ---- Paths: edit these ----
BAM_DIR=/mnt/data/emil/candidaemia/bam   # Kinnex genome BAMs, <sample>.bam (same as ase_haplotag.sh)
BED=/mnt/data/emil/candidaemia/hg38.HouseKeepingGenes.bed
OUT=/mnt/data/emil/candidaemia/gene_body_coverage
JOBS=12                                  # geneBody_coverage.py is single-threaded: one core per sample
CONDA_ENV=rseqc                          # conda create -n rseqc -c bioconda -c conda-forge rseqc

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"
command -v geneBody_coverage.py >/dev/null || { echo "geneBody_coverage.py not found in conda env $CONDA_ENV" >&2; exit 1; }
[ -f "$BED" ] || { echo "Gene file $BED not found: copy hg38.HouseKeepingGenes.bed to the server first" >&2; exit 1; }
mkdir -p "$OUT"

run_one() {
  local s=$1 bam=$BAM_DIR/$1.bam
  [ -f "$OUT/$s.geneBodyCoverage.txt" ] && { echo "$s done already, skipped"; return 0; }
  [ -f "$bam.bai" ] || [ -f "$BAM_DIR/$s.bai" ] || samtools index "$bam"
  geneBody_coverage.py -r "$BED" -i "$bam" -o "$OUT/$s" > "$OUT/$s.log" 2>&1
  echo "$s finished"
}
export -f run_one
export BAM_DIR BED OUT

if [ $# -ge 1 ]; then
  run_one "$1"
else
  ls "$BAM_DIR"/*.bam | xargs -n1 basename | sed 's/\.bam$//' | parallel -j "$JOBS" run_one {}
fi
