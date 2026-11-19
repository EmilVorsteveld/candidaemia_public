#!/usr/bin/env bash
#
# Gene fusions with pbfusion (PacBio) on the Kinnex genome BAMs. GitHub issue #8.
# Downstream analysis: 06_pbfusion.Rmd.
#
# Install once:
#   conda create -n pbfusion -c bioconda pbfusion samtools parallel
#   conda activate pbfusion; pbfusion discover --help    # check that the annotation option is --gtf
#
# Usage:
#   bash pbfusion.sh             # all BAMs in $BAM_DIR, $JOBS at a time
#   bash pbfusion.sh ALEX1_D0    # one sample
#
# Input: the pbmm2-aligned FLNC genome BAMs (<sample>.bam), the same ones used by isocall and the ASE scripts.
# Every event is kept (--min-fusion-quality LOW, so readthroughs and low-support calls too); filtering is done in R.
# Output: $OUT/<sample>.breakpoints.groups.bed per sample, and $OUT/pbfusion_all_samples.tsv with all samples
# (first column: sample). Copy pbfusion_all_samples.tsv to data/kinnex/pbfusion/ in the project.
set -euo pipefail

# ---- Paths: edit these ----
BAM_DIR=/mnt/data/emil/candidaemia/bam                       # Kinnex genome BAMs, <sample>.bam
GTF=/mnt/data/emil/candidaemia/gencode.v49.annotation.gtf.gz # GENCODE v49, as used by isocall (must match hg38)
OUT=/mnt/data/emil/candidaemia/pbfusion
JOBS=4                                                        # samples at a time
THREADS=16                                                    # threads per sample
CONDA_ENV=pbfusion

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"
for tool in pbfusion samtools parallel; do command -v $tool >/dev/null || { echo "$tool not found in conda env $CONDA_ENV" >&2; exit 1; }; done
[ -f "$GTF" ] || { echo "GTF $GTF not found" >&2; exit 1; }
mkdir -p "$OUT"

# ---- Annotation cache, once (built with the installed pbfusion: v0.6 changed the cache format) ----
CACHE=$OUT/gencode.v49.gtf.bin
if [ ! -f "$CACHE" ]; then
  if [[ "$GTF" == *.gz ]]; then
    gunzip -c "$GTF" > "$OUT/gencode.v49.gtf"
    pbfusion gff-cache --gtf "$OUT/gencode.v49.gtf" --gtf-out "$CACHE"
    rm "$OUT/gencode.v49.gtf"
  else
    pbfusion gff-cache --gtf "$GTF" --gtf-out "$CACHE"
  fi
fi

run_one() {
  local s=$1 bam=$BAM_DIR/$1.bam
  [ -f "$OUT/$s.breakpoints.groups.bed" ] && { echo "$s done already, skipped"; return 0; }
  [ -f "$bam" ] || { echo "$s: $bam not found" >&2; return 1; }
  [ -f "$bam.bai" ] || [ -f "$BAM_DIR/$s.bai" ] || samtools index "$bam"
  samtools view -H "$bam" | grep -q '^@PG.*pbmm2' || echo "$s: warning, no pbmm2 in the BAM header (pbfusion expects pbmm2 ISOSEQ alignments)" >&2
  pbfusion discover --gtf "$CACHE" --output-prefix "$OUT/$s" --threads "$THREADS" --min-fusion-quality LOW "$bam" > "$OUT/$s.log" 2>&1
  echo "$s finished"
}
export -f run_one
export BAM_DIR OUT CACHE THREADS

if [ $# -ge 1 ]; then
  run_one "$1"
  exit 0
fi
ls "$BAM_DIR"/*.bam | xargs -n1 basename | sed 's/\.bam$//' | parallel -j "$JOBS" run_one {}

# ---- All samples in one table: column header once, sample name as the first column ----
ALL=$OUT/pbfusion_all_samples.tsv
first=$(ls "$OUT"/*.breakpoints.groups.bed | head -1)
grep '^#' "$first" | tail -1 | sed 's/^#*/#sample\t/' > "$ALL"
for f in "$OUT"/*.breakpoints.groups.bed; do
  s=$(basename "$f" .breakpoints.groups.bed)
  grep -v '^#' "$f" | sed "s/^/$s\t/" >> "$ALL" || true   # grep exits 1 when a sample has no events
done
echo "Done: $(grep -vc '^#' "$ALL") events in $ALL"
