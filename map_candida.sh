#!/usr/bin/env bash
# Map unmapped reads from Kinnex BAMs to a combined Candida reference, in parallel.
# Usage: ./map_unmapped_to_candida.sh [dna|rna]     (default: dna)
set -euo pipefail

MODE="${1:-dna}"
IN=/mnt/data/emil/candidaemia/bam
OUT=/mnt/data/emil/candidaemia/bam_candida
REF=/mnt/data/emil/candidaemia/candida_refs/candida_species_ref.fna

# 64 cores: 16 samples at once x ~4 threads each
JOBS=16
MM_THREADS=3     # minimap2 threads per sample
ST_THREADS=1     # extra samtools threads per sample

case "$MODE" in
  dna) MM_PRESET="-ax map-hifi" ;;
  rna) MM_PRESET="-ax splice:hq -uf" ;;
  *) echo "mode must be dna or rna" >&2; exit 1 ;;
esac

mkdir -p "$OUT"/logs
MMI="$OUT/candida_species_ref.mmi"
[ -s "$MMI" ] || minimap2 -t 16 -d "$MMI" "$REF"

process() {
  bam="$1"
  s=$(basename "$bam" .bam)
  final="$OUT/${s}.candida.bam"
  [ -s "${final}.bai" ] && { echo "skip $s (done)"; return 0; }

  # unmapped reads -> minimap2 -> sorted BAM, streamed (no intermediate FASTQ)
  samtools fastq -@ "$ST_THREADS" -f 4 -F 0x900 "$bam" 2> "$OUT/logs/${s}.fastq.log" \
    | minimap2 $MM_PRESET -t "$MM_THREADS" "$MMI" - 2> "$OUT/logs/${s}.minimap2.log" \
    | samtools sort -@ "$ST_THREADS" -T "$OUT/${s}.tmp" -o "${final}.part" -

  mv "${final}.part" "$final"
  samtools index -@ "$ST_THREADS" "$final"
  samtools idxstats "$final" > "$OUT/${s}.idxstats.txt"
  echo "done $s"
}
export -f process
export OUT MMI MM_PRESET MM_THREADS ST_THREADS

ls "$IN"/*.bam | xargs -P "$JOBS" -I{} bash -c 'process "$@"' _ {}

# combine per-sample idxstats into one table: sample, contig, length, mapped
for f in "$OUT"/*.idxstats.txt; do
  s=$(basename "$f" .idxstats.txt)
  awk -v s="$s" -v OFS='\t' '$1!="*"{print s,$1,$2,$3}' "$f"
done > "$OUT/all_samples.idxstats.tsv"

echo "All samples finished. Output in $OUT"