#!/usr/bin/env bash
# Per-sample, per-species read counts and fractions from the Candida-mapped BAMs.
# Usage: bash quantify_candida.sh [min_mapq]     (default 20)
set -euo pipefail

IN=/mnt/data/emil/candidaemia/bam                    # original BAMs
OUT=/mnt/data/emil/candidaemia/bam_candida           # mapped BAMs from the mapping script
REF=/mnt/data/emil/candidaemia/candida_refs/candida_species_ref.fna
MINQ="${1:-20}"
JOBS=16

# 1. contig -> species table, taken from the .fna headers (first word = contig,
#    next two words = species, e.g. "Candida albicans"). Check this file by eye!
MAP="$OUT/contig2species.tsv"
grep '^>' "$REF" | sed 's/^>//' | awk -F' ' -v OFS='\t' '{print $1, $2" "$3}' > "$MAP"

count() {
  s="$1"
  bam="$OUT/${s}.candida.bam"

  # reads that went into the Candida mapping (each unmapped read once)
  n_unmapped=$(samtools view -@ 2 -c -F 0x900 "$bam")
  # all reads in the original BAM (each read once)
  n_total=$(samtools view -@ 2 -c -F 0x900 "$IN/${s}.bam")

  # primary, mapped, MAPQ >= MINQ -> count per species
  samtools view -F 0x904 -q "$MINQ" "$bam" \
    | awk -F'\t' -v OFS='\t' -v s="$s" -v u="$n_unmapped" -v t="$n_total" '
        FNR==NR { sp[$1]=$2; next }          # first file: contig -> species
        { k = ($3 in sp) ? sp[$3] : $3; n[k]++; tot++ }
        END {
          for (k in n)
            print s, k, n[k], tot, u, t, n[k]/tot, n[k]/u, n[k]/t
          if (tot==0) print s, "none", 0, 0, u, t, 0, 0, 0
        }' "$MAP" - > "$OUT/${s}.species_counts.tsv"
}
export -f count
export IN OUT MAP MINQ

ls "$OUT"/*.candida.bam | xargs -n1 basename | sed 's/\.candida\.bam$//' \
  | xargs -P "$JOBS" -I{} bash -c 'count "$@"' _ {}

{
  printf "sample\tspecies\treads\tcandida_reads\tunmapped_input_reads\ttotal_reads\tfrac_of_candida\tfrac_of_unmapped\tfrac_of_total\n"
  cat "$OUT"/*.species_counts.tsv
} > "$OUT/species_fractions.tsv"

echo "Wrote $OUT/species_fractions.tsv (MAPQ >= $MINQ)"
