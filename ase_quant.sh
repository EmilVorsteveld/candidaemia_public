#!/usr/bin/env bash
#
# Allele-specific expression, step 2: re-quantify each haplotype's reads with oarfish against the filtered isoform set.
# Input: step 1 (ase_haplotag.sh) read -> haplotype lists, plus the transcriptome BAMs used for the main oarfish quant.
# Output: $OUT/quant/<RNA>.H1.quant and <RNA>.H2.quant, same format as data/kinnex/oarfish_out*/<SAMPLE>.quant.
# Downstream analysis: 08_allele_specific.Rmd.
#
# Usage:
#   bash ase_quant.sh ALEX1_D0      # one sample
#   bash ase_quant.sh all 8         # every sample that has a haplotag list, 8 at a time (GNU parallel)
set -euo pipefail

# ---- Paths: edit these ----
TX_BAM_DIR=/mnt/data/emil/candidaemia/bam_tx_2   # the BAMs behind the oarfish run used in 00_kinnex_prep (oarfish_out_v2 <- bam_tx_2)
OUT=/mnt/data/emil/candidaemia/ase
THREADS=8

for tool in samtools oarfish; do command -v $tool >/dev/null || { echo "$tool not found" >&2; exit 1; }; done
mkdir -p "$OUT/reads" "$OUT/quant" logs

if [ "${1:-}" = "all" ]; then
  SELF=$(realpath "$0")
  ls "$OUT"/haplotag/*.haplotag.tsv.gz | xargs -n1 basename | sed 's/\.haplotag\.tsv\.gz$//' |
    parallel -j "${2:-8}" --joblog logs/ase_quant.joblog "bash $SELF {} > logs/ase_quant_{}.log 2>&1"
  exit 0
fi

RNA=${1:?give a Kinnex sample id, e.g. ALEX1_D0}
TAGS=$OUT/haplotag/$RNA.haplotag.tsv.gz
TXBAM=$TX_BAM_DIR/$RNA.tx.bam
[ -f "$TAGS" ] || { echo "no haplotag list for $RNA" >&2; exit 1; }
[ -f "$TXBAM" ] || { echo "missing $TXBAM" >&2; exit 1; }

# ---- Read name lists per haplotype ----
# whatshap list columns: readname, haplotype (H1/H2/none), phaseset, chromosome.
# A read with supplementary alignments can appear more than once; drop reads tagged with both haplotypes.
zcat "$TAGS" | awk '!/^#/ && ($2 == "H1" || $2 == "H2") {print $1"\t"$2}' | sort -u |
  awk '{h[$1] = ($1 in h && h[$1] != $2) ? "conflict" : $2} END {for (r in h) print r"\t"h[r]}' > "$OUT/reads/$RNA.tsv"
for H in H1 H2; do
  awk -v h=$H '$2 == h {print $1}' "$OUT/reads/$RNA.tsv" > "$OUT/reads/$RNA.$H.txt"
done
echo "-- $RNA reads: H1 $(wc -l < "$OUT/reads/$RNA.H1.txt"), H2 $(wc -l < "$OUT/reads/$RNA.H2.txt"), conflicting $(awk '$2 == "conflict"' "$OUT/reads/$RNA.tsv" | wc -l)"

# ---- Subset the transcriptome BAM and run oarfish with the main-run settings ----
for H in H1 H2; do
  samtools view -b -@ 2 -N "$OUT/reads/$RNA.$H.txt" -o "$OUT/quant/$RNA.$H.tx.bam" "$TXBAM"
  oarfish -j "$THREADS" --alignments "$OUT/quant/$RNA.$H.tx.bam" \
    --seq-tech pac-bio-hifi --filter-group no-filters --model-coverage \
    -o "$OUT/quant/$RNA.$H"
  rm "$OUT/quant/$RNA.$H.tx.bam"
done
echo "-- done: $OUT/quant/$RNA.H{1,2}.quant"
