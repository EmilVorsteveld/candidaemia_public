#!/usr/bin/env bash
#SBATCH --job-name=haplotag
#SBATCH --array=2-37
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/haplotag_%A_%a.log
#
# Allele-specific expression, step 1: tag each Kinnex read with its haplotype (H1/H2) using the patient's phased WGS VCF.
# Downstream analysis: 08_allele_specific.Rmd.
#
# Usage:
#   bash ase_haplotag.sh check ALEX1_D0     # input checks for one sample (contigs, index, phasing, read names)
#   bash ase_haplotag.sh run ALEX1_D0       # filter VCF + haplotag + summary for one sample
#   bash ase_haplotag.sh all 12             # all samples on one machine with GNU parallel (12 at a time)
#   sbatch ase_haplotag.sh                  # or: all samples as a SLURM array (line 1 of the sheet is the header)
#
# Sample sheet: sample_match.txt (id <tab> matching_DNA), 36 Kinnex samples. Samples with matching_DNA = NA
# (ATT45, IPPO4, POLY6) have no WGS and are skipped. The file has Windows line endings; these are stripped.
set -euo pipefail

# ---- Paths: edit these ----
REF=/mnt/data/emil/transcriptome_data/hg38.fa # same FASTA as used for the Kinnex genome BAMs and the WGS
BAM_DIR=/mnt/data/emil/candidaemia/bam           # expects $BAM_DIR/<rna_id>.bam
VCF_DIR=/mnt/data/emil/candidaemia/vcf_wgs                     # expects $VCF_DIR/<dna_id>.vcf.gz
OUT=/mnt/data/emil/candidaemia/ase
SHEET=/mnt/data/emil/candidaemia/sample_match.txt

# ---- Environment: bcftools/samtools from base, whatshap stacked on top ----
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate base
conda activate --stack whatshap
for tool in bcftools samtools whatshap; do command -v $tool >/dev/null || { echo "$tool not found" >&2; exit 1; }; done

mkdir -p "$OUT/vcf" "$OUT/haplotag" logs

# ---- All samples on one machine (GNU parallel) ----
MODE=${1:-run}
case "$MODE" in
  check|filter|run|all) ;;
  *) set -- run "$@"; MODE=run ;;   # "bash ase_haplotag.sh ALEX1_D0" = run on that sample
esac
if [ "$MODE" = "all" ]; then
  JOBS=${2:-12}
  SELF=$(realpath "$0")
  # 1) filter each patient's VCF once (D0 and D5-7 share it), 2) haplotag every sample with WGS
  awk 'NR > 1 {sub(/\r$/, ""); if ($2 != "NA" && !seen[$2]++) print $1}' "$SHEET" | parallel -j "$JOBS" "bash $SELF filter {}"
  awk 'NR > 1 {sub(/\r$/, ""); if ($2 != "NA") print $1}' "$SHEET" | parallel -j "$JOBS" --joblog logs/haplotag.joblog "bash $SELF run {} > logs/haplotag_{}.log 2>&1"
  cat "$OUT"/haplotag/*.summary.tsv > "$OUT/haplotag/all_samples.summary.tsv"
  exit 0
fi

# ---- Which sample ----
sheet_line() { awk -v n="$1" 'NR == n {sub(/\r$/, ""); print}' "$SHEET"; }   # handles CRLF and a missing final newline
if [ -n "${SLURM_ARRAY_TASK_ID:-}" ]; then
  read -r RNA DNA <<< "$(sheet_line "$SLURM_ARRAY_TASK_ID")"
else
  RNA=${2:?give a Kinnex sample id, e.g. ALEX1_D0}
  DNA=$(awk -v r="$RNA" '{sub(/\r$/, "")} $1 == r {print $2}' "$SHEET")
  [ -n "$DNA" ] || { echo "$RNA not in $SHEET" >&2; exit 1; }
fi
if [ "$DNA" = "NA" ]; then echo "$RNA has no matching WGS, skipped"; exit 0; fi
BAM=$BAM_DIR/$RNA.bam
VCF=$VCF_DIR/$DNA.vcf.gz
HET=$OUT/vcf/$DNA.het_phased.vcf.gz
echo "== $RNA <-> $DNA"

# ---- Index inputs if needed ----
[ -f "$BAM.bai" ] || [ -f "${BAM%.bam}.bai" ] || [ -f "$BAM.csi" ] || samtools index "$BAM"
[ -f "$VCF.tbi" ] || [ -f "$VCF.csi" ] || bcftools index -t "$VCF"
SAMPLE=$(bcftools query -l "$VCF" | head -1)

# ---- Checks ----
if [ "$MODE" = "check" ]; then
  echo "-- Contigs (BAM vs VCF vs REF): names must match"
  samtools view -H "$BAM" | awk '$1 == "@SQ" {print $2}' | head -3
  bcftools view -h "$VCF" | grep '^##contig' | head -3
  [ -f "$REF.fai" ] || samtools faidx "$REF"
  cut -f1 "$REF.fai" | head -3
  echo "-- BAM contigs missing from REF (should be none):"
  comm -23 <(samtools view -H "$BAM" | awk '$1 == "@SQ" {sub("SN:", "", $2); print $2}' | sort) <(cut -f1 "$REF.fai" | sort) | head
  echo "-- Aligner used for the BAM:"
  samtools view -H "$BAM" | grep '^@PG' | head -2 | cut -c1-200
  echo "-- VCF sample name: $SAMPLE"
  echo "-- PASS heterozygous SNV genotypes (most should be 0|1 or 1|0)"
  bcftools view -f PASS -v snps -g het "$VCF" | bcftools query -f '[%GT]\n' | sort | uniq -c
  echo "-- Phase blocks: number and median span (bp)"
  bcftools query -i 'FMT/PS!="."' -f '%CHROM\t%POS[\t%PS]\n' "$VCF" |
    awk '{k=$1":"$3; if(!(k in s)) s[k]=$2; e[k]=$2} END{for(k in s) print e[k]-s[k]}' | sort -n |
    awk '{a[NR]=$1} END{print "blocks:", NR, " median span:", a[int((NR+1)/2)]}'
  echo "-- First read names (compare with the reads quantified by oarfish)"
  samtools view "$BAM" | head -3 | cut -f1
  exit 0
fi

# ---- Filter WGS VCF to PASS, biallelic, phased heterozygous SNVs (once per patient) ----
if [ ! -f "$HET.tbi" ]; then
  bcftools view -f PASS -m2 -M2 -v snps -i 'GT="het" && FMT/PS!="."' "$VCF" -Oz -o "$HET.tmp.vcf.gz"
  mv "$HET.tmp.vcf.gz" "$HET"
  bcftools index -t "$HET"
fi
echo "-- Phased het SNVs: $(bcftools view -H "$HET" | wc -l)"
[ "$MODE" = "filter" ] && exit 0

# ---- Haplotag ----
whatshap haplotag \
  --reference "$REF" \
  --ignore-read-groups --sample "$SAMPLE" \
  --skip-missing-contigs \
  --output-haplotag-list "$OUT/haplotag/$RNA.haplotag.tsv.gz" \
  -o "$OUT/haplotag/$RNA.haplotagged.bam" \
  "$HET" "$BAM"
samtools index "$OUT/haplotag/$RNA.haplotagged.bam"

# ---- Summary: reads per haplotype (expect ~30-60% tagged, H1 ~ H2) ----
zcat "$OUT/haplotag/$RNA.haplotag.tsv.gz" | awk -v s="$RNA" '
  !/^#/ {n[$2]++; t++}
  END {for (k in n) printf "%s\t%s\t%d\t%.1f%%\n", s, k, n[k], 100*n[k]/t}' | tee "$OUT/haplotag/$RNA.summary.tsv"
