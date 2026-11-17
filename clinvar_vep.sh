#!/usr/bin/env bash
# ClinVar variants annotated with VEP (GTF mode) against GENCODE v49 and against the Kinnex novel isoforms.
# VEP consequences are per transcript, so "GENCODE + Kinnex" = flagged in the GENCODE run OR in the Kinnex run;
# the Kinnex run only needs the novel isoforms, which keeps it fast. Analysis: 05_clinvar_vep.Rmd.
set -euo pipefail

cd "$(dirname "$0")"
SHARED=~/work/shared_data
VEP=~/ensembl-vep/vep
FASTA=$SHARED/hg38/hg38.fa
OUT=data/kinnex/clinvar
mkdir -p "$OUT" "$SHARED/clinvar"

# 1. Latest ClinVar (GRCh38)
if [ ! -f "$SHARED/clinvar/clinvar.vcf.gz" ]; then
  curl -sSf -o "$SHARED/clinvar/clinvar.vcf.gz" https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
  curl -sSf -o "$SHARED/clinvar/clinvar.vcf.gz.tbi" https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
fi
bcftools view -h "$SHARED/clinvar/clinvar.vcf.gz" | grep "^##fileDate"

# 2. Expressed-gene BED and Kinnex novel-isoform GTF
Rscript - <<'EOF'
suppressPackageStartupMessages({library(tidyverse); library(DESeq2)})
out <- here::here("data/kinnex/clinvar")

# GENCODE v49 spans of the genes in the gene-level DESeq2 object
dds <- readRDS(here::here("output/kinnex/dds.rda"))
genes <- data.table::fread(cmd = "awk -F'\\t' '$3 == \"gene\"' ~/work/shared_data/gencode/gencode.v49.annotation.gtf",
                           header = FALSE, sep = "\t", quote = "") %>%
  transmute(chrom = V1, start = V4 - 1L, end = V5, gene_id = sub('.*gene_id "([^"]+)".*', "\\1", V9)) %>%
  subset(gene_id %in% rownames(dds)) %>%
  arrange(chrom, start)
write_tsv(genes, file.path(out, "expressed_genes.bed"), col_names = FALSE)
cat("Expressed genes with GENCODE coordinates:", nrow(genes), "of", nrow(dds), "\n")

# Expressed NIC/NNC isoforms from the ORFannotate GTF (has CDS)
annotated <- readRDS(here::here("data/kinnex/sqanti/annotated.rda"))
dds_tx <- readRDS(here::here("output/kinnex/dds_tx.rda"))
novel_ids <- subset(annotated, structural_category %in% c("novel_in_catalog", "novel_not_in_catalog") &
                      transcript_id %in% rownames(dds_tx))$transcript_id
gtf <- data.table::fread(here::here("data/kinnex/sqanti/filter_out/ORFannotate_annotated.gtf"),
                         header = FALSE, sep = "\t", quote = "")
gtf$tx <- sub('.*transcript_id "([^"]+)".*', "\\1", gtf$V9)
gtf <- gtf[gtf$tx %in% novel_ids & gtf$V3 %in% c("transcript", "exon", "CDS", "five_prime_utr", "three_prime_utr"), ]
cat("Novel expressed isoforms:", length(novel_ids), "- found in ORFannotate GTF:", length(unique(gtf$tx)), "\n")

coding <- unique(gtf$tx[gtf$V3 == "CDS"])

# ORFannotate leaves the CDS frame as "."; VEP needs it: bases to skip = (3 - CDS length so far %% 3) %% 3, in transcript order
cds <- which(gtf$V3 == "CDS")
cds_order <- cds[order(gtf$tx[cds], ifelse(gtf$V7[cds] == "+", gtf$V4[cds], -gtf$V5[cds]))]
cds_len <- gtf$V5[cds_order] - gtf$V4[cds_order] + 1
len_before <- ave(cds_len, gtf$tx[cds_order], FUN = function(x) cumsum(x) - x)
gtf$V8[cds_order] <- as.character((3 - len_before %% 3) %% 3)
gene_names <- setNames(annotated$gene_name, annotated$transcript_id)
biotype <- ifelse(gtf$tx %in% coding, "protein_coding", "lncRNA")
# VEP's CSQ output is |-delimited, so | in IDs is replaced
gtf$V9 <- paste0(gsub("\\|", "_", gtf$V9),
                 ' gene_name "', gene_names[gtf$tx], '"; gene_type "', biotype, '"; transcript_type "', biotype, '";')
gtf$V2 <- "Kinnex"
gtf$tx <- NULL

# VEP needs a parent gene line per gene_id: span of its isoforms
gene_lines <- gtf[gtf$V3 == "transcript", ] %>%
  mutate(gene_id = sub('.*gene_id "([^"]+)".*', "\\1", V9),
         gene_name = sub('.*gene_name "([^"]*)".*', "\\1", V9),
         coding = grepl('transcript_type "protein_coding"', V9)) %>%
  group_by(V1, V7, gene_id) %>%
  summarise(V4 = min(V4), V5 = max(V5), gene_name = dplyr::first(gene_name), coding = any(coding), .groups = "drop") %>%
  transmute(V1, V2 = "Kinnex", V3 = "gene", V4, V5, V6 = ".", V7, V8 = ".",
            V9 = paste0('gene_id "', gene_id, '"; gene_name "', gene_name, '"; gene_type "',
                        ifelse(coding, "protein_coding", "lncRNA"), '";'))
gtf <- bind_rows(gene_lines, gtf)
data.table::fwrite(gtf, file.path(out, "kinnex_novel.gtf"), sep = "\t", col.names = FALSE, quote = FALSE)
EOF

# 3. ClinVar P/LP and B/LB in expressed genes, chr-prefixed
# Review status (>= 1 star) is filtered in 05_clinvar_vep.Rmd: its values contain commas, which bcftools splits on
SEL="INFO/CLNSIG=\"Pathogenic\" || INFO/CLNSIG=\"Likely_pathogenic\" || INFO/CLNSIG=\"Pathogenic/Likely_pathogenic\" || INFO/CLNSIG=\"Benign\" || INFO/CLNSIG=\"Likely_benign\" || INFO/CLNSIG=\"Benign/Likely_benign\""
echo -e "MT chrM\n$(for c in $(seq 1 22) X Y; do echo "$c chr$c"; done)" | tr ' ' '\t' > "$OUT/chr_rename.txt"
bcftools annotate --rename-chrs "$OUT/chr_rename.txt" "$SHARED/clinvar/clinvar.vcf.gz" -Ou |
  bcftools view -T "$OUT/expressed_genes.bed" -Ou |
  bcftools view -i "$SEL" -Ou |
  bcftools norm -d exact -Oz -o "$OUT/clinvar_sel.vcf.gz"
tabix -f -p vcf "$OUT/clinvar_sel.vcf.gz"
echo "Selected variants per class:"
bcftools query -f '%INFO/CLNSIG\n' "$OUT/clinvar_sel.vcf.gz" | sort | uniq -c

# 4. Sorted, bgzipped, indexed GTFs
for name in gencode kinnex_novel; do
  if [ "$name" = gencode ]; then src=$SHARED/gencode/gencode.v49.annotation.gtf; else src=$OUT/kinnex_novel.gtf; fi
  grep -v '^#' "$src" | sort -t$'\t' -k1,1 -k4,4n -k5,5n | bgzip -c > "$OUT/$name.gtf.gz"
  tabix -f -p gff "$OUT/$name.gtf.gz"
done

# 5. VEP (GTF mode), one run per annotation (SKIP_VEP=1 stops after the prep)
[ "${SKIP_VEP:-0}" = 1 ] && exit 0
for name in gencode kinnex_novel; do
  "$VEP" -i "$OUT/clinvar_sel.vcf.gz" --format vcf --gtf "$OUT/$name.gtf.gz" --fasta "$FASTA" \
    --tab --fields "Uploaded_variation,Location,Allele,Gene,Feature,BIOTYPE,Consequence,IMPACT,EXON,INTRON,HGVSc,HGVSp" \
    --hgvs --fork 8 --no_stats --force_overwrite -o "$OUT/vep_$name.tsv"
  gzip -f "$OUT/vep_$name.tsv"
done
