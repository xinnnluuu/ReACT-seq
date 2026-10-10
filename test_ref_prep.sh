#!/usr/bin/env bash
# test_ref_prep.sh — build the ReACT-seq lightweight test references
# (hg38 chr21 subset, Ensembl-style contig names).
#
# The workflow's test_config.yaml points every reference at test_ref/. This
# script (re)builds them from public sources:
#   chr21.fa            samtools faidx on a local GRCh38 fasta, or Ensembl FTP
#   gene_plus.bed       chr21 gene bodies, split by strand (4-col: chr start
#   gene_minus.bed        end gene_id) — consumed by ReACT_seq_stoprate_transcrip.py
#   peaks_chr21.bed     candidate RT-stop regions; built from the workflow's
#                       own IP-vs-Input narrowPeaks once they exist (see step 4)
#   STAR_hg38_chr21     STAR index of chr21 (reduced genomeSAindexNbits)
#   bowtie2_rRNA_hg38   rRNA depletion index (45S U13369.1 + 5S NR_023363.1)
#
# Usage:
#   bash test_ref_prep.sh [GRCh38.fa] [Homo_sapiens.GRCh38.gtf.gz]
# Both arguments optional: omitted files are downloaded from the Ensembl FTP.
# Tools required: samtools, bowtie2-build, STAR, awk, curl.
set -eo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

GENOME=${1:-}
GTF=${2:-}
WORK=${TMPDIR:-/tmp}/react_test_ref.$$
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

# --- 1) chr21 fasta ------------------------------------------------------------
if [ ! -s test_ref/chr21.fa ]; then
  mkdir -p test_ref
  if [ -z "$GENOME" ] && [ ! -s "$WORK/GRCh38.fa" ]; then
    echo "downloading GRCh38 primary assembly (once) ..."
    curl -sS -m 3600 -o "$WORK/GRCh38.fa.gz" \
      'https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
    gunzip "$WORK/GRCh38.fa.gz"; GENOME=$WORK/GRCh38.fa
  fi
  samtools faidx "$GENOME" 21 > test_ref/chr21.fa
  samtools faidx test_ref/chr21.fa
fi
CHR=21

# --- 2) strand-split gene BEDs (chr21) ------------------------------------------
if [ ! -s test_ref/gene_plus.bed ] || [ ! -s test_ref/gene_minus.bed ]; then
  if [ -z "$GTF" ] && [ ! -s "$WORK/grch38.gtf.gz" ]; then
    curl -sS -m 600 -o "$WORK/grch38.gtf.gz" \
      'https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz'
    GTF=$WORK/grch38.gtf.gz
  fi
  for strand in plus minus; do
    if [ "$strand" = plus ]; then s=+; else s=-; fi
    zcat "$GTF" | awk -F'\t' -v OFS='\t' -v c="$CHR" -v s="$s" '
      $1==c && $3=="gene" && $7==s {
        match($9, /gene_id "[^"]+"/)
        print c, $4-1, $5, substr($9, RSTART+9, RLENGTH-10)
      }' | sort -k1,1 -k2,2n > "test_ref/gene_${strand}.bed"
  done
fi

# --- 3) aligner indexes -----------------------------------------------------------
if [ ! -f test_ref/STAR_hg38_chr21/Genome ]; then
  mkdir -p test_ref/STAR_hg38_chr21
  STAR --runMode genomeGenerate --runThreadN 4 \
       --genomeDir test_ref/STAR_hg38_chr21 --genomeFastaFiles test_ref/chr21.fa \
       --genomeSAindexNbases 11 > test_ref/STAR_hg38_chr21/STAR_chr21.log 2>&1
fi
if [ ! -f test_ref/bowtie2_rRNA_hg38/hg38_rRNA.1.bt2 ]; then
  mkdir -p test_ref/bowtie2_rRNA_hg38
  for acc in U13369.1 NR_023363.1; do
    curl -sS -m 120 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc%.1}&rettype=fasta&retmode=text"
  done > test_ref/bowtie2_rRNA_hg38/hg38_rRNA.fa
  bowtie2-build test_ref/bowtie2_rRNA_hg38/hg38_rRNA.fa \
                test_ref/bowtie2_rRNA_hg38/hg38_rRNA \
                2> test_ref/bowtie2_rRNA_hg38/bt2build.log
fi

# --- 4) candidate peak regions ----------------------------------------------------
# find_highstop_hg38.py screens RT-stop sites inside candidate regions; the
# method defines these as IP-enriched (peak-calling) regions. After the
# workflow's peak steps have run, collect the union of IP narrowPeaks:
if [ ! -s test_ref/peaks_chr21.bed ]; then
  if ls test_out/peaks/*_peaks.narrowPeak >/dev/null 2>&1; then
    cat test_out/peaks/*_peaks.narrowPeak | cut -f1-3 | sort -k1,1 -k2,2n \
      > test_ref/peaks_chr21.bed
  else
    echo "NOTE: test_out/peaks/*.narrowPeak not present yet — run the workflow first," >&2
    echo "then re-run this script to build peaks_chr21.bed (or point peaks_bed at a" >&2
    echo "candidate-region track of your own)." >&2
    : > test_ref/peaks_chr21.bed
  fi
fi

echo "=== ReACT-seq test_ref ready: $(date -Is) ==="
ls -la test_ref/
