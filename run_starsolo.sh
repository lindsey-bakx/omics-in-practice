#!/bin/bash
# STARsolo alignment of BRB-seq data (human, GRCh38)
# Run this from: /scratch/Omics_in_practice/lindseybakx/projects/brbseq

STAR_BIN=$HOME/software/STAR-2.7.9a/bin/Linux_x86_64_static/STAR
GENOME_DIR=/scratch/Omics_in_practice/genomes/GRCh38/index/star
R1=fastq/Omics_BRBseq_MQ_UDI_1_R1.fastq.gz
R2=fastq/Omics_BRBseq_MQ_UDI_1_R2.fastq.gz
WHITELIST=barcodes.txt
OUT=results/starsolo/

mkdir -p "$OUT"

nice -n 9 "$STAR_BIN" \
  --runMode alignReads \
  --runThreadN 8 \
  --outBAMsortingThreadN 8 \
  --genomeDir "$GENOME_DIR" \
  --readFilesCommand zcat \
  --readFilesIn "$R2" "$R1" \
  --outFileNamePrefix "$OUT" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMunmapped Within \
  --outSAMmapqUnique 60 \
  --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
  --outFilterMultimapNmax 1 \
  --quantMode GeneCounts \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist "$WHITELIST" \
  --soloCBstart 1 \
  --soloCBlen 14 \
  --soloUMIstart 15 \
  --soloUMIlen 14 \
  --soloBarcodeReadLength 0 \
  --soloStrand Forward \
  --soloUMIdedup NoDedup 1MM_Directional \
  --soloCellFilter None \
  --soloFeatures Gene
