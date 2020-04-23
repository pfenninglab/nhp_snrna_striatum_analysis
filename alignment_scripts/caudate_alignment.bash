#!/bin/bash
cd ~/nhp_snrna_striatum_analysis/data/raw/caudate/L1
STAR --soloCBwhitelist ~/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/3M-february-2018.txt --soloType Droplet --readFilesIn ~/nhp_snrna_striatum_analysis/data/raw/caudate/Caudate_S2_L001_R2_001.fastq.gz ~/nhp_snrna_striatum_analysis/data/raw/caudate/Caudate_S2_L001_R1_001.fastq.gz --readFilesCommand zcat --soloUMIlen 12 --genomeDir ~/nhp_snrna_striatum_analysis/data/raw/rhemac10_star_custom_genome --runThreadN 12 --soloBarcodeReadLength 0

cd ~/nhp_snrna_striatum_analysis/data/raw/caudate/L2
STAR --soloCBwhitelist ~/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/3M-february-2018.txt --soloType Droplet --readFilesIn ~/nhp_snrna_striatum_analysis/data/raw/caudate/Caudate_S2_L002_R2_001.fastq.gz ~/nhp_snrna_striatum_analysis/data/raw/caudate/Caudate_S2_L002_R1_001.fastq.gz --readFilesCommand zcat --soloUMIlen 12 --genomeDir ~/nhp_snrna_striatum_analysis/data/raw/rhemac10_star_custom_genome --runThreadN 12 --soloBarcodeReadLength 0

cd ~/nhp_snrna_striatum_analysis/data/raw/caudate/L3
STAR --soloCBwhitelist ~/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/3M-february-2018.txt --soloType Droplet --readFilesIn ~/nhp_snrna_striatum_analysis/data/raw/caudate/Caudate_S2_L003_R2_001.fastq.gz ~/nhp_snrna_striatum_analysis/data/raw/caudate/Caudate_S2_L003_R1_001.fastq.gz --readFilesCommand zcat --soloUMIlen 12 --genomeDir ~/nhp_snrna_striatum_analysis/data/raw/rhemac10_star_custom_genome --runThreadN 12 --soloBarcodeReadLength 0

cd ~/nhp_snrna_striatum_analysis/data/raw/caudate/L4
STAR --soloCBwhitelist ~/cellranger-3.0.2/cellranger-cs/3.0.2/lib/python/cellranger/barcodes/3M-february-2018.txt --soloType Droplet --readFilesIn ~/nhp_snrna_striatum_analysis/data/raw/caudate/Caudate_S2_L004_R2_001.fastq.gz ~/nhp_snrna_striatum_analysis/data/raw/caudate/Caudate_S2_L004_R1_001.fastq.gz --readFilesCommand zcat --soloUMIlen 12 --genomeDir ~/nhp_snrna_striatum_analysis/data/raw/rhemac10_star_custom_genome --runThreadN 12 --soloBarcodeReadLength 0


