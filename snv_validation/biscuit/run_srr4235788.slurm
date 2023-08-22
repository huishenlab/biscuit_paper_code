#!/bin/bash
#PBS --job-name srr4235788
#PBS -o run_srr4235788.log
#PBS --ntasks=30
#PBS --mem=200000
#PBS --time=72:00:00

WORKDIR=/path/to/working/directory
TOOLS=/path/to/installed_packages
BISCUIT=${TOOLS}/biscuit-release-1.1.0.20220707
FASTQ=${WORKDIR}/trim
REF=${WORKDIR}/refs/biscuit/hg38
DUP=${TOOLS}/dupsifter-1.0.0.20220804

# Change to working directory
cd ${WORKDIR}/snv_validation/biscuit

# biscuit align
${BISCUIT}/biscuit align -@ 30 -b 1 \
    -R '@RG\tLB:ENCODE\tID:SRR4235788\tPL:Illumina\tPU:HiSeq2500\tSM:GM12878' \
    ${REF}/hg38_noContig.fa \
    ${FASTQ}/SRR4235788_1_val_1_500M.fq.gz \
    ${FASTQ}/SRR4235788_2_val_2_500M.fq.gz | \
${DUP}/dupsifter --add-mate-tags ${REF}/hg38_noContig.fa | \
samtools sort -@ 8 -m 5G -o SRR4235788.sorted.markdup.bam -O BAM -
samtools index SRR4235788.sorted.markdup.bam

# biscuit pileup
${BISCUIT}/biscuit pileup -@ 30 \
    ${REF}/hg38_noContig.fa \
    SRR4235788.sorted.markdup.bam \
    -o SRR4235788.pileup.vcf
bgzip -@ 20 SRR4235788.pileup.vcf
tabix -p vcf SRR4235788.pileup.vcf.gz

# biscuit snp vcf2bed
${BISCUIT}/biscuit vcf2bed -t snp \
    SRR4235788.pileup.vcf.gz | \
bgzip -@ 20 > SRR4235788.snp.sorted.bed.gz
tabix -p bed SRR4235788.snp.sorted.bed.gz
