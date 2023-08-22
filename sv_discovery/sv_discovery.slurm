#!/bin/bash
#SBATCH --job-name sv_discovery
#SBATCH -o sv_discovery.log
#SBATCH --ntasks=30
#SBATCH --mem=200000
#SBATCH --time=120:00:00

WORKDIR=/path/to/working/directory
DIRLOC=${WORKDIR}/sv_discovery
REFLOC=${WORKDIR}/refs/biscuit/hg38
FASTQS=${WORKDIR}/data

TOOLS=/path/to/installed_packages
BISCUIT=${TOOLS}/biscuit-release-1.1.0.20220707
DUPSIFT=${TOOLS}/dupsifter-1.0.0.20220804
LUMPY=${TOOLS}/lumpy_sv_git/bin # lumpy version 0.2.13

cd ${DIRLOC}

# Launch biscuit
${BISCUIT}/biscuit align -@ 30 -b 1 -M \
    -R '@RG\tLB:hg38\tID:WGBS_PML_RARa\tPL:Illumina\tPU:hiseq2000\tSM:PML_RARa' \
    ${REFLOC}/hg38_noContig.fa \
    ${FASTQS}/SRR1800202_1.fastq.gz \
    ${FASTQS}/SRR1800202_2.fastq.gz | \
${DUPSIFT}/dupsifter ${REFLOC}/hg38_noContig.fa | \
samtools sort -@ 8 -m 5G -o PML_RARa.sorted.markdup.withdisc_split_clip.hg38.bam -O BAM -

# Index full BAM
samtools index -@ 15 PML_RARa.sorted.markdup.withdisc_split_clip.hg38.bam

# Generate pileup and vcf2bed
${BISCUIT}/biscuit pileup -@ 30 \
    ${REFLOC}/hg38_noContig.fa \
    PML_RARa.sorted.markdup.withdisc_split_clip.hg38.bam | \
bgzip > PML_RARa.pileup.vcf.gz
tabix -p vcf PML_RARa.pileup.vcf.gz

${BISCUIT}/biscuit vcf2bed PML_RARa.pileup.vcf.gz | \
${BISCUIT}/biscuit mergecg ${REFLOC}/hg38_noContig.fa - | \
bgzip > PML_RARa.cg.mergecg.bed.gz
tabix -p bed PML_RARa.cg.mergecg.bed.gz

# Run lumpy for SV calling
${LUMPY}/lumpyexpress \
    -B PML_RARa.sorted.markdup.withdisc_split_clip.hg38.bam \
    -o PML_RARa.lumpy.vcf
gzip PML_RARa.lumpy.vcf

# Run MANTA for SV calling
#
# Create Singularity container from Docker container:
#    singularity pull docker://varishenlab/sv_calling:sv_calling
#    mv sv_calling_sv_calling.sif sv_calling.sif
#
# Make sure to copy sv_calling.sif file and any resource (FASTAs, BEDs, etc) to working directory
#
# Commands that were run (After the first command the rest were run in Singularity)
#     singularity shell sv_calling.sif
#     cd ${WORKDIR}
#     /usr/local/bin/manta/bin/configManta.py \
#        --tumorBam PML_RARa.sorted.markdup.withdisc_split_clip.hg38.bam \
#        --callRegions sv_calling_resources/manta_call_regions_hg38_blklst_exclude_primary_assembly_no_chrM.bed.gz \
#        --runDir manta_sv_calls \
#        --referenceFasta sv_calling_resources/GRCh38.p13.genome.fa
#     cd manta_sv_calls/
#     python runWorkflow.py --quiet
#     exit
