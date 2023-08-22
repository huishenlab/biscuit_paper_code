#!/bin/bash
#SBATCH --job-name downsample
#SBATCH -o downsample_merged_giab_illumina_fastqs.log
#SBATCH --ntasks=8
#SBATCH --time=48:00:00

# Downsample the GIAB Illumina reads to 500M
WORKDIR=/path/to/working/directory
cd ${WORKDIR}/snv_validation/gatk

TOOLS=/path/to/installed_packages
SEQTK=${TOOLS}/seqtk

${SEQTK}/seqtk sample -2 -s100 merged_giab_illumina_R1.fastq.gz 500000000 | pigz -c -p 8 > downsampled_merged_giab_illumina_R1.fastq.gz
${SEQTK}/seqtk sample -2 -s100 merged_giab_illumina_R2.fastq.gz 500000000 | pigz -c -p 8 > downsampled_merged_giab_illumina_R2.fastq.gz
