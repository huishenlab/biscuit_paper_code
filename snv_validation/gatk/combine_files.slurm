#!/bin/bash
#SBATCH --job-name combine_files
#SBATCH -o combine_files.log
#SBATCH --ntasks=8
#SBATCH --mem=75000
#SBATCH --time=12:00:00

WORKDIR=/path/to/working/directory
TOPDIR=${WORKDIR}/snv_validation/gatk

cd ${TOPDIR}

pigz -p 8 -d -c U*_R1_*.gz | pigz -c > merged_giab_illumina_R1.fastq.gz
pigz -p 8 -d -c U*_R2_*.gz | pigz -c > merged_giab_illumina_R2.fastq.gz
