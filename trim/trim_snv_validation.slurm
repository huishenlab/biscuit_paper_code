#!/bin/bash
#SBATCH --job-name trim_snv_validation
#SBATCH -o trim_snv_validation.log
#SBATCH --ntasks=20
#SBATCH --mem=180000
#SBATCH --time=24:00:00

WORKDIR=/path/to/working/directory
cd ${WORKDIR}/trim

trim_galore \
    --illumina \
    --cores 4 \
    --paired \
    ../data/SRR4235788_1.fastq.gz \
    ../data/SRR4235788_2.fastq.gz

trim_galore \
    --illumina \
    --cores 4 \
    --paired \
    ../data/SRR4235789_1.fastq.gz \
    ../data/SRR4235789_2.fastq.gz

zcat SRR4235788_1_val_1.fq.gz | head -n 2000000000 | pigz -c -p 20 > SRR4235788_1_val_1_500M.fq.gz
zcat SRR4235788_2_val_2.fq.gz | head -n 2000000000 | pigz -c -p 20 > SRR4235788_2_val_2_500M.fq.gz

zcat SRR4235789_1_val_1.fq.gz | head -n 2000000000 | pigz -c -p 20 > SRR4235789_1_val_1_500M.fq.gz
zcat SRR4235789_2_val_2.fq.gz | head -n 2000000000 | pigz -c -p 20 > SRR4235789_2_val_2_500M.fq.gz

rm SRR4235788_1_val_1.fq.gz SRR4235788_2_val_2.fq.gz
rm SRR4235789_1_val_1.fq.gz SRR4235789_2_val_2.fq.gz
