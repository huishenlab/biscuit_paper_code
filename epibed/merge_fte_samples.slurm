#!/bin/bash
#SBATCH --job-name merge
#SBATCH -o merge_fte_samples.log
#SBATCH --ntasks=40
#SBATCH --mem=300000
#SBATCH --time=48:00:00

WORKDIR=/path/to/working/directory
cd ${WORKDIR}/epibed

samtools merge \
    -@ 20 \
    -o Aswift.bam \
    ${DATA}/FtubeAswift.sorted.markdup.bam \
    ${DATA}/FtubeAswiftRep2.sorted.markdup.bam

samtools merge \
    -@ 20 \
    -o Aneb.bam \
    ${DATA}/FtubeAneb.sorted.markdup.bam \
    ${DATA}/FtubeAnebRep2.sorted.markdup.bam

samtools index -@ 20 Aswift.bam
samtools index -@ 20 Aneb.bam
