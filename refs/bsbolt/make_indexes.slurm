#!/bin/bash
#SBATCH --job-name indexing
#SBATCH -o make_indexes.log
#SBATCH --ntasks=4
#SBATCH --mem=50000
#SBATCH --time=48:00:00

# Useful variables
WORKDIR=/path/to/working/directory
REFDIRS=${WORKDIR}/refs/bsbolt

TOOLS=/path/to/installed_packages
BSBOLT=${TOOLS}/miniconda3/bin

# HG38
cd ${REFDIRS}/hg38
${BSBOLT}/bsbolt Index -G hg38_noContig.fa -DB ${REFDIRS}/hg38

# MM10
cd ${REFDIRS}/mm10
${BSBOLT}/bsbolt Index -G mm10_noContig.fa -DB ${REFDIRS}/mm10

# Z11
cd ${REFDIRS}/z11
${BSBOLT}/bsbolt Index -G danRer11.fa -DB ${REFDIRS}/z11
