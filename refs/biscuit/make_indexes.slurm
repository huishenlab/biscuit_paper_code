#!/bin/bash
#SBATCH --job-name indexing
#SBATCH -o make_indexes.log
#SBATCH --ntasks=4
#SBATCH --mem=50000
#SBATCH --time=48:00:00

# Useful variables
WORKDIR=/path/to/working/directory
REFDIRS=${WORKDIR}/refs/biscuit

TOOLS=/path/to/installed_packages
BISCUIT=${TOOLS}/biscuit-release-1.2.1.20230601

# HG38
cd ${REFDIRS}/hg38
${BISCUIT}/biscuit index hg38_noContig.fa

# MM10
cd ${REFDIRS}/mm10
${BISCUIT}/biscuit index mm10_noContig.fa

# Z11
cd ${REFDIRS}/z11
${BISCUIT}/biscuit index danRer11.fa
