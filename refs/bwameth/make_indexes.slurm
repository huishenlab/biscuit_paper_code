#!/bin/bash
#SBATCH --job-name indexing
#SBATCH -o make_indexes.log
#SBATCH --ntasks=4
#SBATCH --mem=50000
#SBATCH --time=48:00:00

# Useful variables
WORKDIR=/path/to/working/directory
REFDIRS=${WORKDIR}/refs/bwameth

TOOLS=/path/to/installed_packages
BWAMETH=${TOOLS}/bwa_meth_git # bwa-meth version 0.2.6

# HG38
cd ${REFDIRS}/hg38
${BWAMETH}/bwameth.py index hg38_noContig.fa

# MM10
cd ${REFDIRS}/mm10
${BWAMETH}/bwameth.py index mm10_noContig.fa

# Z11
cd ${REFDIRS}/z11
${BWAMETH}/bwameth.py index danRer11.fa
