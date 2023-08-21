#!/bin/bash
#SBATCH --job-name indexing
#SBATCH -o make_indexes.log
#SBATCH --ntasks=10
#SBATCH --mem=50000
#SBATCH --time=48:00:00

# Useful variables
WORKDIR=/path/to/working/directory
REFDIRS=${WORKDIR}/refs/gembs

TOOLS=/path/to/installed_packages
GEMBS=${TOOLS}/gemBS_rs_git

cd ${REFDIRS}

# HG38
${GEMBS}/gemBS prepare --config hg38.conf --text-metadata hg38.conf
${GEMBS}/gemBS index

# MM10
${GEMBS}/gemBS prepare --config mm10.conf --text-metadata mm10.conf
${GEMBS}/gemBS index

# Z11
${GEMBS}/gemBS prepare --config z11.conf --text-metadata z11.conf
${GEMBS}/gemBS index
