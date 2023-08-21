#!/bin/bash
#SBATCH --job-name indexing
#SBATCH -o make_indexes.log
#SBATCH --ntasks=4
#SBATCH --mem=100000
#SBATCH --time=48:00:00

# Useful variables
WORKDIR=/path/to/working/directory
REFDIRS=${WORKDIR}/refs/biscuit

module load bbc2/bismark/bismark-0.24.0

# HG38
cd ${REFDIRS}/hg38
bismark_genome_preparation .

# MM10
cd ${REFDIRS}/mm10
bismark_genome_preparation .

# Z11
cd ${REFDIRS}/z11
bismark_genome_preparation .
