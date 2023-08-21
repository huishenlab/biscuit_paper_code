#!/bin/bash
#SBATCH --job-name download_data
#SBATCH -o retrieve_datasets.log
#SBATCH --ntasks=20
#SBATCH --mem=100000
#SBATCH --time=2-00:00:00

set -euo pipefail

# Modules to load
bbc2/sratoolkit/sratoolkit-3.0.5

# Working directory
WORKDIR=/path/to/working/directory
cd ${WORKDIR}

# TCGA samples
#------------------------------------------------------------------------------
# NOT SHOWN BUT AVAILABLE THROUGH LEGACY GDC SITE
# Samples used:
#     TCGA_COAD_A00R
#     TCGA_LUSC_2600
#------------------------------------------------------------------------------


# EM-seq
#------------------------------------------------------------------------------
# Download SRA files
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra8/SRR/012770/SRR13076816
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra67/SRR/012770/SRR13076817

# Turn SRA files into FASTQ files
fastq-dump --split-3 SRR13076816
fastq-dump --split-3 SRR13076817

# Compress newly generated FASTQ files
pigz -p 20 SRR13076816_1.fastq
pigz -p 20 SRR13076816_2.fastq
pigz -p 20 SRR13076817_1.fastq
pigz -p 20 SRR13076817_2.fastq

# Clean up SRA files
rm -f SRR13076816 SRR13076817
#------------------------------------------------------------------------------

# TruSeq Methyl Capture EPIC
#------------------------------------------------------------------------------
# Download SRA files
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra62/SRR/012745/SRR13051114
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra39/SRR/012745/SRR13051115
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra47/SRR/012745/SRR13051116
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra38/SRR/012745/SRR13051117
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra65/SRR/012745/SRR13051118
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra14/SRR/012745/SRR13051120
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra74/SRR/012745/SRR13051121
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra9/SRR/012745/SRR13051122
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra52/SRR/012745/SRR13051123
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra72/SRR/012745/SRR13051124

# Turn SRA files into FASTQ files
fastq-dump --split-3 SRR13051114
fastq-dump --split-3 SRR13051115
fastq-dump --split-3 SRR13051116
fastq-dump --split-3 SRR13051117
fastq-dump --split-3 SRR13051118
fastq-dump --split-3 SRR13051120
fastq-dump --split-3 SRR13051121
fastq-dump --split-3 SRR13051122
fastq-dump --split-3 SRR13051123
fastq-dump --split-3 SRR13051124

# Compress newly generated FASTQ files
pigz -p 20 SRR13051114_1.fastq
pigz -p 20 SRR13051114_2.fastq
pigz -p 20 SRR13051115_1.fastq
pigz -p 20 SRR13051115_2.fastq
pigz -p 20 SRR13051116_1.fastq
pigz -p 20 SRR13051116_2.fastq
pigz -p 20 SRR13051117_1.fastq
pigz -p 20 SRR13051117_2.fastq
pigz -p 20 SRR13051118_1.fastq
pigz -p 20 SRR13051118_2.fastq
pigz -p 20 SRR13051120_1.fastq
pigz -p 20 SRR13051120_2.fastq
pigz -p 20 SRR13051121_1.fastq
pigz -p 20 SRR13051121_2.fastq
pigz -p 20 SRR13051122_1.fastq
pigz -p 20 SRR13051122_2.fastq
pigz -p 20 SRR13051123_1.fastq
pigz -p 20 SRR13051123_2.fastq
pigz -p 20 SRR13051124_1.fastq
pigz -p 20 SRR13051124_2.fastq

# Clean up SRA files
rm -f SRR13051114 SRR13051115 SRR13051116 SRR13051117 SRR13051118
rm -f SRR13051120 SRR13051121 SRR13051122 SRR13051123 SRR13051124
#------------------------------------------------------------------------------

# GM12878
#------------------------------------------------------------------------------
# Download SRA files
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR4235788/SRR4235788.1
wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-9/SRR4235789/SRR4235789.1

# Turn SRA files into FASTQ files
fastq-dump --split-3 SRR4235788.1
fastq-dump --split-3 SRR4235789.1

# Compress newly generated FASTQ files
pigz -p 20 SRR4235788.1_1.fastq
pigz -p 20 SRR4235788.1_2.fastq
pigz -p 20 SRR4235789.1_1.fastq
pigz -p 20 SRR4235789.1_2.fastq

# Clean up SRA files
rm -f SRR4235788.1 SRR4235789.1
#------------------------------------------------------------------------------

# Mouse samples
#------------------------------------------------------------------------------
# Download SRA files
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra57/SRR/013166/SRR13482506
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra48/SRR/013166/SRR13482508

# Turn SRA files into FASTQ files
fastq-dump --split-3 SRR13482506
fastq-dump --split-3 SRR13482508

# Compress newly generated FASTQ files
pigz -p 20 SRR13482506_1.fastq
pigz -p 20 SRR13482506_2.fastq
pigz -p 20 SRR13482508_1.fastq
pigz -p 20 SRR13482508_2.fastq

# Clean up SRA files
rm -f SRR13482506 SRR13482508
#------------------------------------------------------------------------------

# Zebrafish samples
#------------------------------------------------------------------------------
# Download SRA files (EM-seq)
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra75/SRR/011342/SRR11614917
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra57/SRR/011342/SRR11614918
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra51/SRR/011342/SRR11614919
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra67/SRR/011342/SRR11614920

# Turn SRA files into FASTQ files (EM-seq)
fastq-dump --split-3 SRR11614917
fastq-dump --split-3 SRR11614918
fastq-dump --split-3 SRR11614919
fastq-dump --split-3 SRR11614920

# Compress newly generated FASTQ files (EM-seq)
pigz -c -p 20 \
    SRR11614917_1.fastq \
    SRR11614918_1.fastq \
    SRR11614919_1.fastq \
    SRR11614920_1.fastq > \
    SRR11614917_SRR11614920_1.fastq.gz
pigz -c -p 20 \
    SRR11614917_2.fastq \
    SRR11614918_2.fastq \
    SRR11614919_2.fastq \
    SRR11614920_2.fastq > \
    SRR11614917_SRR11614920_2.fastq.gz

# Clean up SRA files (EM-seq)
rm -f SRR11614917 SRR11614917_{1,2}.fastq
rm -f SRR11614918 SRR11614918_{1,2}.fastq
rm -f SRR11614919 SRR11614919_{1,2}.fastq
rm -f SRR11614920 SRR11614920_{1,2}.fastq

# Download SRA files (WGBS set 1)
wget https://sra-download.ncbi.nlm.nih.gov/traces/era7/ERR/ERR5005/ERR5005148
wget https://sra-download.ncbi.nlm.nih.gov/traces/era2/ERR/ERR5005/ERR5005149
wget https://sra-download.ncbi.nlm.nih.gov/traces/era2/ERR/ERR5005/ERR5005150
wget https://sra-download.ncbi.nlm.nih.gov/traces/era2/ERR/ERR5005/ERR5005151

# Turn SRA files into FASTQ files (WGBS set 1)
fastq-dump --split-3 ERR5005148
fastq-dump --split-3 ERR5005149
fastq-dump --split-3 ERR5005150
fastq-dump --split-3 ERR5005151

# Compress newly generated FASTQ files (WGBS set 1)
pigz -c -p 20 \
    ERR5005148_1.fastq \
    ERR5005149_1.fastq \
    ERR5005150_1.fastq \
    ERR5005151_1.fastq > \
    ERR5005148_ERR5005151_1.fastq.gz
pigz -c -p 20 \
    ERR5005148_2.fastq \
    ERR5005149_2.fastq \
    ERR5005150_2.fastq \
    ERR5005151_2.fastq > \
    ERR5005148_ERR5005151_2.fastq.gz

# Clean up SRA files (WGBS set 1)
rm -f ERR5005148 ERR5005148_{1,2}.fastq
rm -f ERR5005149 ERR5005149_{1,2}.fastq
rm -f ERR5005150 ERR5005150_{1,2}.fastq
rm -f ERR5005151 ERR5005151_{1,2}.fastq

# Download SRA files (WGBS set 2)
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra57/SRR/012497/SRR12797058
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra1/SRR/012497/SRR12797059
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra60/SRR/012497/SRR12797060
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra76/SRR/012497/SRR12797061
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra33/SRR/012497/SRR12797062
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra51/SRR/012497/SRR12797063
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra6/SRR/012497/SRR12797064
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra45/SRR/012497/SRR12797065

# Turn SRA files into FASTQ files (WGBS set 2)
fastq-dump --split-3 SRR12797058
fastq-dump --split-3 SRR12797059
fastq-dump --split-3 SRR12797060
fastq-dump --split-3 SRR12797061
fastq-dump --split-3 SRR12797062
fastq-dump --split-3 SRR12797063
fastq-dump --split-3 SRR12797064
fastq-dump --split-3 SRR12797065

# Compress newly generated FASTQ files (WGBS set 2)
pigz -c -p 20 \
    SRR12797058_1.fastq \
    SRR12797059_1.fastq \
    SRR12797060_1.fastq \
    SRR12797061_1.fastq \
    SRR12797062_1.fastq \
    SRR12797063_1.fastq \
    SRR12797064_1.fastq \
    SRR12797065_1.fastq > \
    SRR12797058_SRR12797065_1.fastq.gz
pigz -c -p 20 \
    SRR12797058_2.fastq \
    SRR12797059_2.fastq \
    SRR12797060_2.fastq \
    SRR12797061_2.fastq \
    SRR12797062_2.fastq \
    SRR12797063_2.fastq \
    SRR12797064_2.fastq \
    SRR12797065_2.fastq > \
    SRR12797058_SRR12797065_2.fastq.gz

# Clean up SRA files (WGBS set 2)
rm -f SRR12797058 SRR12797058_{1,2}.fastq
rm -f SRR12797059 SRR12797059_{1,2}.fastq
rm -f SRR12797060 SRR12797060_{1,2}.fastq
rm -f SRR12797061 SRR12797061_{1,2}.fastq
rm -f SRR12797062 SRR12797062_{1,2}.fastq
rm -f SRR12797063 SRR12797063_{1,2}.fastq
rm -f SRR12797064 SRR12797064_{1,2}.fastq
rm -f SRR12797065 SRR12797065_{1,2}.fastq
#------------------------------------------------------------------------------
