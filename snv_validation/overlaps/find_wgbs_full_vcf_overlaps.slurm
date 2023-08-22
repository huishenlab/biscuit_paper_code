#!/bin/bash
#SBATCH --job-name wgbs_overlaps
#SBATCH -o find_wgbs_full_vcf_overlaps.log
#SBATCH --ntasks=10
#SBATCH --mem=50000
#SBATCH --time=12:00:00

# Create directory location variable
WORKDIR=/path/to/working/directory

#Change to working directory
cd ${WORKDIR}/snv_validation/overlaps/wgbs_overlaps

bcftools isec \
    ${WORKLOC}/biscuit/SRR4235788.pileup.vcf.gz \
    ${WORKLOC}/biscuit/SRR4235789.pileup.vcf.gz \
    -p full_vcf_overlap

bgzip -@ 10 full_vcf_overlap/0000.vcf
bgzip -@ 10 full_vcf_overlap/0001.vcf
bgzip -@ 10 full_vcf_overlap/0002.vcf
bgzip -@ 10 full_vcf_overlap/0003.vcf
tabix -p vcf full_vcf_overlap/0000.vcf.gz
tabix -p vcf full_vcf_overlap/0001.vcf.gz
tabix -p vcf full_vcf_overlap/0002.vcf.gz
tabix -p vcf full_vcf_overlap/0003.vcf.gz
