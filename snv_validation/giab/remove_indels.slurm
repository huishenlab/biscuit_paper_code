#!/bin/bash
#SBATCH --job-name remove_indels
#SBATCH -o remove_indels.log
#SBATCH --time=5:00:00
#SBATCH --ntasks=5

#Change to working directory
WORKDIR=/path/to/working/directory
DIRLOC=${WORKDIR}/snv_validation/giab

cd ${DIRLOC}

# Genome-in-a-Bottle downloaded from NIST and comes several directory levels deep in the downloaded data
VCFDIR=${DIRLOC}/ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh38/supplementaryFiles

vcftools \
    --gzvcf ${VCFDIR}/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf.vcf.gz \
    --remove-indels \
    --recode \
    --recode-INFO-all \
    --stdout > HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_indels_removed.vcf
bgzip HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_indels_removed.vcf
tabix -p vcf HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_indels_removed.vcf.gz
