#!/bin/bash
#SBATCH -N make_epibed
#SBATCH -o make_epibed.log
#SBATCH --ntasks=40
#SBATCH --mem=150000
#SBATCH --time=4:00:00

WORKDIR=/path/to/working/directory
DATA=${WORKDIR}/epibed
FAST=${WORKDIR}/refs/biscuit/hg38

cd ${WORKDIR}

echo -e "chr15\t24954913\t24955009\tchr15_24000000" > locations.bed

reg="chr15:24954913-24955009"
tag="chr15_24954913_24955009"

biscuit pileup \
    -g ${reg} \
    -@ 40 \
    ${FAST}/hg38_noContig.fa \
    ${DATA}/Aneb.bam | \
biscuit vcf2bed -t snp -k 1 - > Aneb.${tag}.snp.bed

biscuit pileup \
    -g ${reg} \
    -@ 40 \
    ${FAST}/hg38_noContig.fa \
    ${DATA}/Aswift.bam | \
biscuit vcf2bed -t snp -k 1 - > Aswift.${tag}.snp.bed

# If there are SNPs in region, then make the epiBEDs
if [ `cat Aneb.${tag}.snp.bed | head -n1 | wc -l` -eq "1" ]; then
    biscuit epiread \
        -g ${reg} \
        -@ 40 \
        -B Aneb.${tag}.snp.bed \
        ${FAST}/hg38_noContig.fa \
        ${DATA}/Aneb.bam | \
    sort -k1,1 -k2,2n > Aneb.${tag}.epibed

    if [ `cat Aneb.${tag}.epibed | head -n1 | wc -l` -eq "1" ]; then
        bgzip Aneb.${tag}.epibed
    else
        rm Aneb.${tag}.epibed
    fi
else
    rm Aneb.${tag}.snp.bed
fi

if [ `cat Aswift.${tag}.snp.bed | head -n1 | wc -l` -eq "1" ]; then
    biscuit epiread \
        -g ${reg} \
        -@ 40 \
        -B Aswift.${tag}.snp.bed \
        ${FAST}/hg38_noContig.fa \
        ${DATA}/Aswift.bam | \
    sort -k1,1 -k2,2n > Aswift.${tag}.epibed

    if [ `cat Aswift.${tag}.epibed | head -n1 | wc -l` -eq "1" ]; then
        bgzip Aswift.${tag}.epibed
    else
        rm Aswift.${tag}.epibed
    fi
else
    rm Aswift.${tag}.snp.bed
fi

# Create a combined epiBED if possible. If only one epiBED exists, the rename
if [[ -f Aneb.${tag}.epibed.gz && -f Aswift.${tag}.epibed.gz ]]; then
    zcat Aneb.${tag}.epibed.gz Aswift.${tag}.epibed.gz | \
    sort -k1,1 -k2,2n | \
    uniq | \
    bgzip > Asample.${tag}.epibed.gz

    rm Aneb.${tag}.epibed.gz Aswift.${tag}.epibed.gz
elif [[ -f Aneb.${tag}.epibed.gz && ! -f Aswift.${tag}.epibed.gz ]]; then
    mv Aneb.${tag}.epibed.gz Asample.${tag}.epibed.gz
elif [[ ! -f Aneb.${tag}.epibed.gz && -f Aswift.${tag}.epibed.gz ]]; then
    mv Aswift.${tag}.epibed.gz Asample.${tag}.epibed.gz
else
    exit
fi

# Index if we have an epiBED to process
tabix -p bed Asample.${tag}.epibed.gz

# Clean up
if [ -f Aneb.${tag}.snp.bed ]; then
    rm Aneb.${tag}.snp.bed
fi

if [ -f Aswift.${tag}.snp.bed ]; then
    rm Aswift.${tag}.snp.bed
fi
